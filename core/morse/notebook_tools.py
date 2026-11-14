"""Notebook forms and input bundles for the MORSE Python API.

Scorers are prepared independently. This module neither installs model runtimes
nor compiles installation IDs into experiment configuration. Widgets are an
optional dependency, imported only when the editor is constructed.
"""
from __future__ import annotations

import copy
import getpass
import json
from pathlib import Path
import shutil
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import yaml

from morse.codon import normalize_dna, translate_dna
from morse.construct import PART_ORDER, STOP_CODONS
from morse.input_bundle import input_references as _input_references, relative_input_name as _safe_input_name
from morse.run_spec import load_run_spec
from morse.scorer_manifest import load_scorer_manifest, request_parameter_specs
from morse.visualization import render_sequence_scope_svg

EXPECTED_BACKEND_API_VERSION = 6
TERMINAL_JOB_STATUSES = {"succeeded", "failed", "canceled", "interrupted"}
_OPERATORS = {part: "choose_one" for part in PART_ORDER} | {"cds": "synonymous_codons"}


class GeneratedRunConfiguration(dict):
    """Ordinary YAML mapping plus source files retained only in Python."""

    def __init__(self, configuration, *, input_files=None):
        super().__init__(configuration)
        self.input_files = dict(input_files or {})


yaml.SafeDumper.add_representer(
    GeneratedRunConfiguration,
    lambda dumper, value: dumper.represent_dict(dict(value)),
)


def _fasta_text(records):
    return "".join(f">{name}\n{sequence.replace('T', 'U')}\n" for name, sequence in records.items())


def export_configuration_bundle(path: str | Path) -> Path:
    """Download a YAML and its declared sequence/library inputs together."""
    path = Path(path).expanduser().resolve()
    payload = MorseBackendClient._job_payload(configuration=read_configuration(path))
    target = path.with_suffix(".zip")
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        bundle.writestr(path.name, yaml.safe_dump(payload["configuration"], sort_keys=False))
        for name, content in payload.get("input_files", {}).items():
            bundle.writestr(name, content)
    return target


def design_from_rna(sequence_text, coordinates, *, name, five_prime_utrs="", three_prime_utrs="",
                    poly_a_lengths=(), mutable_parts=("cds",), allow_upstream_start_codons=False,
                    start_codons=("AUG",), stop_codons=("UAA",), schema_version=5):
    records = parse_fasta_text(sequence_text)
    if len(records) != 1:
        raise ValueError("Provide exactly one RNA sequence")
    sequence = next(iter(records.values()))
    order = ("five_prime_utr", "cds", "three_prime_utr", "poly_a")
    if not isinstance(coordinates, dict) or "cds" not in coordinates or set(coordinates) - set(order):
        raise ValueError("Coordinates require cds and may include UTRs and poly_a")
    chunks, cursor = {}, 0
    for part in order:
        if part not in coordinates:
            chunks[part] = ""
            continue
        bounds = coordinates[part]
        if not isinstance(bounds, list) or len(bounds) != 2 or any(type(v) is not int for v in bounds):
            raise ValueError(f"{part}: use [start, end] integer coordinates")
        start, end = bounds
        if start != cursor or not start < end <= len(sequence):
            raise ValueError("Part coordinates must cover the RNA consecutively without gaps")
        chunks[part], cursor = sequence[start:end], end
    if cursor != len(sequence) or set(chunks["poly_a"]) - {"A"}:
        raise ValueError("Coordinates must cover the RNA; poly_a must contain only A")
    lengths = list(poly_a_lengths) or [len(chunks["poly_a"])]
    if schema_version not in {4, 5}:
        raise ValueError("Notebook sequence forms support schema 4 or 5")
    builder = design_from_cds_body if schema_version == 5 else design_from_full_cds
    codons = {"start_codons": start_codons, "stop_codons": stop_codons} if schema_version == 5 else {}
    design = builder(chunks["cds"], name=name, **codons,
        five_prime_utrs=five_prime_utrs or chunks["five_prime_utr"],
        three_prime_utrs=three_prime_utrs or chunks["three_prime_utr"],
        poly_a_lengths=lengths, mutable_parts=mutable_parts,
        allow_upstream_start_codons=allow_upstream_start_codons)
    design["reference_parts"] = chunks
    return design


def parse_fasta_text(text: str, *, prefix: str = "sequence", allow_empty: bool = False) -> dict[str, str]:
    """Read named choices from FASTA, or one unlabelled DNA/RNA sequence."""
    records: dict[str, str] = {}
    name, chunks = prefix, []

    def finish() -> None:
        if chunks:
            if name in records:
                raise ValueError(f"Duplicate sequence name: {name}")
            records[name] = normalize_dna("".join(chunks))

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            name, chunks = line[1:].strip() or f"{prefix}_{len(records)+1}", []
        else:
            chunks.append(line)
    finish()
    if not records:
        if allow_empty:
            return {prefix: ""}
        raise ValueError("Provide a DNA/RNA sequence or FASTA")
    return records


def design_from_full_cds(
    full_cds: str, *, name: str, five_prime_utrs: str = "", three_prime_utrs: str = "",
    poly_a_lengths: Sequence[int] = (0,), mutable_parts: Sequence[str] = ("cds",),
    allow_upstream_start_codons: bool = False,
) -> dict[str, Any]:
    records = parse_fasta_text(full_cds, prefix="cds")
    if len(records) != 1:
        raise ValueError("Provide exactly one full CDS, including its ATG start")
    coding = next(iter(records.values()))
    if not coding.startswith("ATG") or len(coding) % 3 or len(coding) < 6:
        raise ValueError("Full CDS must start with ATG and contain at least two codons")
    stop = coding[-3:] if coding[-3:] in STOP_CODONS else ""
    body = coding[3:-3] if coding[-3:] in STOP_CODONS else coding[3:]
    design = build_design_payload(
        name=name, five_prime_utrs=parse_fasta_text(five_prime_utrs, prefix="five", allow_empty=True),
        start_codons=["ATG"], cds=body, stop_codons=[stop],
        three_prime_utrs=parse_fasta_text(three_prime_utrs, prefix="three", allow_empty=True),
        poly_a_lengths=poly_a_lengths, mutable_parts=mutable_parts,
        allow_upstream_start_codons=allow_upstream_start_codons,
    )
    design["schema_version"] = 4
    return design


def design_from_cds_body(
    cds_body: str, *, name: str, five_prime_utrs: str = "", three_prime_utrs: str = "",
    start_codons: Sequence[str] = ("AUG",), stop_codons: Sequence[str] = ("UAA",),
    poly_a_lengths: Sequence[int] = (0,), mutable_parts: Sequence[str] = ("cds",),
    allow_upstream_start_codons: bool = False,
) -> dict[str, Any]:
    """Read one schema-5 CDS body, keeping its first codon intact."""
    records = parse_fasta_text(cds_body, prefix="cds_body")
    if len(records) != 1:
        raise ValueError("Provide exactly one CDS body without start or stop codons")
    design = build_design_payload(
        name=name, five_prime_utrs=parse_fasta_text(five_prime_utrs, prefix="five", allow_empty=True),
        start_codons=start_codons, cds=next(iter(records.values())), stop_codons=stop_codons,
        three_prime_utrs=parse_fasta_text(three_prime_utrs, prefix="three", allow_empty=True),
        poly_a_lengths=poly_a_lengths, mutable_parts=mutable_parts,
        allow_upstream_start_codons=allow_upstream_start_codons,
    )
    if any(not value for value in design["stop_codons"]):
        raise ValueError("Schema 5 requires a terminal stop codon")
    design["schema_version"] = 5
    return design


def build_design_payload(
    *, name: str, five_prime_utrs: Mapping[str, str], start_codons: Sequence[str], cds: str,
    stop_codons: Sequence[str] | None, three_prime_utrs: Mapping[str, str],
    poly_a_lengths: Sequence[int], mutable_parts: Sequence[str] = PART_ORDER,
    allow_upstream_start_codons: bool = False,
) -> dict[str, Any]:
    body = normalize_dna(cds)
    if not body or len(body) % 3 or "*" in translate_dna(body):
        raise ValueError("CDS body must be in frame and contain no stop codons")
    starts = [normalize_dna(v) for v in start_codons]
    stops = [normalize_dna(v) for v in (stop_codons or ())]
    if not starts or len(starts) != len(set(starts)) or any(len(v) != 3 or v in STOP_CODONS for v in starts):
        raise ValueError("Provide non-stop start codons of length three")
    if not stops or len(stops) != len(set(stops)) or any(v not in STOP_CODONS | {""} for v in stops):
        raise ValueError("Provide TAA, TAG, TGA or an absent stop codon")
    lengths = [int(v) for v in poly_a_lengths]
    if not lengths or any(v < 0 for v in lengths):
        raise ValueError("poly(A) lengths must be nonnegative integers")
    mutable = list(mutable_parts)
    if set(mutable) - set(PART_ORDER) or len(set(mutable)) != len(mutable):
        raise ValueError("Choose distinct known mutable parts")
    return {
        "name": str(name), "five_prime_utrs": {str(k): normalize_dna(v) for k, v in five_prime_utrs.items()},
        "start_codons": starts, "cds": body, "stop_codons": stops,
        "three_prime_utrs": {str(k): normalize_dna(v) for k, v in three_prime_utrs.items()},
        "poly_a_lengths": lengths, "mutable_parts": mutable,
        "allow_upstream_start_codons": bool(allow_upstream_start_codons),
    }


def build_notebook_run_configuration(
    design: Mapping[str, Any], *, run_name: str, optimizer_config: Mapping[str, Any],
    scorers: Sequence[Mapping[str, Any]] | None = None,
    manufacturability: Mapping[str, Any] | None = None, initialization: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a YAML mapping and its small FASTA input bundle."""
    version = design.get("schema_version", 5)
    if version not in {4, 5}:
        raise ValueError("Notebook sequence forms support schema 4 or 5")
    if version == 4 and design["start_codons"][0] != "ATG":
        raise ValueError("The reference CDS must start with ATG")
    mutable = set(design["mutable_parts"])
    five, three = dict(design["five_prime_utrs"]), dict(design["three_prime_utrs"])
    chunks = {
        "five_prime_utr": next(iter(five.values())),
        "cds": (design["cds"] if version == 5 else
                design["start_codons"][0] + design["cds"] + design["stop_codons"][0]),
        "three_prime_utr": next(iter(three.values())),
        "poly_a": "A" * design["poly_a_lengths"][0],
    }
    chunks.update(design.get("reference_parts", {}))
    coordinates, changes, files = {}, {}, {}
    cursor = 0
    for part, sequence in chunks.items():
        if sequence:
            coordinates[part] = [cursor, cursor + len(sequence)]
        cursor += len(sequence)
    files["sequence.fasta"] = _fasta_text({str(run_name): "".join(chunks.values())})
    if "cds" in mutable:
        changes["cds"] = "synonymous"
    for part, choices in (("five_prime_utr", five), ("three_prime_utr", three)):
        if part in mutable and any(choices.values()):
            filename = part + "s.fasta"
            if any(not value for value in choices.values()):
                raise ValueError("UTR libraries cannot mix empty and nonempty choices")
            files[filename] = _fasta_text(choices)
            changes[part] = {"library": filename}
    if "poly_a" in mutable and any(design["poly_a_lengths"]):
        changes["poly_a"] = {"lengths": list(design["poly_a_lengths"])}
    for part, values in (("start_codon", design["start_codons"]), ("stop_codon", design["stop_codons"])):
        if version == 4 and part in mutable and len(values) > 1:
            changes[part] = {"sequences": list(values)}
    optimizer = copy.deepcopy(dict(optimizer_config))
    optimizer["algorithm"] = optimizer.pop("type", optimizer.get("algorithm", "nsga2"))
    optimizer.pop("variables", None)
    configuration = {
        "schema_version": version, "name": str(run_name), "output_dir": "results/" + str(run_name),
        "construct": {"fasta": "sequence.fasta", "coordinates": coordinates, "changes": changes},
        "optimizer": optimizer,
    }
    if version == 5:
        for part in ("start_codon", "stop_codon"):
            values = list(design[part + "s"])
            if part == "stop_codon" and any(not value for value in values):
                raise ValueError("Schema 5 requires a terminal stop codon")
            configuration["construct"][part] = [value.replace("T", "U") for value in values]
    if design.get("allow_upstream_start_codons", False):
        configuration["construct"]["allow_upstream_start_codons"] = True
    if scorers is not None:
        configuration["scorers"] = copy.deepcopy(list(scorers))
    for key, value in (("manufacturability", manufacturability), ("initialization", initialization)):
        if value is not None:
            configuration[key] = copy.deepcopy(dict(value))
    return GeneratedRunConfiguration(configuration, input_files=files)


def write_generated_notebook_config(path: str | Path, configuration: Mapping[str, Any], *, input_files=None):
    """Validate the actual saved YAML through the same loader as morse run."""
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    files = dict(getattr(configuration, "input_files", {}) if input_files is None else input_files)
    referenced = {mapping[key] for mapping, key in _input_references(configuration)}
    files = {name: content for name, content in files.items() if name in referenced}
    text = yaml.safe_dump(dict(configuration), sort_keys=False, allow_unicode=True)
    candidate = path.with_name(path.name + ".validating.yaml")
    try:
        # Validate new inputs before replacing any previously saved input file.
        with tempfile.TemporaryDirectory(prefix=".morse-input-check-", dir=path.parent) as temporary:
            staged = copy.deepcopy(dict(configuration))
            replacements = {}
            for name, content in files.items():
                target = Path(temporary) / _safe_input_name(name)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(content)
                replacements[name] = str(target)
            for mapping, key in _input_references(staged):
                if mapping[key] in replacements:
                    mapping[key] = replacements[mapping[key]]
            staged_path = Path(temporary) / "experiment.yaml"
            staged_path.write_text(yaml.safe_dump(staged, sort_keys=False))
            load_run_spec(staged_path, input_root=path.parent)
        for name, content in files.items():
            target = path.parent / _safe_input_name(name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content)
        candidate.write_text(text)
        load_run_spec(candidate)
        candidate.replace(path)
    finally:
        candidate.unlink(missing_ok=True)
    return load_run_spec(path)


def render_generated_scope_preview(path: str | Path) -> str:
    return render_sequence_scope_svg(load_run_spec(Path(path)))


def read_configuration(path: str | Path) -> dict[str, Any]:
    path = Path(path).expanduser().resolve()
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict) or raw.get("schema_version") not in {3, 4, 5}:
        raise ValueError("Expected one complete schema_version: 3, 4 or 5 YAML configuration")
    files = {}
    for mapping, key in _input_references(raw):
        value = mapping[key]
        source = Path(value).expanduser()
        source = source if source.is_absolute() else path.parent / source
        if source.is_file():
            files[value] = source.read_text()
    return GeneratedRunConfiguration(raw, input_files=files)


def default_cpu_scorer() -> dict[str, Any]:
    return {
        "name": "sequence_metrics", "command": '"${MORSE_PYTHON}" -m morse.resources.simple_scorer',
        "scope": ["start_codon", "cds_body"],
        "metrics": {"cpg_per_kb": {"role": "optimize", "direction": "min"},
                    "gc_distance": {"role": "optimize", "direction": "min"}, "gc_fraction": {"role": "log"}},
    }


def load_manifest_catalog(root: str | Path) -> dict[str, dict[str, Any]]:
    """Inspect metadata only; no scorer module imports or model loading."""
    from morse.schema4 import scope_from_internal
    path = Path(root)
    if path.name == "core":
        path = path.parent
    result = {}
    for manifest_path in sorted((path / "scorers").glob("*/morse-scorer*.yaml")):
        manifest = load_scorer_manifest(manifest_path)
        metrics = {}
        for name, output in manifest["outputs"].items():
            policy = copy.deepcopy(output["default_policy"])
            policy.pop("type", None)
            if policy.get("role") == "objective":
                policy["role"] = "optimize"
            if output["type"] != "float":
                policy["value_type"] = output["type"]
            metrics[name] = policy
        entry = {
            "name": manifest["id"], "command": "REPLACE_WITH_YOUR_SCORER_COMMAND",
            "timeout_seconds": manifest["runtime"]["timeout_seconds"],
            "scope": scope_from_internal(manifest["scope"]["default"]), "metrics": metrics,
        }
        parameters = {
            name: copy.deepcopy(spec["default"])
            for name, spec in request_parameter_specs(manifest).items()
            if "default" in spec
        }
        if parameters:
            entry["parameters"] = parameters
        result[manifest["id"]] = entry
    return result


class ScorerEditor:
    """A prepared command, calculation parameters, scope and metric policies."""
    def __init__(self, scorer: Mapping[str, Any]):
        import ipywidgets as w
        from morse.schema4 import scope_expression
        self.name = w.Text(value=str(scorer["name"]), description="Name")
        self.command = w.Textarea(value=str(scorer["command"]), description="Bash", layout=w.Layout(width="100%", height="70px"))
        from morse.config import scorer_parameters_from_mapping
        parameters = scorer_parameters_from_mapping(scorer)
        self.parameters = w.Textarea(
            value=yaml.safe_dump(parameters, sort_keys=False, allow_unicode=True) if parameters else "",
            description="Parameters", placeholder="parameter_name: value (optional)",
            layout=w.Layout(width="100%", height="80px"),
        )
        scope = scorer["scope"]
        self.scope = w.Text(value=scope_expression(scope), description="Scope", layout=w.Layout(width="100%"))
        self.timeout = w.IntText(value=int(scorer.get("timeout_seconds", 120)), description="Timeout (s)")
        self.enabled = w.Checkbox(value=True, description="Include scorer")
        self.metrics = {}
        rows = [w.HBox([w.Label("Metric", layout=w.Layout(width="190px")),
                        w.Label("Role", layout=w.Layout(width="125px")), w.Label("Direction / bounds")])]
        for name, policy in scorer["metrics"].items():
            role = w.Dropdown(options=["optimize", "threshold", "log"], value=policy["role"], layout=w.Layout(width="125px"))
            direction = w.Dropdown(options=["min", "max"], value=policy.get("direction", "max"), layout=w.Layout(width="85px"))
            lower = w.Text(value=str(policy.get("min", "")), placeholder="min", layout=w.Layout(width="85px"))
            upper = w.Text(value=str(policy.get("max", "")), placeholder="max", layout=w.Layout(width="85px"))
            equals = w.Text(value=yaml.safe_dump(policy["equals"], default_flow_style=True).splitlines()[0] if "equals" in policy else "", placeholder="equals", layout=w.Layout(width="85px"))
            row = (role, direction, lower, upper, equals, policy.get("value_type", "float"))
            self.metrics[name] = row

            def update(change=None, controls=row):
                r, d, lo, hi, eq, _ = controls
                d.disabled = r.value != "optimize"
                lo.disabled = hi.disabled = eq.disabled = r.value != "threshold"
                d.layout.display = "" if r.value == "optimize" else "none"
                for control in (lo, hi, eq):
                    control.layout.display = "" if r.value == "threshold" else "none"
            role.observe(update, names="value")
            update()
            rows.append(w.HBox([w.Label(name, layout=w.Layout(width="190px")), role, direction, lower, upper, equals]))
        self.widget = w.VBox([self.enabled, self.name, self.command, self.parameters, self.scope, self.timeout, *rows])

    def configuration(self) -> dict[str, Any]:
        from morse.config import parse_scorer_parameters
        from morse.schema4 import scope_to_internal
        if not self.command.value.strip() or self.command.value == "REPLACE_WITH_YOUR_SCORER_COMMAND":
            raise ValueError(f"Set the prepared Bash command for {self.name.value}")
        metrics = {}
        for name, (role, direction, lower, upper, equals, value_type) in self.metrics.items():
            metric = {"role": role.value}
            if role.value == "optimize":
                metric["direction"] = direction.value
            elif role.value == "threshold":
                for key, control in (("min", lower), ("max", upper), ("equals", equals)):
                    if control.value.strip():
                        metric[key] = yaml.safe_load(control.value)
            if value_type != "float":
                metric["value_type"] = value_type
            metrics[name] = metric
        scope_to_internal(self.scope.value)
        result = {"name": self.name.value, "command": self.command.value, "timeout_seconds": self.timeout.value,
                  "scope": self.scope.value.strip(), "metrics": metrics}
        parameters = parse_scorer_parameters(yaml.safe_load(self.parameters.value) if self.parameters.value.strip() else {})
        if parameters:
            result["parameters"] = parameters
        return result


class ExperimentEditor:
    """Small widget form that produces an editable, downloadable ordinary YAML."""
    def __init__(self, repository: str | Path):
        import ipywidgets as w
        self.name = w.Text(value="morse_notebook_run", description="Run name")
        self.input_mode = w.Dropdown(options=["CDS", "RNA + coordinates"], value="CDS", description="Input")
        self.cds = w.Textarea(value=">example\nGCUGCUGGUGCUGCU", description="CDS body", layout=w.Layout(width="100%", height="100px"))
        self.cds_upload = w.FileUpload(accept=".fa,.fasta,.fna,.txt", multiple=False, description="Upload FASTA")
        self.coordinates = w.Textarea(value="cds: [0, 15]", description="Coordinates", layout=w.Layout(width="100%", height="95px"))
        self.coordinates.layout.display = "none"
        self.start_codons = w.Text(value="AUG", description="Start choices", placeholder="AUG, CUG")
        self.stop_codons = w.Text(value="UAA", description="Stop choices", placeholder="UAA, UAG, UGA")
        self.input_mode.observe(self._input_mode_changed, names="value")
        self.five = w.Textarea(value="", description="5′UTR choices", layout=w.Layout(width="100%", height="80px"))
        self.three = w.Textarea(value="", description="3′UTR choices", layout=w.Layout(width="100%", height="80px"))
        self.poly_a = w.Text(value="", placeholder="Optional, e.g. 70, 100", description="poly(A) lengths")
        self.mutable = w.SelectMultiple(options=["five_prime_utr", "cds", "three_prime_utr", "poly_a"], value=("cds",), description="Change parts")
        self.upstream_starts = w.Checkbox(value=False, description="Allow upstream AUG")
        self.algorithm = w.Dropdown(options=["ga", "nsga2", "nsga3"], value="nsga2", description="Algorithm")
        self.population = w.IntText(value=32, description="Population")
        self.generations = w.IntText(value=20, description="Generations")
        self.seed = w.IntText(value=43, description="Seed")
        self.output_dir = w.Text(value="results/morse_notebook_run", description="Output folder", layout=w.Layout(width="100%"))
        self.advanced = w.Textarea(value="", placeholder="optimizer:\n  initial_population: /path/to/population.json", description="Advanced YAML", layout=w.Layout(width="100%", height="90px"))
        self.catalog = {"sequence_metrics": default_cpu_scorer(), **load_manifest_catalog(repository)}
        self.catalog_choice = w.Dropdown(options=list(self.catalog), description="Add scorer")
        add = w.Button(description="Add scorer")
        self.scorers: list[ScorerEditor] = []
        self.scorer_box = w.VBox()
        self.add_scorer(default_cpu_scorer())
        add.on_click(lambda _: self.add_scorer(self.catalog[self.catalog_choice.value]))
        generate = w.Button(description="Generate YAML", button_style="primary")
        self.yaml = w.Textarea(layout=w.Layout(width="100%", height="400px"), description="Run YAML")
        self.messages = w.Output()
        generate.on_click(self._generate)
        self.cds_upload.observe(self._uploaded, names="value")
        panels = w.Tab(children=[
            w.VBox([w.HTML("Input CDS excludes start and stop. Coordinates refer to input parts (0-based, end-exclusive). Codon choices: comma-separated; first is the reference. UTR fields accept FASTA choices."),
                self.name, self.input_mode, self.cds_upload, self.cds, self.coordinates, self.start_codons, self.stop_codons, self.five, self.three, self.poly_a, self.mutable, self.upstream_starts]),
            w.VBox([w.HTML("Bash reads one JSON batch from stdin then EOF; diagnostics go to stderr."),
                w.HBox([self.catalog_choice, add]), self.scorer_box]),
            w.VBox([w.HBox([self.algorithm, self.population, self.generations, self.seed]), self.output_dir,
                w.HTML("Advanced YAML can add initialization or override optimizer fields."), self.advanced]),
        ])
        for index, title in enumerate(("1 · Design", "2 · Scorers", "3 · Search")):
            panels.set_title(index, title)
        self.widget = w.VBox([panels, generate, self.messages, self.yaml])
        self.yaml.value = yaml.safe_dump(self.configuration(), sort_keys=False, allow_unicode=True)

    def _uploaded(self, change):
        values = change["new"]
        records = values.values() if isinstance(values, dict) else values
        for item in records:
            self.cds.value = bytes(item["content"]).decode("utf-8-sig")

    def _input_mode_changed(self, change):
        whole = change["new"] == "RNA + coordinates"
        self.cds.description = "RNA input" if whole else "CDS body"
        self.coordinates.layout.display = "" if whole else "none"

    def add_scorer(self, scorer: Mapping[str, Any]) -> ScorerEditor:
        editor = ScorerEditor(scorer)
        self.scorers.append(editor)
        self.scorer_box.children = tuple(s.widget for s in self.scorers)
        return editor

    def configuration(self) -> dict[str, Any]:
        lengths = [int(s.strip()) for s in self.poly_a.value.split(",") if s.strip()]
        options = dict(name=self.name.value, five_prime_utrs=self.five.value, three_prime_utrs=self.three.value,
            poly_a_lengths=lengths or ([0] if self.input_mode.value == "CDS" else []),
            mutable_parts=self.mutable.value, allow_upstream_start_codons=self.upstream_starts.value,
            start_codons=[value.strip() for value in self.start_codons.value.split(",") if value.strip()],
            stop_codons=[value.strip() for value in self.stop_codons.value.split(",") if value.strip()])
        design = (design_from_rna(self.cds.value, yaml.safe_load(self.coordinates.value), **options)
                  if self.input_mode.value == "RNA + coordinates" else design_from_cds_body(self.cds.value, **options))
        optimizer = {"type": self.algorithm.value, "population_size": self.population.value,
                     "generations": self.generations.value, "seed": self.seed.value}
        if self.algorithm.value == "nsga3":
            optimizer["reference_directions"] = {"method": "energy", "n_points": self.population.value, "seed": self.seed.value}
        result = build_notebook_run_configuration(design, run_name=self.name.value, optimizer_config=optimizer,
            scorers=[s.configuration() for s in self.scorers if s.enabled.value])
        result["output_dir"] = self.output_dir.value.strip()
        self._input_files = result.input_files
        if self.advanced.value.strip():
            extra = yaml.safe_load(self.advanced.value)
            if not isinstance(extra, dict) or set(extra) - {"initialization", "optimizer", "manufacturability", "construct"}:
                raise ValueError("Advanced YAML accepts initialization, optimizer, construct and manufacturability")
            for key, value in extra.items():
                if key in {"optimizer", "construct"}:
                    result[key].update(value)
                else:
                    result[key] = value
        self._generated_construct = copy.deepcopy(result["construct"])
        return result

    def _generate(self, _):
        with self.messages:
            self.messages.clear_output()
            try:
                self.yaml.value = yaml.safe_dump(self.configuration(), sort_keys=False, allow_unicode=True)
                print("YAML generated. Review or edit it, then save and validate below.")
            except Exception as exc:
                print(f"Configuration error: {exc}")

    def save(self, path: str | Path):
        raw = yaml.safe_load(self.yaml.value)
        if not isinstance(raw, dict):
            raise ValueError("Run YAML must contain a mapping")
        files = self._input_files
        if raw.get("construct") != self._generated_construct:
            parent = Path(path).expanduser().resolve().parent
            files = {name: content for name, content in files.items() if not (parent / name).is_file()}
        return write_generated_notebook_config(path, raw, input_files=files)


class MorseBackendClient:
    """HTTP client for authenticated job execution, progress and artifacts."""
    def __init__(self, base_url: str, token: str, timeout_seconds: int = 60):
        self.base_url, self.token, self.timeout_seconds = base_url.rstrip("/"), token, timeout_seconds

    @staticmethod
    def _job_payload(*, config=None, configuration=None):
        if (config is None) == (configuration is None):
            raise ValueError("Provide exactly one of config or configuration")
        if config is not None:
            return {"config": str(config)}
        raw = copy.deepcopy(dict(configuration))
        available = getattr(configuration, "input_files", {})
        files, replacements = {}, {}
        for mapping, key in _input_references(raw):
            value = mapping[key]
            if value in replacements:
                mapping[key] = replacements[value]
                continue
            content = available.get(value)
            source = Path(value).expanduser()
            if content is None and source.is_absolute() and source.is_file():
                content = source.read_text()
            if content is None:
                continue  # An already prepared path on the execution host.
            suffix = source.suffix.lower()
            name = f"inputs/input_{len(files)}{suffix}"
            files[name] = content
            replacements[value] = name
            mapping[key] = name
        if sum(len(content.encode("utf-8")) for content in files.values()) > 10 * 1024 * 1024:
            raise ValueError("Notebook sequence inputs exceed the 10 MiB upload limit")
        result = {"configuration": raw}
        if files:
            result["input_files"] = files
        return result

    def health(self):
        return self._json_request("GET", "/health", authenticated=False)

    def configs(self):
        return self._json_request("GET", "/v1/configs")["configs"]

    def jobs(self):
        return self._json_request("GET", "/v1/jobs")["jobs"]

    def preflight_job(self, *, config=None, configuration=None):
        return self._json_request("POST", "/v1/jobs/preflight", self._job_payload(config=config, configuration=configuration))

    def submit_job(self, *, config=None, configuration=None):
        return self._json_request("POST", "/v1/jobs", self._job_payload(config=config, configuration=configuration))

    def job(self, job_id):
        return self._json_request("GET", f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}")

    def job_metrics(self, job_id):
        return self._json_request("GET", f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}/metrics")

    def log(self, job_id, tail=80):
        return self._request("GET", f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}/log?tail={int(tail)}").decode()

    def cancel(self, job_id):
        return self._json_request("DELETE", f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}")

    def resume_job(self, job_id):
        return self._json_request("POST", f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}/resume", {})

    def sequence_scope_svg(self, job_id):
        return self._request("GET", f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}/sequence-scope").decode()

    def download_artifacts(self, job_id, target):
        path = Path(target)
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(
            self.base_url + f"/v1/jobs/{urllib.parse.quote(job_id, safe='')}/artifacts",
            headers={"Authorization": f"Bearer {self.token}"}, method="GET")
        try:
            with urllib.request.urlopen(request, timeout=max(300, self.timeout_seconds)) as response, path.open("wb") as output:
                shutil.copyfileobj(response, output)
        except urllib.error.HTTPError as exc:
            raise RuntimeError(f"Backend HTTP {exc.code}: {exc.read().decode(errors='replace')}") from exc
        return path

    def _json_request(self, method, path, payload=None, *, authenticated=True):
        result = json.loads(self._request(method, path, payload, authenticated=authenticated))
        if not isinstance(result, dict):
            raise RuntimeError("Backend returned a non-object JSON response")
        return result

    def _request(self, method, path, payload=None, *, authenticated=True, timeout=None):
        headers = {"Accept": "application/json"}
        if authenticated:
            headers["Authorization"] = f"Bearer {self.token}"
        data = None
        if payload is not None:
            headers["Content-Type"] = "application/json"
            data = json.dumps(payload).encode()
        request = urllib.request.Request(self.base_url + path, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(request, timeout=timeout or self.timeout_seconds) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode(errors="replace")
            raise RuntimeError(f"Backend HTTP {exc.code}: {detail}") from exc


class LocalMorseClient:
    """Use the backend's existing job manager locally, without an HTTP server."""
    def __init__(self, project_root: str | Path, state_root: str | Path):
        from morse.backend import BackendSettings, JobManager, runtime_metadata
        self.project_root = Path(project_root).expanduser().resolve()
        if (self.project_root / "core/configs").is_dir():
            self.project_root /= "core"
        self.state_root = Path(state_root).expanduser().resolve()
        self.state_root.mkdir(parents=True, exist_ok=True)
        settings = BackendSettings.from_values(config_root=self.state_root,
            results_root=self.state_root / "jobs", workdir=self.project_root, token="local-in-process")
        self._manager, self._runtime_metadata = JobManager(settings), runtime_metadata

    def health(self):
        return {"status": "ok", "execution_mode": "local-in-process", "runtime": self._runtime_metadata()}

    def configs(self):
        return self._manager.list_configs()

    def jobs(self):
        return self._manager.list_jobs()

    def preflight_job(self, *, config=None, configuration=None):
        return self._manager.preflight(MorseBackendClient._job_payload(config=config, configuration=configuration))

    def submit_job(self, *, config=None, configuration=None):
        return self._manager.submit(MorseBackendClient._job_payload(config=config, configuration=configuration))

    def job(self, job_id):
        return self._manager.get(job_id)

    def job_metrics(self, job_id):
        return self._manager.metrics(job_id)

    def log(self, job_id, tail=80):
        return self._manager.log_tail(job_id, tail)

    def cancel(self, job_id):
        return self._manager.cancel(job_id)

    def resume_job(self, job_id):
        return self._manager.resume(job_id)

    def sequence_scope_svg(self, job_id):
        return self._manager.sequence_scope_path(job_id).read_text()

    def download_artifacts(self, job_id, target):
        path = Path(target)
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self._manager.artifact_path(job_id), path)
        return path

    def shutdown(self):
        self._manager.shutdown()


def require_backend_compatibility(health: Mapping[str, Any]) -> None:
    actual = health.get("runtime", {}).get("backend_api_version")
    if actual != EXPECTED_BACKEND_API_VERSION:
        raise RuntimeError(f"Backend API {actual!r}; this notebook requires API {EXPECTED_BACKEND_API_VERSION}")


def colab_secret_or_prompt(name: str = "MORSE_BACKEND_TOKEN") -> str:
    try:
        from google.colab import userdata
        token = userdata.get(name)
    except Exception:
        token = None
    return token or getpass.getpass(f"{name}: ")


def monitor_job(client, job_id: str, *, poll_seconds: float = 3, max_seconds: float = 300):
    """Observe existing progress; interrupting monitoring does not cancel the job."""
    from IPython.display import clear_output, display
    import ipywidgets as w
    started = time.monotonic()
    while True:
        record = client.job(job_id)
        clear_output(wait=True)
        try:
            metrics = client.job_metrics(job_id)
            fraction = float(metrics.get("progress_fraction") or 0)
            display(w.FloatProgress(value=min(1, max(0, fraction)), min=0, max=1, description="Progress"))
            display({key: metrics.get(key) for key in ("status", "generation", "total_generations", "evaluated_candidates", "hypervolume", "hard_feasible_fraction")})
            figure = plot_live_metrics(metrics)
            if figure is not None:
                import matplotlib.pyplot as plt
                display(figure)
                plt.close(figure)
        except Exception as exc:
            print(f"Progress is not available yet: {exc}")
        print(client.log(job_id, tail=15))
        if record["status"] in TERMINAL_JOB_STATUSES or time.monotonic() - started >= max_seconds:
            return record
        time.sleep(max(0.1, poll_seconds))


def plot_live_metrics(metrics: Mapping[str, Any]):
    """Plot the optimizer's recorded history without recalculating statistics."""
    import matplotlib.pyplot as plt
    from morse.publication_style import INK, MUTED, palette_color
    history = [row for row in metrics.get("history", []) if isinstance(row, dict)]
    values = [(row["generation"], row["hypervolume"]) for row in history
              if row.get("hypervolume") is not None and row.get("generation") is not None]
    label = "Hypervolume"
    if not values:
        name = next((item.get("name") for row in history for item in row.get("objectives", []) if item.get("best") is not None), None)
        values = [(row["generation"], item["best"]) for row in history for item in row.get("objectives", [])
                  if item.get("name") == name and item.get("best") is not None]
        label = f"Best {name}" if name else ""
    if not values:
        return None
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": MUTED, "text.color": INK}):
        figure, axis = plt.subplots(figsize=(7, 2.6), constrained_layout=True)
        axis.plot(*zip(*values), color=palette_color(0), linewidth=1.5)
        axis.set(xlabel="Generation", ylabel=label)
        return figure


@dataclass(frozen=True)
class ResultArtifacts:
    root: Path
    report: Path | None
    pareto_csv: Path | None
    pareto_fasta: Path | None
    images: tuple[Path, ...]
    summary: Path | None = None


def discover_result_artifacts(root: str | Path) -> ResultArtifacts:
    root = Path(root)
    def one(name):
        return next(iter(sorted(root.rglob(name))), None)
    return ResultArtifacts(root, one("report.md"), one("pareto.csv"), one("pareto.fasta"),
        tuple(sorted(root.rglob("*.png"))), one("summary.json"))


def extract_artifacts(archive_path: str | Path, destination: str | Path) -> ResultArtifacts:
    target = Path(destination).resolve()
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as bundle:
        for item in bundle.infolist():
            if not (target / item.filename).resolve().is_relative_to(target):
                raise ValueError("Artifact ZIP contains an unsafe path")
        bundle.extractall(target)
    return discover_result_artifacts(target)


def display_result_artifacts(artifacts: ResultArtifacts) -> None:
    """Display already generated Core plots/report; do not rerun model scoring."""
    from IPython.display import display, Image, Markdown, FileLink
    if artifacts.report:
        # Relative image URLs in report.md belong to its result directory. Show
        # those files explicitly below so Colab does not resolve them as notebook
        # URLs; keep the original report available for download.
        text = "\n".join(line for line in artifacts.report.read_text().splitlines() if not line.lstrip().startswith("!["))
        display(Markdown(text))
        display(FileLink(str(artifacts.report)))
    for path in artifacts.images:
        display(Image(filename=str(path)))
    for path in (artifacts.pareto_csv, artifacts.pareto_fasta, artifacts.summary):
        if path:
            display(FileLink(str(path)))


def download_file_to_browser(path: str | Path):
    path = Path(path)
    try:
        from google.colab import files
    except ImportError:
        from IPython.display import display, FileLink
        display(FileLink(str(path)))
    else:
        files.download(str(path))
    return path
def build_notebook_optimizer_config(
    *,
    population_size: int,
    generations: int,
    seed: int,
    n_offsprings: int,
    reference_points: int,
    crossover_probability: float,
    mutation_probability: str | float | None,
    eliminate_duplicates: bool,
    constraint_method: str,
    initial_epsilon: float,
    relaxation_generations: int,
    epsilon_curve: str,
    early_stopping_enabled: bool,
    early_stopping_min_generations: int,
    early_stopping_patience: int,
    early_stopping_rel_tol: float,
    early_stopping_abs_tol: float,
) -> dict[str, Any]:
    """Build the complete NSGA-III section used by the notebook.

    The notebook deliberately does not inherit algorithm values from a checked-in
    run YAML. Every value that can affect an optimization is represented here.
    """

    mutation: str | float
    if (
        mutation_probability is None
        or str(mutation_probability).strip().lower() == "auto"
    ):
        mutation = "auto"
    else:
        mutation = float(mutation_probability)

    constraint = str(constraint_method).strip().lower()
    if constraint == "adaptive_epsilon":
        constraint_handling = {
            "method": "adaptive_epsilon",
            "normalization": "relative",
            "initial_epsilon": float(initial_epsilon),
            "relaxation_generations": int(relaxation_generations),
            "curve": str(epsilon_curve).strip().lower(),
        }
    elif constraint == "feasibility_first":
        constraint_handling = {
            "method": "feasibility_first",
            "normalization": "none",
            "initial_epsilon": 0.0,
            "relaxation_generations": 0,
            "curve": "linear",
        }
    else:
        raise ValueError(
            "constraint_method must be feasibility_first or adaptive_epsilon"
        )

    return {
        "type": "nsga3",
        "population_size": int(population_size),
        "generations": int(generations),
        "seed": int(seed),
        "n_offsprings": int(n_offsprings),
        "eliminate_duplicates": bool(eliminate_duplicates),
        "crossover": {
            "type": "uniform",
            "probability": float(crossover_probability),
        },
        "mutation": {
            "type": "integer_reset",
            "probability": mutation,
        },
        "reference_directions": {
            "method": "energy",
            "n_points": int(reference_points),
            "seed": int(seed),
        },
        "constraint_handling": constraint_handling,
        "early_stopping": {
            "enabled": bool(early_stopping_enabled),
            "metric": "hypervolume",
            "min_generations": int(early_stopping_min_generations),
            "patience": int(early_stopping_patience),
            "rel_tol": float(early_stopping_rel_tol),
            "abs_tol": float(early_stopping_abs_tol),
        },
    }
