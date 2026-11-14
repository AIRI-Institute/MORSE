"""Create a portable schema-5 experiment from a CDS body FASTA."""
from __future__ import annotations

from pathlib import Path
import tempfile

import yaml

from morse.codon import normalize_dna, translate_dna
from morse.fasta import FastaRecord, read_fasta, write_fasta


def create_experiment(
    cds_path: Path,
    output_dir: Path,
    *,
    five_prime_utrs: list[Path] | None = None,
    three_prime_utrs: list[Path] | None = None,
    seed: int = 43,
    generations: int = 10,
    population_size: int = 16,
) -> Path:
    """Copy a CDS body and UTRs; configure AUG and UAA for assembly."""
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(f"Experiment directory already exists: {output_dir}")
    record = read_fasta(cds_path)
    coding = normalize_dna(record.sequence)
    if not coding or len(coding) % 3:
        raise ValueError("Input CDS body must be non-empty and have a length divisible by three")
    if "*" in translate_dna(coding):
        raise ValueError("Input CDS body must not contain stop codons; supply it without start/stop codons")
    if population_size < 2 or generations < 1:
        raise ValueError("population_size must be >= 2 and generations must be >= 1")

    libraries = {}
    for part, paths in (("five_prime_utr", five_prime_utrs), ("three_prime_utr", three_prime_utrs)):
        libraries[part] = [FastaRecord(f"utr_{index}", read_fasta(path).sequence)
                           for index, path in enumerate(paths or [], start=1)]
    chunks = {}
    if libraries["five_prime_utr"]:
        chunks["five_prime_utr"] = libraries["five_prime_utr"][0].sequence
    chunks["cds"] = coding
    if libraries["three_prime_utr"]:
        chunks["three_prime_utr"] = libraries["three_prime_utr"][0].sequence
    coordinates = {}
    cursor = 0
    for part, sequence in chunks.items():
        coordinates[part] = [cursor, cursor + len(sequence)]
        cursor += len(sequence)
    changes = {"cds": "synonymous"}
    for part, records in libraries.items():
        if len(records) > 1:
            changes[part] = {"library": f"{part}.fasta"}
    config = {
        "schema_version": 5,
        "name": output_dir.name,
        "output_dir": "results",
        "construct": {"fasta": "sequence.fasta", "coordinates": coordinates,
                      "start_codon": ["AUG"], "stop_codon": ["UAA"], "changes": changes},
        "optimizer": {
            "algorithm": "nsga2", "population_size": population_size,
            "generations": generations, "seed": seed,
        },
        "scorers": [{
            "name": "sequence_metrics",
            "command": '"${MORSE_PYTHON}" -m morse.resources.simple_scorer',
            "scope": ["start_codon", "cds_body"],
            "metrics": {
                "cpg_per_kb": {"role": "optimize", "direction": "min"},
                "gc_distance": {"role": "optimize", "direction": "min"},
                "gc_fraction": {"role": "log"},
            },
        }],
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".morse-experiment-", dir=output_dir.parent) as temporary:
        stage = Path(temporary) / "experiment"
        stage.mkdir()
        write_fasta(stage / "sequence.fasta", [FastaRecord(record.header, "".join(chunks.values()).replace("T", "U"))])
        for part, records in libraries.items():
            if len(records) > 1:
                write_fasta(stage / f"{part}.fasta", records)
        (stage / "README.md").write_text(
            f"# {output_dir.name}\n\nRNA input contains the CDS body and any supplied UTRs. "
            "AUG and UAA are added from the config.\n\n"
            "Run with `morse run --config experiment.yaml`.\n"
        )
        (stage / "experiment.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        from morse.run_spec import load_run_spec
        load_run_spec(stage / "experiment.yaml")
        stage.rename(output_dir)
    return output_dir / "experiment.yaml"
