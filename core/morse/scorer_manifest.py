from __future__ import annotations

import copy
import os
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml

from morse.scorer_scope import (
    PART_NAMES,
    normalize_scorer_scope,
    parse_scorer_scope,
    scope_parts,
    scorer_scope_expression,
)


MANIFEST_FILE_NAME = "morse-scorer.yaml"
SCORER_ID_PATTERN = re.compile(r"^[a-z][a-z0-9_-]{1,63}$")
OUTPUT_NAME_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
BUILD_CONTEXT_NAME_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,63}$")
SUPPORTED_OUTPUT_TYPES = {"float", "int", "bool", "string", "structure"}
SUPPORTED_PARAMETER_TYPES = {"float", "int", "bool", "string", "path"}
SUPPORTED_PARAMETER_OWNERS = {"experiment", "backend"}
ENVIRONMENT_NAME_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def load_scorer_manifest(path: str | Path) -> dict[str, Any]:
    manifest_path = Path(path)
    try:
        raw = yaml.safe_load(manifest_path.read_text())
    except OSError as exc:
        raise ValueError(f"Cannot read scorer manifest: {manifest_path}") from exc
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid scorer manifest YAML: {manifest_path}") from exc
    if not isinstance(raw, dict):
        raise ValueError("Scorer manifest must contain a YAML mapping")
    return validate_scorer_manifest(raw)


def validate_scorer_manifest(raw: Mapping[str, Any]) -> dict[str, Any]:
    unknown = sorted(
        set(raw)
        - {
            "schema_version",
            "id",
            "name",
            "version",
            "description",
            "runtime",
            "portable_runtime",
            "portable_parameters",
            "scope",
            "outputs",
            "parameters",
            "self_test",
        }
    )
    if unknown:
        raise ValueError("Unknown scorer manifest fields: " + ", ".join(unknown))
    if int(raw.get("schema_version", 0)) != 2:
        raise ValueError("Scorer manifest must define schema_version: 2")

    scorer_id = str(raw.get("id", "")).strip()
    if not SCORER_ID_PATTERN.fullmatch(scorer_id):
        raise ValueError(
            "Scorer id must start with a lowercase letter and contain only "
            "lowercase letters, numbers, underscore, or hyphen"
        )
    name = str(raw.get("name", "")).strip()
    version = str(raw.get("version", "")).strip()
    description = str(raw.get("description", "")).strip()
    if not name or not version or not description:
        raise ValueError("Scorer manifest requires name, version, and description")

    runtime = _validate_runtime(_mapping(raw.get("runtime"), "runtime"))
    portable_runtime_raw = raw.get("portable_runtime")
    portable_runtime = (
        None
        if portable_runtime_raw is None
        else _validate_runtime(
            _mapping(portable_runtime_raw, "portable_runtime")
        )
    )
    scope = _validate_scope_block(_mapping(raw.get("scope"), "scope"))
    outputs = _validate_outputs(_mapping(raw.get("outputs"), "outputs"))
    parameters = _validate_parameters(raw.get("parameters", {}))
    portable_parameters_raw = raw.get("portable_parameters")
    portable_parameters = (
        None
        if portable_parameters_raw is None
        else _validate_parameters(portable_parameters_raw)
    )
    if portable_parameters is not None and portable_runtime is None:
        raise ValueError(
            "portable_parameters requires portable_runtime"
        )
    self_test = _validate_self_test(
        _mapping(raw.get("self_test"), "self_test"),
        outputs,
        scope,
    )
    normalized = {
        "schema_version": 2,
        "id": scorer_id,
        "name": name,
        "version": version,
        "description": description,
        "runtime": runtime,
        "scope": scope,
        "outputs": outputs,
        "parameters": parameters,
        "self_test": self_test,
    }
    if portable_runtime is not None:
        normalized["portable_runtime"] = portable_runtime
    if portable_parameters is not None:
        normalized["portable_parameters"] = portable_parameters
    return normalized








def scorer_scope_label(manifest: Mapping[str, Any]) -> str:
    scope = _mapping(
        _mapping(manifest.get("scope"), "manifest.scope").get("default"),
        "manifest.scope.default",
    )
    return scorer_scope_expression(scope, "manifest.scope.default")


def _validate_runtime(raw: Mapping[str, Any]) -> dict[str, Any]:
    runtime_type = str(raw.get("type", "")).lower()
    if runtime_type not in {"venv", "container"}:
        raise ValueError("runtime.type must be venv or container")
    allowed_fields = {
        "type",
        "timeout_seconds",
        "env",
        *(
            {"entrypoint", "requirements", "system_site_packages"}
            if runtime_type == "venv"
            else {
                "image",
                "dockerfile",
                "build_contexts",
                "docker_args",
                "gpus",
            }
        ),
    }
    unknown = sorted(set(raw) - allowed_fields)
    if unknown:
        raise ValueError("Unknown runtime fields: " + ", ".join(unknown))
    normalized = dict(raw)
    normalized["type"] = runtime_type
    timeout = int(raw.get("timeout_seconds", 180))
    if timeout < 1:
        raise ValueError("runtime.timeout_seconds must be positive")
    normalized["timeout_seconds"] = timeout
    env = raw.get("env", {})
    if not isinstance(env, Mapping):
        raise ValueError("runtime.env must be a mapping")
    normalized["env"] = {str(name): str(value) for name, value in env.items()}

    if runtime_type == "venv":
        entrypoint = raw.get("entrypoint")
        if (
            not isinstance(entrypoint, list)
            or not entrypoint
            or not all(isinstance(item, str) and item for item in entrypoint)
        ):
            raise ValueError("venv runtime requires a non-empty entrypoint list")
        normalized["entrypoint"] = list(entrypoint)
        requirements = raw.get("requirements")
        if requirements is not None:
            normalized["requirements"] = _safe_relative_path(
                requirements, "runtime.requirements"
            )
        normalized["system_site_packages"] = bool(
            raw.get("system_site_packages", False)
        )
    else:
        image = str(raw.get("image", "")).strip()
        dockerfile = raw.get("dockerfile")
        if not image and dockerfile is None:
            raise ValueError("container runtime requires image and/or dockerfile")
        if image:
            normalized["image"] = image
        if dockerfile is not None:
            normalized["dockerfile"] = _safe_relative_path(
                dockerfile, "runtime.dockerfile"
            )
        build_contexts = raw.get("build_contexts", {})
        if not isinstance(build_contexts, Mapping):
            raise ValueError("runtime.build_contexts must be a mapping")
        normalized_contexts: dict[str, str] = {}
        for raw_name, raw_path in build_contexts.items():
            name = str(raw_name)
            if not BUILD_CONTEXT_NAME_PATTERN.fullmatch(name):
                raise ValueError(f"Invalid Docker build context name: {name!r}")
            if not isinstance(raw_path, str) or not raw_path.strip():
                raise ValueError(
                    f"runtime.build_contexts.{name} must be a non-empty path"
                )
            normalized_contexts[name] = raw_path.strip()
        normalized["build_contexts"] = normalized_contexts
        docker_args = raw.get("docker_args", [])
        if not isinstance(docker_args, list) or not all(
            isinstance(item, str) for item in docker_args
        ):
            raise ValueError("runtime.docker_args must be a list of strings")
        normalized["docker_args"] = list(docker_args)
    return normalized


def _validate_scope_block(raw: Mapping[str, Any]) -> dict[str, Any]:
    unknown = sorted(
        set(raw)
        - {
            "configurable",
            "default",
            "allowed_parts",
            "allow_slices",
            "allow_window",
        }
    )
    if unknown:
        raise ValueError("Unknown scope fields: " + ", ".join(unknown))
    default = normalize_scorer_scope(raw.get("default"), "scope.default")
    allowed_raw = raw.get("allowed_parts", sorted(PART_NAMES))
    if (
        not isinstance(allowed_raw, list)
        or not allowed_raw
        or not all(isinstance(part, str) for part in allowed_raw)
    ):
        raise ValueError("scope.allowed_parts must be a non-empty list")
    allowed_parts = tuple(dict.fromkeys(allowed_raw))
    unknown_parts = sorted(set(allowed_parts) - PART_NAMES)
    if unknown_parts:
        raise ValueError(
            "Unknown scope.allowed_parts: " + ", ".join(unknown_parts)
        )
    normalized = {
        "configurable": bool(raw.get("configurable", False)),
        "default": default,
        "allowed_parts": list(allowed_parts),
        "allow_slices": bool(raw.get("allow_slices", True)),
        "allow_window": bool(raw.get("allow_window", False)),
    }
    if normalized["allow_window"] and set(allowed_parts) != PART_NAMES:
        raise ValueError(
            "scope.allow_window requires every sequence part in allowed_parts"
        )
    _validate_selected_scope(normalized, default, scorer_id="manifest default")
    return normalized


def _validate_selected_scope(
    scope_policy: Mapping[str, Any],
    selected_scope: Mapping[str, Any],
    *,
    scorer_id: str,
) -> None:
    config = parse_scorer_scope(selected_scope, f"scope for {scorer_id}")
    allowed_parts = set(scope_policy.get("allowed_parts", PART_NAMES))
    unsupported_parts = sorted(
        set(scope_parts(selected_scope, f"scope for {scorer_id}"))
        - allowed_parts
    )
    if unsupported_parts:
        raise ValueError(
            f"Scope for scorer {scorer_id} uses unsupported parts: "
            + ", ".join(unsupported_parts)
        )
    if (
        config.scope == "parts"
        and not scope_policy.get("allow_slices", True)
        and any(span.start is not None or span.end is not None for span in config.spans)
    ):
        raise ValueError(f"Scope slices are disabled for scorer {scorer_id}")
    if config.scope == "window" and not scope_policy.get("allow_window", False):
        raise ValueError(f"Window scope is disabled for scorer {scorer_id}")


def _validate_outputs(raw: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    if not raw:
        raise ValueError("Scorer manifest outputs must not be empty")
    outputs: dict[str, dict[str, Any]] = {}
    for name, raw_spec in raw.items():
        output_name = str(name)
        if not OUTPUT_NAME_PATTERN.fullmatch(output_name):
            raise ValueError(f"Invalid scorer output name: {output_name!r}")
        spec = _mapping(raw_spec, f"outputs.{output_name}")
        output_type = str(spec.get("type", "")).lower()
        if output_type not in SUPPORTED_OUTPUT_TYPES:
            raise ValueError(
                f"outputs.{output_name}.type must be one of: "
                + ", ".join(sorted(SUPPORTED_OUTPUT_TYPES))
            )
        description = str(spec.get("description", "")).strip()
        if not description:
            raise ValueError(f"outputs.{output_name}.description is required")
        default_policy = _validate_output_policy(
            spec.get("default_policy", {"role": "log"}),
            output_type=output_type,
            context=f"outputs.{output_name}.default_policy",
        )
        normalized = {
            "type": output_type,
            "description": description,
            "default_policy": default_policy,
        }
        if spec.get("unit") is not None:
            normalized["unit"] = str(spec["unit"])
        if spec.get("value_range") is not None:
            value_range = str(spec["value_range"]).strip()
            if not value_range:
                raise ValueError(
                    f"outputs.{output_name}.value_range must not be empty"
                )
            normalized["value_range"] = value_range
        outputs[output_name] = normalized
    return outputs


def _validate_output_policy(
    raw: Any,
    *,
    output_type: str,
    context: str,
) -> dict[str, Any]:
    spec = _mapping(raw, context)
    role = str(spec.get("role", "log")).lower()
    if role not in {"objective", "threshold", "log"}:
        raise ValueError(f"{context}.role must be objective, threshold, or log")
    if output_type == "structure" and role != "log":
        raise ValueError(f"{context}: structure outputs can only be logged")
    if output_type in {"bool", "string"} and role == "objective":
        raise ValueError(f"{context}: {output_type} outputs cannot be objectives")

    normalized: dict[str, Any] = {"type": output_type, "role": role}
    if role == "objective":
        direction = str(spec.get("direction", "")).lower()
        if direction not in {"min", "max"}:
            raise ValueError(f"{context}.direction must be min or max")
        normalized["direction"] = direction
    elif role == "threshold":
        supplied = [
            name for name in ("min", "max", "equals") if spec.get(name) is not None
        ]
        if not supplied:
            raise ValueError(f"{context} threshold requires min, max, and/or equals")
        if output_type in {"bool", "string"} and supplied != ["equals"]:
            raise ValueError(
                f"{context}: {output_type} thresholds support only equals"
            )
        for name in supplied:
            value = spec[name]
            if name in {"min", "max"}:
                value = float(value)
            normalized[name] = value
    return normalized


def _validate_parameters(raw: Any) -> dict[str, dict[str, Any]]:
    parameters_raw = _mapping(raw, "parameters")
    parameters: dict[str, dict[str, Any]] = {}
    for name, raw_spec in parameters_raw.items():
        spec = _mapping(raw_spec, f"parameters.{name}")
        unknown_fields = sorted(
            set(spec)
            - {
                "type",
                "description",
                "required",
                "default",
                "min",
                "max",
                "choices",
                "env",
                "allowed_root",
                "mount",
                "configurable_by",
                "backend_env",
            }
        )
        if unknown_fields:
            raise ValueError(
                f"Unknown fields for parameters.{name}: "
                + ", ".join(unknown_fields)
            )
        parameter_type = str(spec.get("type", "")).lower()
        if parameter_type not in SUPPORTED_PARAMETER_TYPES:
            raise ValueError(
                f"parameters.{name}.type must be one of: "
                + ", ".join(sorted(SUPPORTED_PARAMETER_TYPES))
            )
        description = str(spec.get("description", "")).strip()
        if not description:
            raise ValueError(f"parameters.{name}.description is required")
        normalized: dict[str, Any] = {
            "type": parameter_type,
            "description": description,
            "required": bool(spec.get("required", False)),
        }
        configurable_by = str(
            spec.get("configurable_by", "experiment")
        ).strip().lower()
        if configurable_by not in SUPPORTED_PARAMETER_OWNERS:
            raise ValueError(
                f"parameters.{name}.configurable_by must be one of: "
                + ", ".join(sorted(SUPPORTED_PARAMETER_OWNERS))
            )
        normalized["configurable_by"] = configurable_by
        backend_env = spec.get("backend_env")
        if configurable_by == "backend":
            if backend_env is None or not ENVIRONMENT_NAME_PATTERN.fullmatch(
                str(backend_env)
            ):
                raise ValueError(
                    f"parameters.{name}.backend_env must be a valid environment "
                    "variable name for backend-owned parameters"
                )
            normalized["backend_env"] = str(backend_env)
        elif backend_env is not None:
            raise ValueError(
                f"parameters.{name}.backend_env is valid only for "
                "backend-owned parameters"
            )
        for key in ("default", "min", "max", "choices", "env"):
            if key in spec:
                normalized[key] = copy.deepcopy(spec[key])
        allowed_root = spec.get("allowed_root")
        if allowed_root is not None:
            if parameter_type != "path":
                raise ValueError(
                    f"parameters.{name}.allowed_root is valid only for path parameters"
                )
            normalized_root = Path(_expand_asset_path(str(allowed_root)))
            if not normalized_root.is_absolute():
                raise ValueError(
                    f"parameters.{name}.allowed_root must be an absolute path"
                )
            normalized["allowed_root"] = str(allowed_root)
        mount = spec.get("mount")
        if mount is not None:
            mount_spec = _mapping(mount, f"parameters.{name}.mount")
            container_path = str(mount_spec.get("container_path", "")).strip()
            if not container_path.startswith("/"):
                raise ValueError(
                    f"parameters.{name}.mount.container_path must be absolute"
                )
            normalized["mount"] = {
                "container_path": container_path,
                "read_only": bool(mount_spec.get("read_only", True)),
            }
        parameters[str(name)] = normalized

    # Validate defaults through the same policy used for notebook values.
    _resolve_parameter_values(parameters, {})
    return parameters


def _resolve_parameter_values(
    specs: Mapping[str, Any],
    raw_values: Any,
) -> dict[str, Any]:
    if not isinstance(raw_values, Mapping):
        raise ValueError("scorer selection parameters must be a mapping")
    unknown = sorted(set(raw_values) - set(specs))
    if unknown:
        raise ValueError("Unknown scorer parameters: " + ", ".join(unknown))
    values: dict[str, Any] = {}
    for name, raw_spec in specs.items():
        spec = _mapping(raw_spec, f"parameters.{name}")
        if name in raw_values:
            raw_value = raw_values[name]
        elif "default" in spec:
            raw_value = spec["default"]
        elif spec.get("required"):
            raise ValueError(f"Required scorer parameter is missing: {name}")
        else:
            continue
        parameter_type = spec["type"]
        if parameter_type == "path":
            raw_value = _expand_asset_path(str(raw_value))
        if parameter_type == "int":
            if isinstance(raw_value, bool):
                raise ValueError(f"Parameter {name} must be an integer")
            value: Any = int(raw_value)
        elif parameter_type == "float":
            if isinstance(raw_value, bool):
                raise ValueError(f"Parameter {name} must be numeric")
            value = float(raw_value)
        elif parameter_type == "bool":
            if isinstance(raw_value, bool):
                value = raw_value
            elif str(raw_value).strip().lower() in {"true", "1", "yes"}:
                value = True
            elif str(raw_value).strip().lower() in {"false", "0", "no"}:
                value = False
            else:
                raise ValueError(f"Parameter {name} must be boolean")
        else:
            value = str(raw_value)
        if spec.get("min") is not None and value < spec["min"]:
            raise ValueError(f"Parameter {name} is below its minimum")
        if spec.get("max") is not None and value > spec["max"]:
            raise ValueError(f"Parameter {name} is above its maximum")
        if spec.get("choices") is not None and value not in spec["choices"]:
            raise ValueError(f"Parameter {name} is not one of its allowed choices")
        if parameter_type == "path" and not str(value).strip():
            raise ValueError(f"Parameter {name} path must not be empty")
        values[str(name)] = value
    return values


def _validate_self_test(
    raw: Mapping[str, Any],
    outputs: Mapping[str, Any],
    scope: Mapping[str, Any],
) -> dict[str, Any]:
    unknown_fields = sorted(
        set(raw) - {"sequence", "sequence_context", "expected_outputs"}
    )
    if unknown_fields:
        raise ValueError(
            "Unknown scorer self_test fields: " + ", ".join(unknown_fields)
        )
    sequence = str(raw.get("sequence", "")).strip()
    if not sequence:
        raise ValueError("self_test.sequence is required")
    sequence_context = _validate_sequence_context(
        _mapping(raw.get("sequence_context"), "self_test.sequence_context"),
        sequence_length=len(sequence),
    )
    disallowed_parts = sorted(
        {
            str(item["part"])
            for item in sequence_context["part_ranges"]
        }
        - set(scope["allowed_parts"])
    )
    if disallowed_parts:
        raise ValueError(
            "self_test.sequence_context uses parts outside scope.allowed_parts: "
            + ", ".join(disallowed_parts)
        )
    expected = raw.get("expected_outputs")
    if not isinstance(expected, list) or not expected:
        raise ValueError("self_test.expected_outputs must be a non-empty list")
    unknown = sorted(set(map(str, expected)) - set(outputs))
    if unknown:
        raise ValueError(
            "self_test references unknown outputs: " + ", ".join(unknown)
        )
    return {
        "sequence": sequence,
        "sequence_context": sequence_context,
        "expected_outputs": [str(name) for name in expected],
    }


def _validate_sequence_context(
    raw: Mapping[str, Any], *, sequence_length: int
) -> dict[str, Any]:
    unknown = sorted(set(raw) - {"part_ranges"})
    if unknown:
        raise ValueError(
            "Unknown self_test.sequence_context fields: " + ", ".join(unknown)
        )
    ranges = raw.get("part_ranges")
    if not isinstance(ranges, list) or not ranges:
        raise ValueError(
            "self_test.sequence_context.part_ranges must be a non-empty list"
        )
    normalized: list[dict[str, Any]] = []
    previous_end = 0
    required_fields = {"part", "start", "end", "part_start", "part_end"}
    for index, item in enumerate(ranges):
        context = f"self_test.sequence_context.part_ranges[{index}]"
        item = _mapping(item, context)
        item_unknown = sorted(set(item) - required_fields)
        missing = sorted(required_fields - set(item))
        if item_unknown or missing:
            details = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if item_unknown:
                details.append("unknown " + ", ".join(item_unknown))
            raise ValueError(f"{context} is invalid: " + "; ".join(details))
        part = str(item["part"])
        if part not in PART_NAMES:
            raise ValueError(f"{context}.part is unknown: {part}")
        try:
            start = int(item["start"])
            end = int(item["end"])
            part_start = int(item["part_start"])
            part_end = int(item["part_end"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{context} coordinates must be integers") from exc
        if (
            start != previous_end
            or not 0 <= start < end <= sequence_length
            or not 0 <= part_start < part_end
            or end - start != part_end - part_start
        ):
            raise ValueError(f"{context} has invalid or inconsistent coordinates")
        normalized.append(
            {
                "part": part,
                "start": start,
                "end": end,
                "part_start": part_start,
                "part_end": part_end,
            }
        )
        previous_end = end
    if previous_end != sequence_length:
        raise ValueError(
            "self_test.sequence_context.part_ranges must cover the routed sequence"
        )
    return {"part_ranges": normalized}


def _safe_relative_path(value: Any, context: str) -> str:
    path = Path(str(value))
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"{context} must be a safe relative path")
    return path.as_posix()




def _mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} must be a mapping")
    return dict(value)


def _expand_asset_path(value: str) -> str:
    """Expand documented path variables without permitting shell evaluation."""
    def replace(match):
        name = match.group(1)
        if name in os.environ:
            return os.environ[name]
        if name == "MORSE_ASSETS_DIR":
            return str(Path.home() / ".cache" / "morse" / "assets")
        raise ValueError(f"Environment variable {name} is required for a scorer path")
    return str(Path(re.sub(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}", replace, value)).expanduser())
