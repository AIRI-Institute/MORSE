"""Publication colours assigned by scorer order, independent of scorer names.

The base colours are copied from the final paper palette. Keeping these plain
Python constants in Core avoids a dependency on the optional paper bundle.
Every output of one scorer uses its owner's colour. New scorers automatically
receive the next colour; no model registration or semantic alias is required.
"""
from __future__ import annotations

from collections.abc import Mapping

INK = "#17212F"
MUTED = "#7B818B"
GRID = "#DEE4EC"
INITIAL = "#9AA4B2"
FIXED_FILL = "#EDF1F5"
MUTABLE_FILL = "#EDE9FE"
MUTABLE_EDGE = "#8B5CF6"

COLORS = (
    "#9870B1", "#337F82", "#659773", "#4395B5", "#507FB1",
    "#8D9355", "#B98145", "#C26D95", "#C76C65", "#85929F", "#7F858D",
)


def palette_color(index: int) -> str:
    """Return the colour for a scorer's zero-based configuration position.

    Beyond the base palette, alternate light and dark variants of each colour
    deterministically. Dash variants can additionally distinguish long legends.
    """
    if index < 0:
        raise ValueError("Scorer colour index must be nonnegative")
    base = COLORS[index % len(COLORS)]
    cycle = index // len(COLORS)
    if not cycle:
        return base
    fraction = min(0.16 * ((cycle + 1) // 2), 0.60)
    target = 255 if cycle % 2 else 0
    channels = [int(base[i:i+2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(value + (target-value)*fraction):02X}" for value in channels)


def palette_linestyle(index: int) -> str:
    return ("-", "--", "-.", ":")[(index // len(COLORS)) % 4]


def run_color_mapping(scorers) -> dict[str, dict[str, str]]:
    """Map scorer names and output names from the run's ordered scorer list.

    Accepts configuration dictionaries or objects with name/outputs attributes.
    Separate namespaces avoid collisions between a binding name and an output.
    """
    names, outputs = {}, {}
    for index, scorer in enumerate(scorers):
        if isinstance(scorer, Mapping):
            name, declared = scorer.get("name", f"scorer_{index}"), scorer.get("metrics", scorer.get("outputs", {}))
            prefix = scorer.get("output_prefix")
            namespace = scorer.get("namespace")
        else:
            name, declared = scorer.name, scorer.outputs
            prefix = getattr(scorer, "output_prefix", None)
            namespace = getattr(scorer, "namespace", None)
        color = palette_color(index)
        names[str(name)] = color
        for output in declared:
            key = output if isinstance(output, str) else getattr(output, "name", getattr(output, "output", ""))
            if key:
                output_name = f"{namespace}.{key}" if namespace else (f"{prefix}_{key}" if prefix else str(key))
                outputs[output_name] = color
    return {"scorers": names, "outputs": outputs}


RC_STYLE = {
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.labelsize": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.titlepad": 10,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.facecolor": "white", "figure.facecolor": "white",
    "savefig.facecolor": "white", "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "grid.color": GRID, "grid.linewidth": 0.5,
    "legend.frameon": False, "legend.fontsize": 8,
    "lines.linewidth": 1.4, "pdf.fonttype": 42, "ps.fonttype": 42,
}
