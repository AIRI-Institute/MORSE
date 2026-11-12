"""Optional notebook display; never stored in an optimizer checkpoint."""
from __future__ import annotations

import html
from pathlib import Path
import time


class ProgressDisplay:
    def __init__(self, output_dir, *, verbose=1, plots=None):
        if type(verbose) is not int or verbose not in {0, 1, 2}:
            raise ValueError("verbose must be 0, 1 or 2")
        if isinstance(plots, str):
            plots = [plots]
        self.plots = tuple(["metrics", "hypervolume"] if plots is None else plots)
        unknown = set(self.plots) - {"metrics", "hypervolume", "utrs"}
        if unknown:
            raise ValueError("Unknown plots: " + ", ".join(sorted(unknown)))
        self.output_dir = Path(output_dir)
        self.verbose = verbose
        self.handle = None
        self.last_payload = None
        self.last_draw = 0.0
        self.shell = None
        if verbose:
            try:
                from IPython import get_ipython
                shell = get_ipython()
                if shell is not None and getattr(shell, "kernel", None) is not None:
                    self.shell = shell
            except ImportError:
                pass

    def _show(self, text):
        if self.shell is None:
            print(text, flush=True)
        else:
            from IPython.display import HTML, display
            block = HTML("<pre>" + html.escape(text) + "</pre>")
            if self.handle is None:
                self.handle = display(block, display_id=True)
            else:
                self.handle.update(block)

    def start(self):
        if self.verbose:
            self._show(f"MORSE · preparing initial scoring\nResults: {self.output_dir}")

    @staticmethod
    def _status(payload):
        total = payload.get("total_generations") or "?"
        n = payload.get("population_size", 0)
        feasible = payload.get("hard_feasible_count")
        if feasible is None:
            feasible = round(payload.get("hard_feasible_fraction", 0) * n)
        elapsed = float(payload.get("elapsed_seconds", 0))
        text = (f"Generation {payload['generation']}/{total} · "
                f"{payload.get('evaluated_candidates', 0)} evaluated · "
                f"{elapsed / 60:.1f} min · feasible {feasible}/{n}")
        hv = payload.get("hypervolume")
        if hv is not None:
            text += f"\nHV {hv:.6g} · best {payload.get('best_hypervolume', hv):.6g}"
        elif payload.get("objectives"):
            item = payload["objectives"][0]
            text += f"\n{item['name']}: best {item.get('best')}"
        early = payload.get("early_stopping", {})
        if early.get("enabled") and early.get("stale_generations") is not None:
            text += f" · early stopping {early['stale_generations']}/{early.get('patience', '?')}"
        return text

    def update(self, payload, *, force=False):
        self.last_payload = payload
        if not self.verbose:
            return
        text = self._status(payload)
        if self.verbose < 2 or self.shell is None:
            self._show(text)
            return
        if not force and time.monotonic() - self.last_draw < 1:
            return
        from IPython.display import HTML, display
        import matplotlib.pyplot as plt
        from morse.results import OptimizationResult

        result = OptimizationResult(self.output_dir)
        figures = []
        try:
            for name in self.plots:
                try:
                    figures.append(result.plot(name))
                except ValueError:
                    # A plot can be unavailable before the first feasible population.
                    continue
            # One display target retains earlier cells and any user output.
            import base64
            import io
            blocks = ["<pre>" + html.escape(text) + "</pre>"]
            for figure in figures:
                buffer = io.BytesIO()
                figure.savefig(buffer, format="png", dpi=110, bbox_inches="tight")
                encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
                blocks.append(f'<img style="max-width:100%" src="data:image/png;base64,{encoded}">')
            block = HTML("\n".join(blocks))
            if self.handle is None:
                self.handle = display(block, display_id=True)
            else:
                self.handle.update(block)
        finally:
            for figure in figures:
                plt.close(figure)
        self.last_draw = time.monotonic()

    def finish(self, message, *, success=True):
        if not self.verbose:
            return
        if success and self.last_payload is not None and self.verbose == 2 and self.shell is not None:
            self.update(self.last_payload, force=True)
            # Keep the final plots in the display; completion is a separate short line.
            from IPython.display import HTML, display
            display(HTML(f"<p>{html.escape(message)} · {html.escape(str(self.output_dir))}</p>"))
        else:
            self._show(f"{message}\nResults: {self.output_dir}")
