# Python API

Use the environment installed from the repository's [README](../README.md#install).
For notebook output, include the notebook extra:

```sh
uv sync --all-packages --extra paper --extra mfe --extra notebook --frozen
```

Select this environment as your notebook kernel. To register it after activating
the virtual environment:

```sh
python -m ipykernel install --user --name morse --display-name MORSE
```

The
[CPU example notebook](../core/notebooks/MORSE_Python_API.ipynb) runs a small
experiment with the bundled sequence-metric scorer.
The [graphical notebook](../core/notebooks/MORSE_Server_Colab.ipynb) provides
forms for creating YAML and calls this API in its active Python kernel.

## Configure and run

```python
import morse

params = morse.load_config("experiment.yaml")
result = morse.optimize(
    params=params,
    verbose=2,
    plots=["metrics", "hypervolume", "utrs"],
)
```

`params` is an ordinary dictionary using the experiment YAML schema. Schemas 3
and 4 remain supported. `load_config()` resolves paths relative to the YAML.
You can edit the returned dictionary before running:

```python
params["optimizer"]["population_size"] = 100
params["optimizer"]["generations"] = 200
```

To supply RNA directly:

```python
result = morse.optimize(
    rna_seq,
    params,
    output_dir="results/my_run",
    verbose=2,
)
```

The RNA replaces `construct.fasta`; coordinates and allowed changes stay in the
configuration. Coordinates are 0-based, end exclusive. A CDS-only input uses
`coordinates: {cds: [0, length]}`. Input can contain RNA U or DNA T.
MORSE does not modify the supplied dictionary.

A dictionary read with `yaml.safe_load()` also works. Set `base_dir` when its
relative input paths belong to a different folder:

```python
result = morse.optimize(rna_seq, params, base_dir="experiments/my_run")
```

Scorer commands and their `parameters` work as in the CLI. Prepare the scorer
environment before starting the run; see [scorer setup](scorers.md).

## Check and score

```python
check = morse.validate(rna_seq, params)
scores = morse.score(rna_seq, params)
scores["outputs"]
```

Validation checks configuration and construct structure. Scoring runs the
configured scorers once and returns outputs, threshold checks and `feasible`.
Omit the RNA argument to read the configured FASTA.

## Progress

| `verbose` | Output |
|---|---|
| `0` | No progress output |
| `1` | Generation, elapsed time, feasibility and optimization progress |
| `2` | Progress and plots updated in one notebook display |

`plots` chooses `metrics`, `hypervolume` and `utrs`. Metric plots show population
mean and standard deviation; UTR plots show variant shares in each population.
The settings affect display only. Results and history are saved in `output_dir`
from the configuration, or the explicit argument. If neither is set, MORSE
creates a new directory under `outputs/`.

## Inspect results

```python
result.summary
result.pareto
result.history
result.output_dir

result.plot("metrics")
result.plot("hypervolume")
result.plot("utrs")
result.report()
```

`pareto` and `history` are lists of dictionaries. For a DataFrame, use
`pandas.DataFrame(result.pareto)`. `plot()` returns a Matplotlib figure;
`report()` displays the saved report in a notebook.

Choose a candidate explicitly:

```python
candidate = result.select("hvc")
candidate.id
candidate.sequence
candidate.parts
candidate.metrics
```

`hvc` selects the largest exclusive hypervolume contribution in the final
feasible Pareto set. For selection by a single metric, use its full name:

```python
candidate = result.select(metric="sequence_metrics.cpg_per_kb")
```

The configured direction determines whether smaller or larger values win.

## Reopen or continue

```python
result = morse.load_run("results/my_run")
result.plot("metrics")
```

Interrupting the optimization cell stops the run. Continue from its latest
saved generation in the same output directory:

```python
result = morse.resume("results/my_run", verbose=2)
```

The saved input, resolved configuration and checkpoint are used for continuation.
