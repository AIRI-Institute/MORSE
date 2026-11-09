# Colab and the optional jobs backend

Colab is a graphical interface for the same schema-4 YAML used by the CLI. It
edits the construct, optimizer, scorer commands and metric roles; exports YAML;
launches and follows jobs; and displays their plots and artifacts. Local mode
uses the same job manager as the optional remote backend.

Scorer scripts, containers and services must already be prepared. Neither the
notebook nor the backend clones scorer repositories, installs packages, builds
images, starts model services or downloads checkpoints.

## Start a backend

Install the monorepo environment as described in the root README. From its root:

```sh
export MORSE_BACKEND_TOKEN=YOUR_RANDOM_TOKEN
uv run morse-backend \
  --config-root core/configs \
  --results-root runs/backend \
  --workdir . \
  --host 127.0.0.1 --port 8000
```

Use a random token and keep the real value outside Git. The bearer token grants
trusted execution: an authenticated user can submit Bash scorer commands that
run as the backend service account. Share it only with users trusted to execute
code on that machine. Static `config` selections are restricted to YAML files
under `--config-root`; full submitted configurations can reference files already
available to that trusted account.

The health endpoint is public; all job/configuration endpoints require
`Authorization: Bearer TOKEN`. Keep the backend on loopback and use an SSH tunnel
or an authenticated HTTPS reverse proxy for remote access.

The service environment supplies variables referenced by scorer commands. Core
provides `MORSE_PYTHON` as its current Python interpreter unless overridden.
Commands run once per batch with one JSON request on stdin followed by EOF and
must return one JSON response. External model-service identity and state are not
verified automatically; declared `provenance_paths` record local dependencies.

## Submit the complete experiment

The notebook uses the same body for preflight and submission:

```python
from morse.notebook_tools import MorseBackendClient, read_configuration

client = MorseBackendClient(base_url, token)
configuration = read_configuration("experiment.yaml")
client.preflight_job(configuration=configuration)
client.submit_job(configuration=configuration)
```

`configuration` is the complete schema-4 document, including `scorers` with Bash
`command` strings, `scope` and `metrics`. There is no separate scorer selection,
installation ID, parameter profile or configuration compiler. Saving the form
writes the RNA FASTA and any UTR libraries beside the YAML. The download ZIP
contains the YAML and these inputs together.

The client sends declared sequence, library and initialization text inputs as
`input_files`, up to 10 MiB in total. The backend stages them in an immutable,
content-addressed input directory. It does not upload scorer scripts, model
weights, credentials or arbitrary provenance files. Input filenames must be
relative, without traversal, and must be referenced by the submitted config.
Prepared server-side paths remain usable when no input attachment is supplied.

Default command working directories and remaining relative server-side inputs
resolve against the backend's `config_root`; use an explicit scorer `workdir` when needed.
For a trusted server-side YAML, submit `{"config": "experiment.yaml"}` instead.
The two forms are mutually exclusive.

The backend validates and stores an immutable resolved YAML snapshot for each
job. A selected file keeps the path meaning of its original directory. Editing
the source YAML after submission does not change a queued job. Preflight returns
counts, scorer commands/metrics, a scope SVG and a configuration fingerprint; it
does not invoke scorer commands or enqueue a run.

Managed jobs store results under the backend's `results_root`, overriding the
YAML `output_dir`. Direct CLI runs use `output_dir` unless overridden explicitly.

## API version 6

| Method | Endpoint | Result |
| --- | --- | --- |
| GET | `/health` | Runtime identity and queue counts; no token required |
| GET | `/v1/configs` | Trusted server-side YAML names |
| POST | `/v1/jobs/preflight` | Validate a complete job payload without running it |
| POST | `/v1/jobs` | Queue an optimization |
| POST | `/v1/jobs/pilot` | Queue a pilot; optional `pilot_samples` |
| GET | `/v1/jobs` | List jobs |
| GET | `/v1/jobs/{id}` | Status, timing and configuration identity |
| GET | `/v1/jobs/{id}/metrics` | Progress and generation history |
| GET | `/v1/jobs/{id}/log?tail=200` | Log tail |
| GET | `/v1/jobs/{id}/sequence-scope` | Scope SVG |
| GET | `/v1/jobs/{id}/pilot` | Completed pilot data |
| GET | `/v1/jobs/{id}/artifacts` | ZIP after a terminal status |
| DELETE | `/v1/jobs/{id}` | Cancel a queued or running job |
| POST | `/v1/jobs/{id}/resume` | Resume a checkpoint; body `{}` |

Jobs retain the existing queue, cancellation, checkpoint/resume and artifact
behavior. Backend restart marks unfinished jobs `interrupted`; resume requires
an existing optimizer checkpoint and the usual Core compatibility checks.
Cancellation stops the managed job and its local scorer clients. The backend
does not stop an external model service or container.

Schema 3 remains readable for saved experiments. New forms generate schema 4:
`cds` is the complete CDS, `cds_body` omits its start and stop, and output names
are qualified as `scorer.metric`.

The former scorer-installation endpoints are absent in API 6. Use a matching
notebook version; its compatibility check reads `/health` before launching work.
