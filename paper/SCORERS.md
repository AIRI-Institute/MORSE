# Paper scorer endpoints and frozen presets

The experiment YAMLs contain plain batch commands. Their example localhost ports
are listed below; replace the URLs when your prepared scorers use other addresses.
No model, process or container is started by `morse run`.

Schema-4 scopes use `cds` for the complete coding region and `cds_body` for the
body without start/stop. Paper scorers retain their exact historical inputs:
`[start_codon, cds_body]` excludes the terminal stop, while
`[five_prime_utr, cds, three_prime_utr]` excludes poly(A). New run outputs are named
`scorer.metric`; paper comparisons map these names to the frozen columns.

The full historical image IDs and environment values are archived as static
provenance in `provenance/scorer-environments.json`. That file is not a runtime
registry. Ready checkpoint identities are in `assets/registry.json` at the repo root.
Match these weights, input preprocessing, CPU/GPU placement, thread limits and
inner model batch sizes when comparing a fresh run with the paper.

| Port | Preset | Scorer | Device | Threads | Inner batch |
|---:|---|---|---|---:|---|
| 8701 | `f2-viennarna-transcript-body` | viennarna | cpu | 1 | — |
| 8702 | `f2-csi-human` | csi | cpu | 1 | — |
| 8703 | `f2-orthrus-translation-efficiency` | orthrus-te-full-fit | cuda | 1 | 8 |
| 8704 | `f2-ribodecode-hek293t` | ribodecode | cuda | 1 | 8 |
| 8705 | `f2-sequence-feasibility` | manufacturability | cpu | 1 | — |
| 8706 | `f3a-viennarna-transcript-body` | viennarna | cpu | 1 | — |
| 8707 | `f3a-csi-human` | csi | cpu | 1 | — |
| 8708 | `f3b-viennarna-transcript-body` | viennarna | cpu | 1 | — |
| 8709 | `f3b-viennarna-isolated-cds` | viennarna | cpu | 1 | — |
| 8710 | `f3b-cai-lineardesign-frequencies` | cai | cpu | 1 | — |
| 8711 | `f3b-orthrus-translation-efficiency` | orthrus-te-full-fit | cuda | 1 | 2 |
| 8712 | `f3b-sequence-feasibility` | manufacturability | cpu | 1 | — |
| 8713 | `f3c-ribodecode-hek293t` | ribodecode | cpu | 1 | 64 |
| 8714 | `f3c-sequence-feasibility` | manufacturability | cpu | 1 | — |
| 8715 | `f4-ribonn` | ribonn | cuda | 2 | 16 |
| 8716 | `f4-orthrus-te` | orthrus-te-held-out | cuda | 2 | 4 (adapter default) |
| 8717 | `f4-ribodecode-hek293t` | ribodecode | cpu | 2 | 256 |

The model batch size above is internal to its adapter. Core sends its reference
construct alone and then the other 63 initial candidates. An HTTP adapter must
forward each request unchanged, preserving sequence/context order and request
boundaries; it must not merge these requests. Later generations keep the same
optimizer evaluation and sampling rules.

If your ready scorer already provides HTTP JSON inference, call it directly. For
a persistent JSONL scorer, the optional utility holds one explicitly supplied
command open and forwards POST `/score` requests:

```bash
morse scorers serve --command 'docker exec -i EXISTING_CONTAINER python -u /app/server.py' --port 8703
```

Configure the frozen environment inside that container beforehand. The utility
does not build images, create containers, manage profiles or restart backends.
`GET /health` reports whether its backend process is still running. Ctrl-C stops
only the process group created by that utility. A supplied Bash command can also
call a local interpreter instead of Docker.

CAI/CSI/profile inputs and native upstream starts remain bundled. This runtime
change does not alter source sequences, model checkpoints, metrics, hard
constraints, candidate-selection rules or archived numeric observations.
