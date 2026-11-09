# Custom scorer

`scorer.py` is a dependency-free scorer script. MORSE invokes its Bash command
once per batch, writes one JSON request to stdin, then sends EOF. The script
returns a matching `request_id` and one output object per sequence, in order.
Send diagnostics to stderr; reserve stdout for the single JSON response.

```sh
printf '%s\n' '{"request_id":1,"sequences":["ATGGCTCGT"]}' | python examples/custom_scorer/scorer.py
```

Replace `score()` with your inference function. This example calculates CpG count
and distance from 50% GC as transparent demonstration objectives; neither is a
validated predictor of biological performance. Core sends sequences already
projected to the configured `scope` and optional `sequence_contexts` with part
coordinates. Keep input order and batch length unchanged. Return
`{"request_id":1,"error":"message"}` if inference fails. All model files, extra
source modules and auxiliary tables must be listed in the scorer's
`provenance_paths` so that changing them invalidates resume and cached scores.
See [the researcher guide](../../docs/researcher-guide.md) for configuration.

A schema-4 scorer section uses a command string and metrics:

```yaml
scorers:
  - name: sequence_metrics
    command: '"${MORSE_PYTHON}" scorer.py'
    provenance_paths: [scorer.py]
    scope: [start_codon, cds_body]
    metrics:
      cpg_per_kb: {role: optimize, direction: min}
      gc_distance: {role: optimize, direction: min}
      gc_fraction: {role: log}
```

Save it in the same experiment folder as the script. `MORSE_PYTHON` defaults to
MORSE's interpreter. Prepare dependencies and checkpoints before the run. A
prepared container or HTTP service can replace the local client command; keep
the same JSON contract and preserve `sequence_contexts` when forwarding it.
