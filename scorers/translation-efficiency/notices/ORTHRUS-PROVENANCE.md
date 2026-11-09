# Orthrus upstream provenance

The bundled `quietflamingo/orthrus-large-6-track` backbone is frozen at revision
`c17ac0db4611c95248001d6fb4bcd95e3293cec4`. Its pinned model card declares MIT:
https://huggingface.co/quietflamingo/orthrus-large-6-track/blob/c17ac0db4611c95248001d6fb4bcd95e3293cec4/README.md

The exact upstream project MIT notice is retained as `ORTHRUS-LICENSE.md`:
https://github.com/bowang-lab/Orthrus/blob/26e894ee33c01e0fd31066bf056cf337c6f11700/LICENSE.md

Notice SHA-256: `253b638b6fe54f99093f54801d52e17a2f0fbc8c5c46a290d606ede9642cb5ac`.

Bundled backbone file SHA-256 values:

- `config.json`: `a8b513b3d73356e7ac8a59516ebbf51bc03f635e6329ad32350c289f84c8a288`
- `orthrus_hf.py`: `35fb82d090003c9bf91be1820d5cc0e22c1092d28d9805951d3f284d728cdbe4`
- `pytorch_model.bin`: `a188c974493679d9eb5d18b5ecc3148f3634cd785e20f8303d6f313f62f00c08`

The registered Python source adds postponed annotations to the upstream file;
its inference architecture and arithmetic are preserved. Each image's retained
Dockerfile records and checks its distinct ready task-checkpoint digest too.
