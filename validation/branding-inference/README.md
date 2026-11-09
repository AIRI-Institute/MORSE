# Inference parity after the MORSE namespace rename

Real A100 checks on 2026-11-05 used the renamed source from isolated snapshot
`09809d98f3f641c53897a8991fdf2dce7949d9ba`. Its tracked files remained unchanged
throughout validation and match the current source checkout byte for byte.

| Adapter | Frozen constructs | Maximum absolute difference | Comparison |
|---|---:|---:|---|
| Full-fit Orthrus TE | 3 | 0 | Native image adapter vs renamed source overlay |
| Orthrus half-life | 3 | 0 | Native image adapter vs renamed source overlay |
| RiboNN | 3 | 0 | Native image adapter vs renamed source overlay |
| RiNALMo | 3 | 0 | Native image adapter vs renamed source overlay |
| Held-out Orthrus TE | 3 | 0 | Renamed adapter vs frozen paper observations |

The first four checks use the same immutable Docker image, ready assets and
input batches on both sides. Held-out TE uses its own verified checkpoint and
adapter in the compatible existing runtime, with a reference-alone request
followed by the candidates. Fixed numerical tolerances were retained; all
observed differences are zero.

Images were reused, not rebuilt. These checks validate renamed inference code,
including the Orthrus head classes and runtime environment names; they do not
repeat complete searches or establish parity for newly built release images.
Only the new validation containers were created and removed. Existing server
checkouts, environments, images and foreign jobs were not changed.

Each receipt records image identity, fixture identity, outputs, actual source
hashes and verified model/backbone/ensemble digests. `source-snapshot.json`
describes the private inference-source snapshot rather than a full monorepo
commit. Full technical logs remain under ignored
`validation/private/branding-inference/`.
