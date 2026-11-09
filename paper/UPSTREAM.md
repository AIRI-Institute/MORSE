# Native upstream design regeneration

These optional commands invoke separately installed authors' software. The
monorepository includes only launch and measurement adapters, not copied upstream
source or predictor-training code. Upstream software and ready weights retain
their own terms; consult `licenses/` and the upstream license before obtaining it.

## LinearDesign

Obtain the authors' LinearDesign repository separately at commit
`f0126ca89a8b853088b4bccfd2cc8c378d3678be`. The runner requires this exact Git
revision, invokes `make` if its executable is absent, and runs
`bin/LinearDesign_2D 3 0 codon_usage_freq_table_human.csv` with the bundled MLANA
protein. The output is checked for preserved protein sequence and saved with its
full stdout, stderr and receipt.

```bash
morse paper upstream lineardesign --source /path/to/LinearDesign \
  --output results/lineardesign-plan
morse paper upstream lineardesign --source /path/to/LinearDesign \
  --output results/lineardesign-native --execute
```

The public continuation recipes use the exact frozen lambda=3 starting CDS. A
new upstream design can be passed to an independent experiment by setting the CDS
source and generating a matched initial population. Do not label a continuation
from a changed starting CDS as an exact historical replay.

## RiboDecode

The frozen native generation recipe uses the ready RiboDecode native optimization
environment recorded in `inputs/native_rd.json`, HEK293T, 10 optimization epochs,
zero MFE weight, alpha=beta=100 and 40,000 final-epoch samples. It seeds Python,
NumPy and PyTorch in the process executing the authors' native sequence optimizer,
and enables deterministic algorithms. The predictor weights are not trained here.

```bash
morse paper upstream ribodecode --seed 7 --device cuda:0 \
  --output results/ribodecode-plan
morse paper upstream ribodecode --seed 7 --device cuda:0 \
  --output results/ribodecode-native --execute
```

Pass `--image` for a separately installed compatible native generation image.
The exact upstream optimizer-source SHA256 is checked before generation. The
wrapper imports the external authors' native optimization entry module; that
module is not included in this repository. Generation requires CUDA because the
frozen upstream code uses it directly. `--pilot` performs one epoch with 400
samples and explicitly labels the receipt diagnostic.

After native generation, selection uses the highest native score among distinct
protein-preserving CDSs before feasibility screening. Lexicographic CDS order
breaks score ties; the stop is standardized to UAA. The selected record includes
the complete native-output hash. A new native start should only replace its
corresponding continuation input when its sequence and matched 64-member local
population have been explicitly validated.

Repeat the command for seeds 7, 43, 101, 137 and 211. Native generation costs are
separate from the MORSE continuation ceiling and early-stopping costs.
