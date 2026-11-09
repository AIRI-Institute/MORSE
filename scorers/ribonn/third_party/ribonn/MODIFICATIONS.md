# Inference-only extraction

Upstream revision: `dc896e1d1f78bf6b4e9bb7525c7852fcf25532e8`.

- Kept `ConvBlock` and `RiboNN` convolutional architecture, parameter names,
  state-dict loading and forward arithmetic.
- Replaced the Lightning base class with `torch.nn.Module`; copied constructor
  hyperparameters directly. No architecture or output normalization changed.
- Removed fitting/evaluation steps, losses, optimizers, metrics, data modules,
  cross-validation launchers, unused Saluki/MeanTE networks and stochastic
  augmentation. Nonzero `max_shift` is rejected explicitly.
- Retained only `extract_config` in the metadata helper.
- Replaced upstream usage documentation with this extraction notice. Original
  LICENSE and MODEL WEIGHTS LICENSE notices are unchanged.

The adapter uses FP32 with TF32 disabled and validates all 50 selected model
hashes. Scientific parity with the historical frozen inference must be checked
using the same ready weights; removing the fitting framework is not a license
change or a grant of additional rights.
