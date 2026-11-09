from pathlib import Path
import hashlib
import json
import os
import sys
import numpy as np
import pandas as pd
import torch
from encoding import encode

sys.path.insert(0, str(Path(__file__).resolve().parent / 'third_party/ribonn'))
from src.model import RiboNN
from src.utils.helpers import extract_config


class Predictor:
    def __init__(self):
        self.directory = Path(os.environ.get('RIBONN_MODEL_DIR', '/models'))
        self.device = os.environ.get('RIBONN_DEVICE', 'cuda' if torch.cuda.is_available() else 'cpu')
        self.batch_size = int(os.environ.get('RIBONN_BATCH_SIZE', '16'))
        if self.batch_size < 1:
            raise ValueError('Batch size must be positive')
        torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', '2')))
        # FP32 reference inference; no AMP, cropping, or stochastic shifts.
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        self.manifest = json.loads((self.directory / 'ensemble.json').read_text())
        if self.manifest['n_models'] != 50 or self.manifest['n_targets'] != 78:
            raise ValueError('Expected the official human 50-model / 78-output ensemble')
        runs = pd.read_csv(self.directory / 'human/runs.csv')
        if hashlib.sha256((self.directory / 'human/runs.csv').read_bytes()).hexdigest() != self.manifest['runs_csv_sha256']:
            raise ValueError('Run metadata checksum mismatch')
        self.config = extract_config(runs, runs.run_id.iloc[0])
        expected = {'pad_5_prime': True, 'label_codons': True, 'max_shift': 0,
                    'split_utr5_cds_utr3_channels': False, 'label_3rd_nt_of_codons': False,
                    'label_utr5': False, 'label_utr3': False, 'label_splice_sites': False,
                    'label_up_probs': False, 'num_targets': 78, 'len_after_conv': 9}
        if any(self.config[k] != v for k, v in expected.items()):
            raise ValueError('Unexpected upstream input architecture')
        self.models = []
        for item in self.manifest['models']:
            path = self.directory / item['path']
            if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
                raise ValueError('Checkpoint checksum mismatch: ' + item['run_id'])
            model = RiboNN(**self.config)
            model.load_state_dict(torch.load(path, map_location='cpu', weights_only=True), strict=True)
            self.models.append(model.eval().to(self.device))

    @torch.inference_mode()
    def predict(self, sequences, sequence_contexts):
        outputs = []
        for offset in range(0, len(sequences), self.batch_size):
            x = torch.from_numpy(encode(sequences[offset:offset + self.batch_size],
                                       sequence_contexts[offset:offset + self.batch_size])).to(self.device)
            # Preserve the authors' order: average 5 models within each fold,
            # then 10 folds, then the 78 cellular TE outputs.
            per_model = torch.stack([model(x) for model in self.models])
            per_fold = per_model.reshape(10, 5, len(x), 78).mean(1)
            predictions = per_fold.mean(0).mean(-1)
            spread = per_model.mean(-1).std(0, unbiased=True)
            if not torch.isfinite(predictions).all():
                raise ValueError('Nonfinite RiboNN output')
            outputs.extend({'ribonn_mean_te': float(p), 'ribonn_model_sd': float(s),
                            'ribonn_model_count': 50} for p, s in zip(predictions.cpu(), spread.cpu()))
        return outputs

