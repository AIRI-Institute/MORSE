"""Seed the native RiboDecode entry module in the process that actually optimizes."""
import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import random
import runpy
import sys
import time

import numpy as np
import torch

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, required=True)
parser.add_argument('--pilot', action='store_true')
parser.add_argument('--config', type=Path, required=True)
args = parser.parse_args()
cfg = json.loads(args.config.read_text())
package = importlib.metadata.distribution('RiboDecode')
source = package.locate_file('RiboDecode/train_optim_sim.py')
assert hashlib.sha256(source.read_bytes()).hexdigest() == 'bd0150e06fe99ee64b594a8c9d3fd3e03cb04a163133746f75d689f39501eb90'
assert os.environ['PYTHONHASHSEED'] == str(args.seed)
random.seed(args.seed)
np.random.seed(args.seed)
torch.manual_seed(args.seed)
torch.cuda.manual_seed_all(args.seed)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True
torch.use_deterministic_algorithms(True)
torch.set_num_threads(2)
# Native code uses cpu_count()-1 only for deterministic folding workers.
# Bound that pool to the assigned four-core lane without changing scoring.
multiprocessing.cpu_count = lambda: 4
epochs = 1 if args.pilot else cfg['epochs']
sys.argv = [str(source), '--cds', 'MLANA', '--cds_seq', cfg['initial_sequence'],
    '--env', cfg['environment'], '--mfe_weight', str(cfg['mfe_weight']), '--optim_epoch', str(epochs),
    '--alpha', str(cfg['alpha']), '--beta', str(cfg['beta'])]
original_loads = json.loads
def loads(value, *a, **kw):
    data = original_loads(value, *a, **kw)
    if isinstance(data, dict) and 'training' in data and 'optim_num' in data['training']:
        training = data['training']
        assert training['optim_num'] == 40000 and training['optim_batchsize'] == 20
        assert training['num_workers'] == 8 and training['sampler_type'] == 'dist-optim'
        if args.pilot:
            training['optim_num'] = 400
    return data
json.loads = loads
started = time.monotonic()
runpy.run_path(str(source), run_name='__main__')
json.loads = original_loads
output = Path(f'results_natural/HEK293T_optim_mfe_dist-optim/{epochs-1}/samples/optim_results.txt')
assert len(output.read_text().splitlines()) == (400 if args.pilot else 40000)
Path('seed_receipt.json').write_text(json.dumps({'seed': args.seed, 'pilot': args.pilot, 'epochs': epochs,
    'final_pool_rows': 400 if args.pilot else 40000, 'elapsed_seconds': time.monotonic()-started,
    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'torch': torch.__version__,
    'deterministic_algorithms': torch.are_deterministic_algorithms_enabled(),
    'final_pool_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}, indent=2) + '\n')
