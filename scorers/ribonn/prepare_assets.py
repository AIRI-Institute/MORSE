"""Verify the official archive and freeze the standard top-5-per-fold ensemble."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import urllib.request
import zipfile

RECORD = '17258709'
MD5 = 'dd088ee2feffcda6ce26124a471a140f'
URL = f'https://zenodo.org/api/records/{RECORD}/files/weights.zip/content'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', required=True, type=Path)
    args = parser.parse_args()
    root = args.directory
    root.mkdir(parents=True, exist_ok=True)
    archive = root / 'weights.zip'
    if not archive.exists():
        temporary = root / 'weights.zip.part'
        urllib.request.urlretrieve(URL, temporary)
        temporary.rename(archive)
    if hashlib.md5(archive.read_bytes()).hexdigest() != MD5:
        raise ValueError('Archive checksum mismatch')
    with zipfile.ZipFile(archive) as z:
        for item in z.infolist():
            destination = (root / item.filename).resolve()
            if not destination.is_relative_to(root.resolve()):
                raise ValueError('Unsafe archive path')
        z.extractall(root)
    with (root / 'human/runs.csv').open() as f:
        rows = list(csv.DictReader(f))
    selected = []
    for fold in range(10):
        subset = [r for r in rows if int(r['params.test_fold']) == fold]
        subset.sort(key=lambda r: -float(r['metrics.val_r2']))
        if len(subset) != 9:
            raise ValueError('Expected exactly nine released model runs per fold')
        for row in subset[:5]:
            weight = root / 'human' / row['run_id'] / 'state_dict.pth'
            selected.append({'run_id': row['run_id'], 'test_fold': fold,
                             'validation_R2': float(row['metrics.val_r2']),
                             'path': str(weight.relative_to(root)),
                             'sha256': hashlib.sha256(weight.read_bytes()).hexdigest()})
    manifest = {'record': RECORD, 'species': 'human', 'selection': 'top 5 validation R2 per outer fold',
                'models': selected, 'n_models': 50, 'n_targets': 78,
                'weights_archive_md5': MD5, 'weights_archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
                'runs_csv_sha256': hashlib.sha256((root / 'human/runs.csv').read_bytes()).hexdigest()}
    (root / 'ensemble.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({k: v for k, v in manifest.items() if k != 'models'}, indent=2))


if __name__ == '__main__':
    main()

