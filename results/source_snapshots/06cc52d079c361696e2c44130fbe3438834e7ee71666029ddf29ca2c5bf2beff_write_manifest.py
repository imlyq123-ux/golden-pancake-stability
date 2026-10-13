"""Write results/artifact_manifest.json: SHA256 of every project file (sources, data, figures, raw outputs).

Build products and environments (.venv*, .wheels, .deps, caches, compiled tools) are excluded.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXCLUDE_DIRS = {'.venv', '.venv-repro', '.wheels', '.deps', '__pycache__', '.pytest_cache'}
EXCLUDE_FILES = {'tools/taylor_mpfr', 'tools/gbs_mpfr', 'results/artifact_manifest.json'}


def main():
    files = {}
    for path in sorted(ROOT.rglob('*')):
        rel = path.relative_to(ROOT)
        if not path.is_file() or EXCLUDE_DIRS & set(rel.parts) or str(rel) in EXCLUDE_FILES:
            continue
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1 << 20), b''):
                digest.update(block)
        files[str(rel)] = digest.hexdigest()
    out = ROOT / 'results' / 'artifact_manifest.json'
    out.write_text(json.dumps({'algorithm': 'sha256', 'count': len(files), 'files': files}, indent=1))
    print(f'{len(files)} files hashed -> {out}')


if __name__ == '__main__':
    main()
