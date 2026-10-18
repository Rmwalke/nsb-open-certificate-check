"""Verify this package's file integrity; this is not a mathematical proof check."""
from pathlib import Path
import hashlib
import json
import sys

base = Path(__file__).resolve().parent
manifest = json.loads((base / 'MANIFEST.json').read_text(encoding='utf-8'))
errors = []
for row in manifest['files']:
    path = base / row['path']
    if not path.is_file():
        errors.append('missing: ' + row['path'])
        continue
    data = path.read_bytes()
    if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
        errors.append('changed: ' + row['path'])
expected = {r['path'] for r in manifest['files']} | {'MANIFEST.json'}
extras = sorted(str(p.relative_to(base)) for p in base.rglob('*')
                if p.is_file() and '__pycache__' not in p.parts
                and str(p.relative_to(base)) not in expected)
print(json.dumps({'status': 'PASS' if not errors and not extras else 'FAIL',
                  'payload_files': len(manifest['files']),
                  'errors': errors, 'extra_files': extras}, indent=2))
sys.exit(bool(errors or extras))
