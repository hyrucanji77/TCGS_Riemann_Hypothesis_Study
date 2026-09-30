#!/usr/bin/env python3
"""Repository delivery utility: import approved source and package verified builds.

This is build bookkeeping, not a scientific diagnostic or a proof checker.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()

def load(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))

def save(path, data):
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')

def require(ok, message):
    if not ok:
        raise SystemExit(message)

def import_source():
    manifest_path = ROOT / '.upload_source/ready.json'
    if not manifest_path.exists():
        require((ROOT / 'main.tex').is_file(), 'No manuscript source or import manifest.')
        return
    manifest = load('.upload_source/ready.json')
    require(manifest['format'] == 'line-edits-v1', 'Unknown import format.')
    base_commit = manifest['base_commit']
    require(re.fullmatch(r'[0-9a-f]{40}', base_commit), 'Invalid base commit.')
    base = subprocess.check_output(['git', 'show', base_commit + ':main.tex'], cwd=ROOT)
    require(hashlib.sha256(base).hexdigest() == manifest['base_sha256'], 'Base source checksum mismatch.')
    lines = base.decode('utf-8').splitlines(keepends=True)
    result, cursor = [], 0
    for part in manifest['parts']:
        start, end = part['start'], part['end']
        path = Path(part['path'])
        require(path.parent == Path('.upload_source') and path.suffix == '.txt', 'Invalid part path.')
        require(0 <= cursor <= start <= end <= len(lines), 'Invalid source range.')
        require(digest(path) == part['sha256'], 'Source part checksum mismatch: ' + str(path))
        result.extend(lines[cursor:start])
        result.append((ROOT / path).read_bytes().decode('utf-8'))
        cursor = end
    result.extend(lines[cursor:])
    data = ''.join(result).encode('utf-8')
    require(hashlib.sha256(data).hexdigest() == manifest['sha256'], 'Final source checksum mismatch.')
    target = ROOT / 'main.tex'
    require(not target.exists() or target.read_bytes() in (base, data), 'Refusing to overwrite a different source.')
    target.write_bytes(data)
    record = dict(manifest)
    record['status'] = 'PASS'
    record['bytes'] = len(data)
    record['statement'] = 'main.tex is byte-identical to the approved review-amended source. PDF compilation is separate.'
    save('verification/github_source_import.json', record)
    print('Exact reviewed source imported:', manifest['sha256'])

def verify_inputs():
    config = load('release_integrity_20260930.json')
    require(digest('main.tex') == config['manuscript_sha256'], 'Approved manuscript checksum mismatch.')
    for name, expected in config['source_sha256'].items():
        require(digest(name) == expected, 'Source changed: ' + name)
    # Reconstruct only the two new outputs, then verify against their independently
    # recorded archive checksums before treating them as reference files.
    for dps in (40, 60):
        name = f'screw_kernel_{dps}.json'
        if not (ROOT / name).exists():
            subprocess.run([sys.executable, '-B', 'screw_kernel_diagnostics.py', '--dps', str(dps), '--output', name], cwd=ROOT, check=True, timeout=600)
    for name, expected in config['reference_outputs_sha256'].items():
        require(digest(name) == expected, 'Reference output checksum mismatch: ' + name)
    print('Approved source, 12 scripts/drivers, and 15 reference outputs verified.')

def package():
    config = load('release_integrity_20260930.json')
    require(digest('main.tex') == config['manuscript_sha256'], 'Approved manuscript checksum mismatch.')
    for mapping in ('source_sha256', 'reference_outputs_sha256'):
        for name, expected in config[mapping].items():
            require(digest(name) == expected, 'Release input changed: ' + name)
    for path in ('verification/github_diagnostics_report.json', 'verification/github_screw_diagnostics_report.json'):
        report = load(path)
        require(report['status'] == 'PASS' and report['manuscript_sha256'] == digest('main.tex'), 'Missing matching verification: ' + path)
    for line in (ROOT / 'verification/github_build_sha256.txt').read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        require(digest(name.strip()) == expected, 'Build checksum mismatch: ' + name)
    require(digest('main.pdf') == digest('riemann_hypothesis_TCGS_SEQUENTION.pdf'), 'PDF aliases differ.')
    info = (ROOT / 'verification/github_pdf_info.txt').read_text()
    match = re.search(r'^Pages:\s+(\d+)', info, re.M)
    require(match and int(match.group(1)) == config['expected_pages'], 'Unexpected PDF page count.')
    files = set()
    for pattern in ('*.tex', '*.pdf', '*.py', '*.json', 'requirements.txt', 'README.md', 'LICENSE', 'CITATION.cff', 'tools/*.py'):
        files.update(path for path in ROOT.glob(pattern) if path.is_file())
    files.discard(ROOT / 'source_manifest.json')
    manifest = {
        'release': config['release'],
        'main_document': 'main.tex',
        'pdf_provenance': 'Recompiled on GitHub from the exact reviewed source; not claimed byte-identical to the local PDF.',
        'approved_source_sha256': config['manuscript_sha256'],
        'approved_local_pdf_sha256': config['reference_pdf_sha256'],
        'pages': config['expected_pages'],
        'scientific_diagnostic_scripts': len(config['scientific_scripts']),
        'verification_drivers': len(config['verification_drivers']),
        'mathematical_status': 'RH is not proved; the source-derived arithmetic identity and required nonlocal limit remain unconstructed.',
        'excluded_material': 'Historical drafts, private correspondence, third-party source-PDF archive, and infographic are not part of this paper-and-code bundle.',
        'files': [{'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': digest(path)} for path in sorted(files)]
    }
    save('source_manifest.json', manifest)
    files.add(ROOT / 'source_manifest.json')
    files.update(path for path in (ROOT / 'verification').glob('github_*') if path.is_file())
    files.discard(ROOT / 'verification/github_package_record.json')
    target = ROOT / 'TCGS_Riemann_Hypothesis_Study_Overleaf.zip'
    with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):
            archive.write(path, path.relative_to(ROOT).as_posix())
    record = {'status': 'PASS', 'release': config['release'], 'manuscript_sha256': digest('main.tex'),
              'pdf_sha256': digest('main.pdf'), 'bundle_sha256': digest(target),
              'bundle_files': len(files), 'pages': config['expected_pages'], 'citation_sha256': digest('CITATION.cff')}
    save('verification/github_package_record.json', record)
    print(json.dumps(record, indent=2))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('import-source', 'verify-inputs', 'package'))
    args = parser.parse_args()
    {'import-source': import_source, 'verify-inputs': verify_inputs, 'package': package}[args.command]()

if __name__ == '__main__':
    main()
