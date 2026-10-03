from pathlib import Path
import re
import subprocess
from common import ROOT, table

git = ['git', '-c', f'safe.directory={ROOT.as_posix()}']
objects = subprocess.check_output(git + ['rev-list', '--objects', '--all'], cwd=ROOT, text=True).splitlines()
batch = subprocess.run(git + ['cat-file', '--batch-check=%(objectname) %(objecttype) %(objectsize)'], input='\n'.join(line.split(' ', 1)[0] for line in objects)+'\n', text=True, capture_output=True, check=True, cwd=ROOT).stdout.splitlines()
checks = []
def record(name, ok, detail):
    checks.append(dict(check=name, status='PASS' if ok else 'FAIL', detail=detail))

blobs = [(row.split()[0], int(row.split()[2])) for row in batch if row.split()[1] == 'blob']
record('historical_blob_size', all(size <= 50_000_000 for _, size in blobs), f'{len(blobs)} historical blobs; largest {max(size for _, size in blobs)} bytes')
paths = [line.split(' ', 1)[1] for line in objects if ' ' in line]
prohibited = [path for path in paths if re.search(r'(^|/)(\.venv|\.cache|__pycache__|node_modules|weights)(/|$)|^data/(raw|processed|external|generated|remote)/|\.(fastq|fq)(\.gz)?$|\.(pt|pth|ckpt|safetensors)$', path)]
record('historical_data_exclusions', not prohibited, f'{len(prohibited)} prohibited environment, sequencing or weight paths')
patterns = [rb'gh[pousr]_[A-Za-z0-9]{20,}', rb'github_pat_[A-Za-z0-9_]{30,}', rb'AKIA[A-Z0-9]{16}', rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----', rb'xox[baprs]-[A-Za-z0-9-]{20,}']
credential_hits = []
for oid, size in blobs:
    if size > 50_000_000:
        continue
    content = subprocess.check_output(git + ['cat-file', 'blob', oid], cwd=ROOT)
    if any(re.search(pattern, content) for pattern in patterns):
        credential_hits.append(oid)
record('credential_pattern_scan', not credential_hits, f'{len(credential_hits)} recognized credential or private-key patterns in historical blobs')
messages = subprocess.check_output(git + ['log', '--format=%s'], cwd=ROOT, text=True).splitlines()
record('commit_subjects', all(re.fullmatch(r'[a-z]+(?:_[a-z]+){1,2}', message) for message in messages), f'{len(messages)} subjects checked against required naming convention')
tracked = subprocess.check_output(git + ['ls-files', '-z'], cwd=ROOT).decode().split('\0')
sizes = [(name, (ROOT/name).stat().st_size) for name in tracked if name and (ROOT/name).is_file()]
record('current_file_size', all(size <= 50_000_000 for _, size in sizes), f'{len(sizes)} tracked files; largest {max(size for _, size in sizes)} bytes')
metadata = __import__('json').loads((ROOT/'docs/github_metadata.json').read_text())
record('github_topics', len(metadata['topics']) <= 20 and all(re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', topic) for topic in metadata['topics']), f'{len(metadata["topics"])} lowercase scientific topics')
table(ROOT/'results/publication_preflight.tsv', checks)
for check in checks:
    print(check['status']+': '+check['check']+': '+check['detail'])
if any(check['status'] != 'PASS' for check in checks):
    raise SystemExit(1)
