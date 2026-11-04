#!/bin/bash
# Build the pristine upstream Ribo-TISH v0.2.8 reference install under work/upstream.
# Same pinned dependencies as the fork env (requirements.lock); only the code differs.
set -euo pipefail
cd "$(dirname "$0")/.."
[ -d work/upstream ] || git worktree add --detach work/upstream v0.2.8
[ -x work/upstream/.venv/bin/python ] || uv venv work/upstream/.venv --python 3.11
uv pip install --python work/upstream/.venv/bin/python -r requirements.lock work/upstream
work/upstream/.venv/bin/ribotish --version
