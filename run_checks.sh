#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")" && pwd)"
cd "$repo_root"

python scripts/01_preprocess.py --force
pytest -q
python -m py_compile scripts/*.py scripts/figures/*.py
Rscript scripts/check_r_environment.R
Rscript -e 'for (f in list.files("scripts", "[.]R$", recursive=TRUE, full.names=TRUE)) parse(f)'
