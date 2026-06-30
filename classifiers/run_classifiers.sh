#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec conda run -n gmconda_py3923 python "${SCRIPT_DIR}/run_classifiers.py" "$@"
