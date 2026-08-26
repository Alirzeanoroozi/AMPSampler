#!/usr/bin/env bash
# Create one conda environment per classifier under classifiers/envs/*.yml
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/envs"

for yml in "${ENV_DIR}"/*.yml; do
  name="$(basename "${yml}" .yml)"
  env_name="ampsampler-${name}"
  if conda env list | awk '{print $1}' | grep -qx "${env_name}"; then
    echo "Conda env ${env_name} already exists, skipping create"
  else
    echo "Creating conda env ${env_name} from ${yml} ..."
    conda env create -f "${yml}" -n "${env_name}"
  fi
done

# Install HydrAMP package in editable mode inside its env.
HYDRAMP_DIR="${SCRIPT_DIR}/hydramp"
if [[ -d "${HYDRAMP_DIR}/amp/inference" ]]; then
  echo "Installing HydrAMP package into ampsampler-hydramp ..."
  # Dependencies are already installed from envs/hydramp.yml; avoid legacy TF pins in setup.py.
  conda run -n ampsampler-hydramp pip install -e "${HYDRAMP_DIR}" --no-deps -q
else
  echo "WARNING: HydrAMP source missing. Run download_models.sh first." >&2
fi

echo "All classifier environments are ready."
echo "  ampsampler-ascan2"
echo "  ampsampler-macrel"
echo "  ampsampler-hydramp"
echo "  ampsampler-apex"
echo "  ampsampler-toxinpred"
