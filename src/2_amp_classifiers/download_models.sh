#!/usr/bin/env bash
# Download or refresh classifier model assets into classifiers/_downloads and link
# them into each classifier directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOADS="${SCRIPT_DIR}/_downloads"
GIT_BIN="/opt/ohpc/pub/apps/git/2.9.5/bin"
export PATH="${GIT_BIN}:${PATH}"

mkdir -p "${DOWNLOADS}"

clone_or_pull() {
  local url="$1"
  local dest="$2"
  if [[ -d "${dest}/.git" ]]; then
    echo "Updating ${dest} ..."
    git -C "${dest}" pull --ff-only
  else
    echo "Cloning ${url} -> ${dest} ..."
    git clone --depth 1 "${url}" "${dest}"
  fi
}

# AMP Scanner v2 models
clone_or_pull "https://github.com/dan-veltri/amp-scanner-v2.git" "${DOWNLOADS}/amp-scanner-v2"
mkdir -p "${SCRIPT_DIR}/amp-scanner-v2/trained-models"
for model in OriginalPaper_081917_FULL_MODEL.h5 020419_FULL_MODEL.h5; do
  src="${DOWNLOADS}/amp-scanner-v2/trained-models/${model}"
  dst="${SCRIPT_DIR}/amp-scanner-v2/trained-models/${model}"
  if [[ -s "${src}" ]]; then
    cp -f "${src}" "${dst}"
    echo "Installed ${model}"
  else
    echo "WARNING: missing ${src}" >&2
  fi
done

# APEX pathogen models
clone_or_pull "https://gitlab.com/machine-biology-group-public/apex-pathogen.git" "${DOWNLOADS}/apex-pathogen"
mkdir -p "${SCRIPT_DIR}/apex/APEX_pathogen_models"
rm -rf "${SCRIPT_DIR}/apex/APEX_pathogen_models"
ln -sfn "${DOWNLOADS}/apex-pathogen/APEX_pathogen_models" "${SCRIPT_DIR}/apex/APEX_pathogen_models"
for f in APEX_predict.py APEX_models.py utils.py aaindex1.csv; do
  cp -f "${DOWNLOADS}/apex-pathogen/${f}" "${SCRIPT_DIR}/apex/${f}"
done
echo "Installed APEX models and scripts"

# HydrAMP full source + classifier models
clone_or_pull "https://github.com/szczurek-lab/hydramp.git" "${DOWNLOADS}/hydramp-full"
HYDRAMP_SRC="${DOWNLOADS}/hydramp-full"
HYDRAMP_DST="${SCRIPT_DIR}/hydramp"

if [[ ! -d "${HYDRAMP_SRC}/models/amp_classifier/layers" ]]; then
  echo "Downloading HydrAMP models from Google Drive (first run only, ~500 MB) ..."
  conda run -n gmconda_py3923 python -m pip install -q gdown 2>/dev/null || true
  (cd "${HYDRAMP_SRC}" && python3 data_setup/download_data.py)
  unzip -o "${HYDRAMP_SRC}/downloaded_data_zips/models.zip" -d "${HYDRAMP_SRC}"
fi

# Sync HydrAMP Python package and models into the project copy.
rsync -a --delete \
  --exclude '.git' \
  --exclude 'downloaded_data_zips' \
  --exclude 'genbio' \
  "${HYDRAMP_SRC}/amp/" "${HYDRAMP_DST}/amp/"
rsync -a "${HYDRAMP_SRC}/models/" "${HYDRAMP_DST}/models/"
bash "${SCRIPT_DIR}/apply_hydramp_patches.sh"
echo "Installed HydrAMP source and classifier models"

echo "Model assets ready."
