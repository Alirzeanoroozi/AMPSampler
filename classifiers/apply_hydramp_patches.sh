#!/usr/bin/env bash
# Re-apply AMPSampler compatibility patches to HydrAMP after syncing upstream source.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYDRAMP_DIR="${SCRIPT_DIR}/hydramp"
export HYDRAMP_DIR

cat > "${HYDRAMP_DIR}/amp/models/model_garden.py" <<'EOF'
from amp.models.discriminators import amp_classifier_noCONV
from amp.models.discriminators import veltri_amp_classifier

MODEL_GAREDN = {
    'VeltriAMPClassifier': veltri_amp_classifier.VeltriAMPClassifier,
    'NoConvAMPClassifier': amp_classifier_noCONV.NoConvAMPClassifier,
}


def ensure_model_registered(model_type):
    """Import generation models only when needed (avoids legacy Keras deps for classifiers)."""
    if model_type in MODEL_GAREDN:
        return
    from amp.models.decoders import amp_expanded_decoder
    from amp.models.encoders import amp_expanded_encoder
    from amp.models.master import master

    MODEL_GAREDN.update({
        'AMPExpandedDecoder': amp_expanded_decoder.AMPDecoder,
        'AMPExpandedEncoder': amp_expanded_encoder.AMPEncoder,
        'MasterAMPTrainer': master.MasterAMPTrainer,
    })
EOF

cat > "${HYDRAMP_DIR}/amp/inference/__init__.py" <<'EOF'
"""Classifier-only package init (skip HydrAMPGenerator import)."""

__all__ = []
EOF

python3 <<'PY'
import os
from pathlib import Path

hydramp = Path(os.environ["HYDRAMP_DIR"])

path = hydramp / "amp/utils/basic_model_serializer.py"
text = path.read_text()
if "ensure_model_registered" not in text:
    text = text.replace(
        "        model_class = model_garden.MODEL_GAREDN[model_config['type']]",
        "        model_garden.ensure_model_registered(model_config['type'])\n"
        "        model_class = model_garden.MODEL_GAREDN[model_config['type']]",
    )
if "from amp.models.master.master import MasterAMPTrainer" in text.split("def load_master_model_components")[0]:
    text = text.replace("from amp.models.master.master import MasterAMPTrainer\n", "")
    text = text.replace(
        "def load_master_model_components(model_path: str, return_master=False, softmax=False) -> Tuple[models.Model, ...]:\n    serializer = BasicModelSerializer()",
        "def load_master_model_components(model_path: str, return_master=False, softmax=False) -> Tuple[models.Model, ...]:\n"
        "    from amp.models.master.master import MasterAMPTrainer\n\n"
        "    serializer = BasicModelSerializer()",
    )
path.write_text(text)

seq_path = hydramp / "amp/data_utils/sequence.py"
seq_text = seq_path.read_text()
if "tensorflow.keras.utils" not in seq_text:
    seq_text = seq_text.replace(
        "def pad(x, max_length: int = 25) -> np.ndarray:\n"
        "    return preprocessing.sequence.pad_sequences(\n"
        "        x,\n"
        "        maxlen=max_length,\n"
        "        padding='post',\n"
        "        value=0.0\n"
        "    )",
        "def pad(x, max_length: int = 25) -> np.ndarray:\n"
        "    try:\n"
        "        return preprocessing.sequence.pad_sequences(\n"
        "            x,\n"
        "            maxlen=max_length,\n"
        "            padding='post',\n"
        "            value=0.0,\n"
        "        )\n"
        "    except AttributeError:\n"
        "        from tensorflow.keras.utils import pad_sequences\n\n"
        "        return pad_sequences(x, maxlen=max_length, padding='post', value=0.0)",
    )
    seq_path.write_text(seq_text)
PY

echo "Applied HydrAMP compatibility patches."
