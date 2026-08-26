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
