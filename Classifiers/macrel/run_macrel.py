import numpy as np
from macrel_features import get_sequence_features
import gzip
import onnxruntime as rt

model1 = "models/AMP.onnx.gz"
model2 = "models/Hemo.onnx.gz"

with gzip.open(model1, 'rb') as f:
    model1 = rt.InferenceSession(f.read(), providers=["CPUExecutionProvider"])

with gzip.open(model2, 'rb') as f:
    model2 = rt.InferenceSession(f.read(), providers=["CPUExecutionProvider"])

def predict(features):
    [amp_prob] = model1.run(['output_probability'], {'input_features': features.astype(np.float32)})
    [hemo_prob] = model2.run(['output_probability'], {'input_features': features.astype(np.float32)})
    return amp_prob[0]['AMP'], hemo_prob[0]['Hemo']

seq = "MKTYTKPTLTKKGKLSAITAGPTGNGTSPV"
def run_macrel(seq):
    features = np.expand_dims(get_sequence_features(seq), axis=0)
    amp_prob, hemo_prob = predict(features)
    return amp_prob, hemo_prob

print(run_macrel(seq))