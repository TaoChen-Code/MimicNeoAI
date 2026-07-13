# Binding Prediction Tests

Run the dependency-light contract and golden-fixture suite from the repository root:

```bash
PYTHONPATH=. python -m unittest discover -v \
  mimicneoai/functions/binding_prediction/tests
```

The suite does not invoke external binding predictors. It covers mutation peptide
window construction, WT/MT handling, HLA-II pairing, normalized predictor fields,
wide-table summary metrics, resume behavior, and the oversized-task scale gate.
