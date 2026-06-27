__all__ = [
    "InferenceConfig",
    "run_inference",
    "export_model_to_onnx",
]


def __getattr__(name):
    if name in __all__:
        from mimicneoai.immunogenicity_prediction.core import (
            InferenceConfig,
            export_model_to_onnx,
            run_inference,
        )

        values = {
            "InferenceConfig": InferenceConfig,
            "run_inference": run_inference,
            "export_model_to_onnx": export_model_to_onnx,
        }
        return values[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
