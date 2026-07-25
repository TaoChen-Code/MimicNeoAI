__all__ = [
    "InferenceConfig",
    "prepare_inference_features",
    "run_inference",
    "export_model_to_onnx",
    "default_model_paths",
    "run_default_inference",
]


def __getattr__(name):
    if name in __all__:
        from mimicneoai.immunogenicity_prediction.core import (
            InferenceConfig,
            export_model_to_onnx,
            prepare_inference_features,
            run_inference,
        )
        from mimicneoai.immunogenicity_prediction.default_models import (
            default_model_paths,
            run_default_inference,
        )

        values = {
            "InferenceConfig": InferenceConfig,
            "prepare_inference_features": prepare_inference_features,
            "run_inference": run_inference,
            "export_model_to_onnx": export_model_to_onnx,
            "default_model_paths": default_model_paths,
            "run_default_inference": run_default_inference,
        }
        return values[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
