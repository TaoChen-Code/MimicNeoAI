# coding=utf-8
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import yaml

from mimicneoai.immunogenicity_prediction.core import (
    InferenceConfig,
    export_model_to_onnx,
    run_inference,
)


def _load_yaml(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data


def _build_cfg_dict(raw: Dict[str, Any]) -> Dict[str, Any]:
    path_cfg = raw.get("path", {})
    args_cfg = raw.get("args", {})
    io_cfg = raw.get("io", {})

    cfg = {
        "input_csv": path_cfg.get("input_csv"),
        "output_csv": path_cfg.get("output_csv"),
        "model_path": path_cfg.get("model_path"),
        "hla_fasta": path_cfg.get("hla_fasta"),
        "peptide_col": io_cfg.get("peptide_col", "peptide"),
        "hla_col": io_cfg.get("hla_col", "hla"),
        "score_col": io_cfg.get("score_col", "immunogenicity_score"),
        "batch_size": int(args_cfg.get("batch_size", 512)),
        "num_processes": int(args_cfg.get("num_processes", 1)),
        "device": str(args_cfg.get("device", "auto")),
        "verbose": bool(args_cfg.get("verbose", True)),
        "export_onnx": path_cfg.get("export_onnx", ""),
        "export_onnx_only": bool(args_cfg.get("export_onnx_only", False)),
        "onnx_opset": int(args_cfg.get("onnx_opset", 17)),
    }
    return cfg


def run_from_config(config_path: str) -> int:
    cfg_raw = _load_yaml(config_path)
    cfg = _build_cfg_dict(cfg_raw)

    if not cfg["model_path"]:
        raise ValueError("path.model_path is required in config.")

    if cfg["export_onnx"]:
        export_path = Path(cfg["export_onnx"])
        export_path.parent.mkdir(parents=True, exist_ok=True)
        export_model_to_onnx(
            model_path=cfg["model_path"],
            onnx_path=str(export_path),
            device_name=cfg["device"],
            opset_version=cfg["onnx_opset"],
        )
        print(f"ONNX exported: {export_path}")
        if cfg["export_onnx_only"]:
            return 0

    for req in ["input_csv", "output_csv", "hla_fasta"]:
        if not cfg[req]:
            raise ValueError(
                f"path.{req} is required for inference mode. "
                "If exporting only ONNX, set args.export_onnx_only: true."
            )

    df = pd.read_csv(cfg["input_csv"])

    infer_cfg = InferenceConfig(
        model_path=cfg["model_path"],
        hla_fasta_path=cfg["hla_fasta"],
        peptide_col=cfg["peptide_col"],
        hla_col=cfg["hla_col"],
        output_score_col=cfg["score_col"],
        batch_size=cfg["batch_size"],
        device=cfg["device"],
        num_processes=cfg["num_processes"],
        verbose=cfg["verbose"],
    )

    out_df = run_inference(df, infer_cfg)

    output_path = Path(cfg["output_csv"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, index=False)

    print(f"Done. Saved: {output_path}")
    print(f"Rows: {len(out_df)}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Immunogenicity prediction pipeline")
    parser.add_argument("-c", "--configure", required=True, help="Path to configuration YAML")
    args = parser.parse_args(argv)
    return run_from_config(args.configure)


if __name__ == "__main__":
    raise SystemExit(main())
