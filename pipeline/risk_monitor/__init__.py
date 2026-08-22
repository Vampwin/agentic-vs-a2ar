"""
risk_monitor — Calibrated learned risk monitor for the hybrid KDSS audit architecture.

Public interface
----------------
PolicyAction, evaluate_policy  — frozen three-region policy (policy.py)
InferenceResult, predict_risk   — fail-closed inference (inference.py)
InferenceError                  — raised on any unsafe inference condition
ModelRegistry, RegistryError    — trusted artifact registry (registry.py)

Training and calibration (train.py, calibrate.py) are separate from inference so
that the acting agent cannot invoke training code during an evaluation run.
"""

from .policy import PolicyAction, PolicyResult, evaluate_policy
from .inference import InferenceResult, InferenceError, predict_risk
from .registry import ModelRegistry, RegistryError, RegistryEntry

__all__ = [
    "PolicyAction",
    "PolicyResult",
    "evaluate_policy",
    "InferenceResult",
    "InferenceError",
    "predict_risk",
    "ModelRegistry",
    "RegistryError",
    "RegistryEntry",
]
