"""
F02 HASH_MISMATCH — G0 fault, severity: high

Injection: Modify one byte of a registered input file after its hash has been
logged in the artifact registry.  The G0 hash-integrity check should catch the
mismatch between the stored hash and the current file content.

Parameterisation (15 instances): byte_offset (0-indexed, relative to file start)
× flip_bit_position (0–7).  Offsets sample every ~6.7% of the file; bit positions
cycle 0–7.  All values are structural—no assay data involved.

No assay data used: only a byte position and bit-flip mask are touched.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fault_injection.base import FaultInjector, FaultSpec, InjectionRecord

# 15 instances: 5 offset fractions × 3 bit positions
_OFFSETS = [0.1, 0.2, 0.4, 0.6, 0.8]   # fraction of file length
_BITS = [0, 3, 7]
_PARAM_GRID: list[dict[str, Any]] = [
    {"offset_fraction": o, "flip_bit": b}
    for o in _OFFSETS
    for b in _BITS
]


class F02HashMismatch(FaultInjector):
    """Flip one bit in a registered input file to break its stored hash."""

    @property
    def spec(self) -> FaultSpec:
        return FaultSpec(
            family_id="F02",
            family_name="HASH_MISMATCH",
            gate="G0",
            severity="high",
            description=(
                "Flip one bit at a deterministic byte offset in a registered "
                "workflow input file after its SHA-256 has been logged.  "
                "Expected detection: G0 hash-integrity re-check."
            ),
            param_grid=_PARAM_GRID,
        )

    def verify_no_assay_leakage(self) -> None:
        for p in _PARAM_GRID:
            assert "camp" not in str(p).lower()
            assert "assay" not in str(p).lower()

    def inject(
        self,
        workspace: Path,
        params: dict[str, Any],
        seed: int,
        instance_index: int,
    ) -> InjectionRecord:
        # Expect a target file called 'registered_input.bin' or any first file found
        candidates = list(workspace.glob("registered_input.*"))
        if not candidates:
            candidates = [p for p in workspace.iterdir() if p.is_file()]
        if not candidates:
            raise FileNotFoundError("F02: no file to corrupt in workspace")

        target = sorted(candidates)[0]
        data = bytearray(target.read_bytes())
        if not data:
            raise ValueError(f"F02: {target} is empty")

        byte_idx = max(0, int(len(data) * params["offset_fraction"]) - 1)
        byte_idx = min(byte_idx, len(data) - 1)
        flip_mask = 1 << params["flip_bit"]
        original_byte = data[byte_idx]
        data[byte_idx] ^= flip_mask
        target.write_bytes(bytes(data))

        return InjectionRecord(
            instance_id=self.instance_id(instance_index),
            family_id="F02",
            family_name="HASH_MISMATCH",
            gate="G0",
            severity="high",
            params=params,
            affected_items=[str(target.name)],
            ground_truth=(
                f"Flipped bit {params['flip_bit']} at byte {byte_idx} "
                f"(offset_fraction={params['offset_fraction']}) in {target.name}. "
                f"Original byte: 0x{original_byte:02x}"
            ),
            containment_criterion="G0 contract SHA-256 re-verification of registered inputs",
        )


F02HashMismatch().verify_no_assay_leakage()
