"""ADR-AIIM-420: the F07 redockings (the box centre shifted by the offset along the axis), one after the other with the frozen script."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import ext2

OUTD = ext2.OUT / "poses_f07"
if __name__ == "__main__":
    for i, p in enumerate(ext2.F07_GRID):
        out = OUTD / f"F07_{i:02d}_{p['axis']}{p['offset_angstrom']}.sdf"
        if out.exists():
            continue
        t0 = time.time()
        ext2.redock(out, ext2.f07_center(p), ext2.SEED, f"f07_{i}")
        print(i, p, "done", round(time.time() - t0, 1), "s", flush=True)
    print("all F07 redockings done", flush=True)
