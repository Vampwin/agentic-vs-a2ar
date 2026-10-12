"""ADR-AIIM-420: the clean redocks of the G3 anchor (seeds SEED+1 to SEED+8), run in parallel with the frozen redocking script."""
import sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import ext2

OUTD = ext2.OUT / "poses_clean"

def one(k):
    t0 = time.time()
    ext2.redock(OUTD / f"redock_seed{ext2.SEED + k}.sdf", ext2.CENTER, ext2.SEED + k, f"c{k}")
    return k, round(time.time() - t0, 1)

if __name__ == "__main__":
    with ThreadPoolExecutor(8) as ex:
        for k, s in ex.map(one, range(1, 9)):
            print("seed +", k, "done", s, "s", flush=True)
