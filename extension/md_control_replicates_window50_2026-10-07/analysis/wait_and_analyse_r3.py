"""Wait for control replicate 3 (ADR-AIIM-324), then analyse it and grade R1 to R6.

  python wait_and_analyse_r3.py [max_minutes]

Polls every 5 minutes for prod_r3.gro. Liveness is checked by PROCESS, not by log
timestamp (a lesson of 2026-10-06, when a dead job was reported alive for 11 hours):
if no mdrun and no driver exist and the .gro is absent, it stops at once and says so.
When the .gro appears it runs the unchanged analysis pipeline on replicate 3 and then
the summary script, whose output is printed.
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(r"D:\Agentic VS Workflow")
REP = ROOT / "20_new_runs" / "md_control_replicates_2026-10-06"
ANA = ROOT / "20_new_runs" / "md_control_replicates_window50_2026-10-07" / "analysis"
GRO = REP / "prod_r3.gro"
WSL = "A2A-MD"


def wsl(cmd: str, timeout: int = 3000) -> str:
    p = subprocess.run(["wsl", "-d", WSL, "--", "bash", "-lc", cmd], capture_output=True,
                       text=True, timeout=timeout, env={**__import__("os").environ, "MSYS_NO_PATHCONV": "1"})
    return (p.stdout or "") + (p.stderr or "")


def alive() -> str:
    out = wsl("pgrep -af 'gmx mdrun|run_control_replicates' | grep -v pgrep | cut -c1-80")
    return out.strip()


def last_step() -> str:
    try:
        lines = (REP / "prod_r3.log").read_text(errors="ignore").splitlines()
        for i in range(len(lines) - 1, 0, -1):
            if lines[i].split()[:2] == ["Step", "Time"]:
                return lines[i + 1].strip()
    except OSError:
        pass
    return "?"


def main() -> int:
    max_min = float(sys.argv[1]) if len(sys.argv) > 1 else 110.0
    t_end = time.time() + max_min * 60
    while not GRO.is_file():
        if time.time() > t_end:
            print(f"[{time.strftime('%F %T')}] replicate 3 not finished after {max_min:g} min; "
                  f"last step {last_step()}. Relaunch to keep waiting.")
            return 2
        a = alive()
        if not a and not GRO.is_file():
            print(f"[{time.strftime('%F %T')}] NO mdrun and NO driver, and no prod_r3.gro. "
                  f"The run died. last step {last_step()}")
            return 3
        time.sleep(300)
    time.sleep(90)                                   # let the files close
    print(f"[{time.strftime('%F %T')}] prod_r3.gro present; last step {last_step()}")
    out = wsl("bash '/mnt/d/Agentic VS Workflow/20_new_runs/md_control_replicates_window50_2026-10-07/"
              "analysis/run_control_analysis.sh' '/mnt/d/Agentic VS Workflow/20_new_runs/"
              "md_control_replicates_2026-10-06' prod_r3 control_r3")
    print(out[-1500:])
    p = subprocess.run([sys.executable, str(ANA / "md_control_replicates_window50.py")],
                       capture_output=True, text=True, encoding="utf-8",
                       env={**__import__("os").environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"})
    print(p.stdout[-6000:], p.stderr[-1500:])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
