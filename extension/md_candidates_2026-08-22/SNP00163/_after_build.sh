#!/bin/bash
# ADR-AIIM-331: when the SNP00163 build finishes, run md_09 (frozen), rewrite the production
# mdp files to 50 ns, then start the retry driver. Timeout 5 h on the build.
set -u
D="/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163"
export PYTHONUTF8=1 PYTHONPATH="D:/Agentic VS Workflow/10_new_system/scripts" MSYS_NO_PATHCONV=1
for i in $(seq 1 300); do
  grep -q "build complete" "$D/build_system.log" 2>/dev/null && break
  grep -q "FAILED\|FATAL" "$D/build_system.log" 2>/dev/null && { echo "BUILD ERROR"; tail -5 "$D/build_system.log"; exit 2; }
  sleep 60
done
grep -q "build complete" "$D/build_system.log" || { echo "BUILD TIMEOUT"; exit 3; }
echo "[$(date '+%F %T')] build complete"; tail -3 "$D/build_system.log"
python "D:/Agentic VS Workflow/10_new_system/scripts/md_09_setup_candidate_runs.py" SNP00163 2>&1 | tee "$D/md_09.out" | tail -8
for f in "$D"/04_run/prod_r*.mdp; do
  sed -i 's/^nsteps .*/nsteps                  = 25000000   ; 50 ns (ADR-AIIM-331; md_09 writes 100 ns)/' "$f"
done
grep -H "^nsteps" "$D"/04_run/prod_r*.mdp
wsl.exe -d A2A-MD -- bash -c 'cd "/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163"; setsid nohup bash run_snp00163.sh > driver_snp00163.log 2>&1 < /dev/null & sleep 3; pgrep -af run_snp00163 | head -2'
cat > "$D/_launch_stopper163.sh" <<'EOS'
#!/bin/bash
cd "/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163" || exit 1
setsid nohup bash stop_after_r2.sh SNP00163 run_snp00163.sh > /dev/null 2>&1 < /dev/null &
sleep 2; pgrep -af stop_after_r2 | cut -c1-100
EOS
sed -i 's/\r$//' "$D/_launch_stopper163.sh"
wsl.exe -d A2A-MD -- bash "/mnt/d/Agentic VS Workflow/20_new_runs/md_candidates_2026-08-22/SNP00163/_launch_stopper163.sh"
echo "[$(date '+%F %T')] AFTER_BUILD DONE"
