# Package lists of the three environments, exported on 2026-10-11 (ADR-AIIM-401)

These files list what was installed in the environments on the day of export. They were taken after the runs reported in the
article and are not a reconstruction of the environments as they stood when each run was made; the versions named in
Supplementary Section 4 (read at the time from each run's own logs) agree with them for the packages checked below.

| File | Environment | How it was made |
|---|---|---|
| `windows_python312_pip_freeze.txt`, `windows_python_version.txt` | Windows, Python 3.12.10: evaluators, standardisation, modelling, tests, analysis scripts | `python -m pip freeze` |
| `diffdock_env.yml`, `diffdock_pip_freeze.txt`, `diffdock_python.txt` | WSL Ubuntu 22.04, conda environment `diffdock`: conformer generation, AutoDock Vina, DiffDock | `conda env export --no-builds`, `pip freeze` |
| `md_ambertools26_env.yml`, `md_gromacs_version.txt` | WSL, simulation host: AmberTools 26 environment and the GROMACS build | `conda env export --no-builds -p <env>`, `gmx --version` (selected lines) |

The `prefix:` line of the two conda exports (a path in the exporting user's home directory) was removed and the `name:` line of the simulation export, which was a path, was shortened to the name of the environment; nothing else was edited.

Checked against Supplementary Section 4: RDKit 2026.03.3, scikit-learn 1.9.0, pandas 3.0.3 and NumPy 2.4.6 (Windows); RDKit
2025.09.5, AutoDock Vina 1.2.7 and PyTorch 2.7.0+cu128 (docking); GROMACS 2026.3 with the CUDA runtime 13.10 (simulation).

To rebuild an environment from a list: `pip install -r windows_python312_pip_freeze.txt` in Python 3.12, or
`conda env create -f diffdock_env.yml`. This was not tested on a second machine.
