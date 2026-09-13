"""Where a file lives, in whichever of the two layouts this code is running in.

The scripts in this project were written against the development tree and name
their inputs by its paths -- `20_new_runs/<run>/<file>`, `02_frozen_inputs/...`,
`30_manuscript/...`. The public deposit ships the same evidence under a
different arrangement: code under `pipeline/`, artifacts regrouped by subject
under `verification/`, and figures under `figures/`. A script that hardcodes the
first layout cannot run against the second, which is what a reader who clones
the repository actually has.

Until 2026-09-13 that was resolved by not resolving it: every script carried

    for cand in (Path(r"D:\\Agentic VS Workflow"), Path("/mnt/d/...")):

and exited if neither existed. On the author's machine it found the development
tree even when invoked from inside a copy of the deposit, so the scripts read --
and wrote -- the author's working files rather than the copy. On anyone else's
machine they exited. Both behaviours were wrong and the second was the visible
one.

**The mapping is not guessed.** `build_deposit.py` records the development-tree
`source` of every file it copies, and `MANIFEST.json` publishes that alongside
the SHA-256. This module reads the manifest and inverts it, so a deposit
resolves its own paths from its own published record. Where the manifest does
not name a file, `resolve()` says so rather than returning a path that does not
exist.

Usage, in place of `ROOT / "20_new_runs/<run directory>/<file>"`:

    from kdss_paths import resolve
    resolve("20_new_runs/<run directory>/<file>")

In the development tree that returns exactly the path it did before. In a
deposit it returns `verification/<group>/<name>`. `project_root()` remains
available for code that only needs the root itself.
"""
from __future__ import annotations

import functools
import json
from pathlib import Path

__all__ = ["project_root", "resolve", "under", "dev_dir", "layout",
           "wsl_root", "exists", "PathNotInDeposit"]

# Files that identify each layout. The development tree is identified by a
# project-control document that is deliberately not deposited; the deposit by
# its own manifest together with the two directories it reorganises into.
_DEV_MARKER = "00_project_control/08_DECISION_LOG.md"
_DEPOSIT_MARKERS = ("MANIFEST.json", "pipeline", "verification")


class PathNotInDeposit(FileNotFoundError):
    """A development-tree path with no counterpart in this deposit.

    Raised rather than returning a path that does not exist, because a silent
    miss is how an analysis ends up reporting a number from the wrong file.
    """


def _is_dev(p: Path) -> bool:
    return (p / _DEV_MARKER).is_file()


def _is_deposit(p: Path) -> bool:
    return all((p / m).exists() for m in _DEPOSIT_MARKERS)


@functools.lru_cache(maxsize=1)
def _root_and_layout() -> tuple[Path, str]:
    """Walk up from this file. Nearest enclosing layout wins.

    Walking from `__file__` rather than from the working directory is what lets
    a copy of the deposit resolve to itself: the module that is imported is the
    one inside that copy.
    """
    here = Path(__file__).resolve()
    for cand in (here.parent, *here.parents):
        if _is_dev(cand):
            return cand, "development"
        if _is_deposit(cand):
            return cand, "deposit"
    raise SystemExit(
        "kdss_paths: no project root above " + str(here) + ". Expected either "
        f"a development tree (with {_DEV_MARKER}) or a deposit (with "
        f"{', '.join(_DEPOSIT_MARKERS)})."
    )


def project_root() -> Path:
    """The root of whichever layout this code is running inside."""
    return _root_and_layout()[0]


def layout() -> str:
    """``"development"`` or ``"deposit"``."""
    return _root_and_layout()[1]


@functools.lru_cache(maxsize=1)
def _reverse_index() -> dict[str, str]:
    """development-tree path -> path inside this deposit.

    Two sources, both written by `build_deposit.py`. `MANIFEST.json` records
    the `source` of every file it copied, which covers anything taken straight
    from the development tree. Most of the evidence, though, reaches the
    deposit through the S10 supplement bundle, so its recorded source is a path
    under `30_manuscript/supplement/` rather than the run directory the
    analysis scripts name -- `PATH_MAP.json` carries that second hop, derived
    from `build_supplement_s10.origin_map()`.

    `PATH_MAP.json` wins where both speak, because it is the end of the chain.
    """
    root = project_root()
    index: dict[str, str] = {}

    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    entries = manifest["files"] if isinstance(manifest, dict) else manifest
    for e in entries:
        src = e.get("source")
        if src:
            index[src.replace("\\", "/")] = e["path"].replace("\\", "/")

    pm = root / "PATH_MAP.json"
    if pm.is_file():
        payload = json.loads(pm.read_text(encoding="utf-8"))
        for k, v in (payload.get("map", payload) or {}).items():
            index[k.replace("\\", "/")] = v.replace("\\", "/")
    return index


def resolve(rel: str) -> Path:
    """Locate a file named by its development-tree path.

    In a development tree this is `project_root() / rel` and nothing else
    happens. In a deposit the manifest decides, and a path the manifest does
    not carry raises `PathNotInDeposit` rather than resolving to nothing.
    """
    root, lay = _root_and_layout()
    rel = str(rel).replace("\\", "/").lstrip("/")
    if lay == "development":
        return root / rel

    index = _reverse_index()
    if rel in index:
        return root / index[rel]

    # A directory, named as the prefix of files the manifest does carry.
    prefix = rel.rstrip("/") + "/"
    mapped = {v for k, v in index.items() if k.startswith(prefix)}
    if mapped:
        parents = {str(Path(m).parent) for m in mapped}
        if len(parents) == 1:
            return root / parents.pop()
        raise PathNotInDeposit(
            f"{rel!r} is spread across {len(parents)} directories in this "
            f"deposit ({', '.join(sorted(parents))}); resolve its files "
            f"individually rather than the directory."
        )

    raise PathNotInDeposit(
        f"{rel!r} is not in this deposit. MANIFEST.json lists what is, and "
        f"README.md lists what is deliberately not -- held-out fault "
        f"definitions, molecular-dynamics trajectories, and four inputs that "
        f"are cited because they are already archived elsewhere."
    )


def exists(rel: str) -> bool:
    """Whether `rel` resolves to something that is here."""
    try:
        return resolve(rel).exists()
    except PathNotInDeposit:
        return False


def under(prefix: str, *parts) -> Path:
    """`resolve()` for a path assembled from pieces.

    The scripts commonly bind a directory once -- `RUNS = ROOT / "20_new_runs"`
    -- and then index into it. A deposit has no `20_new_runs` directory to bind,
    so the binding has to be deferred until the whole relative path is known.
    `under("20_new_runs", name, "scores.csv")` does that.
    """
    tail = "/".join(str(p).replace("\\", "/").strip("/") for p in parts if str(p))
    return resolve(f"{prefix.rstrip('/')}/{tail}" if tail else prefix)


def wsl_root() -> str:
    """The project root as WSL sees it.

    The docking and molecular-dynamics steps shell out to WSL and `cd` into the
    project before running. They hardcoded `/mnt/d/Agentic VS Workflow`, which
    is the author's drive letter and directory name. Deriving it keeps those
    commands correct wherever the tree is checked out, on the understanding
    that WSL sees Windows drives under `/mnt/<letter>`.

    These steps need Vina, DiffDock or GROMACS and cannot run from a deposit in
    any case; this is about not being wrong on another author's machine.
    """
    r = project_root().resolve()
    drive = r.drive.rstrip(":").lower()
    if drive:
        return "/mnt/" + drive + r.as_posix()[len(r.drive):]
    return r.as_posix()


def dev_dir(name: str) -> Path:
    """Bind a development-tree directory without requiring it to exist here.

    A deposit has no `20_new_runs/`; its contents are regrouped by subject
    under `verification/`, so only whole file paths can be mapped. Scripts
    that bind the directory at import time and then index into it with
    `under()` should still import. Those that genuinely walk the directory --
    `build_evidence_index.py`, `build_supplement_s10.py` -- are development
    tools and are expected to fail here, which they now do at the point of
    use, with a path that names the reason.
    """
    root, lay = _root_and_layout()
    if lay == "development":
        return root / name
    return root / "_not_in_this_deposit" / name
