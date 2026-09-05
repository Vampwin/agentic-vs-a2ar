"""build_evidence_index.py -- what the manuscript rests on, and what it does not.

Four lists already describe overlapping parts of this question and none of them
answers it. `sources()` in g8_number_traceability.py knows which artifacts carry
a manuscript number. `ITEMS` in build_supplement_s10.py knows which files reach
the supplement, and carries a `supports` string naming the section. `TREE` in
build_deposit.py knows what is published. `RUN_INDEX.md` carries a human label
per directory. They are hand-maintained, they overlap, and they disagree.

This index derives the used set from the three machine-readable ones rather than
restating them, so it cannot drift from the code -- and reports where they
already contradict each other. It points at the run directories with SHA-256; it
copies nothing, because 20_new_runs/ is immutable under ADR-AIIM-092 and a
fourth copy would be a fourth thing to diverge.

What it cannot derive is the other half: results produced and deliberately not
used. That is `evidence_not_used.py`, hand-maintained with a required reason.
The two halves are held together by a partition assertion over the directories
actually on disk, so a run cannot be silently absent from both.

Usage:
    python build_evidence_index.py --report-only   # worklist, nothing fatal
    python build_evidence_index.py                 # write the index, fatal on E*
    python build_evidence_index.py --check         # regenerate and diff
"""
from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import re
import sys
from pathlib import Path


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
RUNS = ROOT / "20_new_runs"
SCRIPTS = ROOT / "10_new_system/scripts"
RUN_INDEX = RUNS / "RUN_INDEX.md"
DECISION_LOG = ROOT / "00_project_control/08_DECISION_LOG.md"
OUT_MD = ROOT / "00_project_control/28_EVIDENCE_INDEX.md"
OUT_JSON = ROOT / "00_project_control/28_EVIDENCE_INDEX.json"

LABELS = {"CORE", "USED", "INPUT", "PENDING", "ADMIN", "V4-ONLY", "DO-NOT-REUSE"}
NOT_V6 = {"V4-ONLY", "DO-NOT-REUSE"}
# PENDING is exempt from the used/not-used partition: the run is executing or its
# result has not been adjudicated against its pre-registered criteria, and
# forcing it into either list would be a claim about an outcome nobody has. The
# exemption is reported as a standing warning (W8) so it cannot be forgotten.
EXEMPT = {"PENDING"}

# The fourth witness. sources(), ITEMS and TREE describe what is *published*.
# A run consumed by a producer script whose output is published is used, and
# none of the three sees it -- which is how ad_descriptor_sweep_2026-08-15 came
# within one edit of being recorded as superseded while build_ad_comparison.py
# was reading it to build a panel of Figure 4. See ADR-AIIM-142.
#
# This is a regex over source text and it over-reports: a directory a script
# only *writes* to is counted as referenced. That is the safe direction for a
# check whose failure mode is claiming something is unused.
PRODUCER_GLOBS = ("*.py", "*.sh")
PRODUCER_SKIP = {"g8_number_traceability.py", "build_supplement_s10.py",
                 "build_deposit.py", "build_evidence_index.py"}
RUN_DIR_RE = re.compile(r"20_new_runs[/\\]([A-Za-z0-9_.-]+_20\d\d-\d\d-\d\d)")
BARE_NAME_RE = re.compile(r"""["']([A-Za-z0-9_.-]+_20\d\d-\d\d-\d\d)["']""")


def parse_producers(scripts: Path, on_disk: set[str]) -> dict[str, list[str]]:
    """run directory -> producer scripts that name it.

    Two shapes are matched. `20_new_runs/<dir>` is the common one. The other is
    a bare directory name held in a constant and joined later --
    `RUN1 = "snp_screen_repro_run1_2026-08-15"` and then `RUNS / RUN1` -- which
    the path pattern cannot see. `verify_gate_reproducibility.py` is written
    that way and was reported as referencing nothing on the day it was added.
    Bare names over-match by construction, so they are intersected with what is
    actually on disk before being counted.
    """
    hits: dict[str, list[str]] = {}
    for pat in PRODUCER_GLOBS:
        for f in sorted(scripts.glob(pat)):
            if f.name in PRODUCER_SKIP:
                continue
            text = f.read_text(encoding="utf-8", errors="ignore")
            found = set(RUN_DIR_RE.findall(text))
            found |= set(BARE_NAME_RE.findall(text)) & on_disk
            for d in sorted(found):
                hits.setdefault(d, []).append(f.name)
    return hits


def load_module(name: str, path: Path):
    """Import by path without running main(). build_deposit mutates sys.path."""
    saved = list(sys.path)
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    finally:
        sys.path[:] = saved


# ---------------------------------------------------------------- sources()
def _path_from_expr(node: ast.AST) -> str | None:
    """RUNS / "a/b.json"  ->  "a/b.json".  Returns None if not that shape."""
    if not isinstance(node, ast.BinOp) or not isinstance(node.op, ast.Div):
        return None
    left, right = node.left, node.right
    parts: list[str] = []
    while isinstance(left, ast.BinOp) and isinstance(left.op, ast.Div):
        if isinstance(left.right, ast.Constant) and isinstance(left.right.value, str):
            parts.insert(0, left.right.value)
        left = left.left
    if not (isinstance(left, ast.Name) and left.id in {"RUNS", "ROOT"}):
        return None
    prefix = "" if left.id == "RUNS" else "@ROOT@/"
    if isinstance(right, ast.Constant) and isinstance(right.value, str):
        parts.append(right.value)
    elif isinstance(right, ast.JoinedStr):
        head = right.values[0] if right.values else None
        if not (isinstance(head, ast.Constant) and isinstance(head.value, str)):
            return None
        # keep only the literal prefix: f"md_candidates_x/{c}/build.json"
        parts.append(head.value.split("{")[0])
    else:
        return None
    return prefix + "/".join(p.strip("/") for p in parts if p)


def parse_sources(path: Path) -> tuple[dict[str, list[str]], list[str]]:
    """key -> [paths relative to 20_new_runs]. Second value: unresolved exprs."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    fns = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    if "sources" not in fns:
        raise SystemExit("sources() not found -- extractor is stale")

    out: dict[str, list[str]] = {}
    unresolved: list[str] = []

    def collect(expr: ast.AST) -> list[str]:
        p = _path_from_expr(expr)
        if p:
            return [p]
        found: list[str] = []
        for sub in ast.walk(expr):
            if sub is expr:
                continue
            q = _path_from_expr(sub)
            if q:
                found.append(q)
        return found

    # Locals bound before the return -- e8 = pd.read_csv(RUNS / "..."), and so
    # on. The return dict then refers to them by name, so without this the
    # extractor silently loses every key that is not an inline expression. It
    # lost five, including the sweep eight Section 3.5 numbers depend on, and
    # the index then reported that run as published-but-unused. Resolving the
    # names is not optional politeness; it is the difference between an index
    # and a misleading one.
    local: dict[str, list[str]] = {}
    for node in fns["sources"].body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1                 and isinstance(node.targets[0], ast.Name):
            found = collect(node.value)
            if found:
                local[node.targets[0].id] = found

    ret = next((n for n in ast.walk(fns["sources"]) if isinstance(n, ast.Return)), None)
    if ret is None or not isinstance(ret.value, ast.Dict):
        raise SystemExit("sources() no longer returns a dict literal")
    for k, v in zip(ret.value.keys, ret.value.values):
        if not (isinstance(k, ast.Constant) and isinstance(k.value, str)):
            continue
        paths = collect(v)
        if not paths:
            # a bare name, or a name inside a call: resolve through the locals
            for sub in ast.walk(v):
                if isinstance(sub, ast.Name) and sub.id in local:
                    paths.extend(local[sub.id])
            if not paths and any(isinstance(sub, ast.Name) and sub.id in fns
                                 for sub in ast.walk(v)):
                paths = []            # helper call; handled after the loop
        if not paths:
            unresolved.append(f"sources()[{k.value!r}]")
        out[k.value] = paths

    # sources()["md"] is a call to _md_system(), which reads its paths
    # directly. Resolve it here and clear the flag the loop raised, rather
    # than leaving a warning that is already answered two lines later.
    if "_md_system" in fns:
        md = collect(fns["_md_system"])
        if md:
            out.setdefault("md", []).extend(md)
            unresolved = [u for u in unresolved if u != "sources()['md']"]

    # A key whose value is a call to a helper defined in the same file, where
    # that helper reads no run directory at all, resolves to zero paths and that
    # is the correct answer rather than a failure to parse. sources()["ntests"]
    # calls _collect_test_count(), which runs pytest's collector and touches
    # 20_new_runs not at all. The extractor is deliberately strict -- it warns
    # rather than silently skipping a shape it does not know -- so the
    # resolution has to be explicit: recurse into the helper, and only clear the
    # warning if the recursion genuinely finds nothing.
    for key, node in zip(ret.value.keys, ret.value.values):
        if not (isinstance(key, ast.Constant) and isinstance(key.value, str)):
            continue
        if out.get(key.value):
            continue
        called = [sub.id for sub in ast.walk(node)
                  if isinstance(sub, ast.Name) and sub.id in fns]
        if len(called) != 1:
            continue
        # Recurse into the helper and record whatever it reads. Paths outside
        # 20_new_runs come back @ROOT@-prefixed and are skipped downstream, so a
        # helper that touches only the repository -- _collect_test_count() runs
        # pytest from ROOT/10_new_system -- resolves to no run directory, which
        # is the right answer and not a parse failure. The warning exists for
        # shapes the extractor cannot read, and this one it can.
        out[key.value] = collect(fns[called[0]])
        unresolved = [u for u in unresolved if u != f"sources()[{key.value!r}]"]
    return out, unresolved


# ------------------------------------------------------------- RUN_INDEX.md
def parse_run_index(path: Path) -> dict[str, dict]:
    labels: dict[str, dict] = {}
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip().startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 2:
            continue
        dirs = re.findall(r"`([^`]+)`", cells[0])
        label = cells[1].replace("*", "").strip().upper()
        if not dirs or label not in LABELS:
            continue
        for d in dirs:
            labels[d.strip()] = {"label": label, "row_line": n}
    return labels


def prose_mentions(path: Path) -> set[str]:
    return set(re.findall(r"`([a-zA-Z0-9_.\-]+_20\d\d-\d\d-\d\d)`",
                          path.read_text(encoding="utf-8")))


# --------------------------------------------------------------------- main
# --- emitting the index ------------------------------------------------------

# Above this, a file is pointed at but not hashed. 20_new_runs/ holds 6.7 GB in
# 2,362 files and twelve of them account for most of it -- MD trajectories and
# pooled pose archives. Hashing those on every build would make the check too
# slow to run, and the check that never runs is the one that catches nothing.
HASH_SIZE_CAP = 50 * 1024 * 1024

# Rewritten on every G8 run and every figure build, by design (ADR-AIIM-124
# exempts them). They are hashed so the index is complete, and exempted from the
# change detection in --check, which would otherwise fail after every G8 run and
# train its reader to ignore it.
DERIVED_PROVENANCE = {
    "g8_number_traceability.json",
    "figure_ad_multispace_provenance.json",
}

# Exempt by exact path, not by name. `SESSION_LOG_*.md` as a pattern would also
# unfreeze the eight historical logs under the WP directories, which are finished
# records and must stay frozen. This one is the running log of the live
# remediation: it is appended once per session by design, so a hash-change here
# is the file working, not the directory being edited. Same class as ADR-124.
APPEND_ONLY = {
    "manuscript_remediation_2026-08-25/SESSION_LOG_2026-08-25.md",
}


def hash_run(run_dir: Path, dep) -> tuple[list[dict], dict[str, int]]:
    """Point at every file in a run directory, with a hash where practical."""
    files, stats = [], {"hashed": 0, "size": 0, "suffix": 0}
    for p in sorted(run_dir.rglob("*")):
        if not p.is_file():
            continue
        if set(p.relative_to(run_dir).parts[:-1]) & dep.SKIP_DIRS:
            continue
        rel = p.relative_to(run_dir).as_posix()
        size = p.stat().st_size
        if p.suffix in dep.SKIP_SUFFIX:
            stats["suffix"] += 1
            files.append({"path": rel, "bytes": size, "sha256": None,
                          "not_hashed": "suffix excluded from the deposit"})
            continue
        if size > HASH_SIZE_CAP:
            stats["size"] += 1
            files.append({"path": rel, "bytes": size, "sha256": None,
                          "not_hashed": f"over {HASH_SIZE_CAP // (1024 * 1024)} MB"})
            continue
        stats["hashed"] += 1
        entry = {"path": rel, "bytes": size, "sha256": dep.sha256(p)}
        if p.name in DERIVED_PROVENANCE:
            entry["derived_provenance"] = True
        if f"{run_dir.name}/{rel}" in APPEND_ONLY:
            entry["append_only"] = True
        files.append(entry)
    return files, stats


def build_index(on_disk, labels, used, not_used, exempt, E, W) -> dict:
    dep = load_module("build_deposit", SCRIPTS / "build_deposit.py")
    runs, tot = {}, {"files_hashed": 0, "files_skipped_size": 0,
                     "files_skipped_suffix": 0}
    for name in on_disk:
        files, stats = hash_run(RUNS / name, dep)
        tot["files_hashed"] += stats["hashed"]
        tot["files_skipped_size"] += stats["size"]
        tot["files_skipped_suffix"] += stats["suffix"]
        u = used.get(name, {})
        if name in exempt:
            status = "pending"
        elif name in not_used and not_used[name].get("scope",
                                                     "whole_run") == "whole_run":
            status = "not_used"
        else:
            status = "used"
        runs[name] = {
            "label": labels.get(name, {}).get("label"),
            "status": status,
            "referenced_by": {
                "g8_sources": bool(u.get("g8_sources")),
                "s10_items": bool(u.get("s10_items")),
                "deposit_tree": bool(u.get("deposit_tree")),
                "producers": bool(u.get("producers")),
            },
            "sources_keys": sorted(set(u.get("sources_keys", []))),
            "supports": sorted(set(u.get("supports", []))),
            "producer_scripts": sorted(set(u.get("producer_scripts", []))),
            "not_used": not_used.get(name),
            "n_files": len(files),
            "files": files,
        }
    tot["runs"] = len(runs)
    tot["used"] = sum(1 for r in runs.values() if r["status"] == "used")
    tot["not_used"] = sum(1 for r in runs.values() if r["status"] == "not_used")
    tot["pending"] = sum(1 for r in runs.values() if r["status"] == "pending")
    return {
        "generated_by": "10_new_system/scripts/build_evidence_index.py",
        "authority": "ADR-AIIM-142",
        "what_this_is": "A pointer index. It copies nothing. Every run "
                        "directory under 20_new_runs/ appears exactly once, "
                        "classified as used, not used, or pending, with the "
                        "classification derived from four machineries rather "
                        "than retyped -- sources() in g8_number_traceability.py, "
                        "ITEMS in build_supplement_s10.py, TREE in "
                        "build_deposit.py, and a scan of the producer scripts.",
        "totals": tot,
        "warnings": [f"{c}: {m}" for c, m in W],
        "runs": runs,
    }


def render_markdown(index: dict) -> str:
    t = index["totals"]
    L = ["# Evidence index — every run this project produced, and whether the "
         "manuscript uses it", "",
         f"**Generated by `{index['generated_by']}` under {index['authority']}. "
         f"Do not edit by hand.**", "",
         index["what_this_is"], "",
         f"{t['runs']} run directories: **{t['used']} used**, "
         f"**{t['not_used']} not used**, **{t['pending']} pending**. "
         f"{t['files_hashed']} files hashed; {t['files_skipped_size']} are "
         f"pointed at but not hashed for size, {t['files_skipped_suffix']} by "
         f"suffix.", "",
         "## Used", "",
         "| Run | Label | g8 | S10 | deposit | producers | files |",
         "|---|---|:-:|:-:|:-:|:-:|---:|"]
    m = lambda b: "•" if b else ""
    for name, r in sorted(index["runs"].items()):
        if r["status"] != "used":
            continue
        k = r["referenced_by"]
        L.append(f"| `{name}` | {r['label']} | {m(k['g8_sources'])} | "
                 f"{m(k['s10_items'])} | {m(k['deposit_tree'])} | "
                 f"{m(k['producers'])} | {r['n_files']} |")

    L += ["", "### Which manuscript numbers depend on which run", "",
          "`sources()` keys are the check names `g8_number_traceability.py` "
          "already uses. They are the project's existing vocabulary for "
          "manuscript claims, not a second one invented here.", "",
          "| Run | sources() keys |", "|---|---|"]
    for name, r in sorted(index["runs"].items()):
        if r["sources_keys"]:
            L.append(f"| `{name}` | {', '.join('`' + k + '`' for k in r['sources_keys'])} |")

    L += ["", "## Pending", "",
          "Executing, or the result is not yet adjudicated against its "
          "pre-registered criteria. Exempt from the partition until it is.", ""]
    for name, r in sorted(index["runs"].items()):
        if r["status"] == "pending":
            L.append(f"- `{name}` — {r['n_files']} files")

    L += ["", "## Not used", "",
          "Produced and deliberately left out. Every entry carries a kind, a "
          "prose reason and the deciding ADR, and the reason is printed here in "
          "full: code can force a reason to exist and an ADR to resolve, but it "
          "cannot force the reason to be true. Publication is the only check on "
          "that, so nothing here is summarised.", ""]
    for name, r in sorted(index["runs"].items()):
        nu = r["not_used"]
        if not nu or r["status"] != "not_used":
            continue
        L.append(f"### `{name}`")
        L.append("")
        L.append(f"**{nu['kind']}**"
                 + (f" → `{nu['superseded_by']}`" if nu.get("superseded_by") else "")
                 + f" · {', '.join(nu.get('adr', []))} · decided {nu.get('decided','')}")
        L.append("")
        L.append(nu["reason"])
        L.append("")

    L += ["## Reported, but the material is withheld", "",
          "These are used — a number from each is in the manuscript — and their "
          "raw material is not published. The distinction is recorded rather "
          "than left to a reader to notice.", ""]
    for name, r in sorted(index["runs"].items()):
        nu = r["not_used"]
        if not nu or nu.get("scope") != "artifacts_only":
            continue
        L.append(f"### `{name}`")
        L.append("")
        L.append(f"**{nu['kind']}** · excluded from "
                 f"{', '.join(nu.get('excluded_from', []))} · "
                 f"{', '.join(nu.get('adr', []))}")
        L.append("")
        L.append(nu["reason"])
        if nu.get("published_anyway"):
            L.append("")
            L.append("Published anyway: "
                     + ", ".join(f"`{p}`" for p in nu["published_anyway"]))
        L.append("")

    if index["warnings"]:
        L += ["## Open warnings", "",
              "Not fatal, and not hidden.", ""]
        L += [f"- {w}" for w in index["warnings"]]
        L.append("")
    return "\n".join(L) + "\n"


def run_check(index: dict) -> int:
    """Compare against the committed index. Hash drift is the thing to catch."""
    if not OUT_JSON.is_file():
        print("no committed index to check against; run without --check first")
        return 1
    old = json.loads(OUT_JSON.read_text(encoding="utf-8"))
    problems = []
    for name, r in index["runs"].items():
        o = old["runs"].get(name)
        if o is None:
            problems.append(f"{name}: new run directory, not in the committed index")
            continue
        if o["status"] != r["status"]:
            problems.append(f"{name}: status {o['status']} -> {r['status']}")
        if o["label"] != r["label"]:
            problems.append(f"{name}: label {o['label']} -> {r['label']}")
        oldf = {f["path"]: f for f in o["files"]}
        for f in r["files"]:
            g = oldf.pop(f["path"], None)
            if g is None:
                problems.append(f"{name}/{f['path']}: new file")
            elif (g.get("sha256") != f.get("sha256")
                  and not f.get("derived_provenance")
                  and not f.get("append_only")):
                problems.append(f"{name}/{f['path']}: content changed — "
                                f"20_new_runs is immutable (ADR-AIIM-092)")
        for path in oldf:
            problems.append(f"{name}/{path}: file has disappeared")
    for name in old["runs"]:
        if name not in index["runs"]:
            problems.append(f"{name}: run directory has disappeared")

    if problems:
        print(f"\n--- CHECK FAILED ({len(problems)}) " + "-" * 40)
        for p in problems:
            print(f"    {p}")
        return 1
    print("\ncheck clean: every run, label, status and hash matches the "
          "committed index")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-only", action="store_true",
                    help="print the worklist; downgrade every fatal to a warning "
                         "and write nothing")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    nu_mod = load_module("evidence_not_used", SCRIPTS / "evidence_not_used.py")
    nu_mod.self_test()
    NOT_USED = nu_mod.NOT_USED

    s10 = load_module("build_supplement_s10", SCRIPTS / "build_supplement_s10.py")
    dep = load_module("build_deposit", SCRIPTS / "build_deposit.py")
    src_keys, unresolved = parse_sources(SCRIPTS / "g8_number_traceability.py")

    # run directory -> which machinery references it
    used: dict[str, dict] = {}

    def touch(run: str, key: str, detail: str = "") -> None:
        e = used.setdefault(run, {"g8_sources": False, "s10_items": False,
                                  "deposit_tree": False, "producers": False,
                                  "sources_keys": [], "supports": [],
                                  "producer_scripts": []})
        e[key] = True
        if key == "g8_sources" and detail:
            e["sources_keys"].append(detail)
        if key == "s10_items" and detail:
            e["supports"].append(detail)
        if key == "producers" and detail:
            e["producer_scripts"].append(detail)

    for k, paths in src_keys.items():
        for p in paths:
            if p.startswith("@ROOT@/"):
                continue
            touch(p.split("/")[0], "g8_sources", k)

    for rel, _sub, supports in s10.ITEMS:
        resolved = (RUNS / rel).resolve()
        try:
            inside = resolved.relative_to(RUNS.resolve())
        except ValueError:
            continue                      # ../30_manuscript/figures/... escapes
        touch(inside.parts[0], "s10_items", supports)

    for src_rel, _dst, _desc in dep.TREE:
        if src_rel.startswith("20_new_runs/"):
            touch(src_rel.split("/", 1)[1].split("/")[0], "deposit_tree")

    on_disk = sorted(p.name for p in RUNS.iterdir() if p.is_dir())

    for run, scripts_naming in parse_producers(SCRIPTS, set(on_disk)).items():
        for s in scripts_naming:
            touch(run, "producers", s)

    labels = parse_run_index(RUN_INDEX)
    mentioned = prose_mentions(RUN_INDEX)
    adr_ids = set(re.findall(r"^#{2,3} (ADR-AIIM-\d+)",
                             DECISION_LOG.read_text(encoding="utf-8"), re.M))

    E: list[tuple[str, str]] = []
    W: list[tuple[str, str]] = []

    exempt = {r for r, l in labels.items() if l["label"] in EXEMPT}

    for run in sorted(used):
        if run in exempt:
            continue
        if run not in labels:
            how = "prose mention only" if run in mentioned else "not named at all"
            E.append(("E1", f"{run}: read by the machinery but no status row in "
                            f"RUN_INDEX.md ({how})"))
        elif labels[run]["label"] in NOT_V6 and (used[run]["g8_sources"]
                                                 or used[run]["s10_items"]):
            E.append(("E2", f"{run}: labelled {labels[run]['label']} yet feeds a v6 "
                            f"number or the supplement "
                            f"(RUN_INDEX.md:{labels[run]['row_line']})"))

    for run in sorted(exempt):
        W.append(("W8", f"{run}: PENDING -- exempt from the partition until its "
                        f"result is adjudicated (RUN_INDEX.md:"
                        f"{labels[run]['row_line']})"))

    for run in on_disk:
        if run in exempt:
            continue
        if run not in used and run not in NOT_USED:
            E.append(("E3", f"{run}: on disk, in neither the derived-used set nor "
                            f"evidence_not_used.py"))
    for run in NOT_USED:
        if run not in on_disk:
            E.append(("E6", f"{run}: declared not-used but not on disk"))
        for a in NOT_USED[run].get("adr", []):
            if a not in adr_ids:
                E.append(("E5", f"{run}: cites {a}, which has no heading in the "
                                f"decision log"))
    missing = [f"{run}/{p.split('/', 1)[1] if '/' in p else ''}"
               for run, e in used.items() for p in []]
    for k, paths in src_keys.items():
        for p in paths:
            if p.startswith("@ROOT@/"):
                continue
            if not (RUNS / p).exists():
                E.append(("E7", f"sources()[{k!r}] -> {p} does not exist"))

    for run, lab in labels.items():
        if run not in used and lab["label"] in {"CORE", "USED"}:
            W.append(("W1", f"{run}: labelled {lab['label']} but no machinery "
                            f"references it"))
    for run, e in used.items():
        if e["g8_sources"] and not (e["s10_items"] or e["deposit_tree"]):
            W.append(("W2", f"{run}: carries a manuscript number but is published "
                            f"in neither the supplement nor the deposit"))
        if not e["g8_sources"] and (e["s10_items"] or e["deposit_tree"]):
            W.append(("W7", f"{run}: published but no reported number depends on it"))
    for u in unresolved:
        W.append(("W5", f"{u}: path expression not recognised by the extractor"))
    # README.md and RUN_INDEX.md belong there; stray files do not. The two
    # AD-sweep logs are the exception, and are named rather than tolerated: they
    # are the stdout of the runs whose directories sit beside them, written
    # there by how those runs were invoked. They carry no claim, `.log` is in
    # build_deposit's SKIP_SUFFIX, and moving them into the run directories
    # would modify artifacts ADR-AIIM-092 declares immutable. A warning that
    # will never be resolved teaches its reader to ignore warnings, so this one
    # is closed with a reason instead of left standing.
    TOP_LEVEL_ALLOWED = {
        "README.md", "RUN_INDEX.md",
        "ad_sweep_expanded_2026-08-21.log", "ad_sweep_selfsim_2026-08-25.log",
    }
    for f in RUNS.iterdir():
        if f.is_file() and f.name not in TOP_LEVEL_ALLOWED:
            W.append(("W4", f"{f.name}: loose file at the top of 20_new_runs/"))

    print("=" * 74)
    print("EVIDENCE INDEX -- report only, nothing written")
    print("=" * 74)
    print(f"run directories on disk       : {len(on_disk)}")
    print(f"referenced by some machinery  : {len(used)}")
    print(f"  by g8 sources()             : {sum(1 for e in used.values() if e['g8_sources'])}")
    print(f"  by S10 ITEMS                : {sum(1 for e in used.values() if e['s10_items'])}")
    print(f"  by deposit TREE             : {sum(1 for e in used.values() if e['deposit_tree'])}")
    print(f"  by producer scripts         : {sum(1 for e in used.values() if e['producers'])}")
    print(f"labelled in RUN_INDEX.md      : {len(labels)}")
    print(f"declared in evidence_not_used : {len(NOT_USED)}")
    print()
    for code in ("E1", "E2", "E3", "E5", "E6", "E7"):
        rows = [m for c, m in E if c == code]
        if rows:
            print(f"--- {code}  ({len(rows)}) " + "-" * (56 - len(code)))
            for m in rows:
                print(f"    {m}")
    for code in ("W1", "W2", "W4", "W5", "W7", "W8"):
        rows = [m for c, m in W if c == code]
        if rows:
            print(f"--- {code}  ({len(rows)}) " + "-" * (56 - len(code)))
            for m in rows:
                print(f"    {m}")
    print()
    print(f"{len(E)} would-be fatal, {len(W)} warnings")

    if args.report_only:
        return 0

    if E:
        print("\nrefusing to write an index over an unresolved contradiction.")
        print("fix the fatals above, or re-run with --report-only.")
        return 1

    index = build_index(on_disk, labels, used, NOT_USED, exempt, E, W)

    if args.check:
        return run_check(index)

    OUT_JSON.write_text(json.dumps(index, indent=2), encoding="utf-8")
    OUT_MD.write_text(render_markdown(index), encoding="utf-8")
    tot = index["totals"]
    print(f"\nwritten -> {OUT_MD.relative_to(ROOT).as_posix()}")
    print(f"           {OUT_JSON.relative_to(ROOT).as_posix()}")
    print(f"           {tot['runs']} runs, {tot['files_hashed']} files hashed, "
          f"{tot['files_skipped_size']} skipped for size, "
          f"{tot['files_skipped_suffix']} by suffix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
