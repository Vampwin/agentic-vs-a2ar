"""analyse_branch_coverage.py -- what the 75 % coverage figure hides.

The external reproducibility audit of 13 September 2026 (ADR-AIIM-189) shipped a
full `coverage.py` 7.16.0 branch-aware report, `coverage_core.json`, and the
manuscript quotes three aggregates from it: 75 % over the core packages, and
29 % and 25 % for the contract engine and the orchestrator.

Those two module figures are *statement* percentages, and they are doing
something the sentence does not say. This script reads the same immutable
artifact and reports the branch figures underneath them.

Nothing is recomputed and no test is run: every number below is read out of the
audit artifact, which is immutable under ADR-AIIM-092. The point of the script
is that the audit already contained this and it was reported as two percentages.

Three questions, answered from the artifact:

  1. Which core modules have *zero* executed branches?
  2. For the two load-bearing ones, what do the uncovered branches decide?
  3. Is there any test that imports them at all?

Question 3 is answered by parsing the import graph of the test suite with `ast`,
not by grep, so that `from ... import` and aliased forms are caught too.

Output -> 20_new_runs/coverage_branch_analysis_2026-09-21/
"""
import ast
import io
import json
import os
import sys
from datetime import datetime, timezone

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
SYS = os.path.dirname(HERE)
ROOT = os.path.dirname(SYS)

AUDIT = os.path.join(ROOT, "20_new_runs", "external_repro_audit_2026-09-13",
                     "coverage_core.json")
TESTS = os.path.join(SYS, "tests")
OUT = os.path.join(ROOT, "20_new_runs", "coverage_branch_analysis_2026-09-21")

# The two modules the manuscript names. Keys as coverage.py recorded them.
FOCUS = {
    "kdss_audit\\contract_engine.py": "kdss_audit/contract_engine.py",
    "kdss_audit\\orchestrator.py": "kdss_audit/orchestrator.py",
}

# A branch is fail-closed-relevant if its guarded lines mention any of these.
FAILCLOSED = ("HALT", "ESCALATE", "_halted", "_escalation_pending", "raise ",
              "AUTO_REPAIR")


def load_audit():
    with io.open(AUDIT, encoding="utf-8") as f:
        return json.load(f)


def module_rows(files):
    rows = []
    for key, v in files.items():
        s = v["summary"]
        if s["num_statements"] == 0:
            continue
        rows.append({
            "module": key.replace("\\", "/"),
            "statements": s["num_statements"],
            "percent_statements": round(s["percent_covered"], 1),
            "branches": s["num_branches"],
            "covered_branches": s["covered_branches"],
            "missing_branches": s["missing_branches"],
            "percent_branches": (round(s["percent_branches_covered"], 1)
                                 if s["num_branches"] else None),
        })
    rows.sort(key=lambda r: (r["percent_branches"] if r["percent_branches"]
                             is not None else 999, r["percent_statements"]))
    return rows


def source_lines(rel):
    p = os.path.join(SYS, rel)
    with io.open(p, encoding="utf-8") as f:
        return f.read().split("\n")


def classify(rel, missing):
    """Describe each uncovered branch by the line that decides it."""
    src = source_lines(rel)
    seen, out = set(), []
    for pair in missing:
        ln = pair[0]
        if ln in seen or ln < 1 or ln > len(src):
            continue
        seen.add(ln)
        decide = src[ln - 1].strip()
        # look at the guarded body: up to 6 lines, stopping at a dedent
        body, indent = [], len(src[ln - 1]) - len(src[ln - 1].lstrip())
        for i in range(ln, min(len(src), ln + 6)):
            line = src[i]
            if line.strip() and (len(line) - len(line.lstrip())) <= indent:
                break
            if line.strip():
                body.append(line.strip())
        blob = decide + " " + " ".join(body)
        out.append({
            "line": ln,
            "decides": decide,
            "guards": body[:3],
            "fail_closed_relevant": any(t in blob for t in FAILCLOSED),
        })
    return sorted(out, key=lambda d: d["line"])


def scan_imports(root):
    """Every module name imported anywhere under `root`, via AST.

    Returns (names, n_files, by_name) so a caller can ask not only whether a
    module is imported but from where.
    """
    names, files, by_name = set(), 0, {}
    for dirpath, _dirs, fnames in os.walk(root):
        if "__pycache__" in dirpath:
            continue
        for fn in fnames:
            if not fn.endswith(".py"):
                continue
            files += 1
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, ROOT).replace("\\", "/")
            try:
                tree = ast.parse(io.open(p, encoding="utf-8",
                                         errors="replace").read(), filename=p)
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                got = []
                if isinstance(node, ast.Import):
                    got = [a.name for a in node.names]
                elif isinstance(node, ast.ImportFrom):
                    # Resolve `from .contract_engine import X` against the
                    # importing file's own package, or the relative import is
                    # invisible and a caller concludes nothing reaches it.
                    base = node.module or ""
                    if node.level:
                        pkg = os.path.dirname(os.path.relpath(p, SYS))
                        parts = pkg.replace("\\", "/").split("/")
                        parts = parts[:len(parts) - (node.level - 1)] if \
                            node.level > 1 else parts
                        prefix = ".".join([x for x in parts if x])
                        base = f"{prefix}.{base}" if base else prefix
                    if base:
                        got = [base] + [f"{base}.{a.name}" for a in node.names]
                for g in got:
                    names.add(g)
                    by_name.setdefault(g, set()).add(rel)
    return names, files, {k: sorted(v) for k, v in by_name.items()}


def test_imports():
    """Every module name imported anywhere under tests/, via AST."""
    names, files = set(), 0
    for dirpath, _dirs, fnames in os.walk(TESTS):
        if "__pycache__" in dirpath:
            continue
        for fn in fnames:
            if not fn.endswith(".py"):
                continue
            files += 1
            p = os.path.join(dirpath, fn)
            try:
                tree = ast.parse(io.open(p, encoding="utf-8").read(), filename=p)
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for a in node.names:
                        names.add(a.name)
                elif isinstance(node, ast.ImportFrom):
                    if node.module:
                        names.add(node.module)
                        for a in node.names:
                            names.add(f"{node.module}.{a.name}")
    return names, files


def main():
    audit = load_audit()
    rows = module_rows(audit["files"])
    imported, n_test_files = test_imports()

    # Where, in the whole system, is each module actually reached from? The
    # question the coverage number raises is not "is it tested" alone but
    # "is the untested thing on the path the study itself ran".
    all_names, n_sys_files, by_name = scan_imports(SYS)

    def importers(rel):
        dotted = rel[:-3].replace("/", ".")
        hits = set()
        for name, files in by_name.items():
            if name == dotted or name.startswith(dotted + "."):
                hits |= set(files)
        out = {"tests": [], "benchmark": [], "production": [], "other": []}
        for f in sorted(hits):
            if "/tests/" in f:
                out["tests"].append(f)
            elif "/fault_injection/" in f:
                out["benchmark"].append(f)
            elif "/scripts/" in f:
                out["production"].append(f)
            else:
                out["other"].append(f)
        return out

    zero = [r for r in rows if r["branches"] > 0 and r["covered_branches"] == 0]

    focus = {}
    for key, rel in FOCUS.items():
        v = audit["files"][key]
        s = v["summary"]
        br = classify(rel, v["missing_branches"])
        hits = [n for n in imported
                if n == rel[:-3].replace("/", ".") or n.startswith(
                    rel[:-3].replace("/", ".") + ".")]
        focus[rel] = {
            "percent_statements_reported_in_manuscript": round(s["percent_covered"], 1),
            "num_branches": s["num_branches"],
            "covered_branches": s["covered_branches"],
            "percent_branches": round(s["percent_branches_covered"], 1),
            "executed_branches": v["executed_branches"],
            "executed_lines": v["executed_lines"],
            "n_uncovered_decision_points": len(br),
            "n_fail_closed_relevant": sum(1 for b in br
                                          if b["fail_closed_relevant"]),
            "imported_by_any_test": bool(hits),
            "importing_names": sorted(hits),
            "reached_from": importers(rel),
            "uncovered_decision_points": br,
        }

    total_br = sum(r["branches"] for r in rows)
    total_cov = sum(r["covered_branches"] for r in rows)

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_artifact": "20_new_runs/external_repro_audit_2026-09-13/coverage_core.json",
        "source_is_immutable": True,
        "recomputed_anything": False,
        "coverage_tool": f"coverage.py {audit['meta']['version']}",
        "branch_coverage_enabled": audit["meta"]["branch_coverage"],
        "n_modules_with_statements": len(rows),
        "core_branches_total": total_br,
        "core_branches_covered": total_cov,
        "core_percent_branches": round(100.0 * total_cov / total_br, 1) if total_br else None,
        "n_test_files_scanned": n_test_files,
        "n_system_files_scanned": n_sys_files,
        "modules_with_zero_executed_branches": [
            {k: r[k] for k in ("module", "statements", "percent_statements",
                               "branches", "covered_branches")} for r in zero],
        "focus_modules": focus,
        "all_modules": rows,
    }

    os.makedirs(OUT, exist_ok=True)
    with io.open(os.path.join(OUT, "branch_coverage_analysis.json"), "w",
                 encoding="utf-8", newline="\n") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
        f.write("\n")

    # ------------------------------------------------------------- console
    print("=" * 74)
    print("BRANCH COVERAGE UNDERNEATH THE AUDITED PERCENTAGES")
    print("=" * 74)
    print(f"source     : {summary['source_artifact']}")
    print(f"tool       : {summary['coverage_tool']}  branch={summary['branch_coverage_enabled']}")
    print(f"core branch coverage: {total_cov} / {total_br} = "
          f"{summary['core_percent_branches']} %")
    print()
    print(f"Modules with statements but ZERO executed branches: {len(zero)}")
    for r in zero:
        print(f"  {r['module']:48} stmt {r['percent_statements']:5.1f} %  "
              f"branches {r['covered_branches']}/{r['branches']}")
    print()
    for rel, d in focus.items():
        print("-" * 74)
        print(rel)
        print(f"  manuscript quotes statement coverage : "
              f"{d['percent_statements_reported_in_manuscript']} %")
        print(f"  branch coverage                      : "
              f"{d['covered_branches']}/{d['num_branches']} = {d['percent_branches']} %")
        print(f"  imported by any of {n_test_files} test files      : "
              f"{d['imported_by_any_test']}")
        rf = d["reached_from"]
        print(f"  reached from tests                   : {rf['tests'] or 'NONE'}")
        print(f"  reached from benchmark harness       : {rf['benchmark'] or 'NONE'}")
        print(f"  reached from production scripts      : {rf['production'] or 'NONE'}")
        print(f"  uncovered decision points            : "
              f"{d['n_uncovered_decision_points']} "
              f"({d['n_fail_closed_relevant']} fail-closed-relevant)")
        for b in d["uncovered_decision_points"]:
            flag = "FAIL-CLOSED" if b["fail_closed_relevant"] else "           "
            print(f"    L{b['line']:<4} [{flag}] {b['decides'][:52]}")
    print()
    print(f"written -> {OUT}")


if __name__ == "__main__":
    main()
