"""ADR-AIIM-412: functions shared by the jobs. The registered censored-EM is mirrored pass by pass (same forest, same projection, same
tolerance rule); the ridge-Tobit model and the simulation are new. Nothing in 10_new_system is modified."""
from __future__ import annotations

import importlib.util
import os
import sys
import time
from collections import Counter
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(r"D:\Agentic VS Workflow")
RUNS = ROOT / "20_new_runs"
OUT = RUNS / "v9_censored_sensitivity_2026-10-12"
SCRIPTS = ROOT / "10_new_system" / "scripts"
sys.path.insert(0, str(SCRIPTS))
RS = 42
TOP_ABS, MAJORITY, SHORTLIST = 20, 5, 10        # pool_comparison.py, as imported by em_survivors.py
LAMBDAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)

_cache = {}


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def em_module():
    if "em" not in _cache:
        _cache["em"] = _load("_em", RUNS / "pki_censored_em_2026-10-05" / "analysis" / "pki_censored_em.py")
    return _cache["em"]


def data():
    if "data" not in _cache:
        z = np.load(OUT / "cache.npz", allow_pickle=True)
        _cache["data"] = {k: z[k] for k in z.files}
    return _cache["data"]


def groups():
    """The five distinct splits: the first member of each group of conventions that share a split, and the whole group."""
    import json
    g = json.loads((RUNS / "em_survivors_2026-10-06" / "em_survivors.json").read_text("utf-8"))["distinct_split_groups"]
    return g


def split_for(name):
    em = em_module()
    cfg = next(c for c in em.CELLS if em.cell_name(c) == name)
    d = data()
    return em.split(list(d["scaf"]), len(d["y0"]), **cfg)


# ---------------------------------------------------------------------------------------------- forest and the registered EM
def pipe(n_jobs=1):
    return Pipeline([("var", VarianceThreshold(threshold=0.01)),
                     ("rf", RandomForestRegressor(random_state=RS, n_jobs=n_jobs))])


def start_values(y0, kind, tr, init, seed=42):
    y = y0.copy()
    up, lo = tr[kind[tr] == 1], tr[kind[tr] == -1]
    pts = tr[kind[tr] == 0]
    if init == "bound":
        pass
    elif init == "inside1":
        y[up] = y0[up] - 1.0
        y[lo] = y0[lo] + 1.0
    elif init == "random":
        rng = np.random.default_rng(seed)
        y[up] = y0[up] - rng.uniform(0, 1.5, len(up))
        y[lo] = y0[lo] + rng.uniform(0, 1.5, len(lo))
    elif init == "mean":
        m = float(y0[pts].mean())
        y[up] = np.minimum(m, y0[up])
        y[lo] = np.maximum(m, y0[lo])
    else:
        raise ValueError(init)
    return y


def em_passes(X, y0, kind, tr, XQ, init="bound", max_iter=10, tol=0.01, n_jobs=1):
    """The registered censored-EM (em.censored_em), pass by pass. Per pass: the model fitted on the current targets, its predictions
    on XQ, then the projection into the interval and the change of the imputed targets. As in the registered code, when the change
    falls under the tolerance the model is refitted on the updated targets and that model is the result."""
    y = start_values(y0, kind, tr, init)
    up, lo = tr[kind[tr] == 1], tr[kind[tr] == -1]
    recs, final = [], None
    for it in range(1, max_iter + 1):
        mdl = pipe(n_jobs).fit(X[tr], y[tr])
        pq = mdl.predict(XQ)
        new = y.copy()
        new[up] = np.minimum(mdl.predict(X[up]), y0[up])
        new[lo] = np.maximum(mdl.predict(X[lo]), y0[lo])
        d = np.abs(new[tr] - y[tr])
        recs.append({"pass": it, "max_change": float(d.max()), "mean_change": float(d.mean()),
                     "n_moved_gt_0.01": int((d > 0.01).sum()), "pred": pq})
        y = new
        final = pq
        if d.max() < tol:
            mdl = pipe(n_jobs).fit(X[tr], y[tr])
            final = mdl.predict(XQ)
            break
    return {"recs": recs, "final_pred": final, "n_pass": len(recs), "stopped_by_tol": bool(recs[-1]["max_change"] < tol),
            "final_change": recs[-1]["max_change"], "y_imputed": y}


# ---------------------------------------------------------------------------------------------- the selection rule (em_survivors.py)
def select(preds_by_conv, uids, murcko):
    """preds_by_conv: eight arrays in the order of em.CELLS. Returns the survivors in order and the shortlist."""
    cnt = Counter()
    for p in preds_by_conv:
        top = np.argsort(-p)[:TOP_ABS]
        cnt.update([uids[i] for i in top])
    mean = dict(zip(uids, np.vstack(preds_by_conv).mean(axis=0)))
    surv = {k: v for k, v in cnt.items() if v >= MAJORITY}
    ranked = sorted(surv, key=lambda c: (-surv[c], -mean[c]))
    fw = dict(zip(uids, murcko))
    shortlist, seen = [], {}
    for c in ranked:
        if fw[c] in seen:
            continue
        seen[fw[c]] = c
        shortlist.append(c)
        if len(shortlist) == SHORTLIST:
            break
    return ranked, shortlist


def expand(pred_by_group):
    """pred_by_group: {first_member_name: array}; returns the eight arrays in the order of em.CELLS."""
    em = em_module()
    gs = groups()
    out = []
    for c in em.CELLS:
        nm = em.cell_name(c)
        g = next(g for g in gs if nm in g)
        out.append(pred_by_group[g[0]])
    return out


# ---------------------------------------------------------------------------------------------- ridge-Tobit
def _nll_grad(theta, X, y, kind, lam):
    n, d = X.shape
    w, b, s = theta[:d], theta[d], theta[d + 1]
    sig = np.exp(s)
    mu = X @ w + b
    r = (y - mu) / sig
    P, U, L = kind == 0, kind == 1, kind == -1
    dmu = np.zeros(n)
    nll = 0.5 * np.sum(r[P] ** 2) + P.sum() * s
    ds = np.sum(1.0 - r[P] ** 2)
    dmu[P] = -r[P] / sig
    if U.any():                                   # y <= c: probability Phi(z), z = (c - mu) / sigma
        z = r[U]
        lg = log_ndtr(z)
        ratio = np.exp(-0.5 * z * z - 0.5 * np.log(2 * np.pi) - lg)
        nll += -lg.sum()
        dmu[U] = ratio / sig
        ds += np.sum(ratio * z)
    if L.any():                                   # y >= c: probability Phi(-z)
        z = r[L]
        q = -z
        lg = log_ndtr(q)
        ratio = np.exp(-0.5 * q * q - 0.5 * np.log(2 * np.pi) - lg)
        nll += -lg.sum()
        dmu[L] = -ratio / sig
        ds += -np.sum(ratio * z)
    nll_pen = nll + 0.5 * lam * float(w @ w)
    g = np.empty_like(theta)
    g[:d] = X.T @ dmu + lam * w
    g[d] = dmu.sum()
    g[d + 1] = ds
    return nll_pen, g


def tobit_nll(theta, X, y, kind):
    """Unpenalised negative log-likelihood (used for the held-out score)."""
    return _nll_grad(theta, X, y, kind, 0.0)[0]


def tobit_fit(X, y, kind, lam, maxiter=500):
    d = X.shape[1]
    pts = kind == 0
    theta0 = np.zeros(d + 2)
    theta0[d] = float(y[pts].mean())
    theta0[d + 1] = float(np.log(y[pts].std()))
    res = minimize(_nll_grad, theta0, args=(X, y, kind, lam), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxfun": 3 * maxiter})
    return res.x, res


class RidgeTobit:
    """Pipeline-like wrapper: VarianceThreshold(0.01) and StandardScaler fitted on the training rows, penalty chosen by 3-fold GroupKFold
    on the training scaffolds (held-out censored likelihood), then one fit on all training rows."""

    def __init__(self, lambdas=LAMBDAS, maxiter=500):
        self.lambdas, self.maxiter = lambdas, maxiter

    def _prep_fit(self, Xtr):
        self.var = VarianceThreshold(threshold=0.01).fit(Xtr)
        Z = self.var.transform(Xtr)
        self.sc = StandardScaler().fit(Z)
        return self.sc.transform(Z)

    def _prep(self, X):
        return self.sc.transform(self.var.transform(X))

    def fit(self, Xtr, ytr, ktr, group):
        Z = self._prep_fit(Xtr)
        scores = {}
        for lam in self.lambdas:
            tot = 0.0
            for a, bb in GroupKFold(n_splits=3).split(Z, ytr, group):
                th, _ = tobit_fit(Z[a], ytr[a], ktr[a], lam, self.maxiter)
                tot += tobit_nll(th, Z[bb], ytr[bb], ktr[bb])
            scores[lam] = tot
        self.cv = scores
        self.lam = min(scores, key=scores.get)
        self.theta, self.res = tobit_fit(Z, ytr, ktr, self.lam, self.maxiter)
        self.d = Z.shape[1]
        return self

    def predict(self, X):
        Z = self._prep(X)
        return Z @ self.theta[:self.d] + self.theta[self.d]

    @property
    def sigma(self):
        return float(np.exp(self.theta[self.d + 1]))


def check_gradient():
    """Finite-difference check of the analytic gradient on a small random problem (all three record types)."""
    rng = np.random.default_rng(0)
    n, d = 60, 7
    X = rng.normal(size=(n, d))
    y = X @ rng.normal(size=d) + rng.normal(size=n) * 0.5 + 6
    kind = np.zeros(n, dtype=int)
    kind[:10] = 1
    kind[10:18] = -1
    y[:10] = y[:10] - 0.3
    theta = np.concatenate([rng.normal(size=d) * 0.3, [6.0, -0.4]])
    f0, g = _nll_grad(theta, X, y, kind, 3.0)
    num = np.zeros_like(theta)
    for i in range(len(theta)):
        e = np.zeros_like(theta)
        e[i] = 1e-6
        num[i] = (_nll_grad(theta + e, X, y, kind, 3.0)[0] - _nll_grad(theta - e, X, y, kind, 3.0)[0]) / 2e-6
    return float(np.max(np.abs(num - g)) / max(1.0, np.max(np.abs(g))))


# ---------------------------------------------------------------------------------------------- scores
def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def sp(a, b):
    return float(spearmanr(a, b).statistic)
