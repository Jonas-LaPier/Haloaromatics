"""Multiple-linear-regression QSAR exploration with train/test splits.

For each level, every combination of 1-3 descriptors (config.QSAR_DESCRIPTORS with values for
all compounds that have rate constants) is fitted as ln kobs = b0 + sum_i b_i x_i and judged on
data it was not fitted on:

    repeated stratified splits  MLR_SPLITS random splits into training and test sets that keep
                                MLR_TEST_PER_GROUP compounds of each group in the test set
                                (fixed seed MLR_SEED). Pooled over all splits:
                                  RMSE_test = sqrt(mean of (y - y_hat)^2 over every test prediction)
                                  Q2_F1     = 1 - sum (y - y_hat)^2 / sum (y - mean(y_train))^2
                                  R2_train  = mean training-set R2
                                  RMSE_test_sd = spread of the per-split test RMSE
    transfer split              train on CORR_FIT_GROUP (bromobenzenes), predict all other
                                compounds: RMSE_transfer
    all-data fit                coefficients, standardised coefficients, R2, adj R2, LOO Q2
    collinearity                max |r| between the descriptors of a model and max VIF;
                                models with max |r| > MLR_MAX_R are flagged as collinear

Models are ranked by RMSE_test. A mean-only model (predicts the training mean) is the baseline.
"""
from __future__ import annotations

import itertools
import math

import numpy as np

import config as C


def _ols(X, y):
    A = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return b


def _pred(b, X):
    return b[0] + X @ b[1:]


def splits(groups, n_splits, per_group, seed):
    """List of boolean test masks, each holding per_group[g] compounds of group g."""
    rng = np.random.default_rng(seed)
    idx = {g: np.flatnonzero(groups == g) for g in per_group}
    out = []
    for _ in range(n_splits):
        test = np.zeros(len(groups), bool)
        for g, k in per_group.items():
            if len(idx[g]) > k:
                test[rng.choice(idx[g], size=k, replace=False)] = True
        out.append(test)
    return out


def evaluate(X, y, groups, masks, fit_group):
    """Scores of the model y ~ X (X: n x p)."""
    sq, sq0, r2s, per = [], [], [], []
    for test in masks:
        tr = ~test
        b = _ols(X[tr], y[tr])
        e = y[test] - _pred(b, X[test])
        sq += list(e ** 2)
        per.append(math.sqrt(float((e ** 2).mean())))
        sq0 += list((y[test] - y[tr].mean()) ** 2)
        yt = y[tr]
        r2s.append(1 - ((yt - _pred(b, X[tr])) ** 2).sum() / ((yt - yt.mean()) ** 2).sum())
    b = _ols(X, y)
    res = y - _pred(b, X)
    sst = ((y - y.mean()) ** 2).sum()
    n, p = X.shape
    r2 = 1 - (res ** 2).sum() / sst
    press = 0.0
    for i in range(n):
        m = np.arange(n) != i
        press += float((y[i] - _pred(_ols(X[m], y[m]), X[i:i + 1])[0]) ** 2)
    fg = groups == fit_group
    transfer = None
    if fg.sum() > p + 1 and (~fg).any():
        bt = _ols(X[fg], y[fg])
        transfer = math.sqrt(float(((y[~fg] - _pred(bt, X[~fg])) ** 2).mean()))
    sd = X.std(axis=0, ddof=1)
    return {"RMSE_test": math.sqrt(np.mean(sq)), "RMSE_test_sd": float(np.std(per, ddof=1)),
            "Q2_F1": 1 - sum(sq) / sum(sq0),
            "R2_train": float(np.mean(r2s)), "R2": float(r2),
            "adj_R2": float(1 - (1 - r2) * (n - 1) / (n - p - 1)) if n - p - 1 > 0 else None,
            "Q2_LOO": 1 - press / sst, "RMSE": math.sqrt(float((res ** 2).mean())),
            "RMSE_transfer": transfer, "coef": [float(v) for v in b],
            "std_coef": [float(v) for v in b[1:] * sd / y.std(ddof=1)], "n": n}


def collinearity(X):
    if X.shape[1] < 2:
        return 0.0, 1.0
    r = np.corrcoef(X, rowvar=False)
    maxr = float(np.max(np.abs(r[np.triu_indices_from(r, 1)])))
    try:
        vif = float(np.max(np.diag(np.linalg.inv(r))))
    except np.linalg.LinAlgError:
        vif = float("inf")
    return maxr, vif


def explore(data_rows, levels=None, max_terms=None, descriptors=None):
    """data_rows: correlations.run() data rows. Returns (model rows, info per level)."""
    levels = levels or list(C.LEVELS)
    max_terms = max_terms or C.MLR_MAX_TERMS
    descriptors = descriptors or [d for d, _, _ in C.QSAR_DESCRIPTORS]
    out, info = [], {}
    for lv in levels:
        rows = [r for r in data_rows if r["level"] == lv]
        cand = [d for d in descriptors
                if all(r.get(d) not in (None, "") for r in rows)
                and len({round(float(r[d]), 8) for r in rows}) > 2]
        y = np.array([float(r["ln_kobs"]) for r in rows])
        groups = np.array([r["group"] for r in rows])
        masks = splits(groups, C.MLR_SPLITS, C.MLR_TEST_PER_GROUP, C.MLR_SEED)
        # baseline: predict the training mean
        sq, sq0 = [], []
        for t in masks:
            sq += list((y[t] - y[~t].mean()) ** 2)
        base = math.sqrt(np.mean(sq))
        info[lv] = {"candidates": cand, "n": len(rows), "RMSE_test_mean_model": base,
                    "excluded": [d for d in descriptors if d not in cand]}
        for k in range(1, max_terms + 1):
            for combo in itertools.combinations(cand, k):
                X = np.array([[float(r[d]) for d in combo] for r in rows])
                maxr, vif = collinearity(X)
                s = evaluate(X, y, groups, masks, C.CORR_FIT_GROUP)
                out.append({"level": lv, "k": k, "descriptors": list(combo), "max_abs_r": maxr,
                            "max_VIF": vif, "collinear": maxr > C.MLR_MAX_R, **s})
    out.sort(key=lambda m: m["RMSE_test"])
    return out, info


def predictions(data_rows, level, combo):
    """All-data fit of one model: [(row, observed, predicted)]."""
    rows = [r for r in data_rows if r["level"] == level]
    X = np.array([[float(r[d]) for d in combo] for r in rows])
    y = np.array([float(r["ln_kobs"]) for r in rows])
    b = _ols(X, y)
    return [(r, float(yo), float(yp)) for r, yo, yp in zip(rows, y, _pred(b, X))]


def oos_predictions(data_rows, level, combo):
    """Out-of-sample prediction of each compound: mean of its predictions over the splits in
    which it was in the test set. Returns [(row, observed, mean prediction, times tested)]."""
    rows = [r for r in data_rows if r["level"] == level]
    X = np.array([[float(r[d]) for d in combo] for r in rows])
    y = np.array([float(r["ln_kobs"]) for r in rows])
    groups = np.array([r["group"] for r in rows])
    acc = [[] for _ in rows]
    for test in splits(groups, C.MLR_SPLITS, C.MLR_TEST_PER_GROUP, C.MLR_SEED):
        b = _ols(X[~test], y[~test])
        for i in np.flatnonzero(test):
            acc[i].append(float(_pred(b, X[i:i + 1])[0]))
    return [(r, float(yo), float(np.mean(a)) if a else None, len(a)) for r, yo, a in zip(rows, y, acc)]
