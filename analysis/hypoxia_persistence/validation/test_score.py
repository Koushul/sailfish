import numpy as np

from score import perm_auroc, visium_hex_moran


def test_hex_moran_positive_on_smooth_field():
    rows, cols, vals = [], [], []
    for r in range(8):
        for c in range(r % 2, 16, 2):
            rows.append(r)
            cols.append(c)
            vals.append(float(r + c))
    mi = visium_hex_moran(vals, rows, cols)
    assert mi > 0.5


def test_hex_moran_near_zero_on_noise():
    rng = np.random.default_rng(0)
    rows, cols = [], []
    for r in range(10):
        for c in range(r % 2, 20, 2):
            rows.append(r)
            cols.append(c)
    vals = rng.normal(size=len(rows))
    mi = visium_hex_moran(vals, rows, cols)
    assert abs(mi) < 0.25


def test_perm_auroc_perfect_separator():
    y = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    s = np.array([0.0, 0.1, 0.2, 0.3, 0.8, 0.9, 1.0, 1.1])
    obs, p = perm_auroc(y, s, n=50, seed=0)
    assert obs == 1.0
    assert p < 0.05
