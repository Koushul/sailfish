import numpy as np

from kinetics import classify, size_normalize


def test_size_normalize_uses_shared_library():
    counts = np.array([[1.0, 3.0], [2.0, 2.0]])
    lib = np.array([10.0, 20.0])
    out = size_normalize(counts, lib)
    assert np.allclose(out[0], counts[0] * (15.0 / 10.0))
    assert np.allclose(out[1], counts[1] * (15.0 / 20.0))


def test_classify_requires_hypoxic_spliced_for_reverting():
    sample = np.array(["E15S", "E15S", "E15S", "E14S"])
    theta = np.array([0.1, 0.95, 0.1, 0.9])
    v = np.array([2.0, 2.0, 0.0, 2.0])
    mem = np.zeros(4)
    out = classify(sample, theta, v, mem, v_cut=1.0, t_low=0.3, t_high=0.7, mem_cut=1.0)
    assert out[0] == "reverting"
    assert out[1] == "reverted"
    assert out[2] == "persistent"
    assert out[3] == "never_hypoxic"


def test_classify_requires_not_fully_hypoxic_for_inducing():
    sample = np.array(["E15S", "E15S"])
    theta = np.array([0.05, 0.8])
    v = np.array([-2.0, -2.0])
    mem = np.zeros(2)
    out = classify(sample, theta, v, mem, v_cut=1.0, t_low=0.3, t_high=0.7, mem_cut=1.0)
    assert out[0] == "persistent"
    assert out[1] == "inducing"
