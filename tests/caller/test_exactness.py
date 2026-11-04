"""Seeded checks that the speedups give bit-identical results to the upstream code paths.

Each test compares a fast helper with an inline copy of the upstream v0.2.8 code (or with the
unchanged upstream method). Runs in seconds; no fixtures needed.
"""

import random
import struct

import numpy as np
import pytest

from ribotish.zbio import interval as I
from ribotish.zbio import ribo as R
from ribotish.zbio import stat as S


def _bits(v):
    if v is None:
        return None
    return type(v).__name__, struct.pack("d", v)


def _upstream_combination_log(n, k):
    n, k = int(n), int(k)
    if n < 0:
        return None
    if k > n or k < 0:
        return None
    return S.logsum(n) - S.logsum(k) - S.logsum(n - k)


def _upstream_count_dict(arr):
    data = {}
    for d in arr:
        if d not in data:
            data[d] = 0
        data[d] += 1
    return data


def _upstream_frame_lists(cnts, start, stop):
    n = (stop - start) // 3
    inarr = [0] * n
    outarr = 2 * n * [0]
    for i in range(n):
        p = start + i * 3
        inarr[i] = cnts[p]
        outarr[i * 2] = cnts[p + 1]
        outarr[i * 2 + 1] = cnts[p + 2]
    return (inarr, outarr)


def test_combination_log_bits():
    rng = random.Random(1)
    ns = list(range(0, 400)) + [65534, 65535, 65536, 70000]
    ns += [rng.randint(0, 5000) for _ in range(300)]
    for n in ns:
        for k in (-1, 0, 1, n // 3, n // 2, n - 1, n, n + 1, rng.randint(-2, n + 2)):
            for f in (int, float, np.int64):
                assert _bits(S.combination_log(f(n), f(k))) == _bits(
                    _upstream_combination_log(f(n), f(k))
                ), (n, k, f)


def test_count_dict_items_and_types():
    rng = random.Random(2)
    pool = [0, 0, 0, 1, 2, 3, 7, 1.0, True, np.int64(1), 40]
    for _ in range(5000):
        arr = [rng.choice(pool) for _ in range(rng.randint(0, 60))]
        a, b = _upstream_count_dict(arr), S.countDict(arr)
        assert type(b) is dict
        assert list(a.items()) == list(b.items())
        assert [type(k) for k in a] == [type(k) for k in b]
        assert [type(v) for v in a.values()] == [type(v) for v in b.values()]


def test_frame_test_lists_and_pvalues():
    rng = random.Random(3)
    captured = []
    orig = R.rstest

    def spy(x, y, delta=1e-4, show=False):
        captured.append((list(x), list(y)))
        return orig(x, y, delta=delta, show=show)

    class T:
        pass

    R.rstest = spy
    try:
        for _ in range(3000):
            L = rng.randint(1, 60)
            cnts = [
                rng.choice([0] * 8 + [1, 1, 2, 3, rng.randint(0, 30)]) for _ in range(3 * L + 6)
            ]
            t = T()
            t.cnts = cnts
            st = rng.randint(0, 5)
            for s1, s2 in ((st, st + 3 * rng.randint(0, L)), (st, len(cnts) + 3), (-3, st + 9)):
                try:
                    exp = _upstream_frame_lists(cnts, s1, s2)
                except IndexError:
                    with pytest.raises(IndexError):
                        R.Ribo.frame_test(t, s1, s2)
                    continue
                captured.clear()
                p = R.Ribo.frame_test(t, s1, s2)
                if captured:
                    assert captured[0] == tuple(exp)
                    assert _bits(p) == _bits(R._rstest(exp[0], exp[1], 1e-4, False))
    finally:
        R.rstest = orig


def test_rstest_cache_matches_direct():
    rng = random.Random(4)
    R._rstestCache.clear()
    R._rstestCacheItems = 0
    for _ in range(3000):
        n = rng.randint(1, 40)
        x = [rng.choice([0] * 6 + [1, 2, 3, rng.randint(0, 20)]) for _ in range(n)]
        y = [rng.choice([0] * 8 + [1, 2, rng.randint(0, 10)]) for _ in range(2 * n)]
        direct = R._rstest(x, y, 1e-4, False)
        assert _bits(R.rstest(x, y)) == _bits(direct)
        assert _bits(R.rstest(x, y)) == _bits(direct)  # second call is a cache hit
        xs = x[:]
        rng.shuffle(xs)  # same count dict, possibly different key order
        assert _bits(R.rstest(xs, y)) == _bits(R._rstest(xs, y, 1e-4, False))
    big_x = list(range(200))
    big_y = list(range(100, 400))
    before = len(R._rstestCache)
    assert _bits(R.rstest(big_x, big_y)) == _bits(R._rstest(big_x, big_y, 1e-4, False))
    assert len(R._rstestCache) == before  # key above _rstestKeyMax is not cached


def _rand_lst(rng, wellformed):
    out = []
    p = rng.randint(0, 20)
    for _ in range(rng.randint(1, 5)):
        if wellformed:
            s = p + rng.choice([0, 0, 1, 2, 5, 30])
            e = s + rng.randint(1, 40)
        else:
            s = rng.randint(0, 300)
            e = s + rng.randint(-3, 40)
        out.append((s, e))
        p = e
    return out


def _upstream_compat(blocks, other, mis):
    try:
        return I.Interval(itvs=blocks).is_compatible(other, mis=mis)
    except Exception as e:  # noqa: BLE001 - compare exception types too
        return type(e)


def _fast_compat(blocks, other, mis):
    try:
        return I.blocks_compatible(blocks, other, mis=mis)
    except Exception as e:  # noqa: BLE001
        return type(e)


def test_blocks_compatible_normalized_transcripts():
    rng = random.Random(5)
    seen = set()
    for _ in range(30000):
        blocks = tuple(_rand_lst(rng, rng.random() < 0.85))
        t = I.Interval(itvs=_rand_lst(rng, rng.random() < 0.8))
        if len(t) == 0:
            continue
        mis = rng.choice([0, 1, 2, 3, 5])
        a = _upstream_compat(blocks, t, mis)
        assert _fast_compat(blocks, t, mis) == a, (blocks, t.lst, mis)
        seen.add(a)
    assert True in seen and False in seen


def test_blocks_compatible_unnormalized_transcripts():
    """other.lst appended to after check() (overlapping/unsorted) must fall back exactly."""
    rng = random.Random(6)
    for _ in range(20000):
        blocks = tuple(_rand_lst(rng, True))
        t = I.Interval(itvs=_rand_lst(rng, True))
        t.lst += _rand_lst(rng, rng.random() < 0.5)  # no check(): may overlap or be unsorted
        mis = rng.choice([0, 1, 2, 5])
        assert _fast_compat(blocks, t, mis) == _upstream_compat(blocks, t, mis), (blocks, t.lst)
    # concrete case: a measure-only answer double counts (10, 20) and returns True; upstream
    # re-normalizes to the union (0, 60) + (70, 100), so 10 bases differ and the answer is False
    t = I.Interval(itvs=[(0, 60), (70, 100)])
    t.lst.insert(1, (10, 20))
    assert _upstream_compat(((0, 100),), t, 2) is False
    assert _fast_compat(((0, 100),), t, 2) is False
