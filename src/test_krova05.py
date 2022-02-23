import pytest
import math

import krovak05

krovak = krovak05.Transformation()


@pytest.mark.parametrize("B,L,expected_result", [
    (50, 14, 45.632),
    (50, 15, 44.438),
    (50.5, 15, 43.746),
])
def test_cr2005_undulation(B, L, expected_result):
    assert pytest.approx(krovak.interpolate_undulation(B, L),
                         abs=0.001) == expected_result


@pytest.mark.parametrize("B,L,H, expected_result", [
    (50, 14, 100, (5774041.357, 6048448.769, 54.368)),
    (50, 15, 100, (5703011.867, 6058147.236, 55.562)),
    (50.5, 15, 100, (5695856.617, 6002995.770, 56.254)),
])
def test_check_etrs_jtsk05(B, L, H, expected_result):
    assert pytest.approx(krovak.etrs_jtsk05(B, L, H),
                         abs=0.001) == expected_result


@pytest.mark.parametrize("B,L,H, expected_result", [
    (50, 14, 100, (774041.354, 1048448.752, 54.368)),
    (50, 15, 100, (703011.898, 1058147.296, 55.562)),
    (50.5, 15, 100, (695856.537, 1002995.859, 56.254)),
])
def test_check_etrs_jtsk(B, L, H, expected_result):
    assert pytest.approx(krovak.etrs_jtsk(B, L, H),
                         abs=0.002) == expected_result


@pytest.mark.parametrize("Y,X,H, expected_result", [
    (774041.354, 1048448.752, 54.368, (50, 14)),
    (703011.898, 1058147.296, 55.562, (50, 15)),
    (695856.537, 1002995.859, 56.254, (50.5, 15)),
])
def test_check_jtsk_etrs(Y, X, H, expected_result):

    B, L, _ = krovak.jtsk_etrs(Y, X, H)

    B_1_cm = (0.01 / 6378000) * (180 / math.pi)
    L_1_cm = (
        0.01 / (6378000 * math.cos(expected_result[0]))) * (180 / math.pi)

    assert pytest.approx(B, abs=B_1_cm) == expected_result[0]
    assert pytest.approx(L, abs=L_1_cm) == expected_result[1]
