import sys
import pytest

sys.path.insert(0, '../src/')

from transformation import Krovak05

krovak = Krovak05()


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
