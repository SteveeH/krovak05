import math
from krovak05 import Transformation


k = Transformation()

A_GRS80 = 6378137
E2_GRS80 = 0.006694380022901
E_GRS80 = math.sqrt(E2_GRS80)

A_BESSEL = 6377397.155
E2_BESSEL = 0.00667437223062
E_BESSEL = math.sqrt(E2_BESSEL)


RHO = math.pi / 180
S_0 = 78.50 * RHO
Fi_0 = 49.50 * RHO
n = math.sin(S_0)

N_0 = (A_BESSEL * math.sqrt(1 - E2_BESSEL)) / \
    (1 - E2_BESSEL * math.sin(Fi_0)**2)
Ro0 = 0.9999 * N_0 * (1 / math.tan(S_0))

a_ = (90 - (59 + 42/60 + 42.69689/3600))*RHO
alfa = math.sqrt(1 + (E2_BESSEL*math.cos(Fi_0)**4)/(1-E2_BESSEL))


def reverse(Y, X, H):

    # zpetne zavedeni korekci tabulky
    dY, dX = k.interpolate_dydx(Y, X)

    Y05 = Y + dY
    X05 = X + dX

    # zpetne zavedeni bikubicke dotransformacnich clenu

    dYb, dXb = k.bicub_dotr(Y05, X05)

    Y_ = Y05 + dYb
    X_ = X05 + dXb

    # rovina na plast kuzele
    Ro = math.sqrt(X_**2 + Y_**2)
    Epsilon = math.atan(Y_ / X_)

    # plast kuzele na kouli (kartograficke souradnice) -> Š,D

    D = Epsilon / math.sin(S_0)
    S = 2 * (math.atan(((Ro0 / Ro)**(1 / n)) *
             math.tan(S_0 / 2 + math.pi / 4)) - math.pi / 4)

    # kartograficke souradnice -> zemepisne souradnice -> U,V

    U = math.asin( math.cos(a_)*math.sin(S) - math.sin(a_)*math.cos(S)*math.cos(D) )
    dV = math.asin((math.cos(S)*math.sin(D))/math.cos(U))

    # koule na Besseluv elipsoid -> B, L

    L = (24 + 50/60)*RHO - dV/alfa

    B0
    


    print(dY, dX)

    return None, None, None


if __name__ == "__main__":

    B, L, H = reverse(695856.537, 1002995.859, 56.254)
