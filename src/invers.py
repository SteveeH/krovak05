import math
from .krovak05 import Transformation


def reverse(Y, X, H):
    krovak = Transformation()

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

    a_ = (90 - (59 + 42 / 60 + 42.69689 / 3600)) * RHO
    alfa = math.sqrt(1 + ((E2_BESSEL * (math.cos(Fi_0)**4)) / (1 - E2_BESSEL)))
    U_0 = math.asin(math.sin(Fi_0) / alfa)
    k_ = math.tan(U_0 / 2 + math.pi / 4) * ((1 / math.tan(Fi_0 / 2 + math.pi / 4))**alfa) * \
        (((1 + E_BESSEL * math.sin(Fi_0)) /
          (1 - E_BESSEL * math.sin(Fi_0)))**(alfa * E_BESSEL / 2))

    # zpetne zavedeni korekci tabulky
    dY, dX = krovak.interpolate_dydx(Y, X)

    Y05 = Y + dY
    X05 = X + dX

    # zpetne zavedeni bikubicke dotransformacnich clenu

    dYb, dXb = krovak.bicub_dotr(Y05, X05)

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

    U = math.asin(math.cos(a_) * math.sin(S) -
                  math.sin(a_) * math.cos(S) * math.cos(D))
    dV = math.asin((math.cos(S) * math.sin(D)) / math.cos(U))

    # koule na Besseluv elipsoid -> B_bessel, L

    L_bessel = (24 + 50 / 60) * RHO - dV / alfa

    B_0 = 0
    B_i = U

    counter = 0

    while abs(B_0 - B_i) > 1e-15:
        counter += 1
        B_0 = B_i
        B_i = 2 * (math.atan((k_**(-1 / alfa)) * (math.tan(U / 2 + math.pi / 4)**(1 / alfa)) *
                   (((1 + E_BESSEL * math.sin(B_0)) / (1 - E_BESSEL * math.sin(B_0)))**(E_BESSEL / 2))) - math.pi / 4)

    B_bessel = B_i

    # print(f"{counter=}")

    # prevoda na pravouhle souradnice

    HH = H + krovak.interpolate_undulation(B_bessel / RHO, L_bessel / RHO)

    NN = A_BESSEL / math.sqrt(1 - E2_BESSEL * math.sin(B_bessel)**2)

    X_bessel = (NN + HH) * math.cos(B_bessel) * math.cos(L_bessel)
    Y_bessel = (NN + HH) * math.cos(B_bessel) * math.sin(L_bessel)
    Z_bessel = (NN * (1 - E2_BESSEL) + HH) * math.sin(B_bessel)

    # prevod na elipsoid GRS80 helmertova transformace
    # parametry:
    r_ = 206264.806
    pp1 = 572.213
    pp2 = 85.334
    pp3 = 461.940
    pp4 = 1 + 3.5378 * 1e-6
    pp5 = -5.24836073 / r_
    pp6 = -1.52899176 / r_
    pp7 = -4.97316164 / r_

    X_grs = pp1 + pp4 * (X_bessel + pp5 * Y_bessel - pp6 * Z_bessel)
    Y_grs = pp2 + pp4 * (-pp5 * X_bessel + Y_bessel + pp7 * Z_bessel)
    Z_grs = pp3 + pp4 * (pp6 * X_bessel - pp7 * Y_bessel + Z_bessel)

    print(X_grs, Y_grs, Z_grs)

    # prevod na zemepisne souradnice

    grs_dist = math.sqrt(X_grs**2 + Y_grs**2)
    B_grs_0 = 1
    B_grs_i = math.atan((Z_grs / grs_dist) * (1 + E2_GRS80 / (1 - E2_GRS80)))

    while (B_grs_0 - B_grs_i) > 1e-15:

        B_grs_0 = B_grs_i

        NN_i = A_GRS80 / math.sqrt(1 - E2_GRS80 * math.sin(B_grs_0)**2)
        HH_e = grs_dist / math.cos(B_grs_0) - NN_i
        B_grs_i = math.atan((Z_grs / grs_dist) *
                            ((1 - (NN_i * E2_GRS80) / (NN_i + HH_e))**(-1)))

    L_grs = math.atan(Y_grs / X_grs) / RHO
    B_grs = B_grs_i / RHO

    return B_grs, L_grs, HH


if __name__ == "__main__":

    B, L, H = reverse(695856.537, 1002995.859, 56.254)

    print(B, L, H)

    print("Odchylka B => {} mm".format((B - 50.5) * RHO * 6378000 * 1000))
    print("Odchylka L => {} mm".format(
        (L - 15) * RHO * 6378000 * math.cos(50.5 * RHO) * 1000))
