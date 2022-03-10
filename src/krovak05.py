import os
import math
import linecache


def read_specific_line(file_path, lines_idx: list):
    """ Selecting specific lines from a file
    Args:
        file_path (path): ..
        lines_idx (list): index of search rows

    Returns:
        output_lines (list[str]): list of selected lines
    """
    output_lines = []

    for line in lines_idx:
        output_lines.append(linecache.getline(file_path, int(line)))

    return output_lines


class Transformation:

    POSSIBLE_DIFF_TABLE_NAMES = ["table_yx_3_v1710",
                                 "table_yx_3_v1202", "table_yx_3_v1005"]

    def __init__(self, diff_table_name: str = "table_yx_3_v1710"):

        self.diff_table_name = diff_table_name
        self._heigh_table_path = os.path.join(os.path.dirname(__file__),
                                              "kvazigeoids", "CR-2005_v1005.csv")

    @property
    def diff_table_name(self):
        return self._diff_table_name

    @diff_table_name.setter
    def diff_table_name(self, name):

        if name in self.POSSIBLE_DIFF_TABLE_NAMES:
            self._diff_table_name = name
        else:
            print(
                "DIFF table named {} doesn't exist. Default value was set!".format(name))
            self._diff_table_name = self.POSSIBLE_DIFF_TABLE_NAMES[0]

        self._diff_table_path = os.path.join(os.path.dirname(__file__),
                                             "difference_tables", self._diff_table_name + ".csv")

        print(self._diff_table_path)

    def get_available_diff_tables(self):
        """Return available differential tables

        Usage for measuring:

        table_yx_3_v1710 - (DEFAULT) obtained by measurement from 1.1.2018

        table_yx_3_v1202 - obtained by measurement between 1.7.2012-30.6.2018

        table_yx_3_v1005 - obtained by measurement between 2.1.2011-31.12.2012


        Returns:
            (list): list of differential table names 
        """

        return self.POSSIBLE_DIFF_TABLE_NAMES

    def interpolate_dydx(self, Y, X):
        """Interpolation of difference values dy,dx from difference table grid for conversion
        from S-JTSK/05 to S-JTSK coordinate system.

        Y_jtsk = Y_jstk05 - 5_000_000 + dy

        X_jtsk = X_jstk05 - 5_000_000 + dx

        Example: 

        Krovak05.interpolate_dydx(750000, 1050000) --> (0.072, -0.037) 

        Args:
            Y (float): Y S-JTSK coordinate
            X (float): X S-JTSK coordinate

        Returns:
            dY (float): dY correction [m]
            dX (float): dX correction [m]
        """

        DIFF_X = 2000
        DIFF_Y = 2000

        Xmin = (X - X % DIFF_X)
        Ymin = (Y - Y % DIFF_Y)

        nXmin = (Xmin - 930000) / DIFF_X
        line = nXmin * 241 + 1 + (Ymin - 428000) / DIFF_Y

        lines_data = read_specific_line(
            self._diff_table_path, [line, line + 1, line + 241, line + 242])

        LeftDown, RightDown, LeftUp, RightUp = [
            [float(data) for data in lines_data[idx].strip().split(",")] for idx in range(0, 4)]

        # DY interpolation
        nYUp = LeftUp[2] + ((RightUp[2] - LeftUp[2]) / DIFF_Y) * (Y - Ymin)
        nYDown = LeftDown[2] + \
            ((RightDown[2] - LeftDown[2]) / DIFF_Y) * (Y - Ymin)
        dY1 = nYDown + ((nYUp - nYDown) / DIFF_X) * (X - Xmin)

        nYLeft = LeftDown[2] + \
            ((LeftUp[2] - LeftDown[2]) / DIFF_X) * (X - Xmin)
        nYRight = RightDown[2] + \
            ((RightUp[2] - RightDown[2]) / DIFF_X) * (X - Xmin)
        dY2 = nYLeft + ((nYRight - nYLeft) / DIFF_Y) * (Y - Ymin)

        # DX interpolation
        nXUp = LeftUp[3] + ((RightUp[3] - LeftUp[3]) / DIFF_Y) * (Y - Ymin)
        nXDown = LeftDown[3] + \
            ((RightDown[3] - LeftDown[3]) / DIFF_Y) * (Y - Ymin)
        dX1 = nXDown + ((nXUp - nXDown) / DIFF_X) * (X - Xmin)

        nXLeft = LeftDown[3] + \
            ((LeftUp[3] - LeftDown[3]) / DIFF_X) * (X - Xmin)
        nXRight = RightDown[3] + \
            ((RightUp[3] - RightDown[3]) / DIFF_X) * (X - Xmin)
        dX2 = nXLeft + ((nXRight - nXLeft) / DIFF_Y) * (Y - Ymin)

        return round((dY1 + dY2) / 2, 5), round((dX1 + dX2) / 2, 5)

    def interpolate_undulation(self, B, L):
        """Interpolation of kvasigeoid undulation (N) from CR2005 data set for calculation of Bpv height:
        Hbpv = Hel - N

        example:

        Krovak05.interpolate_undulation(50,15) -> -0.267 m

        Args:
            B (float): ETRS89 latitude
            L (float): ETRS89 longitude

        Returns:
            N (float): undulation [m]
        """
        DIFF_B = 0.01666666666666666666
        DIFF_L = 0.0250000

        Bmin = math.floor((B - B % DIFF_B) * 100000) / 100000
        Lmin = math.floor((L - L % DIFF_L) * 100000) / 100000

        nBmin = round((Bmin - 48.30) / DIFF_B)
        line = int(nBmin * 306 + 1 + (Lmin - 11.70) / DIFF_L)

        lines_data = read_specific_line(
            self._heigh_table_path, [line, line + 1, line + 306, line + 307])

        # closest points data format B_i,L_i,DH_i ==> LeftUp,RightUp,LeftDown,RightDown corners
        LeftDown, RightDown, LeftUp, RightUp = [
            [float(data) for data in lines_data[idx].split(",")] for idx in range(0, 4)]

        nUp = LeftUp[2] + ((RightUp[2] - LeftUp[2]) / DIFF_L) * (L - Lmin)
        nDown = LeftDown[2] + \
            ((RightDown[2] - LeftDown[2]) / DIFF_L) * (L - Lmin)
        N = nDown + ((nUp - nDown) / DIFF_B) * (B - Bmin)

        return round(N, 5)

    def bicub_dotr(self, Y, X):
        """Calculation of bicubic dotransformation correction increments
        to mitigate the S-JTSK and S-JTSK/05 mismatch

        Args:
            Y (float): Y S-JTSK coordinate [m]
            X (float): X S-JTSK coordinate [m]

        Returns:
            dY (float): dY correction [m]
            dX (float): dX correction [m]
        """

        Yred = Y - 654_000
        Xred = X - 1_089_000

        A1 = 0.2946529277 * 10 ** -1
        A2 = 0.2515965696 * 10 ** -1
        A3 = 0.1193845912 * 10 ** -6
        A4 = -0.4668270147 * 10 ** -6
        A5 = 0.9233980362 * 10 ** -11
        A6 = 0.1523735715 * 10 ** -11
        A7 = 0.1696780024 * 10 ** -17
        A8 = 0.4408314235 * 10 ** -17
        A9 = -0.8331083518 * 10 ** -23
        A10 = -0.3689471323 * 10 ** -23

        dY = A2 + A3 * Yred + A4 * Xred + A5 * (2 * Xred * Yred) + A6 * (Xred ** 2 - Yred ** 2) + A8 * (Xred * (Xred ** 2 - 3 * Yred ** 2)) + A7 * (
            Yred * (3 * Xred ** 2 - Yred ** 2)) - A10 * (4 * Yred * Xred * (Xred ** 2 - Yred ** 2)) + A9 * (Xred ** 4 + Yred ** 4 - 6 * Xred ** 2 * Yred ** 2)

        dX = A1 + A3 * Xred - A4 * Yred + A5 * (Xred ** 2 - Yred ** 2) - A6 * (2 * Xred * Yred) + A7 * (Xred * (Xred ** 2 - 3 * Yred ** 2)) - A8 * (
            Yred * (3 * Xred ** 2 - Yred ** 2)) + A9 * (4 * Yred * Xred * (Xred ** 2 - Yred ** 2)) + A10 * (Xred ** 4 + Yred ** 4 - 6 * Xred ** 2 * Yred ** 2)

        return dY, dX

    def etrs_jtsk05(self, B, L, H):
        """Transformation sferic coordinates B, L, H to S-JTSK/05 

        Args:
            B (float): ETRS89 latitude [deg]
            L (float): ETRS89 longitude [deg]
            H (float): Elipsoidal height [m]

        Returns:
            Y (float): y-coordinate S-JTSK/05 [m]
            X (float): x-coordinate S-JTSK/05 [m]
            H (float): height Baltic after adjustment [m]
        """

        # Input transformation parameters --->
        RHO = math.pi / 180

        A_GRS80 = 6378137
        E2_GRS80 = 0.006694380022901
        E_GRS80 = math.sqrt(E2_GRS80)

        A_BESSEL = 6377397.155
        E2_BESSEL = 0.00667437223062
        E_BESSEL = math.sqrt(E2_BESSEL)

        k1 = 0.9999
        F0 = 49.50 * RHO
        S0 = 78.50 * RHO
        Uq = (59 + 42 / 60 + 42.69689 / 3600) * RHO

        alfa = math.sqrt(
            1 + (E2_BESSEL * math.pow(math.cos(F0), 4)) / (1 - E2_BESSEL))
        U0 = math.asin(math.sin(F0) / alfa)
        gf0 = ((1 + E_BESSEL * math.sin(F0)) /
               (1 - E_BESSEL * math.sin(F0))) ** ((alfa * E_BESSEL) / 2)
        k = math.tan(U0 / 2 + math.pi / 4) * \
            ((1 / math.tan((F0 / 2 + math.pi / 4)))**alfa) * gf0
        n = math.sin(S0)
        N0 = (A_BESSEL * math.sqrt(1 - E2_BESSEL)) / \
            (1 - E2_BESSEL * math.sin(F0) ** 2)
        ro0 = k1 * N0 * (1 / math.tan(S0))
        aa_ = math.pi / 2 - Uq

        # <---

        # Parameters of Helmert 3D transformation GRS80 --> Bessel
        # shifts
        p1 = -572.203
        p2 = -85.328
        p3 = -461.934
        # scale
        p4 = (1 + (-3.5393 * 1e-6))
        # rotation angles
        rRho = 206264.806
        p5 = 5.24832714 / rRho
        p6 = 1.52900087 / rRho
        p7 = 4.97311727 / rRho

        # Calculate Bpv height
        H_bpv = H - self.interpolate_undulation(B, L)

        # Convert degrees to radians
        B_ = B * RHO
        L_ = L * RHO

        # Conversion of ellipsoidal to rectangular coordinates
        N_etrs = A_GRS80 / math.sqrt(1 - E2_GRS80 * math.sin(B_)**2)

        X_etrs = (N_etrs + H) * math.cos(B_) * math.cos(L_)
        Y_etrs = (N_etrs + H) * math.cos(B_) * math.sin(L_)
        Z_etrs = (N_etrs * (1 - E2_GRS80) + H) * math.sin(B_)

        # Helmert transformation from GRS80 ellipsoid -> Bessel ellipsoid

        X1 = p1 + p4 * (X_etrs + p5 * Y_etrs - p6 * Z_etrs)
        Y1 = p2 + p4 * (-p5 * X_etrs + Y_etrs + p7 * Z_etrs)
        Z1 = p3 + p4 * (p6 * X_etrs - p7 * Y_etrs + Z_etrs)

        # Reverse transfer to ellipsoidal coordinates
        LL = math.atan(Y1 / X1)
        BB_0 = 0
        BB_i = math.atan((Z1 / (math.sqrt(X1**2 + Y1**2))) *
                         (1 + (E2_BESSEL / (1 - E2_BESSEL))))

        while abs(BB_0 - BB_i) > 1e-12:
            BB_0 = BB_i

            N_i = A_BESSEL / math.sqrt(1 - E2_BESSEL * math.sin(BB_i)**2)
            H_i = (math.sqrt(X1**2 + Y1**2) / math.cos(BB_i)) - N_i
            BB_i = math.atan((Z1 / math.sqrt(X1**2 + Y1**2)) *
                             (1 - (N_i * E2_BESSEL) / (N_i + H_i))**(-1))

        BB = BB_i

        # Modified Krovak projection
        gB = ((1 + E_BESSEL * math.sin(BB)) /
              (1 - E_BESSEL * math.sin(BB)))**((alfa * E_BESSEL) / 2)
        U = 2 * (math.atan(k * (math.tan(BB / 2 + math.pi / 4)**alfa)
                           * (gB**(-1))) - math.pi / 4)
        dV = alfa * ((42.5 * RHO - (LL + (17 + 40 / 60) * RHO)))

        # Conversion of geographical coordinates to Cartografical
        S = math.asin(math.cos(aa_) * math.sin(U) +
                      math.sin(aa_) * math.cos(U) * math.cos(dV))
        D = math.asin((math.cos(U) * math.sin(dV)) / math.cos(S))

        Epsilon = n * D

        ro = ro0 * (math.tan(S0 / 2 + math.pi / 4)**n) * \
            (math.tan(S / 2 + math.pi / 4)**(-n))

        Y_ = ro * math.sin(Epsilon)
        X_ = ro * math.cos(Epsilon)

        # Bicubic do-transformation
        deltaY, deltaX = self.bicub_dotr(Y_, X_)

        return Y_ - deltaY + 5_000_000, X_ - deltaX + 5_000_000, H_bpv

    def jtsk_jtsk05(self, Y_jtsk, X_jtsk):
        """Transition from S-JTSK to S-JTSK/05 coordinate system using
        correction tables defined in <diff_table_name>.

        Args:
            Y_jtsk (float): y-coordinate S-JTSK [m]
            X_jtsk (float): x-coordinate S-JTSK [m]

        Returns:
            Y_jtsk05 (float): y-coordinate S-JTSK/05 [m]
            X_jtsk05 (float): x-coordinate S-JTSK/05 [m]
        """

        dy, dx = self.interpolate_dydx(Y_jtsk, X_jtsk)

        return Y_jtsk + 5_000_000 + dy, X_jtsk + 5_000_000 + dx

    def jtsk05_jtsk(self, Y_jtsk05, X_jtsk05):
        """Transition from S-JTSK/05 to S-JTSK coordinate system using
        correction tables defined in <diff_table_name>.

        Args:
            Y_jtsk05 (float): y-coordinate S-JTSK/05 [m]
            X_jtsk05 (float): x-coordinate S-JTSK/05 [m]

        Returns:
            Y_jtsk (float): y-coordinate S-JTSK [m]
            X_jtsk (float): x-coordinate S-JTSK [m]
        """

        Y_ = Y_jtsk05 - 5_000_000
        X_ = X_jtsk05 - 5_000_000

        dy, dx = self.interpolate_dydx(Y_, X_)

        return Y_ - dy, X_ - dx

    def etrs_jtsk(self, B, L, H):
        """Transformation sferic coordinates B, L, H ETRS89 (ETRF2000) to S-JTSK coordinates system 

        Args:
            B (float): ETRS89 latitude [deg]
            L (float): ETRS89 longitude [deg]
            H (float): Elipsoidal height [m]

        Returns:
            Y (float): y-coordinate S-JTSK [m]
            X (float): x-coordinate S-JTSK [m]
            H (float): height Baltic after adjustment [m]
        """

        Y_jtsk05, X_jtsk05, H_bpv = self.etrs_jtsk05(B, L, H)

        Y_jtsk, X_jtsk = self.jtsk05_jtsk(Y_jtsk05, X_jtsk05)

        return Y_jtsk, X_jtsk, H_bpv

    def jtsk_etrs(self, Y, X, H):
        """Transformation S-JTSK coordinates to sferic coordinates B, L, H ETRS89 (ETRF2000)

        **reverse transformation to Transformation.etrs_jtsk**

        Args:
            Y (float): Y S-JTSK [m]
            X (float): X S-JTSK [m]
            H (float): height Bpv [m]

        Returns:
            B (float): ETRS89 latitude [deg]
            L (float): ETRS89 longitude [deg]
            H (float): Elipsoidal height [m]
        """

        A_GRS80 = 6378137
        E2_GRS80 = 0.006694380022901

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
        alfa = math.sqrt(
            1 + ((E2_BESSEL * (math.cos(Fi_0)**4)) / (1 - E2_BESSEL)))
        U_0 = math.asin(math.sin(Fi_0) / alfa)
        k_ = math.tan(U_0 / 2 + math.pi / 4) * ((1 / math.tan(Fi_0 / 2 + math.pi / 4))**alfa) * \
            (((1 + E_BESSEL * math.sin(Fi_0)) /
              (1 - E_BESSEL * math.sin(Fi_0)))**(alfa * E_BESSEL / 2))

        # Re-introduction of S-JTSK05->S-JTSK corrections
        dY, dX = self.interpolate_dydx(Y, X)

        Y05 = Y + dY
        X05 = X + dX

        # Re-introduction of bicubic dotransformation
        dYb, dXb = self.bicub_dotr(Y05, X05)

        Y_ = Y05 + dYb
        X_ = X05 + dXb

        # Conversion from plane to cone shell
        Ro = math.sqrt(X_**2 + Y_**2)
        Epsilon = math.atan(Y_ / X_)

        # Conversion from a cone shell to a sphere (cartographic coordinates) -> S,D

        D = Epsilon / math.sin(S_0)
        S = 2 * (math.atan(((Ro0 / Ro)**(1 / n)) *
                           math.tan(S_0 / 2 + math.pi / 4)) - math.pi / 4)

        # Cartographic coordinates -> geographical coordinates -> U,V

        U = math.asin(math.cos(a_) * math.sin(S) -
                      math.sin(a_) * math.cos(S) * math.cos(D))
        dV = math.asin((math.cos(S) * math.sin(D)) / math.cos(U))

        # Sphere to Besseluv elipsoid -> B_bessel, L_Bessel

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

        # Conversion to rectangular coordinates

        HH = H + self.interpolate_undulation(B_bessel / RHO, L_bessel / RHO)

        NN = A_BESSEL / math.sqrt(1 - E2_BESSEL * math.sin(B_bessel)**2)

        X_bessel = (NN + HH) * math.cos(B_bessel) * math.cos(L_bessel)
        Y_bessel = (NN + HH) * math.cos(B_bessel) * math.sin(L_bessel)
        Z_bessel = (NN * (1 - E2_BESSEL) + HH) * math.sin(B_bessel)

        # Conversion to elipsoid GRS80 - helmertov transformation
        # parameters:
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

        # Conversion to geographical coordinates

        grs_dist = math.sqrt(X_grs**2 + Y_grs**2)
        B_grs_0 = 1
        B_grs_i = math.atan((Z_grs / grs_dist) *
                            (1 + E2_GRS80 / (1 - E2_GRS80)))

        while (B_grs_0 - B_grs_i) > 1e-15:

            B_grs_0 = B_grs_i

            NN_i = A_GRS80 / math.sqrt(1 - E2_GRS80 * math.sin(B_grs_0)**2)
            HH_e = grs_dist / math.cos(B_grs_0) - NN_i
            B_grs_i = math.atan((Z_grs / grs_dist) *
                                ((1 - (NN_i * E2_GRS80) / (NN_i + HH_e))**(-1)))

        B_grs = B_grs_i / RHO
        L_grs = math.atan(Y_grs / X_grs) / RHO
        H_grs = H + self.interpolate_undulation(B_grs, L_grs)

        return B_grs, L_grs, H_grs


if __name__ == "__main__":

    krovak = Transformation()

    # Undulation of kvasigeoid
    undulation = krovak.interpolate_undulation(50, 15)
    print(undulation)

    # Differences between S-JTSK and S-JTSK/05
    dy, dx = krovak.interpolate_dydx(750000, 1050000)
    print(dy, dx)

    # Get list of possible dydx grid data
    grids = krovak.get_available_diff_tables()
    print(grids)

    # Transform ETRS89 (ETRF2000) coordinates to S-JTSK/05
    B_etrs_in = 50
    L_etrs_in = 15
    H_etrs_in = 100

    Y_sjtsk05, X_sjtsk05, H_bpv = krovak.etrs_jtsk05(
        B_etrs_in, L_etrs_in, H_etrs_in)
    print(Y_sjtsk05, X_sjtsk05, H_bpv)

    # Transform ETRS89 (ETRF2000) coordinates to S-JTSK
    Y_sjtsk, X_sjtsk, H_bpv = krovak.etrs_jtsk(B_etrs_in, L_etrs_in, H_etrs_in)
    print(Y_sjtsk, X_sjtsk, H_bpv)

    # Reverse transformation S-JTSK coordinate to ETRS89
    B_etrs_out, L_etrs_out, H_etrs_out = krovak.jtsk_etrs(
        Y_sjtsk, X_sjtsk, H_bpv)
    print(B_etrs_out, L_etrs_out, H_etrs_out)

    print("Differences:")
    print(f"dB = {(B_etrs_in-B_etrs_out)*(math.pi/180)*6378000*1000} mm")
    print(f"dL = {(B_etrs_in-B_etrs_out)*(math.pi/180)*(6378000*math.cos(B_etrs_in*(math.pi/180)))*1000} mm")
    print(f"dH = {(H_etrs_in-H_etrs_out)*1000} mm")
