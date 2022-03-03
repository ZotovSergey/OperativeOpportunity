import numpy as np
from pyorbital import astronomy, dt2np, tlefile


ECC_EPS = 1.0e-6  # Too low for computing further drops.
ECC_LIMIT_LOW = -1.0e-3
ECC_LIMIT_HIGH = 1.0 - ECC_EPS  # Too close to 1
ECC_ALL = 1.0e-4

EPS_COS = 1.5e-12

NR_EPS = 1.0e-12

CK2 = 5.413080e-4
CK4 = 0.62098875e-6
E6A = 1.0e-6
QOMS2T = 1.88027916e-9
S = 1.01222928
S0 = 78.0
XJ3 = -0.253881e-5
XKE = 0.743669161e-1
XKMPER = 6378.135
XMNPDA = 1440.0
#MFACTOR = 7.292115E-5
AE = 1.0
SECDAY = 8.6400E4

F = 1 / 298.257223563  # Earth flattening WGS-84
A = 6378.137  # WGS84 Equatorial radius


SGDP4_ZERO_ECC = 0
SGDP4_DEEP_NORM = 1
SGDP4_NEAR_SIMP = 2
SGDP4_NEAR_NORM = 3

KS = AE * (1.0 + S0 / XKMPER)
A3OVK2 = (-XJ3 / CK2) * AE**3


class OrbitalError(Exception):
    pass


class Orbital():

    """Class for orbital computations.

        The *satellite* parameter is the name of the satellite to work on and is
        used to retreive the right TLE data for internet or from *tle_file* in case
        it is provided.
        """

    def __init__(self, satellite, tle_file=None, line1=None, line2=None):
        satellite = satellite.upper()
        self.satellite_name = satellite
        self.tle = tlefile.read(satellite, tle_file=tle_file,
                                line1=line1, line2=line2)
        self.orbit_elements = OrbitElements(self.tle)
        self._sgdp4 = _SGDP4(self.orbit_elements)

    def get_position(self, utc_time, normalize=True, infinite_orbit=True):
        """Get the cartesian position and velocity from the satellite.
        """

        kep = self._sgdp4.propagate(utc_time, infinite_orbit=infinite_orbit)
        pos, vel = kep2xyz(kep)

        if normalize:
            pos /= XKMPER
            vel /= XKMPER * XMNPDA / SECDAY

        return pos, vel


class _SGDP4(object):

    """Class for the SGDP4 computations.
    """

    def __init__(self, orbit_elements):
        self.mode = None

        perigee = orbit_elements.perigee
        self.eo = orbit_elements.excentricity
        self.xincl = orbit_elements.inclination
        self.xno = orbit_elements.original_mean_motion
        k_2 = CK2
        k_4 = CK4
        k_e = XKE
        self.bstar = orbit_elements.bstar
        self.omegao = orbit_elements.arg_perigee
        self.xmo = orbit_elements.mean_anomaly
        self.xnodeo = orbit_elements.right_ascension
        self.t_0 = orbit_elements.epoch
        self.xn_0 = orbit_elements.mean_motion
        A30 = -XJ3 * AE**3

        if not(0 < self.eo < ECC_LIMIT_HIGH):
            raise OrbitalError('Eccentricity out of range: %e' % self.eo)
        elif not((0.0035 * 2 * np.pi / XMNPDA) < self.xn_0 < (18 * 2 * np.pi / XMNPDA)):
            raise OrbitalError('Mean motion out of range: %e' % self.xn_0)
        elif not(0 < self.xincl < np.pi):
            raise OrbitalError('Inclination out of range: %e' % self.xincl)


        self.cosIO = np.cos(self.xincl)
        self.sinIO = np.sin(self.xincl)
        theta2 = self.cosIO**2
        theta4 = theta2 ** 2
        self.x3thm1 = 3.0 * theta2 - 1.0
        self.x1mth2 = 1.0 - theta2
        self.x7thm1 = 7.0 * theta2 - 1.0

        a1 = (XKE / self.xn_0) ** (2. / 3)
        betao2 = 1.0 - self.eo**2
        betao = np.sqrt(betao2)
        temp0 = 1.5 * CK2 * self.x3thm1 / (betao * betao2)
        del1 = temp0 / (a1**2)
        a0 = a1 * \
            (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)))
        del0 = temp0 / (a0**2)
        self.xnodp = self.xn_0 / (1.0 + del0)
        self.aodp = (a0 / (1.0 - del0))
        self.perigee = (self.aodp * (1.0 - self.eo) - AE) * XKMPER
        self.apogee = (self.aodp * (1.0 + self.eo) - AE) * XKMPER
        self.period = (2 * np.pi * 1440.0 / XMNPDA) / self.xnodp

        if self.period >= 225:
            # Deep-Space model
            self.mode = SGDP4_DEEP_NORM
        elif self.perigee < 220:
            # Near-space, simplified equations
            self.mode = SGDP4_NEAR_SIMP
        else:
            # Near-space, normal equations
            self.mode = SGDP4_NEAR_NORM

        if self.perigee < 156:
            s4 = self.perigee - 78
            if s4 < 20:
                s4 = 20

            qoms24 = ((120 - s4) * (AE / XKMPER))**4
            s4 = (s4 / XKMPER + AE)
        else:
            s4 = KS
            qoms24 = QOMS2T

        pinvsq = 1.0 / (self.aodp**2 * betao2**2)
        tsi = 1.0 / (self.aodp - s4)
        self.eta = self.aodp * self.eo * tsi
        etasq = self.eta**2
        eeta = self.eo * self.eta
        psisq = np.abs(1.0 - etasq)
        coef = qoms24 * tsi**4
        coef_1 = coef / psisq**3.5

        self.c2 = (coef_1 * self.xnodp * (self.aodp *
                                          (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                                          (0.75 * CK2) * tsi / psisq * self.x3thm1 *
                                          (8.0 + 3.0 * etasq * (8.0 + etasq))))

        self.c1 = self.bstar * self.c2

        self.c4 = (2.0 * self.xnodp * coef_1 * self.aodp * betao2 * (self.eta *
                                                                     (2.0 + 0.5 * etasq) + self.eo * (0.5 + 2.0 *
                                                                                                      etasq) - (2.0 * CK2) * tsi / (self.aodp * psisq) * (-3.0 *
                                                                                                                                                          self.x3thm1 * (1.0 - 2.0 * eeta + etasq *
                                                                                                                                                                         (1.5 - 0.5 * eeta)) + 0.75 * self.x1mth2 * (2.0 *
                                                                                                                                                                                                                     etasq - eeta * (1.0 + etasq)) * np.cos(2.0 * self.omegao))))

        self.c5, self.c3, self.omgcof = 0.0, 0.0, 0.0

        if self.mode == SGDP4_NEAR_NORM:
            self.c5 = (2.0 * coef_1 * self.aodp * betao2 *
                       (1.0 + 2.75 * (etasq + eeta) + eeta * etasq))
            if self.eo > ECC_ALL:
                self.c3 = coef * tsi * A3OVK2 * \
                    self.xnodp * AE * self.sinIO / self.eo
            self.omgcof = self.bstar * self.c3 * np.cos(self.omegao)

        temp1 = 3.0 * CK2 * pinvsq * self.xnodp
        temp2 = temp1 * CK2 * pinvsq
        temp3 = 1.25 * CK4 * pinvsq**2 * self.xnodp

        self.xmdot = (self.xnodp + (0.5 * temp1 * betao * self.x3thm1 + 0.0625 *
                                    temp2 * betao * (13.0 - 78.0 * theta2 +
                                                     137.0 * theta4)))

        x1m5th = 1.0 - 5.0 * theta2

        self.omgdot = (-0.5 * temp1 * x1m5th + 0.0625 * temp2 *
                       (7.0 - 114.0 * theta2 + 395.0 * theta4) +
                       temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4))

        xhdot1 = -temp1 * self.cosIO
        self.xnodot = (xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) +
                                 2.0 * temp3 * (3.0 - 7.0 * theta2)) * self.cosIO)

        if self.eo > ECC_ALL:
            self.xmcof = (-(2. / 3) * AE) * coef * self.bstar / eeta
        else:
            self.xmcof = 0.0

        self.xnodcf = 3.5 * betao2 * xhdot1 * self.c1
        self.t2cof = 1.5 * self.c1

        # Check for possible divide-by-zero for X/(1+cos(xincl)) when
        # calculating xlcof */
        temp0 = 1.0 + self.cosIO
        if np.abs(temp0) < EPS_COS:
            temp0 = np.sign(temp0) * EPS_COS

        self.xlcof = 0.125 * A3OVK2 * self.sinIO * \
            (3.0 + 5.0 * self.cosIO) / temp0

        self.aycof = 0.25 * A3OVK2 * self.sinIO

        self.cosXMO = np.cos(self.xmo)
        self.sinXMO = np.sin(self.xmo)
        self.delmo = (1.0 + self.eta * self.cosXMO)**3

        if self.mode == SGDP4_NEAR_NORM:
            c1sq = self.c1**2
            self.d2 = 4.0 * self.aodp * tsi * c1sq
            temp0 = self.d2 * tsi * self.c1 / 3.0
            self.d3 = (17.0 * self.aodp + s4) * temp0
            self.d4 = 0.5 * temp0 * self.aodp * tsi * \
                (221.0 * self.aodp + 31.0 * s4) * self.c1
            self.t3cof = self.d2 + 2.0 * c1sq
            self.t4cof = 0.25 * \
                (3.0 * self.d3 + self.c1 * (12.0 * self.d2 + 10.0 * c1sq))
            self.t5cof = (0.2 * (3.0 * self.d4 + 12.0 * self.c1 * self.d3 + 6.0 * self.d2**2 +
                                 15.0 * c1sq * (2.0 * self.d2 + c1sq)))

        elif self.mode == SGDP4_DEEP_NORM:
            raise NotImplementedError('Deep space calculations not supported')

    def propagate(self, utc_time, infinite_orbit=True):
        kep = {}

        # get the time delta in minutes
        #ts = astronomy._days(utc_time - self.t_0) * XMNPDA
        # print utc_time.shape
        # print self.t_0
        utc_time = dt2np(utc_time)
        ts = (utc_time - self.t_0) / np.timedelta64(1, 'm')

        em = self.eo
        xinc = self.xincl

        xmp = self.xmo + self.xmdot * ts
        xnode = self.xnodeo + ts * (self.xnodot + ts * self.xnodcf)
        omega = self.omegao + self.omgdot * ts

        if self.mode == SGDP4_ZERO_ECC:
            raise NotImplementedError('Mode SGDP4_ZERO_ECC not implemented')
        elif self.mode == SGDP4_NEAR_SIMP:
            raise NotImplementedError('Mode "Near-space, simplified equations"'
                                      ' not implemented')
        elif self.mode == SGDP4_NEAR_NORM:
            delm = self.xmcof * \
                ((1.0 + self.eta * np.cos(xmp))**3 - self.delmo)
            temp0 = ts * self.omgcof + delm
            xmp += temp0
            omega -= temp0
            if infinite_orbit:
                a = self.aodp
                e = em
                xl = xmp + omega + xnode + self.xnodp
            else:
                tempa = 1.0 - \
                    (ts *
                     (self.c1 + ts * (self.d2 + ts * (self.d3 + ts * self.d4))))
                tempe = self.bstar * \
                    (self.c4 * ts + self.c5 * (np.sin(xmp) - self.sinXMO))
                templ = ts * ts * \
                    (self.t2cof + ts *
                     (self.t3cof + ts * (self.t4cof + ts * self.t5cof)))
                a = self.aodp * tempa**2
                e = em - tempe
                xl = xmp + omega + xnode + self.xnodp * templ

        else:
            raise NotImplementedError('Deep space calculations not supported')

        if np.any(a < 1):
            raise Exception('Satellite crased at time %s', utc_time)
        elif np.any(e < ECC_LIMIT_LOW):
            raise ValueError('Satellite modified eccentricity to low: %s < %e'
                             % (str(e[e < ECC_LIMIT_LOW]), ECC_LIMIT_LOW))

        e = np.where(e < ECC_EPS, ECC_EPS, e)
        e = np.where(e > ECC_LIMIT_HIGH, ECC_LIMIT_HIGH, e)

        beta2 = 1.0 - e**2

        # Long period periodics
        sinOMG = np.sin(omega)
        cosOMG = np.cos(omega)

        temp0 = 1.0 / (a * beta2)
        axn = e * cosOMG
        ayn = e * sinOMG + temp0 * self.aycof
        xlt = xl + temp0 * self.xlcof * axn

        elsq = axn**2 + ayn**2

        if np.any(elsq >= 1):
            raise Exception('e**2 >= 1 at %s', utc_time)

        kep['ecc'] = np.sqrt(elsq)
        epw = np.fmod(xlt - xnode, 2 * np.pi)
        # needs a copy in case of an array
        capu = np.array(epw)
        maxnr = kep['ecc']

        ecosE = 0
        esinE = 0
        cosEPW = 0
        sinEPW = 0

        for i in range(10):
            sinEPW = np.sin(epw)
            cosEPW = np.cos(epw)

            ecosE = axn * cosEPW + ayn * sinEPW
            esinE = axn * sinEPW - ayn * cosEPW
            f = capu - epw + esinE
            if np.all(np.abs(f) < NR_EPS):
                break

            df = 1.0 - ecosE

            # 1st order Newton-Raphson correction.
            nr = f / df

            # 2nd order Newton-Raphson correction.
            nr = np.where(np.logical_and(i == 0, np.abs(nr) > 1.25 * maxnr),
                          np.sign(nr) * maxnr,
                          f / (df + 0.5 * esinE * nr))
            epw += nr

        # Short period preliminary quantities
        temp0 = 1.0 - elsq
        betal = np.sqrt(temp0)
        pl = a * temp0
        r = a * (1.0 - ecosE)
        invR = 1.0 / r
        temp2 = a * invR
        temp3 = 1.0 / (1.0 + betal)
        cosu = temp2 * (cosEPW - axn + ayn * esinE * temp3)
        sinu = temp2 * (sinEPW - ayn - axn * esinE * temp3)
        u = np.arctan2(sinu, cosu)
        sin2u = 2.0 * sinu * cosu
        cos2u = 2.0 * cosu**2 - 1.0
        temp0 = 1.0 / pl
        temp1 = CK2 * temp0
        temp2 = temp1 * temp0

        # Update for short term periodics to position terms.

        rk = r * (1.0 - 1.5 * temp2 * betal * self.x3thm1) + \
            0.5 * temp1 * self.x1mth2 * cos2u
        uk = u - 0.25 * temp2 * self.x7thm1 * sin2u
        xnodek = xnode + 1.5 * temp2 * self.cosIO * sin2u
        xinck = xinc + 1.5 * temp2 * self.cosIO * self.sinIO * cos2u

        if np.any(rk < 1):
            raise Exception('Satellite crashed at time %s', utc_time)

        temp0 = np.sqrt(a)
        temp2 = XKE / (a * temp0)
        rdotk = ((XKE * temp0 * esinE * invR - temp2 * temp1 * self.x1mth2 * sin2u) *
                 (XKMPER / AE * XMNPDA / 86400.0))
        rfdotk = ((XKE * np.sqrt(pl) * invR + temp2 * temp1 *
                   (self.x1mth2 * cos2u + 1.5 * self.x3thm1)) *
                  (XKMPER / AE * XMNPDA / 86400.0))

        kep['radius'] = rk * XKMPER / AE
        kep['theta'] = uk
        kep['eqinc'] = xinck
        kep['ascn'] = xnodek
        kep['argp'] = omega
        kep['smjaxs'] = a * XKMPER / AE
        kep['rdotk'] = rdotk
        kep['rfdotk'] = rfdotk

        return kep


class OrbitElements():

    """Class holding the orbital elements.
    """

    def __init__(self, tle):
        self.epoch = tle.epoch
        self.excentricity = tle.excentricity
        self.inclination = np.deg2rad(tle.inclination)
        self.right_ascension = np.deg2rad(tle.right_ascension)
        self.arg_perigee = np.deg2rad(tle.arg_perigee)
        self.mean_anomaly = np.deg2rad(tle.mean_anomaly)

        self.mean_motion = tle.mean_motion * (np.pi * 2 / XMNPDA)
        self.mean_motion_derivative = tle.mean_motion_derivative * \
            np.pi * 2 / XMNPDA ** 2
        self.mean_motion_sec_derivative = tle.mean_motion_sec_derivative * \
            np.pi * 2 / XMNPDA ** 3
        self.bstar = tle.bstar * AE

        n_0 = self.mean_motion
        k_e = XKE
        k_2 = CK2
        i_0 = self.inclination
        e_0 = self.excentricity

        a_1 = (k_e / n_0) ** (2.0 / 3)
        delta_1 = ((3 / 2.0) * (k_2 / a_1**2) * ((3 * np.cos(i_0)**2 - 1) /
                                                 (1 - e_0**2)**(2.0 / 3)))

        a_0 = a_1 * (1 - delta_1 / 3 - delta_1**2 - (134.0 / 81) * delta_1**3)

        delta_0 = ((3 / 2.0) * (k_2 / a_0**2) * ((3 * np.cos(i_0)**2 - 1) /
                                                 (1 - e_0**2)**(2.0 / 3)))

        # original mean motion
        n_0pp = n_0 / (1 + delta_0)
        self.original_mean_motion = n_0pp

        # semi major axis
        a_0pp = a_0 / (1 - delta_0)
        self.semi_major_axis = a_0pp

        self.period = np.pi * 2 / n_0pp

        self.perigee = (a_0pp * (1 - e_0) / AE - AE) * XKMPER

        self.right_ascension_lon = (self.right_ascension
                                    - astronomy.gmst(self.epoch))

        if self.right_ascension_lon > np.pi:
            self.right_ascension_lon -= 2 * np.pi


def kep2xyz(kep):
    sinT = np.sin(kep['theta'])
    cosT = np.cos(kep['theta'])
    sinI = np.sin(kep['eqinc'])
    cosI = np.cos(kep['eqinc'])
    sinS = np.sin(kep['ascn'])
    cosS = np.cos(kep['ascn'])

    xmx = -sinS * cosI
    xmy = cosS * cosI

    ux = xmx * sinT + cosS * cosT
    uy = xmy * sinT + sinS * cosT
    uz = sinI * sinT

    x = kep['radius'] * ux
    y = kep['radius'] * uy
    z = kep['radius'] * uz

    vx = xmx * cosT - cosS * sinT
    vy = xmy * cosT - sinS * sinT
    vz = sinI * cosT

    v_x = kep['rdotk'] * ux + kep['rfdotk'] * vx
    v_y = kep['rdotk'] * uy + kep['rfdotk'] * vy
    v_z = kep['rdotk'] * uz + kep['rfdotk'] * vz

    return np.array((x, y, z)), np.array((v_x, v_y, v_z))