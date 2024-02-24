using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Runtime.InteropServices;
using static System.Math;

namespace AA.Net {
    public static class Saturn{
        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Saturn, julianEphemerisDay);
        }

        public static (TimeSpan rightAscension, double declination) PositionGeocentric(double julianEphemerisDay) {
            return Body.PositionGeocentric(Bodies.Saturn, julianEphemerisDay);
        }

        public static double NextPerihelion(double start) {
            return Body.NextApsis(Bodies.Saturn, Apsis.Perihelion, start);
        }

        public static double NextAphelion(double start) {
            return Body.NextApsis(Bodies.Saturn, Apsis.Aphelion, start);
        }

        public static (double fraction, double positionAngle, double magnitude) Illumination(double julianEphemerisDay) {
            return Body.Illumination(Bodies.Saturn, julianEphemerisDay);
        }

        public static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(double julianEphemerisDay, bool J2000 = false) {
            return Body.OrbitalElements(Bodies.Saturn, julianEphemerisDay, J2000);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Saturn, julianEphemerisDay);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>
        /// <list type="table">
        /// <item><term>majorAxis</term>
        /// <description>Size of the major axis of the outer edge of the outer ring in arcseconds.</description>
        /// </item>
        /// <item><term>minorAxis</term>
        /// <description>Size of the minor axis in arcseconds.</description>
        /// </item>
        /// <item><term>positionAngle</term>
        /// <description>Eastward angle of the northern semiminor axis from the northern celestial pole.</description>
        /// </item>
        /// <item><term>latitudeEarth</term>
        /// <description>Saturnicentric latitude of Earth referred to plane of the ring. When positive, the northern side of the ring is visible.</description>
        /// </item>
        /// <item><term>latitudeSun</term>
        /// <description>Saturnicentric latitude of Sun. When positive, the northern side of the ring is illuminated.</description>
        /// </item>
        /// <item><term>deltaLongitude</term>
        /// <description>Difference between Saturnicentric longitudes of Sun and Earth.</description>
        /// </item>
        /// </list>
        /// </returns>
        /// <remarks>pp. 317-320</remarks>
        public static (double majorAxis, double minorAxis, double positionAngle, double latitudeEarth, double latitudeSun, double deltaLongitude)
        Rings(double julianEphemerisDay) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = T * T;
            double i = 28.075216 - 0.012998 * T + 0.000004 * T2;
            double omega = 169.50847 + 1.394681 * T + 0.000412 * T2;
            var poz = Body.PositionFromLightTime(Bodies.Saturn, julianEphemerisDay);
            double iRad = i.ToRadians();
            double omegaRad = omega.ToRadians();
            double B = Asin(Sin(iRad) * Cos(poz.beta) * Sin(poz.lambda - omegaRad) - Cos(iRad) * Sin(poz.beta));
            double major = 375.35 / poz.delta;
            double minor = major * Sin(Abs(B));
            double N = 113.6655 + 0.8771 * T;
            double lPrime = (poz.l - 0.01759 / poz.r).ToRadians();
            double bPrime = (poz.b - 0.000764 * (Cos(poz.l.ToRadians()) - N) / poz.r).ToRadians();
            double lmo = lPrime - omegaRad;
            double BPrime = Asin(Sin(iRad) * Cos(bPrime) * Sin(lmo) - Cos(iRad) * Sin(bPrime)).ToDegrees();
            double U1 = Atan2(Sin(iRad) * Sin(bPrime) + Cos(iRad) * Cos(bPrime) * Sin(lmo), Cos(bPrime) * Cos(lmo));
            lmo = poz.lambda - omegaRad;
            double U2 = Atan2(Sin(iRad) * Sin(poz.beta) + Cos(iRad) * Cos(poz.beta) * Sin(lmo), Cos(poz.beta) * Cos(lmo));
            double deltaU = Abs(U1 - U2).ToDegrees();
            double lambda0 = omega - 90;
            double beta0 = 90 - i;
            double lambda1 = poz.lambda.ToDegrees().To360() + 0.005693 * Cos(poz.l0.ToRadians() - poz.lambda) / Cos(poz.beta);
            poz.beta = poz.beta.ToDegrees() + 0.005693 * Sin(poz.l0.ToRadians() - poz.lambda) * Sin(poz.beta);
            poz.lambda = lambda1;
            double nutL = Sky.Nutation(julianEphemerisDay).longitude;
            lambda0 += nutL;
            poz.lambda += nutL;
            double epsilon = Sky.ObliquityOfEcliptic(julianEphemerisDay);
            var equ = Transform.ToEquatorial(lambda0, beta0, epsilon);
            double ra0 = equ.rightAscension.ToRadians();
            double d0 = equ.declination.ToRadians();
            equ = Transform.ToEquatorial(poz.lambda, poz.beta, epsilon);
            double ra = equ.rightAscension.ToRadians();
            double d = equ.declination.ToRadians();
            double P = Atan2(Cos(d0) * Sin(ra0 - ra), Sin(d0) * Cos(d) - Cos(d0) * Sin(d) * Cos(ra0 - ra)).ToDegrees();
            return (major, minor, P, B.ToDegrees(), BPrime, deltaU);
        }

        private struct MoonData {
            internal double lambda;
            internal double r;
            internal double gamma;
            internal double omega;
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>For Satellites I - VIII, in order, apparent rectangular coordinates from the center of Saturn (North, West and Far Side positive) in units of Saturn's equatorial radius.</returns>
        /// <remarks>pp. 323-333</remarks>
        public static (double x, double y, double z)[]
        Satellites (double julianEphemerisDay) {
            var poz = Body.PositionFromLightTime(Bodies.Saturn, julianEphemerisDay);
            var bes = Transform.Precession(poz.lambda.ToDegrees(), poz.beta.ToDegrees(), Time.J2000, Time.B1950);
            double lambdaZero = bes.longitude.ToRadians();
            double betaZero = bes.latitude.ToRadians();
            double JDE = julianEphemerisDay - poz.tau;
            double t1 = JDE - 2411093;
            double t2 = t1 / 365.25;
            double t3 = (JDE - 2433282.423) / 365.25 + 1950;
            double t4 = JDE - 2411368;
            double t5 = t4 / 365.25;
            double t6 = JDE - 2415020;
            double t7 = t6 / 36525;
            double t8 = t6 / 365.25;
            double t9 = (JDE - 2442000.5) / 365.25;
            double t10 = JDE - 2409786;
            double t11 = t10 / 36525;
            double W0 = (5.095 * (t3 - 1866.39)).ToRadians();
            double W1 = (74.4 + 32.39 * t2).ToRadians();
            double W2 = (134.3 + 96.62 * t2).ToRadians();
            double W3 = (42 - 0.5118 * t5).ToRadians();
            double W4 = (276.59 + 0.5118 * t5).ToRadians();
            double W5 = (267.2635 + 1222.1136 * t7).ToRadians();
            double W6 = (175.4762 + 1221.5515 * t7).ToRadians();
            double W7 = (2.4891 + 0.002435 * t7).ToRadians();
            double W8 = (113.35 - 0.2597 * t7).ToRadians();
            double s1 = Sin(28.0817.ToRadians());
            double s2 = Sin(168.8112.ToRadians());
            double c1 = Cos(28.0817.ToRadians());
            double c2 = Cos(168.8112.ToRadians());
            double e1 = 0.05589 - 0.000346 * t7;
            MoonData moon;
            (double x, double y, double z)[] result = new (double, double, double)[8];
            var fict = ABC(0, 0, 1);
            double D = Atan2(fict.A, fict.C);
            for (int m = 0; m < 8; m++) {
                double L, p, M, C, e, lambdaPrime, i, omega, a, K = 0;
                moon = new MoonData();
                switch (m) {
                    case 0:
                        L = 127.64 + 381.994497 * t1 - 43.57 * Sin(W0) - 0.72 * Sin(3 * W0) - 0.02144 * Sin(5 * W0);
                        p = 106.1 + 365.549 * t2;
                        M = (L - p).ToRadians();
                        C = 2.18287 * Sin(M) + 0.025988 * Sin(2 * M) + 0.00043 * Sin(3 * M);
                        moon.lambda = L + C;
                        moon.r = 3.06879 / (1 + 0.01905 * Cos(M + C.ToRadians()));
                        moon.gamma = 1.563;
                        moon.omega = 54.5 - 365.072 * t2;
                        K = 20947;
                        break;
                    case 1:
                        L = 200.317 + 262.7319002 * t1 + 0.25667 * Sin(W1) + 0.20883 * Sin(W2);
                        p = 309.107 + 123.44121 * t2;
                        M = (L - p).ToRadians();
                        C = 0.55577 * Sin(M) + 0.00168 * Sin(2 * M);
                        moon.lambda = L + C;
                        moon.r = 3.94118 / (1 + 0.00485 * Cos(M + C.ToRadians()));
                        moon.gamma = 0.0262;
                        moon.omega = 348 - 151.95 * t2;
                        K = 23715;
                        break;
                    case 2:
                        moon.lambda = 285.306 + 190.69791226 * t1 + 2.063 * Sin(W0) + 0.03409 * Sin(3 * W0) + 0.001015 * Sin(5 * W0);
                        moon.r = 4.880998;
                        moon.gamma = 1.0976;
                        moon.omega = 111.33 - 72.2441 * t2;
                        K = 26382;
                        break;
                    case 3:
                        L = 254.712 + 131.53493193 * t1 - 0.0215 * Sin(W1) - 0.01733 * Sin(W2);
                        p = 174.8 + 30.82 * t2;
                        M = (L - p).ToRadians();
                        C = 0.24717 * Sin(M) + 0.00033 * Sin(2 * M);
                        moon.lambda = L + C;
                        moon.r = 6.24871 / (1 + 0.002157 * Cos(M + C.ToRadians()));
                        moon.gamma = 0.0139;
                        moon.omega = 232 - 30.27 * t2;
                        K = 29876;
                        break;
                    case 4:
                        double pPrime = (342.7 + 10.057 * t2).ToRadians();
                        double a1 = 0.000265 * Sin(pPrime) + 0.01 * Sin(W4);
                        double a2 = 0.000265 * Cos(pPrime) + 0.01 * Cos(W4);
                        e = Sqrt(a1 * a1 + a2 * a2);
                        p = Atan2(a1, a2).ToDegrees();
                        double N = (345 - 10.057 * t2).ToRadians();
                        lambdaPrime = 359.244 + 79.6900472 * t1 + 0.086754 * Sin(N);
                        i = 28.0362 + 0.346898 * Cos(N) + 0.0193 * Cos(W3);
                        omega = 168.8034 + 0.736936 * Sin(N) + 0.041 * Sin(W3);
                        a = 8.725924;
                        Subroutine();
                        K = 35313;
                        break;
                    case 5:
                        L = 261.1582 + 22.57697855 * t4 + 0.074025 * Sin(W3);
                        double iPrime = (27.45141 + 0.295999 * Cos(W3)).ToRadians();
                        double omegaPrime = (168.66925 + 0.628808 * Sin(W3)).ToRadians();
                        a1 = Sin(W7) * Sin(omegaPrime - W8);
                        a2 = Cos(W7) * Sin(iPrime) - Sin(W7) * Cos(iPrime) * Cos(omegaPrime - W8);
                        double g0 = 102.8623.ToRadians();
                        double psi = Atan2(a1, a2);
                        double s = Sqrt(a1 * a1 + a2 * a2);
                        double g = W4 - omegaPrime - psi;
                        double varpi = 0;
                        for (int j = 0; j < 3; j++) {
                            varpi = W4 + 0.37515.ToRadians() * (Sin(2 * g) - Sin(2 * g0));
                            g = W4 - omegaPrime - psi;
                        }
                        double ePrime = 0.029092 + 0.00019048 * (Cos(2 * g) - Cos(2 * g0));
                        double q = 2 * (W5 - varpi);
                        double b1 = Sin(iPrime) * Sin(omegaPrime - W8);
                        double b2 = Cos(W7) * Sin(iPrime) * Cos(omegaPrime - W8) - Sin(W7) * Cos(iPrime);
                        double theta = Atan2(b1, b2) + W8;
                        e = ePrime + 0.002778797 * ePrime * Cos(q);
                        p = varpi.ToDegrees() + 0.159215 * Sin(q);
                        double u = 2 * W5 - 2 * theta + psi;
                        double h = 0.9375 * ePrime * ePrime * Sin(q) + 0.1875 * s * s * Sin(2 * (W5 - theta));
                        lambdaPrime = L - 0.254744 * (e1 * Sin(W6) + 0.75 * e1 * e1 * Sin(2 * W6) + h);
                        i = iPrime.ToDegrees() + 0.031843 * s * Cos(u);
                        omega = omegaPrime.ToDegrees() + 0.031843 * s * Sin(u) / Sin(iPrime);
                        a = 20.216193;
                        Subroutine();
                        K = 53800;
                        break;
                    case 6:
                        double eta = (92.39 + 0.5621071 * t6).ToRadians();
                        double zeta = (148.19 - 19.18 * t8).ToRadians();
                        theta = (184.8 - 35.41 * t9).ToRadians();
                        double thetaPrime = theta - 7.5.ToRadians();
                        double aS = (176 + 12.22 * t8).ToRadians();
                        double bS = (8 + 24.44 * t8).ToRadians();
                        double cS = bS + 5.0.ToRadians();
                        varpi = 69.898 - 18.67088 * t8;
                        double phi = (2 * (varpi - W5)).ToRadians();
                        double chi = (94.9 - 2.292 * t8).ToRadians();
                        a = 24.50601 - 0.08686 * Cos(eta) - 0.00166 * Cos(zeta + eta) + 0.00175 * Cos(zeta - eta);
                        e = 0.103458 - 0.004099 * Cos(eta) - 0.000167 * Cos(zeta + eta)
                            + 0.000235 * Cos(zeta - eta) + 0.02303 * Cos(zeta) - 0.00212 * Cos(2 * zeta)
                            + 0.000151 * Cos(3 * zeta) + 0.00013 * Cos(phi);
                        p = varpi + 0.15648 * Sin(chi) - 0.4457 * Sin(eta) - 0.2657 * Sin(zeta + eta)
                            - 0.3573 * Sin(zeta - eta) - 12.872 * Sin(zeta) + 1.668 * Sin(2 * zeta)
                            - 0.2419 * Sin(3 * zeta) - 0.07 * Sin(phi);
                        lambdaPrime = 177.047 + 16.91993829 * t6 + 0.15648 * Sin(chi) + 9.142 * Sin(eta)
                            + 0.007 * Sin(2 * eta) - 0.014 * Sin(3 * eta) + 0.2275 * Sin(zeta + eta)
                            + 0.2112 * Sin(zeta - eta) - 0.26 * Sin(zeta) - 0.0098 * Sin(2 * zeta)
                            - 0.013 * Sin(aS) + 0.017 * Sin(bS) - 0.0303 * Sin(phi);
                        i = 27.3347 + 0.643486 * Cos(chi) + 0.315 * Cos(W3) + 0.018 * Cos(theta) - 0.018 * Cos(cS);
                        omega = 168.6812 + 1.40136 * Cos(chi) + 0.68599 * Sin(W3) - 0.0392 * Sin(cS) + 0.0366 * Sin(thetaPrime);
                        Subroutine();
                        K = 59222;
                        break;
                    case 7:
                        L = 261.1582 + 22.57697855 * t4;
                        double varpiPrime = (91.796 + 0.562 * t7).ToRadians();
                        psi = (4.367 - 0.195 * t7).ToRadians();
                        theta = (146.819 - 3.198 * t7).ToRadians();
                        phi = (60.470 + 1.521 * t7).ToRadians();
                        double PHI = (205.055 - 2.091 * t7).ToRadians();
                        ePrime = 0.028298 + 0.001156 * t11;
                        double varpiZero = (352.91 + 11.71 * t11).ToRadians();
                        double mu = (76.3852 + 4.53795125 * t10).ToRadians();
                        iPrime = (18.4602 - 0.9518 * t11 - 0.072 * t11 * t11 + 0.0054 * t11 * t11 * t11).ToRadians();
                        omegaPrime = (143.198 - 3.919 * t11 + 0.116 * t11 * t11 + 0.008 * t11 * t11 * t11).ToRadians();
                        double l = mu - varpiZero;
                        g = varpiZero - omegaPrime - psi;
                        double g1 = varpiZero - omegaPrime - phi;
                        double lS = W5 - varpiPrime;
                        double gS = varpiPrime - theta;
                        double lT = L.ToRadians() - W4;
                        double gT = W4 - PHI;
                        double u1 = 2 * (l + g - lS - gS);
                        double u2 = l + g1 - lT - gT;
                        double u3 = l + 2 * (g - lS - gS);
                        double u4 = lT + gT - g1;
                        double u5 = 2 * (lS + gS);
                        a = 58.935028 + 0.004638 * Cos(u1) + 0.058222 * Cos(u2);
                        e = ePrime - 0.0014097 * Cos(g1 - gT) + 0.0003733 * Cos(u5 - 2 * g)
                            + 0.000118 * Cos(u3) + 0.0002408 * Cos(l)
                            + 0.0002849 * Cos(l + u2) + 0.000619 * Cos(u4);
                        double w = 0.08077 * Sin(g1 - gT) + 0.02139 * Sin(u5 - 2 * g) - 0.00676 * Sin(u3)
                            + 0.0138 * Sin(l) + 0.01632 * Sin(l + u2) + 0.03547 * Sin(u4);
                        p = varpiZero.ToDegrees() + w / ePrime;
                        lambdaPrime = mu.ToDegrees() - 0.04299 * Sin(u2) - 0.00789 * Sin(u1) - 0.06312 * Sin(lS)
                            - 0.00295 * Sin(2 * lS) - 0.02231 * Sin(u5) + 0.0065 * Sin(u5 + psi);
                        i = iPrime.ToDegrees() + 0.04204 * Cos(u5 + psi) + 0.00235 * Cos(l + g1 + lT + gT + phi) + 0.0036 * Cos(u2 + phi);
                        double wPrime = 0.04204 * Sin(u5 + psi) + 0.00235 * Sin(l + g1 + lT + gT + phi) + 0.00358 * Sin(u2 + phi);
                        omega = omegaPrime.ToDegrees() + wPrime / Sin(iPrime);
                        Subroutine();
                        K = 91820;
                        break;
                }
                moon.lambda = moon.lambda.To360();
                moon.omega = moon.omega.To360();
                double uu = (moon.lambda - moon.omega).ToRadians();
                double ww = (moon.omega - 168.8112).ToRadians();
                moon.gamma = moon.gamma.ToRadians();
                double X = moon.r * (Cos(uu) * Cos(ww) - Sin(uu) * Cos(moon.gamma) * Sin(ww));
                double Y = moon.r * (Sin(uu) * Cos(ww) * Cos(moon.gamma) + Cos(uu) * Sin(ww));
                double Z = moon.r * Sin(uu) * Sin(moon.gamma);
                var abc = ABC(X, Y, Z);
                X = abc.A * Cos(D) - abc.C * Sin(D);
                Y = abc.A * Sin(D) + abc.C * Cos(D);
                Z = abc.B;
                X += Abs(Z) / K * Sqrt(1 - Pow(X / moon.r, 2));
                double W = poz.delta / (poz.delta + Z / 2475);
                X *= W;
                Y *= W;
                result[m] = (X, Y, Z);

                void Subroutine() {
                    double e2 = e * e;
                    double e3 = e2 * e;
                    double e4 = e3 * e;
                    double e5 = e4 * e;
                    lambdaPrime = lambdaPrime.To360();
                    omega = omega.To360();
                    M = (lambdaPrime - p).ToRadians();
                    C = (2 * e - 0.25 * e3 + 0.0520833333 * e5) * Sin(M)
                        + (1.25 * e2 - 0.458333333 * e4) * Sin(2 * M)
                        + (1.083333333 * e3 - 0.671875 * e5) * Sin(3 * M)
                        + 1.072917 * e4 * Sin(4 * M) + 1.142708 * e5 * Sin(5 * M);
                    moon.r = a * (1 - e2) / (1 + e * Cos(M + C));
                    double g = (omega - 168.8112).ToRadians();
                    i = i.ToRadians();
                    double a1 = Sin(i) * Sin(g);
                    double a2 = c1 * Sin(i) * Cos(g) - s1 * Cos(i);
                    moon.gamma = Asin(Sqrt(a1 * a1 + a2 * a2)).ToDegrees();
                    double u = Atan2(a1, a2).ToDegrees();
                    moon.omega = 168.8112 + u;
                    double h = c1 * Sin(i) - s1 * Cos(i) * Cos(g);
                    double psi = Atan2(s1 * Sin(g), h).ToDegrees();
                    moon.lambda = lambdaPrime + C.ToDegrees() + u - g.ToDegrees() - psi;
                }
            }
            return result;
            
            (double A, double B, double C) ABC(double X, double Y, double Z) {
                double B1 = c1 * Y - s1 * Z;
                double C1 = s1 * Y + c1 * Z;
                double A2 = c2 * X - s2 * B1;
                double B2 = s2 * X + c2 * B1;
                double A3 = A2 * Sin(lambdaZero) - B2 * Cos(lambdaZero);
                double B3 = A2 * Cos(lambdaZero) + B2 * Sin(lambdaZero);
                double B4 = B3 * Cos(betaZero) + C1 * Sin(betaZero);
                double C4 = C1 * Cos(betaZero) - B3 * Sin(betaZero);
                return (A3, B4, C4);
            }
        }
    }
}
