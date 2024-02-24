using System;
using System.Collections.Generic;
using static System.Math;
using System.Text;

namespace AA.Net {
    public static class Jupiter {
        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Jupiter, julianEphemerisDay);
        }

        public static (TimeSpan rightAscension, double declination) PositionGeocentric(double julianEphemerisDay) {
            return Body.PositionGeocentric(Bodies.Jupiter, julianEphemerisDay);
        }

        public static double NextPerihelion(double start) {
            return Body.NextApsis(Bodies.Jupiter, Apsis.Perihelion, start);
        }

        public static double NextAphelion(double start) {
            return Body.NextApsis(Bodies.Jupiter, Apsis.Aphelion, start);
        }

        public static (double fraction, double positionAngle, double magnitude) Illumination(double julianEphemerisDay) {
            return Body.Illumination(Bodies.Jupiter, julianEphemerisDay);
        }

        public static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(double julianEphemerisDay, bool J2000 = false) {
            return Body.OrbitalElements(Bodies.Jupiter, julianEphemerisDay, J2000);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Jupiter, julianEphemerisDay);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>For Satellites I - IV, in order, apparent rectangular coordinates from the center of Jupiter (North, West and Far Side positive) in units of Jupiter's equatorial radius.</returns>
        public static (double x, double y, double z)[] Satellites(double julianEphemerisDay) {
            var poz = Body.PositionFromLightTime(Bodies.Jupiter, julianEphemerisDay);
            double t = julianEphemerisDay - 2443000.5 - poz.tau;
            double l1 = (106.07719 + 203.488955790 * t).ToRadians();
            double l2 = (175.73161 + 101.374724735 * t).ToRadians();
            double l3 = (120.55883 + 50.317609207 * t).ToRadians();
            double l4 = (84.44459 + 21.571071177 * t).ToRadians();
            double pi1 = (97.0881 + 0.16138586 * t).ToRadians();
            double pi2 = (154.8663 + 0.04726307 * t).ToRadians();
            double pi3 = (188.184 + 0.00712734 * t).ToRadians();
            double pi4 = (335.2868 + 0.00184 * t).ToRadians();
            double omega1 = (312.3346 - 0.13279386 * t).ToRadians();
            double omega2 = (100.4411 - 0.03263064 * t).ToRadians();
            double omega3 = (119.1942 - 0.00717703 * t).ToRadians();
            double omega4 = (322.6186 - 0.00175934 * t).ToRadians();
            double Gamma = 0.33033 * Sin((163.679 + 0.0010512 * t).ToRadians()) + 0.03439 * Sin((34.486 - 0.0161731 * t).ToRadians());
            double Philambda = (199.6766 + 0.1737919 * t).ToRadians();
            double psi = (316.5182 - 0.00000208 * t).ToRadians();
            double G = (30.23756 + 0.0830925701 * t + Gamma).ToRadians();
            double GPrime = (31.97853 + 0.0334597339 * t).ToRadians();
            double Pi = 13.469942.ToRadians();
            double Sigma1 = (0.47259 * Sin(2 * (l1 - l2))
                        - 0.03478 * Sin(pi3 - pi4)
                        + 0.01081 * Sin(l2 - 2 * l3 + pi3)
                        + 0.00738 * Sin(Philambda)
                        + 0.00713 * Sin(l2 - 2 * l3 + pi2)
                        - 0.00674 * Sin(pi1 + pi3 - 2 * Pi - 2 * G)
                        + 0.00666 * Sin(l2 - 2 * l3 + pi4)
                        + 0.00445 * Sin(l1 - pi3)
                        - 0.00354 * Sin(l1 - l2)
                        - 0.00317 * Sin(2 * psi - 2 * Pi)
                        + 0.00265 * Sin(l1 - pi4)
                        - 0.00186 * Sin(G)
                        + 0.00162 * Sin(pi2 - pi3)
                        + 0.00158 * Sin(4 * (l1 - l2))
                        - 0.00155 * Sin(l1 - l3)
                        - 0.00138 * Sin(psi + omega3 - 2 * Pi - 2 * G)
                        - 0.00115 * Sin(2 * (l1 - 2 * l2 + omega2))
                        + 0.00089 * Sin(pi2 - pi4)
                        + 0.00085 * Sin(l1 + pi3 - 2 * Pi - 2 * G)
                        + 0.00083 * Sin(omega2 - omega3)
                        + 0.00053 * Sin(psi - omega2)).ToRadians();
            double Sigma2 = (1.06476 * Sin(2 * (l2 - l3))
                        + 0.04256 * Sin(l1 - 2 * l2 + pi3)
                        + 0.03581 * Sin(l2 - pi3)
                        + 0.02395 * Sin(l1 - 2 * l2 + pi4)
                        + 0.01984 * Sin(l2 - pi4)
                        - 0.01778 * Sin(Philambda)
                        + 0.01654 * Sin(l2 - pi2)
                        + 0.01334 * Sin(l2 - 2 * l3 + pi2)
                        + 0.01294 * Sin(pi3 - pi4)
                        - 0.01142 * Sin(l2 - l3)
                        - 0.01057 * Sin(G)
                        - 0.00775 * Sin(2 * (psi - Pi))
                        + 0.00524 * Sin(2 * (l1 - l2))
                        - 0.00460 * Sin(l1 - l3)
                        + 0.00316 * Sin(psi - 2 * G + omega3 - 2 * Pi)
                        - 0.00203 * Sin(pi1 + pi3 - 2 * Pi - 2 * G)
                        + 0.00146 * Sin(psi - omega3)
                        - 0.00145 * Sin(2 * G)
                        + 0.00125 * Sin(psi - omega4)
                        - 0.00115 * Sin(l1 - 2 * l3 + pi3)
                        - 0.00094 * Sin(2 * (l2 - omega2))
                        + 0.00086 * Sin(2 * (l1 - 2 * l2 + omega2))
                        - 0.00086 * Sin(5 * GPrime - 2 * G + 52.225.ToRadians())
                        - 0.00078 * Sin(l2 - l4)
                        - 0.00064 * Sin(3 * l3 - 7 * l4 + 4 * pi4)
                        + 0.00064 * Sin(pi1 - pi4)
                        - 0.00063 * Sin(l1 - 2 * l3 + pi4)
                        + 0.00058 * Sin(omega3 - omega4)
                        + 0.00056 * Sin(2 * (psi - Pi - G))
                        + 0.00056 * Sin(2 * (l2 - l4))
                        + 0.00055 * Sin(2 * (l1 - l3))
                        + 0.00052 * Sin(3 * l3 - 7 * l4 + pi3 + 3 * pi4)
                        - 0.00043 * Sin(l1 - pi3)
                        + 0.00041 * Sin(5 * (l2 - l3))
                        + 0.00041 * Sin(pi4 - Pi)
                        + 0.00032 * Sin(omega2 - omega3)
                        + 0.00032 * Sin(2 * (l3 - G - Pi))).ToRadians();
            double Sigma3 = (0.16490 * Sin(l3 - pi3)
                        + 0.09081 * Sin(l3 - pi4)
                        - 0.06907 * Sin(l2 - l3)
                        + 0.03784 * Sin(pi3 - pi4)
                        + 0.01846 * Sin(2 * (l3 - l4))
                        - 0.01340 * Sin(G)
                        - 0.01014 * Sin(2 * (psi - Pi))
                        + 0.00704 * Sin(l2 - 2 * l3 + pi3)
                        - 0.00620 * Sin(l2 - 2 * l3 + pi2)
                        - 0.00541 * Sin(l3 - l4)
                        + 0.00381 * Sin(l2 - 2 * l3 + pi4)
                        + 0.00235 * Sin(psi - omega3)
                        + 0.00198 * Sin(psi - omega4)
                        + 0.00176 * Sin(Philambda)
                        + 0.00130 * Sin(3 * (l3 - l4))
                        + 0.00125 * Sin(l1 - l3)
                        - 0.00119 * Sin(5 * GPrime - 2 * G + 52.225.ToRadians())
                        + 0.00109 * Sin(l1 - l2)
                        - 0.00100 * Sin(3 * l3 - 7 * l4 + 4 * pi4)
                        + 0.00091 * Sin(omega3 - omega4)
                        + 0.00080 * Sin(3 * l3 - 7 * l4 + pi3 + 3 * pi4)
                        - 0.00075 * Sin(2 * l2 - 3 * l3 + pi3)
                        + 0.00072 * Sin(pi1 + pi3 - 2 * Pi - 2 * G)
                        + 0.00069 * Sin(pi4 - Pi)
                        - 0.00058 * Sin(2 * l3 - 3 * l4 + pi4)
                        - 0.00057 * Sin(l3 - 2 * l4 + pi4)
                        + 0.00056 * Sin(l3 + pi3 - 2 * Pi - 2 * G)
                        - 0.00052 * Sin(l2 - 2 * l3 + pi1)
                        - 0.00050 * Sin(pi2 - pi3)
                        + 0.00048 * Sin(l3 - 2 * l4 + pi3)
                        - 0.00045 * Sin(2 * l2 - 3 * l3 + pi4)
                        - 0.00041 * Sin(pi2 - pi4)
                        - 0.00038 * Sin(2 * G)
                        - 0.00037 * Sin(pi3 - pi4 + omega3 - omega4)
                        - 0.00032 * Sin(3 * l3 - 7 * l4 + 2 * pi3 + 2 * pi4)
                        + 0.00030 * Sin(4 * (l3 - l4))
                        + 0.00029 * Sin(l3 + pi4 - 2 * Pi - 2 * G)
                        - 0.00028 * Sin(omega3 + psi - 2 * Pi - 2 * G)
                        + 0.00026 * Sin(l3 - Pi - G)
                        + 0.00024 * Sin(l2 - 3 * l3 + 2 * l4)
                        + 0.00021 * Sin(2 * (l3 - Pi - G))
                        - 0.00021 * Sin(l3 - pi2)
                        + 0.00017 * Sin(2 * (l3 - pi3))).ToRadians();
            double Sigma4 = (0.84287 * Sin(l4 - pi4)
                        + 0.03431 * Sin(pi4 - pi3)
                        - 0.03305 * Sin(2 * (psi - Pi))
                        - 0.03211 * Sin(G)
                        - 0.01862 * Sin(l4 - pi3)
                        + 0.01186 * Sin(psi - omega4)
                        + 0.00623 * Sin(l4 + pi4 - 2 * G - 2 * Pi)
                        + 0.00387 * Sin(2 * (l4 - pi4))
                        - 0.00284 * Sin(5 * GPrime - 2 * G + 52.225.ToRadians())
                        - 0.00234 * Sin(2 * (psi - pi4))
                        - 0.00223 * Sin(l3 - l4)
                        - 0.00208 * Sin(l4 - Pi)
                        + 0.00178 * Sin(psi + omega4 - 2 * pi4)
                        + 0.00134 * Sin(pi4 - Pi)
                        + 0.00125 * Sin(2 * (l4 - G - Pi))
                        - 0.00117 * Sin(2 * G)
                        - 0.00112 * Sin(2 * (l3 - l4))
                        + 0.00107 * Sin(3 * l3 - 7 * l4 + 4 * pi4)
                        + 0.00102 * Sin(l4 - G - Pi)
                        + 0.00096 * Sin(2 * l4 - psi - omega4)
                        + 0.00087 * Sin(2 * (psi - omega4))
                        - 0.00085 * Sin(3 * l3 - 7 * l4 + pi3 + 3 * pi4)
                        + 0.00085 * Sin(l3 - 2 * l4 + pi4)
                        - 0.00081 * Sin(2 * (l4 - psi))
                        + 0.00071 * Sin(l4 + pi4 - 2 * Pi - 3 * G)
                        + 0.00061 * Sin(l1 - l4)
                        - 0.00056 * Sin(psi - omega3)
                        - 0.00054 * Sin(l3 - 2 * l4 + pi3)
                        + 0.00051 * Sin(l2 - l4)
                        + 0.00042 * Sin(2 * (psi - G - Pi))
                        + 0.00039 * Sin(2 * (pi4 - omega4))
                        + 0.00036 * Sin(psi + Pi - pi4 - omega4)
                        + 0.00035 * Sin(2 * GPrime - G + 188.37.ToRadians())
                        - 0.00035 * Sin(l4 - pi4 + 2 * Pi - 2 * psi)
                        - 0.00032 * Sin(l4 + pi4 - 2 * Pi - G)
                        + 0.00030 * Sin(2 * GPrime - 2 * G + 149.15.ToRadians())
                        + 0.00029 * Sin(3 * l3 - 7 * l4 + 2 * pi3 + 2 * pi4)
                        + 0.00028 * Sin(l4 - pi4 + 2 * psi - 2 * Pi)
                        - 0.00028 * Sin(2 * (l4 - omega4))
                        - 0.00027 * Sin(pi3 - pi4 + omega3 - omega4)
                        - 0.00026 * Sin(5 * GPrime - 3 * G + 188.37.ToRadians())
                        + 0.00025 * Sin(omega4 - omega3)
                        - 0.00025 * Sin(l2 - 3 * l3 + 2 * l4)
                        - 0.00023 * Sin(3 * (l3 - l4))
                        + 0.00021 * Sin(2 * l4 - 2 * Pi - 3 * G)
                        - 0.00021 * Sin(2 * l3 - 3 * l4 + pi4)
                        + 0.00019 * Sin(l4 - pi4 - G)
                        - 0.00019 * Sin(2 * l4 - pi3 - pi4)
                        - 0.00018 * Sin(l4 - pi4 + G)
                        - 0.00016 * Sin(l4 + pi3 - 2 * Pi - 2 * G)).ToRadians();
            double[] L = new double[4];
            L[0] = l1 + Sigma1;
            L[1] = l2 + Sigma2;
            L[2] = l3 + Sigma3;
            L[3] = l4 + Sigma4;
            double[] B = new double[4];
            B[0] = Atan(0.0006393 * Sin(L[0] - omega1)
                    + 0.0001825 * Sin(L[0] - omega2)
                    + 0.0000329 * Sin(L[0] - omega3)
                    - 0.0000311 * Sin(L[0] - psi)
                    + 0.0000093 * Sin(L[0] - omega4)
                    + 0.0000075 * Sin(3 * L[0] - 4 * l2 - 1.9927 * Sigma1 + omega2)
                    + 0.0000046 * Sin(L[0] + psi - 2 * Pi - 2 * G));
            B[1] = Atan(0.0081004 * Sin(L[1] - omega2)
                    + 0.0004512 * Sin(L[1] - omega3)
                    - 0.0003284 * Sin(L[1] - psi)
                    + 0.0001160 * Sin(L[1] - omega4)
                    + 0.0000272 * Sin(l1 - 2 * l3 + 1.0146 * Sigma2 + omega2)
                    - 0.0000144 * Sin(L[1] - omega1)
                    + 0.0000143 * Sin(L[1] + psi - 2 * Pi - 2 * G)
                    + 0.0000035 * Sin(L[1] - psi + G)
                    - 0.0000028 * Sin(l1 - 2 * l3 + 1.0146 * Sigma2 + omega3));
            B[2] = Atan(0.0032402 * Sin(L[2] - omega3)
                    - 0.0016911 * Sin(L[2] - psi)
                    + 0.0006847 * Sin(L[2] - omega4)
                    - 0.0002797 * Sin(L[2] - omega2)
                    + 0.0000321 * Sin(L[2] + psi - 2 * Pi - 2 * G)
                    + 0.0000051 * Sin(L[2] - psi + G)
                    - 0.0000045 * Sin(L[2] - psi - G)
                    - 0.0000045 * Sin(L[2] + psi - 2 * Pi)
                    + 0.0000037 * Sin(L[2] + psi - 2 * Pi - 3 * G)
                    + 0.0000030 * Sin(2 * l2 - 3 * L[2] + 4.03 * Sigma3 + omega2)
                    - 0.0000021 * Sin(2 * l2 - 3 * L[2] + 4.03 * Sigma3 + omega3));
            B[3] = Atan(-0.0076579 * Sin(L[3] - psi)
                    + 0.0044134 * Sin(L[3] - omega4)
                    - 0.0005112 * Sin(L[3] - omega3)
                    + 0.0000773 * Sin(L[3] + psi - 2 * Pi - 2 * G)
                    + 0.0000104 * Sin(L[3] - psi + G)
                    - 0.0000102 * Sin(L[3] - psi - G)
                    + 0.0000088 * Sin(L[3] + psi - 2 * Pi - 3 * G)
                    - 0.0000038 * Sin(L[3] + psi - 2 * Pi - G));
            double SigmaR1 = -0.0041339 * Cos(2 * (l1 - l2))
                        - 0.0000387 * Cos(l1 - pi3)
                        - 0.0000214 * Cos(l1 - pi4)
                        + 0.0000170 * Cos(l1 - l2)
                        - 0.0000131 * Cos(4 * (l1 - l2))
                        + 0.0000106 * Cos(l1 - l3)
                        - 0.0000066 * Cos(l1 + pi3 - 2 * Pi - 2 * G);
            double SigmaR2 = 0.0093848 * Cos(l1 - l2)
                        - 0.0003116 * Cos(l2 - pi3)
                        - 0.0001744 * Cos(l2 - pi4) 
                        - 0.0001442 * Cos(l2 - pi2)
                        + 0.0000553 * Cos(l2 - l3)
                        + 0.0000523 * Cos(l1 - l3)
                        - 0.0000290 * Cos(2 * (l1 - l2))
                        + 0.0000164 * Cos(2 * (l2 - omega2))
                        + 0.0000107 * Cos(l1 - 2 * l3 + pi3)
                        - 0.0000102 * Cos(l2 - pi1)
                        - 0.0000091 * Cos(2 * (l1 - l3));
            double SigmaR3 = -0.0014388 * Cos(l3 - pi3)
                        - 0.0007919 * Cos(l3 - pi4)
                        + 0.0006342 * Cos(l2 - l3)
                        - 0.0001761 * Cos(2 * (l3 - l4))
                        + 0.0000294 * Cos(l3 - l4)
                        - 0.0000156 * Cos(3 * (l3 - l4))
                        + 0.0000156 * Cos(l1 - l3)
                        - 0.0000153 * Cos(l1 - l2)
                        + 0.0000070 * Cos(2 * l2 - 3 * l3 + pi3)
                        - 0.0000051 * Cos(l3 + pi3 - 2 * Pi - 2 * G);
            double SigmaR4 = -0.0073546 * Cos(l4 - pi4)
                        + 0.0001621 * Cos(l4 - pi3)
                        + 0.0000974 * Cos(l3 - l4)
                        - 0.0000543 * Cos(l4 + pi4 - 2 * Pi - 2 * G)
                        - 0.0000271 * Cos(2 * (l4 - pi4))
                        + 0.0000182 * Cos(l4 - Pi)
                        + 0.0000177 * Cos(2 * (l3 - l4))
                        - 0.0000167 * Cos(2 * l4 - psi - omega4)
                        + 0.0000167 * Cos(psi - omega4)
                        - 0.0000155 * Cos(2 * (l4 - Pi - G))
                        + 0.0000142 * Cos(2 * (l4 - psi))
                        + 0.0000105 * Cos(l1 - l4)
                        + 0.0000092 * Cos(l2 - l4)
                        - 0.0000089 * Cos(l4 - Pi - G)
                        - 0.0000062 * Cos(l4 + pi4 - 2 * Pi - 3 * G)
                        + 0.0000048 * Cos(2 * (l4 - omega4));
            double[] R = new double[4];
            R[0] = 5.90569 * (1 + SigmaR1);
            R[1] = 9.39657 * (1 + SigmaR2);
            R[2] = 14.98832 * (1 + SigmaR3);
            R[3] = 26.36273 * (1 + SigmaR4);
            double T0 = (julianEphemerisDay - 2433282.423) / 36525;
            double P = (1.3966626 * T0 + 0.0003088 * T0 * T0).ToRadians();
            L[0] += P;
            L[1] += P;
            L[2] += P;
            L[3] += P;
            psi += P;
            double I = (3.120262 + 0.0006 * Time.JulianCentury(julianEphemerisDay, 1900.0)).ToRadians();
            var orb = Body.OrbitalElements(Bodies.Jupiter, julianEphemerisDay);
            double Omega = orb.longitudeAscendingNode.ToRadians();
            double Phi = psi - Omega;
            double i = orb.inclination.ToRadians();
            var fict = ABC(0, 0, 1);
            double D = Atan2(fict.A, fict.C);
            (double x, double y, double z)[] result = new (double, double, double)[4];
            for (int j = 0; j < 4; j++) {
                double x = R[j] * Cos(L[j] - psi) * Cos(B[j]);
                double y = R[j] * Sin(L[j] - psi) * Cos(B[j]);
                double z = R[j] * Sin(B[j]);
                var abc = ABC(x, y, z);
                x = abc.A * Cos(D) - abc.C * Sin(D);
                y = abc.A * Sin(D) + abc.C * Cos(D);
                z = abc.B;
                int K = 0;
                switch(j) {
                    case 0: K = 17295; break;
                    case 1: K = 21819; break;
                    case 2: K = 27558; break;
                    case 3: K = 36548; break;
                }
                x += Abs(z) / K * Sqrt(1 - Pow(x / R[j], 2));
                double w = poz.delta / (poz.delta + z / 2095);
                x *= w;
                y *= w;
                result[j] = (x, y, z);
            }
            return result;

            (double A, double B, double C) ABC(double x, double y, double z) {
                double B1 = y * Cos(I) - z * Sin(I);
                double C1 = y * Sin(I) + z * Cos(I);
                double A2 = x * Cos(Phi) - B1 * Sin(Phi);
                double B2 = x * Sin(Phi) + B1 * Cos(Phi);
                double B3 = B2 * Cos(i) - C1 * Sin(i);
                double C3 = B2 * Sin(i) + C1 * Cos(i);
                double A4 = A2 * Cos(Omega) - B3 * Sin(Omega);
                double B4 = A2 * Sin(Omega) + B3 * Cos(Omega);
                double A5 = A4 * Sin(poz.lambda) - B4 * Cos(poz.lambda);
                double B5 = A4 * Cos(poz.lambda) + B4 * Sin(poz.lambda);
                double B6 = C3 * Sin(poz.beta) + B5 * Cos(poz.beta);
                double C6 = C3 * Cos(poz.beta) - B5 * Sin(poz.beta);
                return (A5, B6, C6);
            }
        }
    }
}
