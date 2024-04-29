using System;
using System.Data;
using static System.Math;

namespace AA.Net {
    public static class Moon {

        public enum Phase {
            NewMoon,
            FirstQuarter,
            FullMoon,
            LastQuarter
        }

        public enum Apsis {
            Perigee,
            Apogee
        }

        /// <remarks>p. 407</remarks>
        public const int Diameter = 3476;

        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay, bool apparent = true) {
            return Body.Position(Bodies.Moon, julianEphemerisDay, apparent);
        }

        /// <remarks>p. 338</remarks>
        internal static double EarthEccentricity(double T) {
            return 1 - 0.002516 * T - 0.0000074 * T * T;
        }

        /// <remarks>p. 338</remarks>
        private static double ArgumentOfLatitude(double T) {
            return (93.2720950 + 483202.0175233 * T - 0.0036539 * T * T - T * T * T / 3526000 + T * T * T * T / 863310000).To360();
        }
        
        /// <remarks>p. 338</remarks>
        internal static (double D, double M, double Mprime, double F) DMMF(double T, bool radians) {
            double T2 = T * T;
            double T3 = T2 * T;
            double T4 = T3 * T;
            double D = (297.8501921 + 445267.1114034 * T - 0.0018819 * T2 + T3 / 545868 - T4 / 113065000).To360();
            double M = (357.5291092 + 35999.0502909 * T - 0.0001536 * T2 + T3 / 24490000).To360();
            double Mprime = (134.9633964 + 477198.8675055 * T + 0.0087414 * T2 + T3 / 69699 - T4 / 14712000).To360();
            double F = ArgumentOfLatitude(T);
            if (radians) {
                D = D.ToRadians();
                M = M.ToRadians();
                Mprime = Mprime.ToRadians();
                F = F.ToRadians();
            }
            return (D, M, Mprime, F);
        }

        /// <param name="distance">Earth-moon distance in km.</param>
        /// <returns>Equatorial horizontal parallax of the moon in degrees.</returns>
        /// <remarks>p. 337</remarks>
        public static double Parallax(double distance) {
            return Asin(6378.14 / distance).ToDegrees();
        }
        
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <param name="apparent">True (default) for apparent node or False for mean node.</param>
        /// <returns>Longitude of the ascending node in degrees.</returns>
        /// <remarks>p. 343</remarks>
        public static double AscendingNode(double julianEphemerisDay, bool apparent = true) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = Pow(T, 2);
            double T3 = Pow(T, 3);
            double T4 = Pow(T, 4);
            double omega = 125.0445479 - 1934.1362891 * T + 0.0020754 * T2 + T3 / 467441 - T4 / 60616000;
            if (apparent) {
                (double D, double M, double Mprime, double F) = DMMF(T, true);
                omega = omega - 1.4979 * Sin(2 * (D - F)) - 0.1500 * Sin(M) - 0.1226 * Sin(2 * D) + 0.1176 * Sin(2 * F) - 0.0801 * Sin(2 * (Mprime - F));
            }
            return omega;
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Longitude of the perigee in degrees.</returns>
        /// <remarks>p. 343</remarks>
        public static double PerigeeLongitude(double julianEphemerisDay) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            return 83.3532465 + 4069.0137287 * T - 0.0103200 * Pow(T, 2) - Pow(T, 3) / 80053 + Pow(T, 4) / 18999000;
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Illuminated fraction of the disk and position angle of the bright limb in degrees.</returns>
        /// <remarks>pp. 345-346</remarks>
        public static (double fraction, double positionAngle)
        Illumination(double julianEphemerisDay) {
            var sun = Body.Position(Bodies.Sun, julianEphemerisDay);
            var moon = Position(julianEphemerisDay);
            double epsilon = Sky.ObliquityOfEcliptic(julianEphemerisDay);
            var sunP = Transform.ToEquatorial(sun.longitude, 0, epsilon);
            double alpha0 = sunP.rightAscension.ToRadians();
            double delta0 = sunP.declination.ToRadians();
            var moonP = Transform.ToEquatorial(moon.longitude, moon.latitude, epsilon);
            double alpha = moonP.rightAscension.ToRadians();
            double delta = moonP.declination.ToRadians();
            double psi = Acos(Sin(delta0) * Sin(delta) + Cos(delta0) * Cos(delta) * Cos(alpha0 - alpha));
            double i = Atan2(sun.distance * Sin(psi), moon.distance / Sun.AUKilometers - sun.distance * Cos(psi));
            double k = (1 + Cos(i)) / 2;
            double chi = Body.BrightLimbPositionAngle(alpha0, delta0, alpha, delta);
            return (k, chi.To360());
        }

        internal const double NewMoon2000 = 2451550.09766;

        /// <param name="k">Lunation from J2000.</param>
        /// <returns>Time in Julian Centuries.</returns>
        /// <remarks>p. 350</remarks>
        internal static double JulianCentury(double k) {
            return k / 1236.85;
        }
        
        /// <param name="julianDay">Julian Day.</param>
        /// <returns>Lunation from first moon of J2000.</returns>
        /// <remarks>p. 350, revised.</remarks>
        internal static double Lunation(double julianDay) {
            return (julianDay.ToFractionalYear() - 2000.018) * 12.3685;
        }

        private static double NextQuarter(double julianDay) {
            return Ceiling(Lunation(julianDay) * 4) / 4;
        }

        /// <param name="start">Julian Day at which to begin the search.</param>
        /// <returns>Next time and phase.</returns>
        public static (double julianEphemerisDay, Phase phase)
        NextPhase(double start) {
            //Find next quarter phase
            double k = NextQuarter(start);
            Phase next = (Phase)((k - Floor(k)) * 4);
            return (NextPhase(next, start), next);
        }
        
        /// <param name="phase">A value from the <see cref="Phase"/> enumeration.</param>
        /// <param name="start">Julian Day at which to begin the search.</param>
        /// <returns>Next time of requested phase.</returns>
        /// <remarks>pp. 349-352</remarks>
        public static double NextPhase(Phase phase, double start) {
            //Round up start to requested quarter
            double k = Lunation(start);
            double frac = (double)phase / 4;
            k = Ceiling(k - frac) + frac;
            //Meeus magic starts here
            double T = JulianCentury(k);
            double T2 = T * T;
            double T3 = T2 * T;
            double T4 = T3 * T;
            double JDE = NewMoon2000 + 29.530588861 * k
                                       +  0.00015437 * T2
                                       -  0.000000150 * T3
                                       +  0.00000000073 * T4;
            double E = EarthEccentricity(T);
            double M = (2.5534 + 29.1053567 * k
                               -  0.0000014 * T2
                               -  0.00000011 * T3)
                               .To360().ToRadians();
            double mPrime = (201.5643 + 385.81693528 * k
                                      +   0.0107582 * T2
                                      +   0.00001238 * T3
                                      -   0.000000058 * T4)
                                      .To360().ToRadians();
            double F = (160.7108 + 390.67050284 * k
                                 -   0.0016118 * T2
                                 -   0.00000227 * T3
                                 +   0.000000011 * T4)
                                 .To360().ToRadians();
            double omega = (124.7746 - 1.56375588 * k
                                     + 0.0020672 * T2
                                     + 0.00000215 * T3)
                                     .To360().ToRadians();
            double A1 = (299.77 + 0.107408 * k - 0.009173 * T2).ToRadians();
            double A2 = (251.88 + 0.016321 * k).ToRadians();
            double A3 = (251.83 + 26.651886 * k).ToRadians();
            double A4 = (349.42 + 36.412478 * k).ToRadians();
            double A5 = (84.66 + 18.206239 * k).ToRadians();
            double A6 = (141.74 + 53.303771 * k).ToRadians();
            double A7 = (207.14 + 2.453732 * k).ToRadians();
            double A8 = (154.84 + 7.306860 * k).ToRadians();
            double A9 = (34.52 + 27.261239 * k).ToRadians();
            double A10 = (207.19 + 0.121824 * k).ToRadians();
            double A11 = (291.34 + 1.844379 * k).ToRadians();
            double A12 = (161.72 + 24.198154 * k).ToRadians();
            double A13 = (239.56 + 25.513099 * k).ToRadians();
            double A14 = (331.55 + 3.592518 * k).ToRadians();
            if (phase == Phase.NewMoon) {
                JDE += - 0.40720 * Sin(mPrime)
                       + 0.17241 * E * Sin(M)
                       + 0.01608 * Sin(2 * mPrime)
                       + 0.01039 * Sin(2 * F)
                       + 0.00739 * E * Sin(mPrime - M)
                       - 0.00514 * E * Sin(mPrime + M)
                       + 0.00208 * E * E * Sin(2 * M)
                       - 0.00111 * Sin(mPrime - 2 * F)
                       - 0.00057 * Sin(mPrime + 2 * F)
                       + 0.00056 * E * Sin(2 * mPrime + M)
                       - 0.00042 * Sin(3 * mPrime)
                       + 0.00042 * E * Sin(M + 2 * F)
                       + 0.00038 * E * Sin(M - 2 * F)
                       - 0.00024 * E * Sin(2 * mPrime - M)
                       - 0.00017 * Sin(omega)
                       - 0.00007 * Sin(mPrime + 2 * M)
                       + 0.00004 * Sin(2 * mPrime - 2 * F)
                       + 0.00004 * Sin(3 * M)
                       + 0.00003 * Sin(mPrime + M - 2 * F)
                       + 0.00003 * Sin(2 * mPrime + 2 * F)
                       - 0.00003 * Sin(mPrime + M + 2 * F)
                       + 0.00003 * Sin(mPrime - M + 2 * F)
                       - 0.00002 * Sin(mPrime - M - 2 * F)
                       - 0.00002 * Sin(3 * mPrime + M)
                       + 0.00002 * Sin(4 * mPrime);
            }
            else if (phase == Phase.FullMoon) {
                JDE += -0.40614 * Sin(mPrime)
                       + 0.17302 * E * Sin(M)
                       + 0.01614 * Sin(2 * mPrime)
                       + 0.01043 * Sin(2 * F)
                       + 0.00734 * E * Sin(mPrime - M)
                       - 0.00515 * E * Sin(mPrime + M)
                       + 0.00209 * E * E * Sin(2 * M)
                       - 0.00111 * Sin(mPrime - 2 * F)
                       - 0.00057 * Sin(mPrime + 2 * F)
                       + 0.00056 * E * Sin(2 * mPrime + M)
                       - 0.00042 * Sin(3 * mPrime)
                       + 0.00042 * E * Sin(M + 2 * F)
                       + 0.00038 * E * Sin(M - 2 * F)
                       - 0.00024 * E * Sin(2 * mPrime - M)
                       - 0.00017 * Sin(omega)
                       - 0.00007 * Sin(mPrime + 2 * M)
                       + 0.00004 * Sin(2 * mPrime - 2 * F)
                       + 0.00004 * Sin(3 * M)
                       + 0.00003 * Sin(mPrime + M - 2 * F)
                       + 0.00003 * Sin(2 * mPrime + 2 * F)
                       - 0.00003 * Sin(mPrime + M + 2 * F)
                       + 0.00003 * Sin(mPrime - M + 2 * F)
                       - 0.00002 * Sin(mPrime - M - 2 * F)
                       - 0.00002 * Sin(3 * mPrime + M)
                       + 0.00002* Sin(4 * mPrime);
            }
            else {
                JDE += -0.62801 * Sin(mPrime)
                       + 0.17172 * E * Sin(M)
                       - 0.01183 * E * Sin(mPrime + M)
                       + 0.00862 * Sin(2 * mPrime)
                       + 0.00804 * Sin(2 * F)
                       + 0.00454 * E * Sin(mPrime - M)
                       + 0.00204 * E * E * Sin(2 * M)
                       - 0.00180 * Sin(mPrime - 2 * F)
                       - 0.00070 * Sin(mPrime + 2 * F)
                       - 0.00040 * Sin(3 * mPrime)
                       - 0.00034 * E * Sin(2 * mPrime - M)
                       + 0.00032 * E * Sin(M + 2 * F)
                       + 0.00032 * E * Sin(M - 2 * F)
                       - 0.00028 * E * E * Sin(mPrime + 2 * M)
                       + 0.00027 * E * Sin(2 * mPrime + M)
                       - 0.00017 * Sin(omega)
                       - 0.00005 * Sin(mPrime - M - 2 * F)
                       + 0.00004 * Sin(2 * mPrime + 2 * F)
                       - 0.00004 * Sin(mPrime + M + 2 * F)
                       + 0.00004 * Sin(mPrime - 2 * M)
                       + 0.00003 * Sin(mPrime + M - 2 * F)
                       + 0.00003 * Sin(3 * M)
                       + 0.00002 * Sin(2 * mPrime - 2 * F)
                       + 0.00002 * Sin(mPrime - M + 2 * F)
                       - 0.00002 * Sin(3 * mPrime + M);
                double W = 0.00306 - 0.00038 * E * Cos(M) + 0.00026 * Cos(mPrime) - 0.00002 * Cos(mPrime - M) + 0.00002 * Cos(mPrime + M) + 0.00002 * Cos(2 * F);
                if (phase == Phase.FirstQuarter) JDE += W;
                else JDE -= W;
            }
            JDE +=   0.000325 * Sin(A1) + 0.000056 * Sin(A8)
                   + 0.000165 * Sin(A2) + 0.000047 * Sin(A9)
                   + 0.000164 * Sin(A3) + 0.000042 * Sin(A10)
                   + 0.000126 * Sin(A4) + 0.000040 * Sin(A11)
                   + 0.000110 * Sin(A5) + 0.000037 * Sin(A12)
                   + 0.000062 * Sin(A6) + 0.000035 * Sin(A13)
                   + 0.000060 * Sin(A7) + 0.000023 * Sin(A14);

            return JDE;
        }

        /// <param name="start">Julian Day at which to begin the search.</param>
        /// <returns>Characteristics of the next eclipse.</returns>
        /// <remarks>pp. 379-383</remarks>
        public static (EclipseType type, double greatestEclipse, double magnitudeUmbral, double magnitudePenumbral, double semidurationPartial, double semidurationTotal, double semidurationPenumbral, bool isAscending, double axisDistance, double radiusUmbra, double radiusPenumbra)
        NextEclipse(double start) {
            double k = Ceiling(Lunation(start) - 0.5) - 0.5;
            var e = Body.FindNextEclipse(k, true);
            while (e.greatestEclipse < start) e = Body.FindNextEclipse(++k, true);
            double p, t, semiDurPartial = 0, semiDurTotal = 0, semiDurPenumbral;
            double n = 60 / (0.5458 + 0.04 * Cos(e.Mprime));
            double gamma2 = e.gamma * e.gamma;
            double magUmbral = (1.0128 - e.u - Abs(e.gamma)) / 0.545;
            double magPenumbral = (1.5573 + e.u - Abs(e.gamma)) / 0.545;
            EclipseType type;
            if (magUmbral > 0) {
                p = 1.0128 - e.u;
                t = 0.4678 - e.u;
                semiDurPartial = n * Sqrt(p * p - gamma2);
                semiDurTotal = n * Sqrt(t * t - gamma2);
                type = magUmbral >= 1 ? EclipseType.Total : EclipseType.Partial;
            }
            else {
                type = EclipseType.Penumbral;
            }
            double h = 1.5573 + e.u;
            semiDurPenumbral = n * Sqrt(h * h - gamma2);
            return (type, e.greatestEclipse, magUmbral, magPenumbral, semiDurPartial, semiDurTotal, semiDurPenumbral, e.isAscending, e.gamma, e.rho, e.sigma);
        }

        /// <param name="start">A Julian Day in dynamical time.</param>
        /// <param name="type">A value from the <see cref="Apsis"/> enumeration.</param>
        /// <returns>Time and parallax angle of the next apsis of requested type.</returns>
        /// <remarks>pp. 355-360</remarks>
        private static (double julianEphemerisDay, double parallax)
        NextApsis(double start, Apsis type) {
            double k = (start.ToFractionalYear() - 1999.97) * 13.2555;
            if (type == Apsis.Perigee) k = Ceiling(k);
            else k = Ceiling(k - 0.5) + 0.5;
            double T = k / 1325.55;
            double T2 = T * T;
            double T3 = T2 * T;
            double T4 = T3 * T;
            double JDE = 2451534.6698
                + 27.55454989 * k
                - 0.0006691 * T2
                - 0.000001098 * T3
                + 0.0000000052 * T4;
            double D = (171.9179
                + 335.9106046 * k
                - 0.0100383 * T2
                - 0.00001156 * T3
                + 0.000000055 * T4).ToRadians();
            double M = (347.3477
                + 27.1577721 * k
                - 0.0008130 * T2
                - 0.0000010 * T3).ToRadians();
            double F = (316.6109
                + 364.5287911 * k
                - 0.0125053 * T2
                - 0.0000148 * T3).ToRadians();
            double sumTerms, parallax;
            if (type == Apsis.Perigee) {
                sumTerms =
                      Sin(2 * D) * -1.6769
                    + Sin(4 * D) * 0.4589
                    + Sin(6 * D) * -0.1856
                    + Sin(8 * D) * 0.0883
                    + Sin(2 * D - M) * (-0.0773 + 0.00019 * T)
                    + Sin(M) * (0.0502 - 0.00013 * T)
                    + Sin(10 * D) * -0.0460
                    + Sin(4 * D - M) * (0.0422 - 0.00011 * T)
                    + Sin(6 * D - M) * -0.0256
                    + Sin(12 * D) * 0.0253
                    + Sin(D) * 0.0237
                    + Sin(8 * D - M) * 0.0162
                    + Sin(14 * D) * -0.0145
                    + Sin(2 * F) * 0.0129
                    + Sin(3 * D) * -0.0112
                    + Sin(10 * D - M) * -0.0104
                    + Sin(16 * D) * 0.0086
                    + Sin(12 * D - M) * 0.0069
                    + Sin(5 * D) * 0.0066
                    + Sin(2 * D + 2 * F) * -0.0053
                    + Sin(18 * D) * -0.0052
                    + Sin(14 * D - M) * -0.0046
                    + Sin(7 * D) * -0.0041
                    + Sin(2 * D + M) * 0.0040
                    + Sin(20 * D) * 0.0032
                    + Sin(D + M) * -0.0032
                    + Sin(16 * D - M) * 0.0031
                    + Sin(4 * D + M) * -0.0029
                    + Sin(9 * D) * 0.0027
                    + Sin(4 * D + 2 * F) * 0.0027
                    + Sin(2 * D - 2 * M) * -0.0027
                    + Sin(4 * D - 2 * M) * 0.0024
                    + Sin(6 * D - 2 * M) * -0.0021
                    + Sin(22 * D) * -0.0021
                    + Sin(18 * D - M) * -0.0021
                    + Sin(6 * D + M) * 0.0019
                    + Sin(11 * D) * -0.0018
                    + Sin(8 * D + M) * -0.0014
                    + Sin(4 * D - 2 * F) * -0.0014
                    + Sin(6 * D + 2 * F) * -0.0014
                    + Sin(3 * D + M) * 0.0014
                    + Sin(5 * D + M) * -0.0014
                    + Sin(13 * D) * 0.0013
                    + Sin(20 * D - M) * 0.0013
                    + Sin(3 * D + 2 * M) * 0.0011
                    + Sin(4 * D + 2 * F - 2 * M) * -0.0011
                    + Sin(D + 2 * M) * -0.0010
                    + Sin(22 * D - M) * -0.0009
                    + Sin(4 * F) * -0.0008
                    + Sin(6 * D - 2 * F) * 0.0008
                    + Sin(2 * D - 2 * F + M) * 0.0008
                    + Sin(2 * M) * 0.0007
                    + Sin(2 * F - M) * 0.0007
                    + Sin(2 * D + 4 * F) * 0.0007
                    + Sin(2 * F - 2 * M) * -0.0006
                    + Sin(2 * D - 2 * F + 2 * M) * -0.0006
                    + Sin(24 * D) * 0.0006
                    + Sin(4 * D - 4 * F) * 0.0005
                    + Sin(2 * D + 2 * M) * 0.0005
                    + Sin(D - M) * -0.0004;
                parallax = 
                      3629.215
                    + 63.224 * Cos(2 * D)
                    - 6.990 * Cos(4 * D)
                    + 2.834 * Cos(2 * D - M)
                    - 0.0071 * T * Cos(2 * D - M)
                    + 1.927 * Cos(6 * D)
                    - 1.263 * Cos(D)
                    - 0.702 * Cos(8 * D)
                    + 0.696 * Cos(M)
                    - 0.0017 * T * Cos(M)
                    - 0.690 * Cos(2 * F)
                    - 0.629 * Cos(4 * D - M)
                    + 0.0016 * T * Cos(4 * D - M)
                    - 0.392 * Cos(2 * D - 2 * F)
                    + 0.297 * Cos(10 * D)
                    + 0.260 * Cos(6 * D - M)
                    + 0.201 * Cos(3 * D)
                    - 0.161 * Cos(2 * D + M)
                    + 0.157 * Cos(D + M)
                    - 0.138 * Cos(12 * D)
                    - 0.127 * Cos(8 * D - M)
                    + 0.104 * Cos(2 * D + 2 * F)
                    + 0.104 * Cos(2 * D - 2 * M)
                    - 0.079 * Cos(5 * D)
                    + 0.068 * Cos(14 * D)
                    + 0.067 * Cos(10 * D - M)
                    + 0.054 * Cos(4 * D + M)
                    - 0.038 * Cos(12 * D - M)
                    - 0.038 * Cos(4 * D - 2 * M)
                    + 0.037 * Cos(7 * D)
                    - 0.037 * Cos(4 * D + 2 * F)
                    - 0.035 * Cos(16 * D)
                    - 0.030 * Cos(3 * D + M)
                    + 0.029 * Cos(D - M)
                    - 0.025 * Cos(6 * D + M)
                    + 0.023 * Cos(2 * M)
                    + 0.023 * Cos(14 * D - M)
                    - 0.023 * Cos(2 * D + 2 * M)
                    + 0.022 * Cos(6 * D - 2 * M)
                    - 0.021 * Cos(2 * D - 2 * F - M)
                    - 0.020 * Cos(9 * D)
                    + 0.019 * Cos(18 * D)
                    + 0.017 * Cos(6 * D + 2 * F)
                    + 0.014 * Cos(2 * F - M)
                    - 0.014 * Cos(16 * D - M)
                    + 0.013 * Cos(4 * D - 2 * F)
                    + 0.012 * Cos(8 * D + M)
                    + 0.011 * Cos(11 * D)
                    + 0.010 * Cos(5 * D + M)
                    - 0.010 * Cos(20 * D);
            }
            else {
                sumTerms =
                Sin(2 * D) * 0.4392
              + Sin(4 * D) * 0.0684
              + Sin(M) * (0.0456 - 0.00011 * T)
              + Sin(2 * D - M) * (0.0426 - 0.00011 * T)
              + Sin(2 * F) * 0.0212
              + Sin(D) * -0.0189
              + Sin(6 * D) * 0.0144
              + Sin(4 * D - M) * 0.0113
              + Sin(2 * D + 2 * F) * 0.0047
              + Sin(D + M) * 0.0036
              + Sin(8 * D) * 0.0035
              + Sin(6 * D - M) * 0.0034
              + Sin(2 * D - 2 * F) * -0.0034
              + Sin(2 * D - 2 * M) * 0.0022
              + Sin(3 * D) * -0.0017
              + Sin(4 * D + 2 * F) * 0.0013
              + Sin(8 * D - M) * 0.0011
              + Sin(4 * D - 2 * M) * 0.0010
              + Sin(10 * D) * 0.0009
              + Sin(3 * D + M) * 0.0007
              + Sin(2 * M) * 0.0006
              + Sin(2 * D + M) * 0.0005
              + Sin(2 * D + 2 * M) * 0.0005
              + Sin(6 * D + 2 * F) * 0.0004
              + Sin(6 * D - 2 * M) * 0.0004
              + Sin(10 * D - M) * 0.0004
              + Sin(5 * D) * -0.0004
              + Sin(4 * D - 2 * F) * -0.0004
              + Sin(2 * F + M) * 0.0003
              + Sin(12 * D) * 0.0003
              + Sin(2 * D + 2 * F - M) * 0.0003
              + Sin(D - M) * -0.0003;
                parallax =
                      3245.251
                    - 9.147 * Cos(2 * D)
                    - 0.841 * Cos(D)
                    + 0.697 * Cos(2 * F)
                    - 0.656 * Cos(M)
                    + 0.0016 * T * Cos(M)
                    + 0.355 * Cos(4 * D)
                    + 0.159 * Cos(2 * D - M)
                    + 0.127 * Cos(D + M)
                    + 0.065 * Cos(4 * D - M)
                    + 0.052 * Cos(6 * D)
                    + 0.043 * Cos(2 * D + M)
                    + 0.031 * Cos(2 * D + 2 * F)
                    - 0.023 * Cos(2 * D - 2 * F)
                    + 0.022 * Cos(2 * D - 2 * M)
                    + 0.019 * Cos(2 * D + 2 * M)
                    - 0.016 * Cos(2 * M)
                    + 0.014 * Cos(6 * D - M)
                    + 0.010 * Cos(8 * D);
            }
            return (JDE + sumTerms, parallax);
        }

        public static (double julianEphemerisDay, double parallax) NextPerigee(double start) {
            return NextApsis(start, Apsis.Perigee);
        }

        public static (double julianEphemerisDay, double parallax) NextApogee(double start) {
            return NextApsis(start, Apsis.Apogee);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Semidiameters in arcseconds.</returns>
        /// <remarks>p. 390</remarks>
        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            double pi = Asin(6378.14 / Position(julianEphemerisDay).distance);
            double s = Asin(0.272481 * Sin(pi)).ToDegrees() * 3600;
            return (s, s);
        }

        public enum NodeType {
            Ascending,
            Descending
        }

        /// <param name="node">A value from the <see cref="NodeType"/> enumeration.</param>
        /// <param name="start">Julian Day at which to begin the search.</param>
        /// <returns>Time of next occurrence of requested node.</returns>
        /// <remarks>pp. 363-364</remarks>
        public static double NextNode(NodeType node, double start) {
            double k = (start.ToFractionalYear() - 2000.05) * 13.4223;
            if (node == NodeType.Ascending) k = Ceiling(k);
            else k = Round(k, MidpointRounding.AwayFromZero) + 0.5;
            double T = k / 1342.23;
            double T2 = T * T;
            double T3 = T2 * T;
            double T4 = T3 * T;
            double D = (183.638 + 331.73735682 * k + 0.0014852 * T2 + 0.00000209 * T3 - 0.00000001 * T4).ToRadians();
            double M = (17.4006 + 26.8203725 * k + 0.0001186 * T2 + 0.00000006 * T3).ToRadians();
            double Mprime = (38.3776 + 355.52747313 * k + 0.0123499 * T2 + 0.000014627 * T3 - 0.000000069 * T4).ToRadians();
            double Omega = (123.9767 - 1.44098956 * k + 0.0020608 * T2 + 0.00000214 * T3 - 0.000000016 * T4).ToRadians();
            double V = (299.75 + 132.85 * T - 0.009173 * T2).ToRadians();
            double P = (Omega + 272.75 - 2.3 * T).ToRadians();
            double E = EarthEccentricity(T);
            double JDE = 2451565.1619 + 27.212220817 * k
                        + 0.0002762 * T2
                        + 0.000000021 * T3
                        - 0.000000000088 * T4
                        - 0.4721 * Sin(Mprime)
                        - 0.1649 * Sin(2 * D)
                        - 0.0868 * Sin(2 * D - Mprime)
                        + 0.0084 * Sin(2 * D + Mprime)
                        - 0.0083 * Sin(2 * D - M) * E
                        - 0.0039 * Sin(2 * D - M - Mprime) * E
                        + 0.0034 * Sin(2 * Mprime)
                        - 0.0031 * Sin(2 * D - 2 * Mprime)
                        + 0.0030 * Sin(2 * D + M) * E
                        + 0.0028 * Sin(M - Mprime) * E
                        + 0.0026 * Sin(M) * E
                        + 0.0025 * Sin(4 * D)
                        + 0.0024 * Sin(D)
                        + 0.0022 * Sin(M + Mprime) * E
                        + 0.0017 * Sin(Omega)
                        + 0.0014 * Sin(4 * D - Mprime)
                        + 0.0005 * Sin(2 * D + M - Mprime) * E
                        + 0.0004 * Sin(2 * D - M + Mprime) * E
                        - 0.0003 * Sin(2 * D - 2 * M) * E
                        + 0.0003 * Sin(4 * D - M) * E
                        + 0.0003 * Sin(V)
                        + 0.0003 * Sin(P);
            return JDE;
        }

        public enum DeclinationType {
            Northern,
            Southern
        }

        /// <param name="direction">A value of the <see cref="DeclinationType"/> enumeration.</param>
        /// <param name="start">Julian Day at which to begin the search.</param>
        /// <returns>Time and value of next maximum declination.</returns>
        /// <remarks>pp. 367-369</remarks>
        public static (double julianEphemerisDay, double declination)
        NextMaximumDeclination (DeclinationType direction, double start) {
            double k = Ceiling((start.ToFractionalYear() - 2000.03) * 13.3686);
            double T = k / 1336.86;
            double T2 = T * T;
            double T3 = T2 * T;
            double D = ((direction == DeclinationType.Northern ? 152.2029 : 345.6676) + 333.0705546 * k - 0.0004214 * T2 + 0.00000011 * T3).ToRadians();
            double M = ((direction == DeclinationType.Northern ? 14.8591 : 1.3951) + 26.9281592 * k - 0.0000355 * T2 - 0.0000001 * T3).ToRadians();
            double Mprime = ((direction == DeclinationType.Northern ? 4.6881 : 186.21) + 356.9562794 * k + 0.0103066 * T2 + 0.00001251 * T3).ToRadians();
            double F = ((direction == DeclinationType.Northern ? 325.8867 : 145.1633) + 1.4467807 * k - 0.002069 * T2 - 0.00000215 * T3).ToRadians();
            double E = EarthEccentricity(T);
            double timeSum = direction == DeclinationType.Northern ?
                + 0.8975 * Cos(F)
                - 0.4726 * Sin(Mprime)
                - 0.1030 * Sin(2 * F)
                - 0.0976 * Sin(2 * D - Mprime)
                - 0.0462 * Cos(Mprime - F)
                - 0.0461 * Cos(Mprime + F)
                - 0.0438 * Sin(2 * D)
                + 0.0162 * Sin(M) * E
                - 0.0157 * Cos(3 * F)
                + 0.0145 * Sin(Mprime + 2 * F)
                + 0.0136 * Cos(2 * D - F)
                - 0.0095 * Cos(2 * D - Mprime - F)
                - 0.0091 * Cos(2 * D - Mprime + F)
                - 0.0089 * Cos(2 * D + F)
                + 0.0075 * Sin(2 * Mprime)
                - 0.0068 * Sin(Mprime - 2 * F)
                + 0.0061 * Cos(2 * Mprime - F)
                - 0.0047 * Sin(Mprime + 3 * F)
                - 0.0043 * Sin(2 * D - M - Mprime) * E
                - 0.0040 * Cos(Mprime - 2 * F)
                - 0.0037 * Sin(2 * D - 2 * Mprime)
                + 0.0031 * Sin(F)
                + 0.0030 * Sin(2 * D + Mprime)
                - 0.0029 * Cos(Mprime + 2 * F)
                - 0.0029 * Sin(2 * D - M) * E
                - 0.0027 * Sin(Mprime + F)
                + 0.0024 * Sin(M - Mprime) * E
                - 0.0021 * Sin(Mprime - 3 * F)
                + 0.0019 * Sin(2 * Mprime + F)
                + 0.0018 * Cos(2 * D - 2 * Mprime - F)
                + 0.0018 * Sin(3 * F)
                + 0.0017 * Cos(Mprime + 3 * F)
                + 0.0017 * Cos(2 * Mprime)
                - 0.0014 * Cos(2 * D - Mprime)
                + 0.0013 * Cos(2 * D + Mprime + F)
                + 0.0013 * Cos(Mprime)
                + 0.0012 * Sin(3 * Mprime + F)
                + 0.0011 * Sin(2 * D - Mprime + F)
                - 0.0011 * Cos(2 * D - 2 * Mprime)
                + 0.0010 * Cos(D + F)
                + 0.0010 * Sin(M + Mprime) * E
                - 0.0009 * Sin(2 * D - 2 * F)
                + 0.0007 * Cos(2 * Mprime + F)
                - 0.0007 * Cos(3 * Mprime + F)
                :
                - 0.8975 * Cos(F)
                - 0.4726 * Sin(Mprime)
                - 0.1030 * Sin(2 * F)
                - 0.0976 * Sin(2 * D - Mprime)
                + 0.0541 * Cos(Mprime - F)
                + 0.0516 * Cos(Mprime + F)
                - 0.0438 * Sin(2 * D)
                + 0.0112 * Sin(M) * E
                + 0.0157 * Cos(3 * F)
                + 0.0023 * Sin(Mprime + 2 * F)
                - 0.0136 * Cos(2 * D - F)
                + 0.0110 * Cos(2 * D - Mprime - F)
                + 0.0091 * Cos(2 * D - Mprime + F)
                + 0.0089 * Cos(2 * D + F)
                + 0.0075 * Sin(2 * Mprime)
                - 0.0030 * Sin(Mprime - 2 * F)
                - 0.0061 * Cos(2 * Mprime - F)
                - 0.0047 * Sin(Mprime + 3 * F)
                - 0.0043 * Sin(2 * D - M - Mprime) * E
                + 0.0040 * Cos(Mprime - 2 * F)
                - 0.0037 * Sin(2 * D - 2 * Mprime)
                - 0.0031 * Sin(F)
                + 0.0030 * Sin(2 * D + Mprime)
                + 0.0029 * Cos(Mprime + 2 * F)
                - 0.0029 * Sin(2 * D - M) * E
                - 0.0027 * Sin(Mprime + F)
                + 0.0024 * Sin(M - Mprime) * E
                - 0.0021 * Sin(Mprime - 3 * F)
                - 0.0015 * Sin(2 * Mprime + F)
                - 0.0006 * Cos(2 * D - 2 * Mprime - F)
                - 0.0018 * Sin(3 * F)
                - 0.0017 * Cos(Mprime + 3 * F)
                + 0.0017 * Cos(2 * Mprime)
                + 0.0014 * Cos(2 * D - Mprime)
                - 0.0013 * Cos(2 * D + Mprime + F)
                - 0.0013 * Cos(Mprime)
                + 0.0012 * Sin(3 * Mprime + F)
                + 0.0011 * Sin(2 * D - Mprime + F)
                + 0.0011 * Cos(2 * D - 2 * Mprime)
                + 0.0010 * Cos(D + F)
                + 0.0010 * Sin(M + Mprime) * E
                - 0.0009 * Sin(2 * D - 2 * F)
                - 0.0007 * Cos(2 * Mprime + F)
                - 0.0007 * Cos(3 * Mprime + F);
            double decSum = direction == DeclinationType.Northern ?
                + 5.1093 * Sin(F)
                + 0.2658 * Cos(2 * F)
                + 0.1448 * Sin(2 * D - F)
                - 0.0322 * Sin(3 * F)
                + 0.0133 * Cos(2 * D - 2 * F)
                + 0.0125 * Cos(2 * D)
                - 0.0124 * Sin(Mprime - F)
                - 0.0101 * Sin(Mprime + 2 * F)
                + 0.0097 * Cos(F)
                - 0.0087 * Sin(2 * D + M - F) * E
                + 0.0074 * Sin(Mprime + 3 * F)
                + 0.0067 * Sin(D + F)
                + 0.0063 * Sin(Mprime - 2 * F)
                + 0.0060 * Sin(2 * D - M - F) * E
                - 0.0057 * Sin(2 * D - Mprime - F)
                - 0.0056 * Cos(Mprime + F)
                + 0.0052 * Cos(Mprime + 2 * F)
                + 0.0041 * Cos(2 * Mprime + F)
                - 0.0040 * Cos(Mprime - 3 * F)
                + 0.0038 * Cos(2 * Mprime - F)
                - 0.0034 * Cos(Mprime - 2 * F)
                - 0.0029 * Sin(2 * Mprime)
                + 0.0029 * Sin(3 * Mprime + F)
                - 0.0028 * Cos(2 * D + M - F) * E
                - 0.0028 * Cos(Mprime - F)
                - 0.0023 * Cos(3 * F)
                - 0.0021 * Sin(2 * D + F)
                + 0.0019 * Cos(Mprime + 3 * F)
                + 0.0018 * Cos(D + F)
                + 0.0017 * Sin(2 * Mprime - F)
                + 0.0015 * Cos(3 * Mprime + F)
                + 0.0014 * Cos(2 * D + 2 * Mprime + F)
                - 0.0012 * Sin(2 * D - 2 * Mprime - F)
                - 0.0012 * Cos(2 * Mprime)
                - 0.0010 * Cos(Mprime)
                - 0.0010 * Sin(2 * F)
                + 0.0006 * Sin(Mprime + F)
                :
                - 5.1093 * Sin(F)
                + 0.2658 * Cos(2 * F)
                - 0.1448 * Sin(2 * D - F)
                + 0.0322 * Sin(3 * F)
                + 0.0133 * Cos(2 * D - 2 * F)
                + 0.0125 * Cos(2 * D)
                - 0.0015 * Sin(Mprime - F)
                + 0.0101 * Sin(Mprime + 2 * F)
                - 0.0097 * Cos(F)
                + 0.0087 * Sin(2 * D + M - F) * E
                + 0.0074 * Sin(Mprime + 3 * F)
                + 0.0067 * Sin(D + F)
                - 0.0063 * Sin(Mprime - 2 * F)
                - 0.0060 * Sin(2 * D - M - F) * E
                + 0.0057 * Sin(2 * D - Mprime - F)
                - 0.0056 * Cos(Mprime + F)
                - 0.0052 * Cos(Mprime + 2 * F)
                - 0.0041 * Cos(2 * Mprime + F)
                - 0.0040 * Cos(Mprime - 3 * F)
                - 0.0038 * Cos(2 * Mprime - F)
                + 0.0034 * Cos(Mprime - 2 * F)
                - 0.0029 * Sin(2 * Mprime)
                + 0.0029 * Sin(3 * Mprime + F)
                + 0.0028 * Cos(2 * D + M - F) * E
                - 0.0028 * Cos(Mprime - F)
                + 0.0023 * Cos(3 * F)
                + 0.0021 * Sin(2 * D + F)
                + 0.0019 * Cos(Mprime + 3 * F)
                + 0.0018 * Cos(D + F)
                - 0.0017 * Sin(2 * Mprime - F)
                + 0.0015 * Cos(3 * Mprime + F)
                + 0.0014 * Cos(2 * D + 2 * Mprime + F)
                + 0.0012 * Sin(2 * D - 2 * Mprime - F)
                - 0.0012 * Cos(2 * Mprime)
                + 0.0010 * Cos(Mprime)
                - 0.0010 * Sin(2 * F)
                + 0.0037 * Sin(Mprime + F);
            double JDE = (direction == DeclinationType.Northern ? 2451562.5897 : 2451548.9289) + 27.321582247 * k + 0.000119804 * T2 - 0.000000141 * T3 + timeSum;
            double delta = 23.6961 - 0.013004 * T + decSum;
            return (JDE, delta);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <param name="observerLongitude">Optional for topocentric calculation.</param>
        /// <param name="observerLatitude">Optional.</param>
        /// <param name="observerElevation">Optional.</param>
        /// <returns>Total libration and position angle of axis of rotation.</returns>
        /// <remarks>pp. 371-375</remarks>
        public static (double longitude, double latitude, double positionAngleAxis)
        Libration(double julianEphemerisDay, double observerLongitude = double.NaN, double observerLatitude = double.NaN, double observerElevation = double.NaN) {
            var poz = Position(julianEphemerisDay);
            double epsilon = Sky.ObliquityOfEcliptic(julianEphemerisDay);
            if (!double.IsNaN(observerLongitude + observerLatitude + observerElevation)) {
                var equ = Transform.ToEquatorial(poz.longitude, poz.latitude, epsilon);
                equ = Transform.ToTopocentric(equ.rightAscension, equ.declination, Parallax(poz.distance), observerLongitude, observerLatitude, observerElevation, julianEphemerisDay);
                var ecl = Transform.ToEcliptical(equ.rightAscension, equ.declination, epsilon);
                poz = (ecl.longitude, ecl.latitude, poz.distance);
            }
            return LibrationSub(poz.longitude, poz.latitude, julianEphemerisDay, epsilon);
        }

        private static (double l, double b, double P)
        LibrationSub(double longitude, double latitude, double julianEphemerisDay, double epsilon = double.NaN) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double I = 1.54242.ToRadians();
            double lambda = longitude.ToRadians();
            double beta = latitude.ToRadians();
            double deltaPsi = Sky.Nutation(julianEphemerisDay).longitude.ToRadians();
            double Omega = AscendingNode(julianEphemerisDay, false).ToRadians();
            double W = lambda - deltaPsi - Omega;
            double A = Atan2(Sin(W) * Cos(beta) * Cos(I) - Sin(beta) * Sin(I), Cos(W) * Cos(beta));
            (double D, double M, double Mprime, double F) = DMMF(T, true);
            double lPrime = A - F;
            double bPrime = Asin(-Sin(W) * Cos(beta) * Sin(I) - Sin(beta) * Cos(I));
            double K1 = (119.75 + 131.849 * T).ToRadians();
            double K2 = (72.56 + 20.186 * T).ToRadians();
            double E = EarthEccentricity(T);
            double rho = -0.02752 * Cos(Mprime)
                        - 0.02245 * Sin(F)
                        + 0.00684 * Cos(Mprime - 2 * F)
                        - 0.00293 * Cos(2 * F)
                        - 0.00085 * Cos(2 * F - 2 * D)
                        - 0.00054 * Cos(Mprime - 2 * D)
                        - 0.00020 * Sin(Mprime + F)
                        - 0.00020 * Cos(Mprime + 2 * F)
                        - 0.00020 * Cos(Mprime - F)
                        + 0.00014 * Cos(Mprime + 2 * F - 2 * D);
            double tau = +0.02520 * E * Sin(M)
                        + 0.00473 * Sin(2 * Mprime - 2 * F)
                        - 0.00467 * Sin(Mprime)
                        + 0.00396 * Sin(K1)
                        + 0.00276 * Sin(Mprime - 2 * D)
                        + 0.00196 * Sin(Omega)
                        - 0.00183 * Cos(Mprime - F)
                        + 0.00115 * Sin(Mprime - 2 * D)
                        - 0.00096 * Sin(Mprime - D)
                        + 0.00046 * Sin(2 * F - 2 * D)
                        - 0.00039 * Sin(Mprime - F)
                        - 0.00032 * Sin(Mprime - M - D)
                        + 0.00027 * Sin(Mprime - M - 2 * D)
                        + 0.00023 * Sin(K2)
                        - 0.00014 * Sin(2 * D)
                        + 0.00014 * Cos(Mprime - 2 * F)
                        - 0.00012 * Sin(Mprime - 2 * F)
                        - 0.00012 * Sin(2 * Mprime)
                        + 0.00011 * Sin(Mprime - 2 * M - 2 * D);
            double sigma = -0.02816 * Sin(Mprime)
                        + 0.02244 * Cos(F)
                        - 0.00682 * Sin(Mprime - 2 * F)
                        - 0.00279 * Sin(2 * F)
                        - 0.00083 * Sin(2 * F - 2 * D)
                        + 0.00069 * Sin(Mprime - 2 * D)
                        + 0.00040 * Cos(Mprime + F)
                        - 0.00025 * Sin(2 * Mprime)
                        - 0.00023 * Sin(Mprime + 2 * F)
                        + 0.00020 * Cos(Mprime - F)
                        + 0.00019 * Sin(Mprime - F)
                        + 0.00013 * Sin(Mprime + 2 * F - 2 * D)
                        - 0.00010 * Cos(Mprime - 3 * F);
            double lDoublePrime = -tau + (rho * Cos(A) + sigma * Sin(A)) * Tan(bPrime);
            double bDoublePrime = sigma * Cos(A) - rho * Sin(A);
            double l = lPrime.ToDegrees().ToPlusMinus360() + lDoublePrime;
            double b = bPrime.ToDegrees() + bDoublePrime;
            double P = 0;
            if (!double.IsNaN(epsilon)) {
                double alpha = Transform.ToEquatorial(lambda.ToDegrees(), beta.ToDegrees(), epsilon).rightAscension.ToRadians();
                epsilon = epsilon.ToRadians();
                rho = rho.ToRadians();
                sigma = sigma.ToRadians();
                double V = Omega + deltaPsi + sigma / Sin(I);
                double X = Sin(I + rho) * Sin(V);
                double Y = Sin(I + rho) * Cos(V) * Cos(epsilon) - Cos(I + rho) * Sin(epsilon);
                double omega = Atan2(X, Y);
                P = Asin(Sqrt(X * X + Y * Y) * Cos(alpha - omega) / Cos(b.ToRadians())).ToDegrees();
            }
            return (l, b, P);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Coordinates of the subsolar point.</returns>
        /// <remarks>p. 376</remarks>
        public static (double longitude, double latitude, double colongitude)
        SubsolarPoint(double julianEphemerisDay) {
            var poz = Position(julianEphemerisDay);
            var sun = Body.Position(Bodies.Sun, julianEphemerisDay, true);
            double distRatio = poz.distance / (sun.distance * Sun.AUKilometers);
            double lambdaH = sun.longitude + 180 + distRatio * 57.296 * Cos(poz.latitude.ToRadians()) * Sin((sun.longitude - poz.longitude).ToRadians());
            double betaH = distRatio * poz.latitude;
            var sub = LibrationSub(lambdaH, betaH, julianEphemerisDay);
            return (sub.l.To360(), sub.b, (90 - sub.l).To360());
        }

        /// <param name="longitude">Selenographic longitude.</param>
        /// <param name="latitude">Selenographic latitude.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Altitude of sun above horizon in degrees.</returns>
        /// <remarks>p. 377</remarks>
        public static double SunAltitude(double longitude, double latitude, double julianEphemerisDay) {
            var ssp = SubsolarPoint(julianEphemerisDay);
            double eta = longitude.ToRadians();
            double theta = latitude.ToRadians();
            double b0 = ssp.latitude.ToRadians();
            double c0 = ssp.colongitude.ToRadians();
            return Asin(Sin(b0) * Sin(theta) + Cos(b0) * Cos(theta) * Sin(c0 + eta)).ToDegrees();
        }

        private static double SunAtHorizon(double longitude, double latitude, double julianEphemerisDay, int directionFactor) {
            double h, t = julianEphemerisDay;
            do {
                h = SunAltitude(longitude, latitude, t);
                t = t + (h / (12.19075 * Cos(latitude.ToRadians()))) * directionFactor;
            }
            while (Abs(h) > 0.001);
            return t;

        }

        /// <param name="longitude">Selenographic longitude.</param>
        /// <param name="latitude">Selenographic latitude.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Time of sunrise.</returns>
        /// <remarks>p. 377</remarks>
        public static double Sunrise(double longitude, double latitude, double julianEphemerisDay) {
            return SunAtHorizon(longitude, latitude, julianEphemerisDay, -1);
        }

        /// <param name="longitude">Selenographic longitude.</param>
        /// <param name="latitude">Selenographic latitude.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Time of sunset.</returns>
        /// <remarks>p. 377</remarks>
        public static double Sunset(double longitude, double latitude, double julianEphemerisDay) {
            return SunAtHorizon(longitude, latitude, julianEphemerisDay, 1);
        }

    }
}