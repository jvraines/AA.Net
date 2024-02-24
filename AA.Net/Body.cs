using System;
using static System.Math;
using System.Collections.Generic;
using System.Reflection;

namespace AA.Net {
    public enum Bodies {
        Mercury,
        Venus,
        Earth,
        Mars,
        Jupiter,
        Saturn,
        Uranus,
        Neptune,
        Pluto,
        Sun,
        Moon
    }

    public enum EclipseType {
        Total,
        Partial,
        Penumbral,
        Annular,
        Hybrid
    }

    internal enum Apsis {
        Perihelion,
        Aphelion
    }

    public static partial class Body {

        /// <param name="body">A value from the <see cref="Bodies"/> enumeration.</param>
        /// <param name="jde">Julian Day in dynamical time.</param>
        /// <param name="apparent">True for apparent or false for geometric position in the case of Sun or Moon only.</param>
        /// <returns>Geocentric position of Sun or Moon, or heliocentric position of other body, for equinox of date.</returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        /// <remarks>pp. 166-167, pp. 337-342</remarks>
        internal static (double longitude, double latitude, double distance)
        Position(Bodies body, double jde, bool apparent = true) {
            switch (body) {
                case Bodies.Sun:
                    var earth = Position(Bodies.Earth, jde);
                    double L = (earth.longitude + 180).To360();
                    double B = -earth.latitude;
                    double T = Time.JulianCentury2000(jde);
                    double lambdaPrime = (L - 1.397 * T - 0.00031 * Pow(T, 2)).ToRadians();
                    L -= Transform.ToDegrees(0, 0, 0.09033);
                    B += Transform.ToDegrees(0, 0, 0.03916) * (Cos(lambdaPrime) - Sin(lambdaPrime));
                    if (apparent) {
                        L += Sky.Nutation(jde).longitude;
                        L -= Transform.ToDegrees(0, 0, 20.4898) / earth.distance;
                    }
                    return (L, B, earth.distance);
                case Bodies.Moon:
                    T = Time.JulianCentury2000(jde);
                    double T2 = Pow(T, 2);
                    double T3 = Pow(T, 3);
                    double T4 = Pow(T, 4);
                    L = (218.3164477 + 481267.88123421 * T - 0.0015786 * T2 + T3 / 538841 - T4 / 65194000).To360();
                    (double D, double M, double Mprime, double F) = Moon.DMMF(T, false);
                    double A1 = (119.75 + 131.849 * T).To360();
                    double A2 = (53.09 + 479264.290 * T).To360();
                    double A3 = (313.45 + 481266.484 * T).To360();
                    double E = Moon.EarthEccentricity(T);
                    int[,] tableA = {
                        {0,0,1,0,6288774,-20905355},
                        {2,0,-1,0,1274027,-3699111},
                        {2,0,0,0,658314,-2955968},
                        {0,0,2,0,213618,-569925},
                        {0,1,0,0,-185116,48888},
                        {0,0,0,2,-114332,-3149},
                        {2,0,-2,0,58793,246158},
                        {2,-1,-1,0,57066,-152138},
                        {2,0,1,0,53322,-170733},
                        {2,-1,0,0,45758,-204586},
                        {0,1,-1,0,-40923,-129620},
                        {1,0,0,0,-34720,108743},
                        {0,1,1,0,-30383,104755},
                        {2,0,0,-2,15327,10321},
                        {0,0,1,2,-12528,1},
                        {0,0,1,-2,10980,79661},
                        {4,0,-1,0,10675,-34782},
                        {0,0,3,0,10034,-23210},
                        {4,0,-2,0,8548,-21636},
                        {2,1,-1,0,-7888,24208},
                        {2,1,0,0,-6766,30824},
                        {1,0,-1,0,-5163,-8379},
                        {1,1,0,0,4987,-16675},
                        {2,-1,1,0,4036,-12831},
                        {2,0,2,0,3994,-10445},
                        {4,0,0,0,3861,-11650},
                        {2,0,-3,0,3665,14403},
                        {0,1,-2,0,-2689,-7003},
                        {2,0,-1,2,-2602,1},
                        {2,-1,-2,0,2390,10056},
                        {1,0,1,0,-2348,6322},
                        {2,-2,0,0,2236,-9884},
                        {0,1,2,0,-2120,5751},
                        {0,2,0,0,-2069,1},
                        {2,-2,-1,0,2048,-4950},
                        {2,0,1,-2,-1773,4130},
                        {2,0,0,2,-1595,1},
                        {4,-1,-1,0,1215,-3958},
                        {0,0,2,2,-1110,1},
                        {3,0,-1,0,-892,3258},
                        {2,1,1,0,-810,2616},
                        {4,-1,-2,0,759,-1897},
                        {0,2,-1,0,-713,-2117},
                        {2,2,-1,0,-700,2354},
                        {2,1,-2,0,691,1},
                        {2,-1,0,-2,596,1},
                        {4,0,1,0,549,-1423},
                        {0,0,4,0,537,-1117},
                        {4,-1,0,0,520,-1571},
                        {1,0,-2,0,-487,-1739},
                        {2,1,0,-2,-399,1},
                        {0,0,2,-2,-381,-4421},
                        {1,1,1,0,351,1},
                        {3,0,-2,0,-340,1},
                        {4,0,-3,0,330,1},
                        {2,-1,2,0,327,1},
                        {0,2,1,0,-323,1165},
                        {1,1,-1,0,299,1},
                        {2,0,3,0,294,1},
                        {2,0,-1,-2,1,8752}
                    };
                    int[,] tableB = {
                        {0,0,0,1,5128122},
                        {0,0,1,1,280602},
                        {0,0,1,-1,277693},
                        {2,0,0,-1,173237},
                        {2,0,-1,1,55413},
                        {2,0,-1,-1,46271},
                        {2,0,0,1,32573},
                        {0,0,2,1,17198},
                        {2,0,1,-1,9266},
                        {0,0,2,-1,8822},
                        {2,-1,0,-1,8216},
                        {2,0,-2,-1,4324},
                        {2,0,1,1,4200},
                        {2,1,0,-1,-3359},
                        {2,-1,-1,1,2463},
                        {2,-1,0,1,2211},
                        {2,-1,-1,-1,2065},
                        {0,1,-1,-1,-1870},
                        {4,0,-1,-1,1828},
                        {0,1,0,1,-1794},
                        {0,0,0,3,-1749},
                        {0,1,-1,1,-1565},
                        {1,0,0,1,-1491},
                        {0,1,1,1,-1475},
                        {0,1,1,-1,-1410},
                        {0,1,0,-1,-1344},
                        {1,0,0,-1,-1335},
                        {0,0,3,1,1107},
                        {4,0,0,-1,1021},
                        {4,0,-1,1,833},
                        {0,0,1,-3,777},
                        {4,0,-2,1,671},
                        {2,0,0,-3,607},
                        {2,0,2,-1,596},
                        {2,-1,1,-1,491},
                        {2,0,-2,1,-451},
                        {0,0,3,-1,439},
                        {2,0,2,1,422},
                        {2,0,-3,-1,421},
                        {2,1,-1,1,-366},
                        {2,1,0,1,-351},
                        {4,0,0,1,331},
                        {2,-1,1,1,315},
                        {2,-2,0,-1,302},
                        {0,0,1,3,-283},
                        {2,1,1,-1,-229},
                        {1,1,0,-1,223},
                        {1,1,0,1,223},
                        {0,1,-2,-1,-220},
                        {2,1,-1,-1,-220},
                        {1,0,1,1,-185},
                        {2,-1,-2,-1,181},
                        {0,1,2,1,-177},
                        {4,0,-2,-1,176},
                        {4,-1,-1,-1,166},
                        {1,0,1,-1,-164},
                        {4,0,1,-1,132},
                        {1,0,-1,-1,-119},
                        {4,-1,0,-1,115},
                        {2,-2,0,1,107}
                    };
                    double sumL = 0;
                    double sumR = 0;
                    for (int row = 0; row <= tableA.GetUpperBound(0); row++) {
                        double sumA = tableA[row, 0] * D + tableA[row, 1] * M + tableA[row, 2] * Mprime + tableA[row, 3] * F;
                        double adjust = tableA[row, 4] == 1 ? 1 : Pow(E, Abs(tableA[row, 1]));
                        sumL += tableA[row, 4] * adjust * Sin(sumA.ToRadians());
                        adjust = tableA[row, 5] == 1 ? 1 : Pow(E, Abs(tableA[row, 1]));
                        sumR += tableA[row, 5] * adjust * Cos(sumA.ToRadians());
                    }
                    double sumB = 0;
                    for (int row = 0; row <= tableB.GetUpperBound(0); row++) {
                        double sumA = tableB[row, 0] * D + tableB[row, 1] * M + tableB[row, 2] * Mprime + tableB[row, 3] * F;
                        double adjust = Pow(E, Abs(tableB[row, 1]));
                        sumB += tableB[row, 4] * adjust * Sin(sumA.ToRadians());
                    }
                    sumL = sumL + 3958 * Sin(A1.ToRadians()) + 1962 * Sin((L - F).ToRadians()) + 318 * Sin(A2.ToRadians());
                    sumB = sumB - 2235 * Sin(L.ToRadians()) + 382 * Sin(A3.ToRadians()) + 175 * Sin((A1 - F).ToRadians()) + 175 * Sin((A1 + F).ToRadians()) + 127 * Sin((L - Mprime).ToRadians()) - 115 * Sin((L + Mprime).ToRadians());
                    double lambda = L + sumL / 1000000;
                    if (apparent) lambda += Sky.Nutation(jde).longitude;
                    double beta = sumB / 1000000;
                    double delta = 385000.56 + sumR / 1000;
                    return (lambda, beta, delta);
                case Bodies.Pluto:
                    if (jde < new DateTime(1885, 1, 1).JulianDay() || jde > new DateTime(2099, 12, 31).JulianDay()) throw new ArgumentOutOfRangeException("Only the years 1885 to 2099 are supported.");
                    int[,] terms = {
                        {0,0,1,-19799805,19850055,-5452852,-14974862,66865439,68951812},
                        {0,0,2,897144,-4954829,3527812,1672790,-11827535,-332538},
                        {0,0,3,611149,1211027,-1050748,327647,1593179,-1438890},
                        {0,0,4,-341243,-189585,178690,-292153,-18444,483220},
                        {0,0,5,129287,-34992,18650,100340,-65977,-85431},
                        {0,0,6,-38164,30893,-30697,-25823,31174,-6032},
                        {0,1,-1,20442,-9987,4878,11248,-5794,22161},
                        {0,1,0,-4063,-5071,226,-64,4601,4032},
                        {0,1,1,-6016,-3336,2030,-836,-1729,234},
                        {0,1,2,-3956,3039,69,-604,-415,702},
                        {0,1,3,-667,3572,-247,-567,239,723},
                        {0,2,-2,1276,501,-57,1,67,-67},
                        {0,2,-1,1152,-917,-122,175,1034,-451},
                        {0,2,0,630,-1277,-49,164,-129,504},
                        {1,-1,0,2571,-459,-197,199,480,-231},
                        {1,-1,1,899,-1449,-25,217,2,-441},
                        {1,0,-3,-1016,1043,589,-248,-3359,265},
                        {1,0,-2,-2343,-1012,-269,711,7856,-7832},
                        {1,0,-1,7042,788,185,193,36,45763},
                        {1,0,0,1199,-338,315,807,8663,8547},
                        {1,0,1,418,-67,-130,-43,-809,-769},
                        {1,0,2,120,-274,5,3,263,-144},
                        {1,0,3,-60,-159,2,17,-126,32},
                        {1,0,4,-82,-29,2,5,-35,-16},
                        {1,1,-3,-36,-29,2,3,-19,-4},
                        {1,1,-2,-40,7,3,1,-15,8},
                        {1,1,-1,-14,22,2,-1,-4,12},
                        {1,1,0,4,13,1,-1,5,6},
                        {1,1,1,5,2,0,-1,3,1},
                        {1,1,3,-1,0,0,0,6,-2},
                        {2,0,-6,2,0,0,-2,2,2},
                        {2,0,-5,-4,5,2,2,-2,-2},
                        {2,0,-4,4,-7,-7,0,14,13},
                        {2,0,-3,14,24,10,-8,-63,13},
                        {2,0,-2,-49,-34,-3,20,136,-236},
                        {2,0,-1,163,-48,6,5,273,1065},
                        {2,0,0,9,-24,14,17,251,149},
                        {2,0,1,-4,1,-2,0,-25,-9},
                        {2,0,2,-3,1,0,0,9,-2},
                        {2,0,3,1,3,0,0,-8,7},
                        {3,0,-2,-3,-1,0,1,2,-10},
                        {3,0,-1,5,-3,0,0,19,35},
                        {3,0,0,0,0,1,0,10,3}
                    };
                    T = Time.JulianCentury2000(jde);
                    double J = 34.35 + 3034.9057 * T;
                    double S = 50.08 + 1222.1138 * T;
                    double P = 238.96 + 144.96 * T;
                    sumL = 0; sumB = 0; sumR = 0;
                    for (int i = 0; i < 43; i++) {
                        double alpha = (terms[i, 0] * J + terms[i, 1] * S + terms[i, 2] * P).ToRadians();
                        double sinAlpha = Sin(alpha);
                        double cosAlpha = Cos(alpha);
                        sumL += terms[i, 3] * sinAlpha + terms[i, 4] * cosAlpha;
                        sumB += terms[i, 5] * sinAlpha + terms[i, 6] * cosAlpha;
                        sumR += terms[i, 7] * sinAlpha + terms[i, 8] * cosAlpha;
                    }
                    double l = 238.958116 + 144.96 * T + sumL / 1000000;
                    double b = -3.908239 + sumB / 1000000;
                    double r = 40.7241346 + sumR / 10000000;
                    return (l, b, r);
                default:
                    return SumPeriodicTerms(PositionTerms[(int)body], jde);
            }
        }

        /// <param name="terms">A subarray of <see cref="PositionTerms"/>.</param>
        /// <param name="jde">Julian Day in dynamical time.</param>
        /// <returns>Heliocentric coordinates and distance for equinox of date.</returns>
        /// <remarks>pp. 218-219</remarks>
        internal static (double longitude, double latitude, double distance)
        SumPeriodicTerms(double[][,] terms, double jde) {
            double tau = Time.JulianCentury2000(jde) / 10;
            double[] coord = new double[3];
            for (int c = 0; c < 3; c++) {
                double[,] data = terms[c];
                List<double> series = new List<double>();
                double sum = 0;
                double lastIndex = 0;
                for (int row = 0; row <= data.GetUpperBound(0); row++) {
                    if (data[row, 0] <= lastIndex) {
                        series.Add(sum);
                        sum = 0;
                    }
                    sum += data[row, 1] * Cos(data[row, 2] + data[row, 3] * tau);
                    lastIndex = data[row, 0];
                }
                series.Add(sum);
                double coordinate = series[0];
                for (int s = 1; s < series.Count; s++) {
                    coordinate += series[s] * Pow(tau, s);
                }
                coord[c] = coordinate /= 100000000;
            }
            return (coord[0].ToDegrees().To360(), coord[1].ToDegrees(), coord[2]);
        }

        /// <param name="planet">A value from the <see cref="Bodies"/> enumeration, excluding Sun, Moon, and Pluto.</param>
        /// <param name="apsis"></param>
        /// <param name="year"></param>
        /// <returns>Julian Day of requested phenomenon.</returns>
        /// <exception cref="ArgumentOutOfRangeException">A planet other than Pluto must be specified.</exception>
        /// <remarks>p. 269</remarks>
        internal static double NextApsis(Bodies planet, Apsis apsis, double year) {
            double k;
            (double, double, double) co;
            switch (planet) {
                case Bodies.Mercury:
                    k = 4.15201 * (year - 2000.12);
                    co = (2451590.257, 87.96934963, 0);
                    break;
                case Bodies.Venus:
                    k = 1.62549 * (year - 2000.53);
                    co = (2451738.233, 224.7008188, -0.0000000327);
                    break;
                case Bodies.Earth:
                    k = 0.99997 * (year - 2000.01);
                    co = (2451547.507, 365.2596358, 0.0000000156);
                    break;
                case Bodies.Mars:
                    k = 0.53166 * (year - 2001.78);
                    co = (2452195.026, 686.9957857, -0.0000001187);
                    break;
                case Bodies.Jupiter:
                    k = 0.08430 * (year - 2011.20);
                    co = (2455636.936, 4332.897065, 0.0001367);
                    break;
                case Bodies.Saturn:
                    k = 0.03393 * (year - 2003.52);
                    co = (2452830.12, 10764.21676, 0.000827);
                    break;
                case Bodies.Uranus:
                    k = 0.01190 * (year - 2051.1);
                    co = (2470213.5, 30694.8767, -0.00541);
                    break;
                case Bodies.Neptune:
                    k = 0.00607 * (year - 2047.5);
                    co = (2468895.1, 60190.33, 0.03429);
                    break;
                default:
                    throw new ArgumentOutOfRangeException("Sun, Moon, and Pluto not supported.");
            }
            k = apsis == Apsis.Perihelion ? Ceiling(k) : Round(k) + 0.5;
            return co.Item1 + co.Item2 * k + co.Item3 * k * k;
        }

        /// <param name="b">A value from the <see cref="Bodies"/> enumeration excepting Sun, Moon, Earth, and Pluto.</param>
        /// <param name="jde">Julian Day in dynamical time.</param>
        /// <returns>Illuminated fraction of disk, position angle of illuminated portion, and apparent magnitude.</returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        /// <remarks>p. 283, p. 285</remarks>
        internal static (double fraction, double positionAngle, double magnitude)
        Illumination(Bodies b, double jde) {
            var pp = Position(b, jde);
            double L = pp.longitude.ToRadians();
            double B = pp.latitude.ToRadians();
            var ep = Position(Bodies.Earth, jde);
            double L0 = ep.longitude.ToRadians();
            double B0 = ep.latitude.ToRadians();
            double x = pp.distance * Cos(B) * Cos(L) - ep.distance * Cos(B0) * Cos(L0);
            double y = pp.distance * Cos(B) * Sin(L) - ep.distance * Cos(B0) * Sin(L0);
            double z = pp.distance * Sin(B) - ep.distance * Sin(B0);
            double delta = Sqrt(x * x + y * y + z * z);
            double i = Acos((pp.distance - ep.distance * Cos(B) * Cos(L - L0)) / delta);
            double k = (1 + Cos(i)) / 2;
            var sp = Position(Bodies.Sun, jde);
            double epsilon = Sky.ObliquityOfEcliptic(jde);
            var spEqu = Transform.ToEquatorial(sp.longitude, sp.latitude, epsilon);
            var ppGeo = PositionGeocentric(b, jde);
            double sunDecSin = Sin(spEqu.declination.ToRadians());
            double sunDecCos = Cos(spEqu.declination.ToRadians());
            double moonDecSin = Sin(ppGeo.declination.ToRadians());
            double moonDecCos = Cos(ppGeo.declination.ToRadians());
            double raDiff = (spEqu.rightAscension - ppGeo.rightAscension).ToRadians();
            double chi = Atan2(sunDecCos * Sin(raDiff), sunDecSin * moonDecCos - sunDecCos * moonDecSin * Cos(raDiff)).ToDegrees().To360();
            double mag;
            double logTerm = 5 * Log10(pp.distance * delta);
            i = i.ToDegrees();
            switch (b) {
                case Bodies.Mercury:
                    mag = 1.16 + logTerm + 0.02838 * (i - 50) + 0.0001023 * (i - 50) * (i - 50);
                    break;
                case Bodies.Venus:
                    mag = -4 + logTerm + 0.01322 * i + 0.0000004247 * i * i * i;
                    break;
                case Bodies.Mars:
                    mag = -1.3 + logTerm + 0.01486 * i;
                    break;
                case Bodies.Jupiter:
                    mag = -8.93 + logTerm;
                    break;
                case Bodies.Saturn:
                    var r = Saturn.Rings(jde);
                    double eB = r.latitudeEarth.ToRadians();
                    mag = -8.68 + logTerm + 0.044 * r.deltaLongitude - 2.6 * Sin(Abs(eB)) + 1.25 * Sin(eB) * Sin(eB);
                    break;
                case Bodies.Uranus:
                    mag = -6.85 + logTerm;
                    break;
                case Bodies.Neptune:
                    mag = -7.05 + logTerm;
                    break;
                default:
                    throw new ArgumentOutOfRangeException($"{b} is not valid for this calculation.");
            }
            return (k, chi, mag);
        }

        static ArgumentOutOfRangeException OopsNotEarth = new ArgumentOutOfRangeException("Cannot evaluate for Earth.");

        /// <param name="body">A value of the <see cref="Bodies"/> enumeration except Earth.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Position and other values for further evaluation. Lambda and beta are returned in radians, in contrast to nearly all methods in this library.</returns>
        /// <remarks>p. 223-224</remarks>
        internal static (double l, double b, double r, double l0, double R, double x, double y, double z, double delta, double tau, double lambda, double beta)
        PositionFromLightTime(Bodies body, double julianEphemerisDay) {
            if (body == Bodies.Earth) throw OopsNotEarth;
            (double l0, double b0, double R) = Position(Bodies.Earth, julianEphemerisDay);
            double l0Rad = l0.ToRadians();
            double b0Rad = b0.ToRadians();
            double l, lRad, b, r, delta, x, y, z;
            double tau, lastTau, jde = julianEphemerisDay;
            Calc();
            do {
                jde = julianEphemerisDay - tau;
                lastTau = tau;
                Calc();
            }
            while (Abs(tau - lastTau) > 0.000005);
            double lambda = Atan2(y, x);
            double beta = Atan2(z, delta);
            return (l, b, r, l0, R, x, y, z, delta, tau, lambda, beta);

            void Calc() {
                (l, b, r) = Position(body, jde);
                double bRad = b.ToRadians();
                lRad = l.ToRadians();
                x = r * Cos(bRad) * Cos(lRad) - R * Cos(l0Rad);
                y = r * Cos(bRad) * Sin(lRad) - R * Sin(l0Rad);
                z = r * Sin(bRad) - R * Sin(b0Rad);
                delta = Sqrt(x * x + y * y + z * z);
                tau = Earth.LightTime(delta);
            }
        }

        /// <param name="b">A planet Mercury-Neptune of the <see cref="Bodies"/>  enumeration.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <param name="J2000">True to evaluate for J2000 or false for equinox of date.</param>
        /// <returns>Orbital elements for requested time and epoch.</returns>
        /// <remarks>pp. 209-215</remarks>
        internal static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(Bodies b, double julianEphemerisDay, bool J2000 = false) {
            if (b > Bodies.Neptune) throw new ArgumentOutOfRangeException("Must be a planet from Mercury to Neptune.");
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = T * T;
            double T3 = T2 * T;
            double[,,] TermsOfDate = new double[,,] {
                {
                    {252.250906,149474.0722491,0.00030350,0.000000018},
                    {0.387098310,0,0,0},
                    {0.20563175,0.000020407,-0.0000000283,-0.00000000018},
                    {7.004986,0.0018215,-0.00001810,0.000000056},
                    {48.330893,1.1861883,0.00017542,0.000000215},
                    {77.456119,1.5564776,0.00029544,0.000000009}
                },
                {
                    {181.979801,58519.2130302,0.00031014,0.000000015},
                    {0.723329820,0,0,0},
                    {0.00677192,-0.000047765,0.0000000981,0.00000000046},
                    {3.394662,0.0010037,-0.00000088,-0.000000007},
                    {76.679920,0.9011206,0.00040618,-0.000000093},
                    {131.563703,1.4022288,-0.00107618,-0.000005678}
                },
                {
                    {100.466457,36000.7698278,0.00030322,0.000000020},
                    {1.000001018,0,0,0},
                    {0.01670863,-0.000042037,-0.0000001267,0.00000000014},
                    {0,0,0,0},
                    {double.NaN,double.NaN,double.NaN,double.NaN},
                    {102.937348,1.7195366,0.00045688,-0.000000018}
                },
                {
                    {355.433000,19141.6964471,0.00031052,0.000000016},
                    {1.523679342,0,0,0},
                    {0.09340065,0.000090484,-0.0000000806,-0.00000000025},
                    {1.849726,-0.0006011,0.00001276,-0.000000007},
                    {49.558093,0.7720959,0.00001557,0.000002267},
                    {336.060234,1.8410449,0.00013477,0.000000536}
                },
                {
                    {34.351519,3036.3027748,0.00022330,0.000000037},
                    {5.202603209,0.0000001913,0,0},
                    {0.04849793,0.000163225,-0.0000004714,-0.00000000201},
                    {1.303267,-0.0054965,0.00000466,-0.000000002},
                    {100.464407,1.0209774,0.00040315,0.000000404},
                    {14.331207,1.6126352,0.00103042,-0.000004464}
                },
                {
                    {50.077444,1223.5110686,0.00051908,-0.000000030},
                    {9.554909192,-0.0000021390,0.000000004,0},
                    {0.05554814,-0.000346641,-0.0000006436,0.00000000340},
                    {2.488879,-0.0037362,-0.00001519,0.000000087},
                    {113.665503,0.8770880,-0.00012176,-0.000002249},
                    {93.057237,1.9637613,0.00083753,0.000004928}
                  },
                {
                    {314.055005,429.8640561,0.00030390,0.000000026},
                    {19.218446062,-0.0000000372,0.00000000098,0},
                    {0.04638122,-0.000027293,0.0000000789,0.00000000024},
                    {0.773197,0.0007744,0.00003749,-0.000000092},
                    {74.005957,0.5211278,0.00133947,0.000018484},
                    {173.005291,1.4863790,0.00021406,0.000000434}
                },
                {
                    {304.348665,219.8833092,0.00030882,0.000000018},
                    {30.110386869,-0.0000001663,0.00000000069,0},
                    {0.00945575,0.000006033,0,-0.00000000005},
                    {1.769953,-0.0093082,-0.00000708,0.000000027},
                    {131.784057,1.1022039,0.00025952,-0.000000637},
                    {48.120276,1.4262957,0.00038434,0.000000020}
                }
            };
            double[,,] TermsJ2000 = new double[,,] {
                {
                    {252.250906,149472.6746358,-0.00000536,0.000000002},
                    {7.004986,-0.0059516,0.00000080,0.000000043},
                    {48.330893,-0.1254227,-0.00008833,-0.000000200},
                    {77.456119,0.1588643,-0.00001342,-0.000000007}
                },
                {
                    {181.979801,58517.8156760,0.00000165,-0.000000002},
                    {3.394662,-0.0008568,-0.00003244,0.000000009},
                    {76.679920,-0.2780134,-0.00014257,-0.000000164},
                    {131.563703,0.0048746,-0.00138467,-0.000005695}
                },
                {
                    {100.466457,35999.3728565,-0.00000568,-0.000000001},
                    {0,0.0130548,-0.00000931,-0.000000034},
                    {174.873176,-0.2410908,0.00004262,0.000000001},
                    {102.937348,0.3225654,0.00014799,-0.000000039}
                },
                {
                    {355.433000,19140.2993039,0.00000262,-0.000000003},
                    {1.849726,-0.0081477,-0.00002255,-0.000000029},
                    {49.558093,-0.2950250,-0.00064048,-0.000001964},
                    {336.060234,0.4439016,-0.00017313,0.000000518}
                },
                {
                    {34.351519,3034.9056606,-0.00008501,0.000000016},
                    {1.303267,-0.0019877,0.00003320,0.000000097},
                    {100.464407,0.1767232,0.00090700,-0.000007272},
                    {14.331207,0.2155209,0.00072211,-0.000004485}
                },
                {
                    {50.077444,1222.1138488,0.00021004,-0.000000046},
                    {2.488879,0.0025514,-0.00004906,0.000000017},
                    {113.665503,-0.2566722,-0.00018399,0.000000480},
                    {93.057237,0.5665415,0.00052850,0.000004912}
                },
                {
                    {314.055005,428.4669983,-0.00000486,0.000000006},
                    {0.773197,-0.0016869,0.00000349,0.000000016},
                    {74.005957,0.0741431,0.00040539,0.000000119},
                    {173.005291,0.0893212,-0.00009470,0.000000414}
                },
                {
                    {304.348665,218.4862002,0.00000059,-0.000000002},
                    {1.769953,0.0002256,0.00000023,0},
                    {131.784057,-0.0061651,-0.00000219,-0.000000078},
                    {48.120276,0.0291866,0.00007610,0}
                }
            };
            double a = Sum(TermsOfDate, 1);
            double e = Sum(TermsOfDate, 2);
            double L, i, omega, pi;
            if (J2000) {
                L = Sum(TermsJ2000, 0).To360();
                i = Sum(TermsJ2000, 1);
                omega = Sum(TermsJ2000, 2);
                pi = Sum(TermsJ2000, 3);
            }
            else {
                L = Sum(TermsOfDate, 0).To360();
                i = Sum(TermsOfDate, 3);
                omega = Sum(TermsOfDate, 4);
                pi = Sum(TermsOfDate, 5);
            }
            return (L, a, e, i, omega, pi);

            double Sum(double[,,] arr, int r) {
                return arr[(int)b, r, 0] + arr[(int)b, r, 1] * T + arr[(int)b, r, 2] * T2 + arr[(int)b, r, 3] * T3;
            }
        }

        /// <param name="b">A value of the <see cref="Bodies"/> enumeration except Earth.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Semidiameter of disk in arcseconds.</returns>
        /// <remarks>pp. 389-391</remarks>
        internal static (double equatorial, double polar) Semidiameter(Bodies b, double julianEphemerisDay) {
            if (b == Bodies.Earth) throw OopsNotEarth;
            var result = Position(b, julianEphemerisDay);
            double delta = result.distance;
            if (b == Bodies.Moon) {
                double pi = Asin(6378.14 / delta);
                double s = Asin(0.272481 * Sin(pi)).ToDegrees() * 3600;
                return (s, s);
            }
            double j = 0, k = 0, CosB = 0;
            switch (b) {
                case Bodies.Sun:
                    j = 959.63;
                    break;
                case Bodies.Mercury:
                    j = 3.36;
                    break;
                case Bodies.Venus:
                    j = 8.41;
                    break;
                case Bodies.Mars:
                    j = 4.68;
                    break;
                case Bodies.Jupiter:
                    j = 98.44;
                    k = 92.06;
                    break;
                case Bodies.Saturn:
                    j = 82.73;
                    k = 73.82;
                    CosB = Cos(Saturn.Rings(julianEphemerisDay).latitudeEarth.ToRadians());
                    break;
                case Bodies.Uranus:
                    j = 35.02;
                    break;
                case Bodies.Neptune:
                    j = 33.5;
                    break;
                case Bodies.Pluto:
                    j = 2.07;
                    break;
            }
            double pol;
            double equ = j / delta;
            if (k == 0) pol = equ;
            else {
                pol = k / delta;
                if (b == Bodies.Saturn) pol = equ * Sqrt(1 - k * CosB * CosB);
            }
            return (equ, pol);
        }

        /// <param name="body">A value of the <see cref="Bodies"/> enumeration except Earth.</param>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Apparent right ascension in HMS and declination in degrees for equinox of date.</returns>
        /// <remarks>pp. 223-225</remarks>
        internal static (TimeSpan rightAscension, double declination) PositionGeocentric(Bodies body, double julianEphemerisDay) {
            if (body == Bodies.Earth) throw OopsNotEarth;
            if (body == Bodies.Sun || body == Bodies.Moon) {
                var p = Position(body, julianEphemerisDay);
                return Transform.ToEquatorial(p.longitude, p.latitude, Sky.ObliquityOfEcliptic(julianEphemerisDay));
            }
            var poz = PositionFromLightTime(body, julianEphemerisDay);
            poz.lambda = poz.lambda.ToDegrees().To360();
            poz.beta = poz.beta.ToDegrees();
            var ab = Sky.Aberration(poz.lambda, poz.beta, julianEphemerisDay);
            poz.lambda += ab.longitude;
            poz.beta += ab.latitude;
            var fk = Transform.ToFK5(poz.lambda, poz.beta, julianEphemerisDay);
            var equ = Transform.ToEquatorial(fk.longitude + Sky.Nutation(julianEphemerisDay).longitude, fk.latitude, Sky.ObliquityOfEcliptic(julianEphemerisDay, true));
            return (equ.rightAscension, equ.declination);
        }

        /// <param name="semimajorAxis">Degrees. Ignored in the case of a parabolic orbit.</param>
        /// <param name="inclination">Degrees.</param>
        /// <param name="ascendingNode">Degrees.</param>
        /// <param name="argumentPerihelion">Degrees.</param>
        /// <param name="timePerihelion">Julian Days.</param>
        /// <param name="distancePerihelion">AU. Ignored in the case of an elliptic orbit.</param>
        /// <returns>Apparent right ascension in HMS, declination in degrees, elongation in degrees, and phase angle in degrees for equinox of J2000.0</returns>
        /// <remarks>pp. 227-230. Arguments must be specified for equinox of J2000.0.</remarks>
        public static (TimeSpan rightAscension, double declination, double elongation, double phaseAngle) PositionGeocentric(double julianEphemerisDay, double semimajorAxis, double eccentricity, double inclination, double ascendingNode, double argumentPerihelion, double timePerihelion, double distancePerihelion) {
            double e = eccentricity;
            if (e < 0.98 && semimajorAxis <= 0) throw new ArgumentOutOfRangeException("semimajorAxis", "Must be greater than zero.");
            double n = 0.9856076686 / (semimajorAxis * Sqrt(semimajorAxis));
            double i = inclination.ToRadians();
            double omega = argumentPerihelion.ToRadians();
            double Omega = ascendingNode.ToRadians();
            double q = distancePerihelion;
            double sinEpsilon = 0.397777156;
            double cosEpsilon = 0.917482062;
            double F = Cos(Omega);
            double G = Sin(Omega) * cosEpsilon;
            double H = Sin(Omega) * sinEpsilon;
            double P = -Sin(Omega) * Cos(i);
            double Q = Cos(Omega) * Cos(i) * cosEpsilon - Sin(i) * sinEpsilon;
            double R = Cos(Omega) * Cos(i) * sinEpsilon + Sin(i) * cosEpsilon;
            double A = Atan2(F, P);
            double B = Atan2(G, Q);
            double C = Atan2(H, R);
            double a = Sqrt(F * F + P * P);
            double b = Sqrt(G * G + Q * Q);
            double c = Sqrt(H * H + R * R);
            var poz = Sun.PositionRectangular(julianEphemerisDay, Time.J2000);
            double M, E, v, r, x, y, z;
            xyz(julianEphemerisDay);
            double xi, eta, zeta, Delta, delta;
            xiEtaZeta();
            double tau = Earth.LightTime(Delta);
            xyz(julianEphemerisDay - tau);
            xiEtaZeta();
            double alpha = Atan2(eta, xi);
            R = Sqrt(poz.x * poz.x + poz.y * poz.y + poz.z * poz.z);
            double psi = Acos((R * R + Delta * Delta - r * r) / (2 * R * Delta));
            double beta = Acos((r * r + Delta * Delta - R * R) / (2 * r * Delta));
            return (alpha.ToDegrees().To360().ToRightAscension(), delta.ToDegrees(), psi.ToDegrees(), beta.ToDegrees());

            void xyz(double t) {
                if (e < 0.98) {
                    M = (t - timePerihelion) * n;
                    E = EccentricAnomaly(M, e).ToRadians();
                    v = 2 * Atan(Sqrt((1 + e) / (1 - e)) * Tan(E / 2));
                    r = semimajorAxis * (1 - e * Cos(E));
                }
                else if (e == 1) {
                    double W = 0.03649116245 / (q * Sqrt(q)) * (t - timePerihelion);
                    double s = 0, lastS;
                    do {
                        lastS = s;
                        s = (2 * s * s * s + W) / (3 * (s * s + 1));
                    }
                    while (Abs(s - lastS) > 0.00000005);
                    v = Atan(s) * 2;
                    r = q * (1 + s * s);
                }
                else {
                    double K = 0.01720209895;
                    int D1 = 10000;
                    double CC = 1 / 3, DD = 1E-9;
                    double Q1 = K * Sqrt((1 + e) / q) / (2 * q);
                    double GG = (1 - e) / (1 + e);
                    double T = t - timePerihelion;
                    if (T == 0) {
                        r = q; v = 0;
                    }
                    else {
                        Exception noConvergence = new Exception("True anomaly/radius vector solution does not converge.");
                        double Q2 = Q1 * T;
                        double S = 2 / (3 * Abs(Q2));
                        S = 2 / Tan(2 * Atan(Pow(Tan(Atan(S) / 2), CC)));
                        if (T < 0) S = -S;
                        double L = 0, S0;
                        do {
                            S0 = S;
                            double Z = 1, Y = S * S, G1 = -Y * S;
                            double Q3 = Q2 + 2 * GG * S * Y / 3;
                            double FF;
                            do {
                                Z++;
                                G1 = -G1 * GG * Y;
                                double Z1 = (Z - (Z + 1) * GG) / (2 * Z + 1);
                                FF = Z1 * G1;
                                Q3 += FF;
                                if (Z > 50 || Abs(FF) > D1) throw noConvergence;
                            }
                            while (Abs(FF) > DD);
                            L++; if (L > 50) throw noConvergence;
                            double S1;
                            do {
                                S1 = S;
                                S = (2 * S * S * S / 3 + Q3) / (S * S + 1);
                            }
                            while (Abs(S - S1) > DD);
                        }
                        while (Abs(S - S0) > DD);
                        v = 2 * Atan(S);
                        r = q * (1 + e) / (1 + e * Cos(v));
                        if (v < 0) v += 2 * PI;
                    }
                }
                x = r * a * Sin(A + omega + v);
                y = r * b * Sin(B + omega + v);
                z = r * c * Sin(C + omega + v);
            }

            void xiEtaZeta() {
                xi = poz.x + x;
                eta = poz.y + y;
                zeta = poz.z + z;
                Delta = Sqrt(xi * xi + eta * eta + zeta * zeta);
                delta = Asin(zeta / Delta);
            }
        }

        /// <param name="meanAnomaly">Mean anomaly in degrees.</param>
        /// <param name="eccentricity">Eccentricity of orbit.</param>
        /// <returns>Eccentric anomaly in degrees.</returns>
        /// <remarks>p. 206</remarks>
        public static double EccentricAnomaly(double meanAnomaly, double eccentricity) {
            double M = meanAnomaly.To360().ToRadians();
            double E = eccentricity;
            double F = Sign(M);
            M = Abs(M) / (2 * PI);
            M = (M - Floor(M)) * 2 * PI * F;
            F = 1;
            if (M > PI) {
                F = -1;
                M = 2 * PI - M;
            }
            double E0 = PI / 2, D = PI / 4;
            for (int j = 1; j < 53; j++) {
                double M1 = E0 - E * Sin(E0);
                E0 += D * Sign(M - M1);
                D /= 2.0;
            }
            return (E0 * F).ToDegrees();
        }

        /// <param name="eccentricity">When equal to 1, orbit is parabolic and <paramref name="distancePerihelion"/> must be specified.</param>
        /// <param name="argumentPerihelion"></param>
        /// <param name="timePerihelion">Julian Day.</param>
        /// <param name="semimajorAxis">When <paramref name="eccentricity"/> is less than 1, either this argument or <paramref name="distancePerihelion"/> must be specified.</param>
        /// <param name="distancePerihelion">When <paramref name="eccentricity"/> is less than 1, either this argument or <paramref name="semimajorAxis"/> must be specified.</param>
        /// <param name="meanMotion"></param>
        /// <returns>Julian Day and distance of each node in AU.</returns>
        /// <exception cref="ArgumentException">The correct combination of arguments has not been supplied.</exception>
        /// <remarks>pp. 275-276</remarks>
        public static (double timeAscending, double distanceAscending, double timeDescending, double distanceDescending)
        Nodes(double eccentricity, double argumentPerihelion, double timePerihelion, double semimajorAxis = -1, double distancePerihelion = -1, double meanMotion = -1) {
            if (eccentricity == 1) {
                if (distancePerihelion == -1) throw new ArgumentException("When eccentricity = 1, distancePerihelion must be specified.");
            }
            else {
                if (semimajorAxis == -1) {
                    if (distancePerihelion == -1) throw new ArgumentException("At least one of semimajorAxis or distancePerihelion must be specified.");
                    semimajorAxis = distancePerihelion / (1 - eccentricity);
                }
                if (meanMotion == -1) meanMotion = 0.9856076686 / (semimajorAxis * Sqrt(semimajorAxis));
            }
            double tA = 0, dA = 0, tD = 0, dD = 0;
            for (int j = 0; j < 2; j++) {
                double v = (j == 0 ? 360 - argumentPerihelion : 180 - argumentPerihelion).ToRadians();
                double t, r;
                if (eccentricity == 1) {
                    double s = Tan(v / 2);
                    t = timePerihelion + 27.403895 * (s * s * s + 3 * s) * distancePerihelion * Sqrt(distancePerihelion);
                    r = distancePerihelion * (1 + s * s);
                }
                else {
                    double E = Atan(Sqrt((1 - eccentricity) / (1 + eccentricity)) * Tan(v / 2)) * 2;
                    double M = (E - eccentricity * Sin(E)).ToDegrees();
                    t = timePerihelion + M / meanMotion;
                    r = semimajorAxis * (1 - eccentricity * Cos(E));
                }
                if (j == 0) {
                    tA = t; dA = r;
                }
                else {
                    tD = t; dD = r;
                }
            }
            return (tA, dA, tD, dD);
        }

        /// <param name="k">Lunation number from J2000.</param>
        /// <param name="isLunar">True for a lunar eclipse or false for a solar eclipse.</param>
        /// <returns>Time of greatest eclipse and other data for further evaluation.</returns>
        /// <remarks>pp. 379-383</remarks>
        internal static (double greatestEclipse, bool isAscending, double gamma, double rho, double sigma, double u, double Mprime)
        FindNextEclipse(double k, bool isLunar) {
            double T, T2, T3, T4, F;
            double gamma, greatestEclipse, Mprime, u, rho, sigma;
            bool isAscending, fail;
            do {
                do {
                    k++;
                    T = Moon.JulianCentury(k);
                    T2 = Pow(T, 2);
                    T3 = Pow(T, 3);
                    T4 = Pow(T, 4);
                    F = (160.7108 + 390.67050284 * k - 0.0016118 * T2 - 0.00000227 * T3 + 0.000000011 * T4).To360();
                }
                while (Abs(Sin(F.ToRadians())) > 0.36);
                double lineup = Moon.NewMoon2000
                    + 29.530588861 * k
                    + 0.00015437 * T2
                    - 0.00000015 * T3
                    + 0.00000000073 * T4;
                double M = (2.5534 + 29.1053567 * k - 0.0000014 * T2 - 0.00000011 * T3).To360().ToRadians();
                Mprime = (201.5643 + 385.81693528 * k + 0.0107582 * T2 + 0.00001238 * T3 - 0.000000058 * T4).To360().ToRadians();
                double omega = (124.7746 - 1.56375588 * k + 0.0020672 * T2 + 0.00000215 * T3).To360().ToRadians();
                double E = Moon.EarthEccentricity(T);
                double F1 = (F - 0.02665 * Sin(omega)).ToRadians();
                double A1 = (299.77 + 0.107408 * k - 0.009173 * T2).ToRadians();
                greatestEclipse = lineup
                    - (isLunar ? 0.4065 : 0.4075) * Sin(Mprime)
                    + (isLunar ? 0.1727 : 0.1721) * E * Sin(M)
                    + 0.0161 * Sin(2 * Mprime)
                    - 0.0097 * Sin(2 * F1)
                    + 0.0073 * E * Sin(Mprime - M)
                    - 0.0050 * E * Sin(Mprime + M)
                    - 0.0023 * Sin(Mprime - 2 * F1)
                    + 0.0021 * E * Sin(2 * M)
                    + 0.0012 * Sin(Mprime + 2 * F1)
                    + 0.0006 * E * Sin(2 * Mprime + M)
                    - 0.0004 * Sin(3 * Mprime)
                    - 0.0003 * E * Sin(M + 2 * F1)
                    + 0.0003 * Sin(A1)
                    - 0.0002 * E * Sin(M - 2 * F1)
                    - 0.0002 * E * Sin(2 * Mprime - M)
                    - 0.0002 * Sin(omega);
                isAscending = F < 90 | F > 270;
                double P =
                      0.2070 * E * Sin(M)
                    + 0.0024 * E * Sin(2 * M)
                    - 0.0392 * Sin(Mprime)
                    + 0.0116 * Sin(2 * Mprime)
                    - 0.0073 * E * Sin(Mprime + M)
                    + 0.0067 * E * Sin(Mprime - M)
                    + 0.0118 * Sin(2 * F1);
                double Q =
                      5.2207
                    - 0.0048 * E * Cos(M)
                    + 0.0020 * E * Cos(2 * M)
                    - 0.3299 * Cos(Mprime)
                    - 0.0060 * E * Cos(Mprime + M)
                    + 0.0041 * E * Cos(Mprime - M);
                double W = Abs(Cos(F1));
                gamma = (P * Cos(F1) + Q * Sin(F1)) * (1 - 0.0048 * W);
                u =
                      0.0059
                    + 0.0046 * E * Cos(M)
                    - 0.0182 * Cos(Mprime)
                    + 0.0004 * Cos(2 * Mprime)
                    - 0.0005 * Cos(M + Mprime);
                rho = 1.2848 + u;
                sigma = 0.7403 - u;
                if (isLunar) fail = 1.5573 + u - Abs(gamma) < 0 && 1.0128 - u - Abs(gamma) < 0;
                else fail = Abs(gamma) > 1.5433 + u;
            }
            while (fail);
            return (greatestEclipse, isAscending, gamma, rho, sigma, u, Mprime);
        }

        /// <param name="alpha0">Geocentric right ascension of Sun.</param>
        /// <param name="delta0">Geocentric declination of Sun.</param>
        /// <param name="alpha">Geocentric right ascension of body.</param>
        /// <param name="delta">Geocentric declination of body.</param>
        /// <returns>Eastward angle of limb midpoint from North Point of disk in degrees.</returns>
        /// <remarks>p. 346</remarks>
        internal static double BrightLimbPositionAngle(double alpha0, double delta0, double alpha, double delta) {
            return Atan2(Cos(delta0) * Sin(alpha0 - alpha), Sin(delta0) * Cos(delta) - Cos(delta0) * Sin(delta) * Cos(alpha0 - alpha)).ToDegrees();
        }
    }
}