using System;
using System.Collections.Generic;
using static System.Math;

namespace AA.Net {

    public static class Sun {

        /// <summary>
        /// Length of astronomical unit (AU) in kilometers. (IAU 2012 Resolution B2)
        /// </summary>
        public const double AUKilometers = 149597870.7;

        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay, bool apparent = true) {
            return Body.Position(Bodies.Sun, julianEphemerisDay, apparent);
        }

        /// <summary>Find the rectangular coordinates at a specified time.</summary>
        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <param name="epoch">Julian Day of the epoch of the reference equinox. Set equal to the first parameter for equinox of date.</param>
        /// <returns>Rectangular geocentric equatorial coordinates in AU.</returns>
        /// <remarks>pp. 171-175</remarks>
        public static (double x, double y, double z) PositionRectangular(double julianEphemerisDay, double epoch) {
            double lambda, beta, epsilon, x, y, z;
            if (epoch == julianEphemerisDay) {
                var poz = Body.Position(Bodies.Sun, julianEphemerisDay, false); 
                lambda = poz.longitude.ToRadians();
                beta = poz.latitude.ToRadians();
                epsilon = Sky.ObliquityOfEcliptic(julianEphemerisDay, false).ToRadians();
                x = poz.distance * Cos(beta) * Cos(lambda);
                y = poz.distance * (Cos(beta) * Sin(lambda) * Cos(epsilon) - Sin(beta) * Sin(epsilon));
                z = poz.distance * (Cos(beta) * Sin(lambda) * Sin(epsilon) + Sin(beta) * Cos(epsilon));
                return (x, y, z);
            }
            double[][,] eterms = new double[][,] {
                new double[,] {
                    {1,175347046,0,0},
                    {2,3341656,4.6692568,6283.0758500},
                    {3,34894,4.62610,12566.15170},
                    {4,3497,2.7441,5753.3849},
                    {5,3418,2.8289,3.5231},
                    {6,3136,3.6277,77713.7715},
                    {7,2676,4.4181,7860.4194},
                    {8,2343,6.1352,3930.2097},
                    {9,1324,0.7425,11506.7698},
                    {10,1273,2.0371,529.6910},
                    {11,1199,1.1096,1577.3435},
                    {12,990,5.233,5884.927},
                    {13,902,2.045,26.298},
                    {14,857,3.508,398.149},
                    {15,780,1.179,5223.694},
                    {16,753,2.533,5507.553},
                    {17,505,4.583,18849.228},
                    {18,492,4.205,775.523},
                    {19,357,2.920,0.067},
                    {20,317,5.849,11790.629},
                    {21,284,1.899,796.298},
                    {22,271,0.315,10977.079},
                    {23,243,0.345,5486.778},
                    {24,206,4.806,2544.314},
                    {25,205,1.869,5573.143},
                    {26,202,2.458,6069.777},
                    {27,156,0.833,213.299},
                    {28,132,3.411,2942.463},
                    {29,126,1.083,20.775},
                    {30,115,0.645,0.980},
                    {31,103,0.636,4694.003},
                    {32,102,0.976,15720.839},
                    {33,102,4.267,7.114},
                    {34,99,6.21,2146.17},
                    {35,98,0.68,155.42},
                    {36,86,5.98,161000.69},
                    {37,85,1.30,6275.96},
                    {38,85,3.67,71430.70},
                    {39,80,1.81,17260.15},
                    {40,79,3.04,12036.46},
                    {41,75,1.76,5088.63},
                    {42,74,3.50,3154.69},
                    {43,74,4.68,801.82},
                    {44,70,0.83,9437.76},
                    {45,62,3.98,8827.39},
                    {46,61,1.82,7084.90},
                    {47,57,2.78,6286.60},
                    {48,56,4.39,14143.50},
                    {49,56,3.47,6279.55},
                    {50,52,0.19,12139.55},
                    {51,52,1.33,1748.02},
                    {52,51,0.28,5856.48},
                    {53,49,0.49,1194.45},
                    {54,41,5.37,8429.24},
                    {55,41,2.40,19651.05},
                    {56,39,6.17,10447.39},
                    {57,37,6.04,10213.29},
                    {58,37,2.57,1059.38},
                    {59,36,1.71,2352.87},
                    {60,36,1.78,6812.77},
                    {61,33,0.59,17789.85},
                    {62,30,0.44,83996.85},
                    {63,30,2.74,1349.87},
                    {64,25,3.16,4690.48},
                    {1,628307584999,0,0},
                    {2,206059,2.678235,6283.075850},
                    {3,4303,2.6351,12566.1517},
                    {4,425,1.590,3.523},
                    {5,119,5.796,26.298},
                    {6,109,2.966,1577.344},
                    {7,93,2.59,18849.23},
                    {8,72,1.14,529.69},
                    {9,68,1.87,398.15},
                    {10,67,4.41,5507.55},
                    {11,59,2.89,5223.69},
                    {12,56,2.17,155.42},
                    {13,45,0.40,796.30},
                    {14,36,0.47,775.52},
                    {15,29,2.65,7.11},
                    {16,21,5.34,0.98},
                    {17,19,1.85,5486.78},
                    {18,19,4.97,213.30},
                    {19,17,2.99,6275.96},
                    {20,16,0.03,2544.31},
                    {21,16,1.43,2146.17},
                    {22,15,1.21,10977.08},
                    {23,12,2.83,1748.02},
                    {24,12,3.26,5088.63},
                    {25,12,5.27,1194.45},
                    {26,12,2.08,4694.00},
                    {27,11,0.77,553.57},
                    {28,10,1.30,6286.60},
                    {29,10,4.24,1549.87},
                    {30,9,2.70,242.73},
                    {31,9,5.64,951.72},
                    {32,8,5.30,2352.87},
                    {33,6,2.65,9437.76},
                    {34,6,4.67,4690.48},
                    {1,8722,1.0725,6283.0758},
                    {2,991,3.1416,0},
                    {3,295,0.437,12566.152},
                    {4,27,0.05,3.52},
                    {5,16,5.19,26.30},
                    {6,16,3.69,155.42},
                    {7,9,0.30,18849.23},
                    {8,9,2.06,77713.77},
                    {9,7,0.83,775.52},
                    {10,5,4.66,1577.34},
                    {11,4,1.03,7.11},
                    {12,4,3.44,5573.14},
                    {13,3,5.14,796.30},
                    {14,3,6.05,5507.55},
                    {15,3,1.19,242.73},
                    {16,3,6.12,529.69},
                    {17,3,0.30,398.15},
                    {18,3,2.28,553.57},
                    {19,2,4.38,5223.69},
                    {20,2,3.75,0.98},
                    {1,289,5.842,6283.076},
                    {2,21,6.05,12566.15},
                    {3,3,5.20,155.42},
                    {4,3,3.14,0},
                    {5,1,4.72,3.52},
                    {6,1,5.97,242.73},
                    {7,1,5.54,18849.23},
                    {1,8,4.14,6283.08},
                    {2,1,3.28,12566.15}
                },
                new double[,] {
                    {1,280,3.199,84334.662},
                    {2,102,5.422,5507.553},
                    {3,80,3.88,5223.69},
                    {4,44,3.70,2352.87},
                    {5,32,4.00,1577.34},
                    {1,227778,3.413766,6283.075850},
                    {2,3806,3.3706,12566.1517},
                    {3,3620,0,0},
                    {4,72,3.33,18849.23},
                    {5,8,3.89,5507.55},
                    {6,8,1.79,5223.69},
                    {7,6,5.20,2352.87},
                    {1,9721,5.1519,6283.07585},
                    {2,233,3.1416,0},
                    {3,134,0.644,12566.152},
                    {4,7,1.07,18849.23},
                    {1,276,0.595,6283.076},
                    {2,17,3.14,0},
                    {3,4,0.12,12566.15},
                    {1,6,2.27,6283.08},
                    {2,1,0,0}
                },
                new double[,] {
                    {1,100013989,0,0},
                    {2,1670700,3.0984635,6283.0758500},
                    {3,13956,3.05525,12566.15170},
                    {4,3084,5.1985,77713.7715},
                    {5,1628,1.1739,5753.3849},
                    {6,1576,2.8469,7860.4194},
                    {7,925,5.453,11506.770},
                    {8,542,4.564,3930.210},
                    {9,472,3.661,5884.927},
                    {10,346,0.964,5507.553},
                    {11,329,5.900,5223.694},
                    {12,307,0.299,5573.143},
                    {13,243,4.273,11790.629},
                    {14,212,5.847,1577.344},
                    {15,186,5.022,10977.079},
                    {16,175,3.012,18849.228},
                    {17,110,5.055,5486.778},
                    {18,98,0.89,6069.78},
                    {19,86,5.69,15720.84},
                    {20,86,1.27,161000.69},
                    {21,65,0.27,17260.15},
                    {22,63,0.92,529.69},
                    {23,57,2.01,83996.85},
                    {24,56,5.24,71430.70},
                    {25,49,3.25,2544.31},
                    {26,47,2.58,775.52},
                    {27,45,5.54,9437.76},
                    {28,43,6.01,6275.96},
                    {29,39,5.36,4694.00},
                    {30,38,2.39,8827.39},
                    {31,37,0.83,19651.05},
                    {32,37,4.90,12139.55},
                    {33,36,1.67,12036.46},
                    {34,35,1.84,2942.46},
                    {35,33,0.24,7084.90},
                    {36,32,0.18,5088.63},
                    {37,32,1.78,398.15},
                    {38,28,1.21,6286.60},
                    {39,28,1.90,6279.55},
                    {40,26,4.59,10447.39},
                    {1,103019,1.107490,6283.075850},
                    {2,1721,1.0644,12566.1517},
                    {3,702,3.142,0},
                    {4,32,1.02,18849.23},
                    {5,31,2.84,5507.55},
                    {6,25,1.32,5223.69},
                    {7,18,1.42,1577.34},
                    {8,10,5.91,10977.08},
                    {9,9,1.42,6275.96},
                    {10,9,0.27,5486.78},
                    {1,4359,5.7846,6283.0758},
                    {2,124,5.579,12566.152},
                    {3,12,3.14,0},
                    {4,9,3.63,77713.77},
                    {5,6,1.87,5573.14},
                    {6,3,5.47,18849.23},
                    {1,145,4.273,6283.076},
                    {2,7,3.92,12566.15},
                    {1,4,2.56,6283.08}
                }
            };
            var earth = Body.SumPeriodicTerms(eterms, julianEphemerisDay);
            lambda = (earth.longitude + 180).To360().ToRadians();
            beta = (-earth.latitude).ToRadians();
            x = earth.distance * Cos(beta) * Cos(lambda);
            y = earth.distance * Cos(beta) * Sin(lambda);
            z = earth.distance * Sin(beta);
            double x0 = x + 0.00000044036 * y - 0.000000190919 * z;
            double y0 = -0.000000479966 * x + 0.917482137087 * y - 0.397776982902 * z;
            double z0 = 0.397776982902 * y + 0.917482137087 * z;
            if (epoch == Time.J2000) return (x0, y0, z0);
            double t = (epoch - Time.J2000) / 36525;
            double t2 = t * t;
            double t3 = t2 * t;
            double zeta = ((2306.2181 * t + 0.30188 * t2 + 0.017998 * t3) / 3600).ToRadians();
            double zed = ((2306.2181 * t + 1.09468 * t2 + 0.018203 * t3) / 3600).ToRadians();
            double theta = ((2004.3109 * t - 0.42665 * t2 - 0.041833 * t3) / 3600).ToRadians();
            double Xx = Cos(zeta) * Cos(zed) * Cos(theta) - Sin(zeta) * Sin(zed);
            double Xy = Sin(zeta) * Cos(zed) + Cos(zeta) * Sin(zed) * Cos(theta);
            double Xz = Cos(zeta) * Sin(theta);
            double Yx = -Cos(zeta) * Sin(zed) - Sin(zeta) * Cos(zed) * Cos(theta);
            double Yy = Cos(zeta) * Cos(zed) - Sin(zeta) * Sin(zed) * Cos(theta);
            double Yz = -Sin(zeta) * Sin(theta);
            double Zx = -Cos(zed) * Sin(theta);
            double Zy = -Sin(zed) * Sin(theta);
            double Zz = Cos(theta);
            x = Xx * x0 + Yx * y0 + Zx * z0;
            y = Xy * x0 + Yy * y0 + Zy * z0;
            z = Xz * x0 + Yz * y0 + Zz * z0;
            return (x, y, z);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Sun, julianEphemerisDay);
        }

        /// <param name="start">Julian Day at which to begin the search.</param>
        /// <returns>Time and characteristics of the next solar eclipse.</returns>
        /// <remarks>pp. 379-382</remarks>
        public static (EclipseType type, double greatestEclipse, double? magnitude, bool isCentral, bool isAscending, double axisDistance, double radiusUmbra, double radiusPenumbra)
        NextEclipse(double start) {
            double k = Ceiling(Moon.Lunation(start) - 0.5) - 1;
            var e = Body.FindNextEclipse(k, false);
            while (e.greatestEclipse < start) e = Body.FindNextEclipse(++k, false);
            double aGamma = Abs(e.gamma);
            bool central;
            double? mag = null;
            EclipseType type;
            if (aGamma > 0.9972) {
                central = false;
                mag = (1.5433 + e.u - aGamma) / (0.5461 + 2 * e.u);
                if (aGamma < 0.9972 + Abs(e.u)) type = mag >= 1 ? EclipseType.Total : EclipseType.Annular;
                else type = EclipseType.Partial;
            }
            else {
                central = true;
                if (e.u < 0) type = EclipseType.Total;
                else if (e.u > 0.0047) type = EclipseType.Annular;
                else {
                    double omega = 0.00464 * Sqrt(1 - e.gamma * e.gamma);
                    if (e.u < omega) type = EclipseType.Hybrid;
                    else type = EclipseType.Annular;
                }
            }
            return (type, e.greatestEclipse, mag, central, e.isAscending, e.gamma, e.u, e.u + 0.5461);
        }

        ///<param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Eastward angle of the axis of rotation from the North Point, heliographic longitude and latitude of the center.</returns>
        /// <remarks>pp. 189-190</remarks>
        public static (double positionAngle, double longitudeCenter, double latitudeCenter)
        Disk (double julianEphemerisDay) {
            double theta = (julianEphemerisDay - 2398220) * 360 / 25.38;
            double I = 7.25.ToRadians();
            double K = (73.6667 + 1.3958333 * (julianEphemerisDay - 2396758) / 36525).ToRadians();
            double lambda =  Body.Position(Bodies.Sun, julianEphemerisDay, false).longitude.ToRadians();
            double epsilon = Sky.ObliquityOfEcliptic(julianEphemerisDay).ToRadians();
            double lambdaPrime = lambda + Sky.Nutation(julianEphemerisDay).longitude.ToRadians();
            double x = Atan(-Cos(lambdaPrime) * Tan(epsilon)).ToDegrees();
            double y = Atan(-Cos(lambda - K) * Tan(I)).ToDegrees();
            double P = x + y;
            double B0 = Asin(Sin(lambda - K) * Sin(I)).ToDegrees();
            double eta = Atan2(-Sin(lambda - K) * Cos(I), -Cos(lambda - K)).ToDegrees();
            double L0 = (eta - theta).To360();
            return (P, B0, L0);
        }

        /// <param name="rotationNumber">Number of desired rotation with Rotation #1 = November 9, 1853.</param>
        /// <returns>Julian Day of rotation start to the nearest hundredth in Universal Time.</returns>
        /// <remarks>p. 191</remarks>
        public static double RotationStart(double rotationNumber) {
            double jde = 2398140.227 + 27.2752316 * rotationNumber;
            double M = (281.96 + 26.882476 * rotationNumber).To360().ToRadians();
            jde += 0.1454 * Sin(M) - 0.0085 * Sin(2 * M) - 0.0141 * Cos(2 * M);
            return Round(jde.JulianUniversalDay(), 2);
        }
    }
}