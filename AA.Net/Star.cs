using System;
using System.Collections.Generic;
using System.Text;
using static System.Math;

namespace AA.Net {
    public static class Star {
        
        /// <param name="magnitude">Array of stellar magnitudes.</param>
        /// <returns>Sum of magnitudes.</returns>
        /// <remarks>p. 393</remarks>
        public static double SumMagnitudes(double[] magnitude) {
            double sum = 0;
            for (int i = 0; i < magnitude.Length; i++) sum += Pow(10, -0.4 * magnitude[i]);
            return -2.5 * Log10(sum);
        }

        /// <param name="magnitude1">Magnitude of first star.</param>
        /// <param name="magnitude2">Magnitude of second star.</param>
        /// <returns>Ratio of apparent luminosities.</returns>
        /// <remarks>p. 395</remarks>
        public static double BrightnessRatio(double magnitude1, double magnitude2) {
            return Pow(10, 0.4 * (magnitude2 - magnitude1));
        }

        /// <param name="brightnessRatio">Ratio of apparent luminosities.</param>
        /// <returns>Equivalent difference in magnitude.</returns>
        /// <remarks>p. 395</remarks>
        public static double MagnitudeDifference(double brightnessRatio) {
            return 2.5 * Log10(brightnessRatio);
        }

        /// <param name="parallax">Parallax in arcseconds.</param>
        /// <returns>Distance in parsecs.</returns>
        /// <remarks>p. 396</remarks>
        public static double DistanceParsecs(double parallax) {
            return 1 / parallax;
        }

        /// <param name="parallax">Parallax in arcseconds.</param>
        /// <returns>Distance in lightyears.</returns>
        /// <remarks>p. 396</remarks>
        public static double DistanceLightYears(double parallax) {
            return 3.2616 / parallax;
        }

        /// <param name="parallax">Parallax in arcseconds.</param>
        /// <param name="apparentMagnitude">Apparent magnitude.</param>
        /// <returns>Absolute magnitude.</returns>
        /// <remarks>p. 396</remarks>
        public static double AbsoluteMagnitude(double parallax, double apparentMagnitude) {
            return apparentMagnitude + 5 + 5 * Log10(parallax);
        }

        /// <param name="year">Time of observation in fractional years.</param>
        /// <param name="period">Period of revolution in mean solar years.</param>
        /// <param name="periastronTime">Time of periastron passage in fractional years.</param>
        /// <param name="eccentricity">Eccentricity of true orbit.</param>
        /// <param name="semimajorAxis">Semimajor axis of orbit in arcseconds.</param>
        /// <param name="inclination">Inclination of plane of true orbit to plane at right angles to line of sight in degrees.</param>
        /// <param name="ascendingNode">Position angle of ascending node in degrees.</param>
        /// <param name="periastronLongitude">Longitude of the primary periastron in degrees.</param>
        /// <returns>Apparent position angle in degrees, angular distance in arcseconds, and eccentricity of apparent orbit.</returns>
        /// <remarks>pp. 397-398</remarks>
        public static (double angle, double distance, double apparentEccentricity) BinaryPosition(double year, double period, double periastronTime, double eccentricity, double semimajorAxis, double inclination, double ascendingNode, double periastronLongitude) {
            double n = 360 / period;
            double M = n * (year - periastronTime);
            double E = Body.EccentricAnomaly(M, eccentricity).ToRadians();
            double r = semimajorAxis * (1 - eccentricity * Cos(E));
            double v = Atan(Sqrt((1 + eccentricity) / (1 - eccentricity)) * Tan(E / 2)) * 2;
            double omega = periastronLongitude.ToRadians();
            double i = inclination.ToRadians();
            double sinVo = Sin(v + omega);
            double cosVo = Cos(v + omega);
            double cosI = Cos(i);
            double theta = (Atan2(sinVo * cosI, cosVo).ToDegrees() + ascendingNode).To360();
            double rho = r * Sqrt(sinVo * sinVo * cosI * cosI + cosVo * cosVo);
            double e2 = eccentricity * eccentricity;
            double cosO = Cos(omega);
            double sinO = Sin(omega);
            double A = (1 - e2 * cosO * cosO) * cosI * cosI;
            double B = e2 * sinO * cosO * cosI;
            double C = 1 - e2 * sinO * sinO;
            double D = Pow(A - C, 2) + 4 * B * B;
            double ePrime = Sqrt(2 * Sqrt(D) / (A + C + Sqrt(D)));
            return (theta, rho, ePrime);
        }
    }
}
