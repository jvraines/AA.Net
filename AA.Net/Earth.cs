using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Reflection;
using System.Runtime.CompilerServices;
using static System.Math;

namespace AA.Net {
    /// <summary>
    /// Constants and methods pertaining to Earth quantities and phenomena.
    /// </summary>
    public static class Earth {
        
        /// <summary>
        /// Astronomical seasons beginning with the March equinox, the June solstice, the September equinox, and the December solstice.
        /// </summary>
        public enum Season {
            March,
            June,
            September,
            December
        }
        
        //IERS updated values 2010
        public const double RadiusEquatorial = 6378.1366;
        public const double RadiusPolar = 6356.752;
        public const double RadiusMean = (RadiusEquatorial + RadiusPolar) / 2;
        public const double Flattening = 1 / 298.25642;
        public static readonly double EccentricityMeridian = Sqrt(2 * Flattening - Flattening * Flattening);
        //NASA updated value 2022
        public const double AngularVelocity = 7.292124E-5;

        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Earth, julianEphemerisDay);
        }

        public static double NextPerihelion(double start) {
            return Body.NextApsis(Bodies.Earth, Apsis.Perihelion, start);
        }

        public static double NextAphelion(double start) {
            return Body.NextApsis(Bodies.Earth, Apsis.Aphelion, start);
        }

        public static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(double julianEphemerisDay, bool J2000 = false) {
            return Body.OrbitalElements(Bodies.Earth, julianEphemerisDay, J2000);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Earth, julianEphemerisDay);
        }

        /// <param name="year">A year between -1000 and +3000.</param>
        /// <param name="season">An <cref>Earth.Season</cref> value.</param>
        /// <returns>Julian Ephemeris Day of the beginning of the season.</returns>
        /// <remarks>p. 178</remarks>
        public static double SeasonStart(int year, Season season) {
            //if (year < -1000 || year > 3000) throw new ArgumentOutOfRangeException("year", "Year cannot be less than -1000 or greater than +3000.");
            double JDE = 0;
            if (year < 1000) {
                double Y = year / 1000d;
                double Y2 = Y * Y;
                double Y3 = Y2 * Y;
                double Y4 = Y3 * Y;
                switch (season) {
                    case Season.March:
                        JDE = 1721139.29189 + 365242.1374 * Y + 0.06134 * Y2 + 0.00111 * Y3 - 0.00071 * Y4;
                        break;
                    case Season.June:
                        JDE = 1721233.24501 + 365241.72562 * Y - 0.05323 * Y2 + 0.00907 * Y3 + 0.00025 * Y4;
                        break;
                    case Season.September:
                        JDE = 1721325.70455 + 365242.49558 * Y - 0.11677 * Y2 - 0.00297 * Y3 + 0.00074 * Y4;
                        break;
                    case Season.December:
                        JDE = 1721414.39987 + 365242.88257 * Y - 0.00769 * Y2 - 0.00933 * Y3 - 0.00006 * Y4;
                        break;
                }
            }
            else {
                double Y = (year - 2000) / 1000d;
                double Y2 = Y * Y;
                double Y3 = Y2 * Y;
                double Y4 = Y3 * Y;
                switch (season) {
                    case Season.March:
                        JDE = 2451623.80984 + 365242.37404 * Y + 0.05169 * Y2 - 0.00411 * Y3 - 0.00057 * Y4;
                        break;
                    case Season.June:
                        JDE = 2451716.56767 + 365241.62603 * Y + 0.00325 * Y2 + 0.00888 * Y3 - 0.0003 * Y4;
                        break;
                    case Season.September:
                        JDE = 2451810.21715 + 365242.01767 * Y - 0.11575 * Y2 + 0.00337 * Y3 + 0.00078 * Y4;
                        break;
                    case Season.December:
                        JDE = 2451900.05952 + 365242.74049 * Y - 0.06223 * Y2 - 0.00823 * Y3 + 0.00032 * Y4;
                        break;
                }
            }
            double lambda, error;
            double target = (int)season * 90;
            do {
                lambda = Sun.Position(JDE).longitude;
                error = target - lambda;
                JDE += 58 * Sin(error.ToRadians());
            } while (Abs(error) > 0.000005);

            return JDE;
        }

        /// <param name="year">The year for which to calculate the date of Easter. Values less than 1583 will be treated as years in the Julian calendar.</param>
        /// <returns>A DateTime representing the date of Easter as midnight in the corresponding calendar.</returns>
        /// <remarks>pp. 67-69</remarks>
        public static DateTime Easter(int year) {
            int n, p, x;
            if (year > 1582) {
                int a = year % 19;
                int b = year / 100;
                int c = year % 100;
                int d = b / 4;
                int e = b % 4;
                int f = (b + 8) / 25;
                int g = (b - f + 1) / 3;
                int h = (19 * a + b - d - g + 15) % 30;
                int i = c / 4;
                int k = c % 4;
                int l = (32 + 2 * e + 2 * i - h - k) % 7;
                int m = (a + 11 * h + 22 * l) / 451;
                x = h + l - 7 * m + 114;
            }
            else {
                int a = year % 4;
                int b = year % 7;
                int c = year % 19;
                int d = (19 * c + 15) % 30;
                int e = (2 * a + 4 * b - d + 34) % 7;
                x = d + e + 114;
            }
            n = x / 31;
            p = x % 31;
            return new DateTime(year, n, p + 1);
        }

        /// <param name="latitude">A latitude in decimal degrees. North is positive and south is negative.</param>
        /// <returns>Radius of the parallel passing through <paramref name="latitude"/> in kilometers.</returns>
        /// <remarks>p. 83</remarks>
        public static double RadiusParallel(double latitude) {
            double latR = latitude.ToRadians();
            double sLat = Sin(latR);
            return RadiusEquatorial * Cos(latR) / Sqrt(1 - EccentricityMeridian * EccentricityMeridian * sLat * sLat);
        }

        /// <param name="latitude">A latitude in decimal degrees. North is positive and south is negative.</param>
        /// <returns>Length of one degree of longitude at <paramref name="latitude"/> in kilometers.</returns>
        /// <remarks>p. 83</remarks>
        public static double LengthDegreeLongitude(double latitude) {
            return PI / 180 * RadiusParallel(latitude);
        }

        /// <param name="latitude">A latitude in decimal degrees. North is positive and south is negative.</param>
        /// <returns>Radius of the meridian at <paramref name="latitude"/> in kilometers.</returns>
        /// <remarks>p. 83</remarks>
        public static double RadiusMeridian(double latitude) {
            double sLat = Sin(latitude.ToRadians());
            return RadiusEquatorial * (1 - EccentricityMeridian * EccentricityMeridian) / Pow(1 - EccentricityMeridian * EccentricityMeridian * sLat * sLat, 1.5);
        }

        /// <param name="latitude">A latitude in decimal degrees. North is positive and south is negative.</param>
        /// <returns>Length of one degree of latitude at <paramref name="latitude"/> in kilometers.</returns>
        /// <remarks>p. 84</remarks>
        public static double LengthDegreeLatitude(double latitude) {
            return PI / 180 * RadiusMeridian(latitude);
        }

        /// <param name="latitude">A latitude in decimal degrees. North is positive and south is negative.</param>
        /// <returns>Linear rotational velocity at <paramref name="latitude"/> in kilometers/second.</returns>
        /// <remarks>p. 84</remarks>
        public static double LinearVelocity(double latitude) {
            return AngularVelocity * RadiusParallel(latitude);
        }

        /// <param name="longitude1">Starting longitude in decimal degrees. West is positive and east is negative.</param>
        /// <param name="latitude1">Starting latitude in decimal degrees. North is positive and south is negative.</param>
        /// <param name="longitude2">Ending longitude in decimal degrees. West is positive and east is negative.</param>
        /// <param name="latitude2">Ending latitude in decimal degrees. North is positive and south is negative.</param>
        /// <returns>Distance between starting and ending points in kilometers.</returns>
        /// <remarks>p. 85</remarks>
        public static double GeodesicDistance(double longitude1, double latitude1, double longitude2, double latitude2) {
            double F = ((latitude1 + latitude2) / 2).ToRadians();
            double G = ((latitude1 - latitude2) / 2).ToRadians();
            double lambda = ((longitude1 - longitude2) / 2).ToRadians();
            double sin2G = Sin(G) * Sin(G);
            double cos2G = Cos(G) * Cos(G);
            double sin2F = Sin(F) * Sin(F);
            double cos2F = Cos(F) * Cos(F);
            double sin2Lambda = Sin(lambda) * Sin(lambda);
            double cos2Lambda = Cos(lambda) * Cos(lambda);
            double S = sin2G * cos2Lambda + cos2F * sin2Lambda;
            double C = cos2G * cos2Lambda + sin2F * sin2Lambda;
            double omega = Atan(Sqrt(S / C));
            double R = Sqrt(S * C) / omega;
            double D = 2 * omega * RadiusEquatorial;
            double H1 = (3 * R - 1) / (2 * C);
            double H2 = (3 * R + 1) / (2 * S);
            return D * (1 + Flattening * H1 * sin2F * cos2G - Flattening * H2 * cos2F * sin2G);
        }

        /// <param name="distance">A distance in AU.</param>
        /// <returns>Time in days for light to travel <paramref name="distance"/>.</returns>
        /// <remarks>p. 224</remarks>
        public static double LightTime(double distance) {
            return distance * 0.0057755183;
        }

        /// <param name="latitude">A latitude in decimal degrees.</param>
        /// <param name="elevation">An elevation above sea level in meters.</param>
        /// <returns>S = ρ * Sin(φ′), C = ρ * Cos(φ′)</returns>
        /// <remarks>p. 82</remarks>
        public static (double S, double C) RhoPhiPrimeQuantities (double latitude, double elevation) {
            double phi = latitude.ToRadians();
            double u = Atan(RadiusPolar / RadiusEquatorial * Tan(phi));
            double S = RadiusPolar / RadiusEquatorial * Sin(u) + elevation / 6378140 * Sin(phi);
            double C = Cos(u) + elevation / 6378140 * Cos(phi);
            return (S, C);
        }
    }
}
