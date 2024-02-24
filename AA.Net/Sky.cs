using System;
using System.Collections.Generic;
using static System.Math;

namespace AA.Net {
    
    /// <summary>
    /// Methods relating to the appearance of bodies and coordinate systems in Earth's sky.
    /// </summary>
    public static class Sky {

        /// <param name="julianDay">A Julian Day in Universal Time.</param>
        /// <param name="apparent">True (default) for apparent time or False for mean time.</param>
        /// <returns>Mean or apparent sidereal time at Greenwich.</returns>
        /// <remarks>p. 88</remarks>
        public static TimeSpan GreenwichSiderealTime(double julianDay, bool apparent = true) {
            double T = Time.JulianCentury2000(julianDay);
            double theta = (280.46061837 + 360.98564736629 * (julianDay - Time.J2000) + 0.000387933 * T * T - T * T * T / 38710000).To360();
            if (apparent) {
                var nutation = Sky.Nutation(julianDay);
                var epsilon = Sky.ObliquityOfEcliptic(julianDay);
                theta += nutation.longitude * Cos(epsilon.ToRadians());
            }
            return TimeSpan.FromHours(theta / 15);
        }

        /// <param name="julianDay">A Julian Day in Universal Time.</param>
        /// <param name="longitude">Longitude in decimal degrees. Positive west and negative east.</param>
        /// <returns>The local sidereal time in hours, minutes, and seconds.</returns>
        /// <remarks>p. 92</remarks>
        public static TimeSpan LocalSiderealTime(double julianDay, double longitude) {
            double gst = GreenwichSiderealTime(julianDay).ToDegrees();
            return (gst - longitude).To360().ToRightAscension();
        }

        /// <param name="time">A local date and time.</param>
        /// <param name="longitude">Longitude in decimal degrees. Positive west and negative east.</param>
        /// <returns>The local sidereal time in hours, minutes, and seconds.</returns>
        public static TimeSpan LocalSiderealTime(DateTime time, double longitude) {
            return LocalSiderealTime(time.ToUniversalTime().JulianDay(), longitude);
        }
        
        /// <param name="greenwichSiderealTime">Sidereal time at Greenwich.</param>
        /// <param name="longitude">Observer's longitude (positive West, negative East) in decimal degrees.</param>
        /// <param name="rightAscension">Right ascension of object.</param>
        /// <returns>Local hour angle of object in decimal degrees.</returns>
        /// <remarks>p. 92. Note reversal of usual geographic longitude sign.</remarks>
        public static double HourAngle(TimeSpan greenwichSiderealTime, double longitude, TimeSpan rightAscension) {
            return (greenwichSiderealTime - longitude.ToRightAscension() - rightAscension).ToDegrees().To360();
        }

        /// <summary>
        /// Find the horizontal parallax at a specified distance.
        /// </summary>
        /// <param name="distance">AU.</param>
        /// <returns>Parallax in degrees.</returns>
        public static double HorizontalParallax(double distance) {
            return Transform.ToDegrees(0, 0, 8.794) / distance;
        }

        /// <param name="latitude">Observer's latitude in decimal degrees.</param>
        /// <param name="declination">Declination of object in decimal degrees.</param>
        /// <param name="hourAngle">Local hour angle of object in decimal degrees.</param>
        /// <returns>Parallactic angle in decimal degrees.</returns>
        /// <remarks>p. 98</remarks>
        public static double ParallacticAngle(double latitude, double declination, double hourAngle) {
            double haRad = hourAngle.ToRadians();
            double decRad = declination.ToRadians();
            return Atan2(Sin(haRad), Tan(latitude.ToRadians()) * Cos(decRad) - Sin(decRad) * Cos(haRad)).ToDegrees();
        }

        /// <remarks>p. 92</remarks>
        public const double J2000ObliquityOfEcliptic = 23.4392911;
        /// <remarks>p. 92</remarks>
        public const double B1950ObliquityOfEcliptic = 23.4457889;

        /// <param name="julianEphemerisDay">Julian Day in dynamical time., corrected from Universal Time by ΔT.</param>
        /// <param name="apparent">True (default) to return apparent obliquity or False to return mean obliquity.</param>
        /// <returns>Mean or true obliquity of the ecliptic in decimal degrees.</returns>
        /// <remarks>p. 147</remarks>
        public static double ObliquityOfEcliptic(double julianEphemerisDay, bool apparent = true) {
            double U = Time.JulianCentury2000(julianEphemerisDay) / 100;
            if (Abs(U) >= 1) throw new ArgumentOutOfRangeException("julianEphemerisDay", "Years outside -8000 to +12000 are not supported.");
            double epsilon = J2000ObliquityOfEcliptic
                             - 4680.93 / 3600 * U
                             - 1.55 / 3600 * Pow(U, 2)
                             + 1999.25 / 3600 * Pow(U, 3)
                             - 51.38 / 3600 * Pow(U, 4)
                             - 249.67 / 3600 * Pow(U, 5)
                             - 39.05 / 3600 * Pow(U, 6)
                             + 7.12 / 3600 * Pow(U, 7)
                             + 27.87 / 3600 * Pow(U, 8)
                             + 5.79 / 3600 * Pow(U, 9)
                             + 2.45 / 3600 * Pow(U, 10);
            return apparent ? epsilon + Sky.Nutation(julianEphemerisDay).obliquity : epsilon;
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time., corrected from Universal Time by ΔT.</param>
        /// <returns>Nutations in longitude and obliquity in decimal degrees.</returns>
        /// <remarks>pp. 143-146</remarks>
        public static (double longitude, double obliquity)
        Nutation(double julianEphemerisDay) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = Pow(T, 2);
            double T3 = Pow(T, 3);
            double D = (297.85036 + 445267.111480 * T - 0.0019142 * T2 + T3 / 189474);
            double M = (357.52772 + 35999.050340 * T - 0.0001603 * T2 - T3 / 300000);
            double Mprime = (134.96298 + 477198.867398 * T + 0.0086972 * T2 + T3 / 56250);
            double F = (93.27191 + 483202.017538 * T - 0.0036825 * T2 + T3 / 327270);
            double omega = (125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000);
            double[,] tableA = {
                {0,0,0,0,1,-171996,-174.2,92025,8.9},
                {-2,0,0,2,2,-13187,-1.6,5736,-3.1},
                {0,0,0,2,2,-2274,-0.2,977,-0.5},
                {0,0,0,0,2,2062,0.2,-895,0.5},
                {0,1,0,0,0,1426,-3.4,54,-0.1},
                {0,0,1,0,0,712,0.1,-7,0},
                {-2,1,0,2,2,-517,1.2,224,-0.6},
                {0,0,0,2,1,-386,-0.4,200,0},
                {0,0,1,2,2,-301,0,129,-0.1},
                {-2,-1,0,2,2,217,-0.5,-95,0.3},
                {-2,0,1,0,0,-158,0,1,0},
                {-2,0,0,2,1,129,0.1,-70,0},
                {0,0,-1,2,2,123,0,-53,0},
                {2,0,0,0,0,63,0,1,0},
                {0,0,1,0,1,63,0.1,-33,0},
                {2,0,-1,2,2,-59,0,26,0},
                {0,0,-1,0,1,-58,-0.1,32,0},
                {0,0,1,2,1,-51,0,27,0},
                {-2,0,2,0,0,48,0,1,0},
                {0,0,-2,2,1,46,0,-24,0},
                {2,0,0,2,2,-38,0,16,0},
                {0,0,2,2,2,-31,0,13,0},
                {0,0,2,0,0,29,0,1,0},
                {-2,0,1,2,2,29,0,-12,0},
                {0,0,0,2,0,26,0,1,0},
                {-2,0,0,2,0,-22,0,1,0},
                {0,0,-1,2,1,21,0,-10,0},
                {0,2,0,0,0,17,-0.1,1,0},
                {2,0,-1,0,1,16,0,-8,0},
                {-2,2,0,2,2,-16,0.1,7,0},
                {0,1,0,0,1,-15,0,9,0},
                {-2,0,1,0,1,-13,0,7,0},
                {0,-1,0,0,1,-12,0,6,0},
                {0,0,2,-2,0,11,0,1,0},
                {2,0,-1,2,1,-10,0,5,0},
                {2,0,1,2,2,-8,0,3,0},
                {0,1,0,2,2,7,0,-3,0},
                {-2,1,1,0,0,-7,0,1,0},
                {0,-1,0,2,2,-7,0,3,0},
                {2,0,0,2,1,-7,0,3,0},
                {2,0,1,0,0,6,0,1,0},
                {-2,0,2,2,2,6,0,-3,0},
                {-2,0,1,2,1,6,0,-3,0},
                {2,0,-2,0,1,-6,0,3,0},
                {2,0,0,0,1,-6,0,3,0},
                {0,-1,1,0,0,5,0,1,0},
                {-2,-1,0,2,1,-5,0,3,0},
                {-2,0,0,0,1,-5,0,3,0},
                {0,0,2,2,1,-5,0,3,0},
                {-2,0,2,0,1,4,0,1,0},
                {-2,1,0,2,1,4,0,1,0},
                {0,0,1,-2,0,4,0,1,0},
                {-1,0,1,0,0,-4,0,1,0},
                {-2,1,0,0,0,-4,0,1,0},
                {1,0,0,0,0,-4,0,1,0},
                {0,0,1,2,0,3,0,1,0},
                {0,0,-2,2,2,-3,0,1,0},
                {-1,-1,1,0,0,-3,0,1,0},
                {0,1,1,0,0,-3,0,1,0},
                {0,-1,1,2,2,-3,0,1,0},
                {2,-1,-1,2,2,-3,0,1,0},
                {0,0,3,2,2,-3,0,1,0},
                {2,-1,0,2,2,-3,0,1,0}
            };
            double deltaPsi = 0;
            double deltaEpsilon = 0;
            for (int r = 0; r <= tableA.GetUpperBound(0); r++) {
                double sumA = (tableA[r, 0] * D + tableA[r, 1] * M + tableA[r, 2] * Mprime + tableA[r, 3] * F + tableA[r, 4] * omega).ToRadians();
                deltaPsi += (tableA[r, 5] + tableA[r, 6] * T) * Sin(sumA);
                deltaEpsilon += (tableA[r, 7] + tableA[r, 8] * T) * Cos(sumA);
            }
            return (Transform.ToDegrees(0, 0, deltaPsi / 10000), Transform.ToDegrees(0, 0, deltaEpsilon / 10000));
        }

        /// <summary>
        /// Find the aberration at a specified position and time.
        /// </summary>
        /// <param name="rightAscension">HMS.</param>
        /// <param name="declination">Degrees</param>
        /// <param name="julianEphemerisDay"></param>
        /// <param name="isFK5">False if coordinates are supplied in FK4 frame.</param>
        /// <returns>Aberration in right ascension and degrees of declination.</returns>
        /// <remarks>pp. 151-152</remarks>
        public static (TimeSpan rightAscension, double declination)
        Aberration(TimeSpan rightAscension, double declination, double julianEphemerisDay, bool isFK5 = true) {
            const double kappa = 20.49552;
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = T * T;
            double L = Body.Position(Bodies.Sun, julianEphemerisDay, false).longitude.ToRadians();
            double e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T2;
            double pi = (102.93735 + 1.71946 * T + 0.00046 * T2).ToRadians();
            double alpha = rightAscension.ToRadians();
            double delta = declination.ToRadians();
            double epsilon = Sky.ObliquityOfEcliptic(julianEphemerisDay).ToRadians();
            double deltaA = -kappa * (Cos(alpha) * Cos(L) * Cos(epsilon) + Sin(alpha) * Sin(L)) / Cos(delta) + (isFK5 ? e * kappa * (Cos(alpha) * Cos(pi) * Cos(epsilon) + Sin(alpha) * Sin(pi)) / Cos(delta) : 0);
            double deltaD = -kappa * (Cos(L) * Cos(epsilon) * (Tan(epsilon) * Cos(delta) - Sin(alpha) * Sin(delta)) + Cos(alpha) * Sin(delta) * Sin(L)) + (isFK5 ?  e * kappa * (Cos(pi) * Cos(epsilon) * (Tan(epsilon) * Cos(delta) - Sin(alpha) * Sin(delta)) + Cos(alpha) * Sin(delta) * Sin(pi)) : 0);
            return ((deltaA / 3600).ToRightAscension(), deltaD / 3600);
        }

        /// <summary>
        /// Find the aberration at a specified position and time.
        /// </summary>
        /// <param name="longitude">Degrees.</param>
        /// <param name="latitude">Degrees.</param>
        /// <param name="julianEphemerisDay"></param>
        /// <param name="isFK5">False if coordinates are supplied in FK4 frame.</param>
        /// <returns>Aberration in degrees.</returns>
        /// <remarks>pp. 151-152</remarks>
        public static (double longitude, double latitude)
        Aberration(double longitude, double latitude, double julianEphemerisDay, bool isFK5 = true) {
            const double kappa = 20.49552;
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = T * T;
            double L = Body.Position(Bodies.Sun, julianEphemerisDay, false).longitude.ToRadians();
            double e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T2;
            double pi = (102.93735 + 1.71946 * T + 0.00046 * T2).ToRadians();
            double lambdaR = longitude.ToRadians(), betaR = latitude.ToRadians();
            double deltaL = (-kappa * Cos(L - lambdaR) + (isFK5 ? e * kappa * Cos(pi - lambdaR) : 0)) / Cos(betaR);
            double deltaB = -kappa * Sin(betaR) * (Sin(L - lambdaR) - (isFK5 ? e * Sin(pi - lambdaR) : 0));
            return (deltaL / 3600, deltaB / 3600);
        }

        /// <param name="obliquity">Apparent obliquity of the ecliptic.</param>
        /// <param name="latitude">Latitude in decimal degrees. Positive north and negative south.</param>
        /// <param name="siderealTime">Local sidereal time.</param>
        /// <returns>Two ecliptic longitudes which intersect the horizon and the angle formed with the horizon.</returns>
        /// <remarks>p. 99</remarks>
        public static (double longitude1, double longitude2, double angle)
        EclipticAndHorizon(double obliquity, double latitude, TimeSpan siderealTime) {
            double epsilon = obliquity.ToRadians();
            double theta = siderealTime.ToRadians();
            double phi = latitude.ToRadians();
            double l = Atan2(-Cos(theta), Sin(epsilon) * Tan(phi) + Cos(epsilon) * Sin(theta)).ToDegrees();
            double I = Acos(Cos(epsilon) * Sin(phi) - Sin(epsilon) * Cos(phi) * Sin(theta));
            return (l.To360(), (l + 180).To360(), I.ToDegrees());
        }

        /// <param name="time">A local date and time.</param>
        /// <param name="longitude">Longitude in decimal degrees. Positive west and negative east.</param>
        /// <param name="latitude">Latitude in decimal degrees. Positive north and negative south.</param>
        /// <returns>Two ecliptic longitudes which intersect the horizon and the angle formed with the horizon.</returns>
        /// <remarks>p. 99</remarks>
        public static (double longitude1, double longitude2, double angle)
        EclipticAndHorizon(DateTime time, double longitude, double latitude) {
            double jd = time.ToUniversalTime().JulianDay();
            double epsilon = Sky.ObliquityOfEcliptic(jd.JulianEphemerisDay());
            TimeSpan lst = LocalSiderealTime(jd, longitude);
            var result = EclipticAndHorizon(epsilon, latitude, lst);
            return (result.longitude1, result.longitude2, result.angle);
        }

        /// <param name="longitude">Ecliptical longitude in decimal degrees.</param>
        /// <param name="latitude">Ecliptical latitude in decimal degrees.</param>
        /// <param name="obliquity">Obliquity of the ecliptic.</param>
        /// <returns>Angle between the north celestial pole and the north pole of the ecliptic; and angle between the ecliptic and a parallel to the celestial equator.</returns>
        /// <remarks>p. 100</remarks>
        public static (double angleToPole, double angleToEquator)
        EclipticAndEquator(double longitude, double latitude, double obliquity) {
            double longR = longitude.ToRadians();
            double latR = latitude.ToRadians();
            double epsilon = obliquity.ToRadians();
            double q = Atan2(Cos(longR) * Tan(epsilon), Sin(latR) * Sin(longR) * Tan(epsilon) - Cos(latR)).ToDegrees();
            double q0 = Atan(-Cos(longR) * Tan(epsilon)).ToDegrees();
            return (q, q0);
        }

        /// <param name="declination">Equatorial declination.</param>
        /// <param name="latitude">Latitude of observing position in decimal degrees. Positive north and negative south.</param>
        /// <returns>Rising and setting angle of object to the horizon.</returns>
        /// <remarks>p. 100</remarks>
        public static double AngleToHorizon(double declination, double latitude) {
            double decR = declination.ToRadians();
            double latR = latitude.ToRadians();
            double B = Tan(decR) * Tan(latR);
            double C = Sqrt(1 - B * B);
            return Atan2(C * Cos(decR), Tan(latR)).ToDegrees();
        }

        /// <param name="altitude">Altitude of object in degrees.</param>
        /// <param name="apparent">True for apparent altitude or false for true altitude.</param>
        /// <returns>Offset of apparent or true altitude in degrees due to atmospheric refraction.</returns>
        /// <remarks>p. 106</remarks>
        public static double Refraction(double altitude, bool apparent = true) {
            if (altitude == 90) return 0;
            if (apparent) {
                double R = 1 / Tan((altitude + 7.31 / (altitude + 4.4)).ToRadians());
                return R - 0.06 * Sin((14.7 * R + 13).ToRadians());
            }
            return 1.02 / Tan((altitude + 10.3 / (altitude + 5.11)).ToRadians());
        }

        /// <param name="altitude">Altitude of object in degrees.</param>
        /// <param name="pressure">Atmospheric pressure in millibars.</param>
        /// <param name="temperature">Atmospheric temperature in degrees Celsius.</param>
        /// <param name="apparent">True for apparent altitude or false for true altitude.</param>
        /// <returns>Offset of apparent or true altitude in degrees due to atmospheric refraction.</returns>
        /// <remarks>pp. 106-107. This method is accurate for yellow wavelengths only.</remarks>
        public static double Refraction(double altitude, double pressure, double temperature, bool apparent = true) {
            return Refraction(altitude, apparent) * (pressure / 1010) * (283 / (273 + temperature));
        }

        /// <param name="longitude1">First body longitude or right ascension in degrees.</param>
        /// <param name="latitude1">First body latitude or declination.</param>
        /// <param name="longitude2">Second body longitude or right ascension in degrees.</param>
        /// <param name="latitude2">Second body latitude or declination.</param>
        /// <returns>Separation between bodies in degrees.</returns>
        /// <remarks>p. 115</remarks>
        public static double AngularSeparation(double longitude1, double latitude1, double longitude2, double latitude2) {
            double lon1 = longitude1.ToRadians();
            double lon2 = longitude2.ToRadians();
            double lat1 = latitude1.ToRadians();
            double lat2 = latitude2.ToRadians();
            double x = Cos(lat1) * Sin(lat2) - Sin(lat1) * Cos(lat2) * Cos(lon2 - lon1);
            double y = Cos(lat2) * Sin(lon2 - lon1);
            double z = Sin(lat1) * Sin(lat2) + Cos(lat1) * Cos(lat2) * Cos(lon2 - lon1);
            return Atan2(Sqrt(x * x + y * y), z).ToDegrees();
        }

        /// <param name="longitude1">First body longitude or right ascension in degrees.</param>
        /// <param name="latitude1">First body latitude or declination.</param>
        /// <param name="longitude2">Second body longitude or right ascension in degrees.</param>
        /// <param name="latitude2">Second body latitude or declination.</param>
        /// <returns>Relative position angle of first body from second body.</returns>
        /// <remarks>p. 116</remarks>
        public static double PositionAngle (double longitude1, double latitude1, double longitude2, double latitude2) {
            double lonD = (longitude1 - longitude2).ToRadians();
            double lat1 = latitude1.ToRadians();
            double lat2 = latitude2.ToRadians();
            return Atan2(Sin(lonD), Cos(lat2) * Tan(lat1) - Sin(lat2) * Cos(lonD)).ToDegrees();
        }

        /// <param name="start">A Julian Day to begin the chronological search.</param>
        /// <param name="body1">First body.</param>
        /// <param name="body2">Second body.</param>
        /// <returns>Julian Day of the next conjunction and separation in declination.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Returned if Sun and/or Earth are specified.</exception>
        /// <exception cref="ArgumentException">Returned if both bodies are the same.</exception>
        /// <remarks>p. 117 et seq.</remarks>
        public static (double conjunction, double degrees)
        NextConjunction(double start, Bodies body1, Bodies body2) {
            if (body1 == Bodies.Earth || body1 == Bodies.Sun || body2 == Bodies.Earth || body2 == Bodies.Sun) throw new ArgumentOutOfRangeException(body1.ToString(), "Not valid for conjunction.");
            if (body1 == body2) throw new ArgumentException("Bodies must be different.");
            double jde = Floor(start) - 1.5;
            double thisRa1 = 0, lastRa1, thisRa2 = 0, lastRa2;
            double thisRaDiff = 0, lastRaDiff;
            double thisDecDiff;
            //prime the pump
            Update();
            //create a queue to hold last 5 dates and deltas
            (double jde, double raDiff, double decDiff)[] diffs = new(double jde, double raDiff, double decDiff)[5];
            int diffPtr = 0;
            //search ahead while sign of difference remains the same
            do {
                Update();
                diffPtr = ++diffPtr % 5;
                diffs[diffPtr] = (jde, thisRaDiff, thisDecDiff);
            }
            while (Sign(thisRaDiff) == Sign(lastRaDiff));
            //search one more day to place pivot in the center of 5
            Update();
            diffPtr = ++diffPtr % 5;
            diffs[diffPtr] = (jde, thisRaDiff, thisDecDiff);

            //find the zero time through 5-value interpolation
            double zeroN = Interpolation.Zero(diffs[(diffPtr + 1) % 5].raDiff, diffs[(diffPtr + 2) % 5].raDiff, diffs[(diffPtr + 3) % 5].raDiff, diffs[(diffPtr + 4) % 5].raDiff, diffs[diffPtr].raDiff);
            //find the declination difference at zero time
            double deltaDec = Interpolation.Interpolate(diffs[(diffPtr + 1) % 5].decDiff, diffs[(diffPtr + 2) % 5].decDiff, diffs[(diffPtr + 3) % 5].decDiff, diffs[(diffPtr + 4) % 5].decDiff, diffs[diffPtr].decDiff, zeroN);
            //return time and difference in declination
            return (diffs[(diffPtr + 3) % 5].jde + zeroN, deltaDec);

            void Update() {
                jde++;
                lastRa1 = thisRa1; lastRa2 = thisRa2;
                lastRaDiff = thisRaDiff;
                var result1 = Body.PositionGeocentric(body1, jde);
                var result2 = Body.PositionGeocentric(body2, jde);
                thisRa1 = result1.rightAscension.TotalSeconds;
                thisRa2 = result2.rightAscension.TotalSeconds;
                //Was equinox crossed?
                if (lastRa1 > 23 * 3600 && thisRa1 < 3600) thisRa1 += 86400;
                else if (thisRa1 > 23 * 3600 && lastRa1 < 3600) thisRa1 -= 86400;
                if (lastRa2 > 23 * 3600 && thisRa2 < 3600) thisRa2 += 86400;
                else if (thisRa2 > 23 * 3600 && lastRa2 < 3600) thisRa2 -= 86400;
                thisRaDiff = thisRa1 - thisRa2;
                thisDecDiff = result1.declination - result2.declination;
            }
        }

        /// <param name="body">A value of the Bodies enum.</param>
        /// <param name="date">A date in Universal Time.</param>
        /// <param name="longitude">Longitude in decimal degrees. Positive west and negative east.</param>
        /// <param name="latitude">Latitude in decimal degrees. Positive north and negative south.</param>
        /// <param name="chronological">True for events in chronological order, or false for all events occurring on the same day.</param>
        /// <returns>Times of rising, transit, and setting for <paramref name="date"/>.</returns>
        /// <remarks>pp. 102-103</remarks>
        public static (DateTime? rise, DateTime transit, DateTime? set)
        RiseTransitSet(Bodies body, DateTime date, double longitude, double latitude, bool chronological = true) {
            DateTime theDate = date.Date.AddDays(chronological ? 1 : 0);
            double deltaT = Time.DeltaT(theDate.Year, theDate.Month);
            List<(double rightAscension, double declination)> position = new List<(double, double)>();
            double lastRaD = 0;
            for (int i = -1; i < 2; i++) {
                double jde = theDate.AddDays(i).JulianEphemerisDay();
                (TimeSpan rightAscension, double declination) p;
                if (body == Bodies.Sun || body == Bodies.Moon) {
                    var thisPosition = Body.Position(body, jde);
                    p = Transform.ToEquatorial(thisPosition.longitude, thisPosition.latitude, Sky.ObliquityOfEcliptic(jde));
                }
                else p = Body.PositionGeocentric(body, jde);
                double raD = p.rightAscension.ToDegrees();
                if (raD < lastRaD) raD += 360;  //algorithm falls apart if we straddle the equinox
                lastRaD = raD;
                position.Add((raD, p.declination));
            }
            double h0 = body == Bodies.Sun ? -0.8333 : body == Bodies.Moon ? 0.125 : -0.5667;
            double theta0 = GreenwichSiderealTime(theDate.JulianDay()).ToDegrees();
            double dec = position[1].declination.ToRadians();
            double lat = latitude.ToRadians();
            double cosine = (Sin(h0.ToRadians()) - Sin(lat) * Sin(dec)) / (Cos(lat) * Cos(dec));
            double H0 = -1;
            if (cosine >= -1 && cosine <= 1) H0 = Acos(cosine).ToDegrees();
            List<double> m = new List<double>();
            m.Add((position[1].rightAscension + longitude - theta0) / 360);
            if (H0 > -1) {
                m.Add(m[0] - H0 / 360);
                m.Add(m[0] + H0 / 360);
            }
            DateTime transit = default;
            DateTime? rise = null;
            DateTime? set = null;
            for (int i = 0; i < (H0 > -1 ? 3 : 1); i++) {
                double thisM = m[i];
                if (!chronological) {
                    if (thisM > 1) thisM--;
                    else if (thisM < 0) thisM++;
                }
                double thisTheta0 = theta0 + 360.985647 * thisM;
                double n = thisM + deltaT / 86400;
                n = Min(1, Max(-1, n)); //clamp in case of equinox straddling in chrono order
                double alpha = Interpolation.Interpolate(position[0].rightAscension, position[1].rightAscension, position[2].rightAscension, n);
                double delta = Interpolation.Interpolate(position[0].declination, position[1].declination, position[2].declination, n);
                double H = (thisTheta0 - longitude - alpha).ToPlusMinus180();
                if (i == 0) transit = theDate.AddDays(thisM - H / 360);
                else {
                    double h = Transform.ToHorizontal(H, latitude, delta).altitude;
                    double days = thisM + (h - h0) / (360 * Cos(delta.ToRadians()) * Cos(latitude.ToRadians()) * Sin(H.ToRadians()));
                    if (i == 1) rise = theDate.AddDays(days);
                    else set = theDate.AddDays(days);
                }
            }
            return (rise, transit, set);
        }

        /// <param name="julianEphemerisDay">A Julian Day in dynamical time.</param>
        /// <returns>Difference between apparent and mean time in minutes.</returns>
        /// <remarks>p. 183</remarks>
        public static double EquationOfTime(double julianEphemerisDay) {
            double tau = Time.JulianCentury2000(julianEphemerisDay) / 10;
            double tau2 = tau * tau;
            double tau3 = tau2 * tau;
            double tau4 = tau3 * tau;
            double tau5 = tau4 * tau;
            double L0 = (280.4664567 + 360007.6982779 * tau
                        + 0.03032028 * tau2 + tau3 / 49931
                        - tau4 / 15300 - tau5 / 2000000).To360();
            double epsilon = ObliquityOfEcliptic(julianEphemerisDay).ToRadians();
            double deltaLambda = Nutation(julianEphemerisDay).longitude;
            double alpha = Body.PositionGeocentric(Bodies.Sun, julianEphemerisDay).rightAscension.ToDegrees();
            double e = (L0 - 0.0057183 - alpha + deltaLambda * Cos(epsilon)) * 4;
            return Abs(e) < 20 ? e : e - 1440 * Sign(e);
        }

        /// <param name="rightAscension">Right ascension for the given epoch.</param>
        /// <param name="declination">Declination for the given epoch.</param>
        /// <param name="fromEpoch">Epoch of coordinates.</param>
        /// <param name="julianEphemerisDay">A Julian Ephemeris Day at which to calculate apparent position.</param>
        /// <param name="properMotionRA">Annual proper motion in seconds of right ascension.</param>
        /// <param name="properMotionDec">Annual proper motion in arcseconds of declination.</param>
        /// <returns>Right ascension and declination corrected for proper motion, precession, nutation, and annual aberration.</returns>
        /// <remarks>pp. 151-152</remarks>
        public static (TimeSpan rightAscension, double declination)
        ApparentPositionStar (TimeSpan rightAscension, double declination, double fromEpoch, double julianEphemerisDay, double properMotionRA, double properMotionDec) {
            var mean = Transform.Precession(rightAscension, declination, fromEpoch, julianEphemerisDay, properMotionRA, properMotionDec);
            double alpha = mean.rightAscension.ToRadians();
            double delta = mean.declination.ToRadians();
            double epsilon = ObliquityOfEcliptic(julianEphemerisDay).ToRadians();
            var nut = Nutation(julianEphemerisDay);
            double deltaPsi = nut.longitude;
            double deltaEpsilon = nut.obliquity;
            TimeSpan deltaAlpha1 = ((Cos(epsilon) + Sin(epsilon) * Sin(alpha) * Tan(delta)) * deltaPsi - Cos(alpha) * Tan(delta) * deltaEpsilon).ToRightAscension();
            double deltaDelta1 = Sin(epsilon) * Cos(alpha) * deltaPsi + Sin(alpha) * deltaEpsilon;
            (TimeSpan deltaAlpha2, double deltaDelta2) = Aberration(rightAscension, declination, julianEphemerisDay);
            return (mean.rightAscension + deltaAlpha1 + deltaAlpha2, mean.declination + deltaDelta1 + deltaDelta2);
        }

        private static ArgumentOutOfRangeException coordinate3Exception = new ArgumentOutOfRangeException("Exactly three coordinates must be supplied.");

        /// <param name="coordinate">An array of exactly three coordinate tuples.</param>
        /// <param name="tolerance">A tolerance value in degrees.</param>
        /// <returns>True if the coordinates lie in a straight line within the tolerance specified.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if the array does not contain exactly three coordinate tuples.</exception>
        /// <remarks>p. 121</remarks>
        public static bool IsStraightLine((TimeSpan rightAscension, double declination)[] coordinate, double tolerance = 0.0001) {
            if (coordinate.Length != 3) throw coordinate3Exception;
            double a1 = coordinate[0].rightAscension.ToRadians();
            double a2 = coordinate[1].rightAscension.ToRadians();
            double a3 = coordinate[2].rightAscension.ToRadians();
            double d1 = coordinate[0].declination.ToRadians();
            double d2 = coordinate[1].declination.ToRadians();
            double d3 = coordinate[2].declination.ToRadians();
            return Tan(d1) * Sin(a2 - a3) + Tan(d2) * Sin(a3 - a1) + Tan(d3) * Sin(a1 - a2) <= tolerance;
        }

        /// <param name="coordinate">An array of exactly three coordinate tuples.</param>
        /// <param name="tolerance">A tolerance value in degrees.</param>
        /// <returns>True if the coordinates lie in a straight line within the tolerance specified.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if the array does not contain exactly three coordinate tuples.</exception>
        /// <remarks>p. 121</remarks>
        public static bool IsStraightLine((double longitude, double latitude)[] coordinate, double tolerance = 0.0001) {
            if (coordinate.Length != 3) throw coordinate3Exception;
            double a1 = coordinate[0].longitude.ToRadians();
            double a2 = coordinate[1].longitude.ToRadians();
            double a3 = coordinate[2].longitude.ToRadians();
            double d1 = coordinate[0].latitude.ToRadians();
            double d2 = coordinate[1].latitude.ToRadians();
            double d3 = coordinate[2].latitude.ToRadians();
            return Tan(d1) * Sin(a2 - a3) + Tan(d2) * Sin(a3 - a1) + Tan(d3) * Sin(a1 - a2) <= tolerance;
        }

        /// <param name="coordinate">An array of exactly three coordinate tuples. The first coordinate is the one whose distance will be calculated from the great circle defined by the second and third coordinates.</param>
        /// <returns>A distance in decimal degrees.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if the array does not contain exactly three coordinate tuples.</exception>
        /// <remarks>p. 124</remarks>
        public static double DegreesFromGreatCircle((TimeSpan rightAscension, double declination)[] coordinate) {
            if (coordinate.Length != 3) throw coordinate3Exception;
            double a1 = coordinate[1].rightAscension.ToRadians();
            double a2 = coordinate[2].rightAscension.ToRadians();
            double a0 = coordinate[0].rightAscension.ToRadians();
            double d1 = coordinate[1].declination.ToRadians();
            double d2 = coordinate[2].declination.ToRadians();
            double d0 = coordinate[0].declination.ToRadians();
            double X1 = Cos(d1) * Cos(a1);
            double Y1 = Cos(d1) * Sin(a1);
            double Z1 = Sin(d1);
            double X2 = Cos(d2) * Cos(a2);
            double Y2 = Cos(d2) * Sin(a2);
            double Z2 = Sin(d2);
            double A = Y1 * Z2 - Z1 * Y2;
            double B = Z1 * X2 - X1 * Z2;
            double C = X1 * Y2 - Y1 * X2;
            double m = Tan(a0);
            double n = Tan(d0) / Cos(a0);
            return Abs(Asin((A + B * m + C * n) / (Sqrt(A * A + B * B + C * C) * Sqrt(1 + m * m + n * n))).ToDegrees());
        }

        /// <param name="coordinate">An array of exactly three coordinate tuples.</param>
        /// <returns>Diameter, in degrees, of the smallest circle which contains all three coordinates.</returns>
        /// <remarks>pp. 127-128</remarks>
        public static double SmallestDiameter ((TimeSpan rightAscension, double declination)[] coordinate) {
            if (coordinate.Length != 3) throw coordinate3Exception;
            double a0 = coordinate[0].rightAscension.ToRadians(); 
            double a1 = coordinate[1].rightAscension.ToRadians();
            double a2 = coordinate[2].rightAscension.ToRadians();
            double d0 = coordinate[0].declination.ToRadians();
            double d1 = coordinate[1].declination.ToRadians();
            double d2 = coordinate[2].declination.ToRadians();
            double[] separation = new double[3];
            separation[0] = AngularDistance(a0, d0, a1, d1);
            separation[1] = AngularDistance(a0, d0, a2, d2);
            separation[2] = AngularDistance(a1, d1, a2, d2);
            Array.Sort(separation);
            double a = separation[2];
            double b = separation[1];
            double c = separation[0];
            return a > Sqrt(b * b + c * c) ? a : 2 * a * b * c / Sqrt((a + b + c) * (a + b - c) * (b + c - a) * (a + c - b));

            double AngularDistance(double A1, double D1, double A2, double D2) {
                return Acos(Sin(D1) * Sin(D2) + Cos(D1) * Cos(D2) * Cos(A1 - A2)).ToDegrees();
            }
        }
    }
}
