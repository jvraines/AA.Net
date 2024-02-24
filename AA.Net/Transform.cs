using System;
using System.ComponentModel;
using System.Reflection;
using static System.Math;

namespace AA.Net {
    public static class Transform {

        /// <param name="degrees">Degrees of arc.</param>
        /// <param name="minutes">Arcminutes.</param>
        /// <param name="seconds">Arcseconds.</param>
        /// <returns>Decimal degrees. If any argument is negative, the result is also negative.</returns>
        public static double ToDegrees(int degrees, int minutes, double seconds) {
            int sign = degrees < 0 || minutes < 0 || seconds < 0 ? -1 : 1;
            return (Abs(degrees) + Abs(minutes) / 60.0 + Abs(seconds) / 3600) * sign;
        }

        /// <param name="rightAscension">Right ascension in H:M:S.</param>
        /// <returns>Decimal degrees.</returns>
        public static double ToDegrees(this TimeSpan rightAscension) {
            return rightAscension.TotalHours * 15;
        }

        /// <param name="degrees">Decimal degrees.</param>
        /// <returns>Degrees of arc, arcminutes, and arcseconds.</returns>
        public static (int degrees, int minutes, double seconds)
        ToDegreesMinutesSeconds(this double degrees) {
            int d = (int)degrees;
            degrees -= d;
            degrees *= 60;
            int m = (int)degrees;
            degrees -= m;
            double s = degrees * 60;
            return (d, m, s);
        }

        /// <returns>Formatted degrees, minutes, and seconds.</returns>
        public static string ToDMSString(this double degrees) {
            var dms = degrees.ToDegreesMinutesSeconds();
            return $"{dms.degrees}°{Abs(dms.minutes)}'{Abs(Round(dms.seconds))}\"";
        }

        /// <returns>Right ascension in H:M:S.</returns>
        public static TimeSpan ToRightAscension(this double degrees) {
            return TimeSpan.FromHours(degrees.To360() / 15);
        }

        /// <returns>Right ascension in H:M:S.</returns>
        public static TimeSpan ToRightAscension(int hours, int minutes, double seconds) {
            return TimeSpan.FromSeconds(hours * 3600 + minutes * 60 + seconds);
        }

        public static double ToRadians(this double degrees) {
            return degrees * PI / 180;
        }

        public static double ToRadians(this TimeSpan hours) {
            return hours.ToDegrees().ToRadians();
        }

        public static double ToDegrees(this double radians) {
            return radians * 180 / PI;
        }

        /// <returns>Angle reduced to range -180 to 180 degrees.</returns>
        public static double ToPlusMinus180(this double degrees) {
            degrees = degrees.To360();
            if (degrees > 180) degrees -= 360;
            return degrees;
        }

       /// <returns>Angle reduced to range -90 to 90 degrees.</returns>
       /// <exception cref="ArgumentOutOfRangeException">The argument cannot be reduced below 180 degrees.</exception>
        public static double ToPlusMinus90(this double degrees) {
            degrees = degrees.To360();
            if (degrees > 180) throw new ArgumentOutOfRangeException("degrees", "Must reduce to an angle from 0 to 180.");
            return 90 - degrees;
        }

        /// <returns>Reduce to a positive angle less than or equal to 360.</returns>
        public static double To360(this double degrees) {
            if ((degrees %= 360) < 0) degrees += 360;
            return degrees;
        }

        /// <returns>Reduce to an angle from -360 to 360 degrees.</returns>
        public static double ToPlusMinus360(this double degrees) {
            return degrees % 360;
        }

        /// <returns>Ecliptical coordinates from equatorial coordinates.</returns>
        /// <remarks>p. 93</remarks>
        public static (double longitude, double latitude)
        ToEcliptical(TimeSpan rightAscension, double declination, double obliquityOfEcliptic) {
            double alpha = rightAscension.ToRadians();
            double delta = declination.ToRadians();
            double epsilon = obliquityOfEcliptic.ToRadians();
            double longitude = Atan2(Sin(alpha) * Cos(epsilon) + Tan(delta) * Sin(epsilon), Cos(alpha));
            double latitude = Asin(Sin(delta) * Cos(epsilon) - Cos(delta) * Sin(epsilon) * Sin(alpha));
            return (longitude.ToDegrees().To360(), latitude.ToDegrees());
        }

        /// <returns>Equatorial coordinates from ecliptical coordinates.</returns>
        /// <remarks>p. 93</remarks>
        public static (TimeSpan rightAscension, double declination)
        ToEquatorial(double longitude, double latitude, double obliquityOfEcliptic) {
            double lambda = longitude.ToRadians();
            double beta = latitude.ToRadians();
            double epsilon = obliquityOfEcliptic.ToRadians();
            double alpha = Atan2(Sin(lambda) * Cos(epsilon) - Tan(beta) * Sin(epsilon), Cos(lambda)).ToDegrees();
            double delta = Asin(Sin(beta) * Cos(epsilon) + Cos(beta) * Sin(epsilon) * Sin(lambda)).ToDegrees();
            return (alpha.ToRightAscension(), delta);
        }

        /// <returns>Spherical coordinates from rectangular coordinates.</returns>
        public static (double theta, double phi, double rho)
        ToSpherical(double x, double y, double z) {
            double rho = Sqrt(x * x + y * y + z * z);
            double theta = Atan2(y, x);
            double phi = Acos(z / rho);
            return (theta.ToDegrees().To360(), phi.ToDegrees().ToPlusMinus90(), rho);
        }

        /// <returns>Rectangular coordinates from spherical coordinates.</returns>
        /// <remarks>p. 172</remarks>
        public static(double x, double y, double z)
        ToRectangular(double theta, double phi, double rho) {
            double lambda = theta.ToRadians();
            double beta = phi.ToRadians();
            double x = rho * Cos(beta) * Cos(lambda);
            double y = rho * Cos(beta) * Sin(lambda);
            double z = rho * Sin(beta);
            return (x, y, z);
        }

        /// <returns>Horizontal coordinates from local hour angle, latitude, and declination.</returns>
        /// <remarks>p. 93. Azimuth is measured westward from the South.</remarks>
        public static (double azimuth, double altitude)
        ToHorizontal(double hourAngle, double latitude, double declination) {
            double H = hourAngle.ToRadians();
            double phi = latitude.ToRadians();
            double delta = declination.ToRadians();
            double azimuth = Atan2(Sin(H), Cos(H) * Sin(phi) - Tan(delta) * Cos(phi));
            double altitude = Sin(phi) * Sin(delta) + Cos(phi) * Cos(delta) * Cos(H);
            return (azimuth.ToDegrees().To360(), altitude.ToDegrees());
        }

        /// <returns>Hour angle and declination from horizontal coordinates and latitude.</returns>
        /// <remarks>p. 94</remarks>
        public static (double hourAngle, double declination)
        ToEquatorialWestward(double azimuth, double altitude, double latitude) {
            double A = azimuth.ToRadians();
            double phi = latitude.ToRadians();
            double h = altitude.ToRadians();
            double hourAngle = Atan2(Sin(A), Cos(A) * Sin(phi) + Tan(h) * Cos(phi));
            double declination = Sin(phi) * Sin(h) - Cos(phi) * Cos(h) * Cos(A);
            return (hourAngle.ToDegrees(), declination.ToDegrees());
        }

        /// <param name="rightAscension">Right ascension in the From epoch.</param>
        /// <param name="declination">Declination in the From epoch.</param>
        /// <param name="fromEpoch">Julian Day of the epoch to convert From.</param>
        /// <param name="toEpoch">Julian Day of the epoch to convert To.</param>
        /// <param name="properMotionRightAscension">Annual proper motion of the body in seconds of right ascension.</param>
        /// <param name="properMotionDeclination">Annual proper motion of the body in arcseconds of declination.</param>
        /// <returns>Right ascension and declination precessed to the To epoch.</returns>
        /// <remarks>pp. 134-135</remarks>
        public static (TimeSpan rightAscension, double declination)
        Precession(TimeSpan rightAscension, double declination, double fromEpoch, double toEpoch, double properMotionRightAscension = 0, double properMotionDeclination = 0) {
            double T = Time.JulianCentury2000(fromEpoch);
            double T2 = T * T;
            double t = (toEpoch - fromEpoch) / 36525;
            double t2 = t * t;
            double t3 = t2 * t;
            double raDegrees = rightAscension.ToDegrees();
            if (properMotionDeclination != 0 || properMotionRightAscension != 0) {
                double tYears = t * 100;
                raDegrees += tYears * ToDegrees(0, 0, properMotionRightAscension * 15);
                declination += tYears * ToDegrees(0, 0, properMotionDeclination);
            }
            double zeta = ToDegrees(0, 0, (2306.2181 + 1.39656 * T - 0.000139 * T2) * t + (0.30188 - 0.000344 * T) * t2 + 0.017998 * t3).ToRadians();
            double z = ToDegrees(0, 0, (2306.2181 + 1.39656 * T - 0.000139 * T2) * t + (1.09468 + 0.000066 * T) * t2 + 0.018203 * t3).ToRadians();
            double theta = ToDegrees(0, 0, (2004.3109 - 0.85330 * T - 0.000217 * T2) * t - (0.42665 + 0.000217 * T) * t2 - 0.041833 * t3).ToRadians();
            double raZeta = raDegrees.ToRadians() + zeta;
            declination = declination.ToRadians();
            double A = Cos(declination) * Sin(raZeta);
            double B = Cos(theta) * Cos(declination) * Cos(raZeta) - Sin(theta) * Sin(declination);
            double C = Sin(theta) * Cos(declination) * Cos(raZeta) + Cos(theta) * Sin(declination);
            TimeSpan ra = (Atan2(A, B) + z).ToDegrees().ToRightAscension();
            double dec = Asin(C).ToDegrees();
            return (ra, dec);
        }

        /// <param name="longitude">Ecliptical longitude in the From epoch.</param>
        /// <param name="latitude">Ecliptical latitude in the From epoch.</param>
        /// <param name="fromEpoch">Julian Day of the epoch to convert From.</param>
        /// <param name="toEpoch">Julian Day of the epoch to convert Fo.</param>
        /// <returns>Longitude and latitude precessed to the To epoch.</returns>
        /// <remarks>pp. 136-137</remarks>
        public static (double longitude, double latitude)
        Precession(double longitude, double latitude, double fromEpoch, double toEpoch) {
            double T = Time.JulianCentury2000(fromEpoch);
            double T2 = T * T;
            double t = (toEpoch - fromEpoch) / 36525;
            double t2 = t * t;
            double t3 = t2 * t;
            double eta = ToDegrees(0, 0, (47.0029 - 0.06603 * T + 0.000598 * T2) * t + (-0.03302 + 0.000598 * T) * t2 + 0.00006 * t3).ToRadians();
            double pi = (174.876384 + ToDegrees(0, 0, 3289.4789 * T + 0.60622 * T2 - (869.8089 + 0.50491 * T) * t + 0.03536 * t2)).ToRadians();
            double rho = ToDegrees(0, 0, (5029.0966 + 2.22226 * T - 0.000042 * T2) * t + (1.11113 - 0.000042 * T) * t2 - 0.000006 * t3).ToRadians();
            double lambda = longitude.ToRadians();
            double beta = latitude.ToRadians();
            double A = Cos(eta) * Cos(beta) * Sin(pi - lambda) - Sin(eta) * Sin(beta);
            double B = Cos(beta) * Cos(pi - lambda);
            double C = Cos(eta) * Sin(beta) + Sin(eta) * Cos(beta) * Sin(pi - lambda);
            double lat = -(Atan2(A, B) - rho - pi).ToDegrees();
            double lon = Asin(C).ToDegrees();
            return (lat, lon);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Coordinates converted to FK5 system.</returns>
        /// <remarks>p. 219</remarks>
        public static (double longitude, double latitude)
        ToFK5(double longitude, double latitude, double julianEphemerisDay) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double T2 = T * T;
            double LPrime = (longitude - 1.397 * T - 0.00031 * T2).ToRadians();
            double L = longitude - ToDegrees(0, 0, 0.09033) + ToDegrees(0, 0, 0.03916) * (Cos(LPrime) + Sin(LPrime)) * Tan(latitude.ToRadians());
            double B = latitude + ToDegrees(0, 0, 0.03916) * (Cos(LPrime) - Sin(LPrime));
            return (L, B);
        }

        /// <param name="observerElevation">In meters.</param>
        /// <param name="parallax">Horizontal parallax in degrees.</param>
        /// <returns>Coordinates corrected for topocentric parallax.</returns>
        /// <remarks>p. 279</remarks>
        public static (TimeSpan rightAscension, double declination)
        ToTopocentric(TimeSpan rightAscension, double declination, double parallax, double observerLongitude, double observerLatitude, double observerElevation, double julianDay) {
            (double S, double C) = Earth.RhoPhiPrimeQuantities(observerLatitude, observerElevation);
            double pi = parallax.ToRadians();
            double H = Sky.HourAngle(Sky.GreenwichSiderealTime(julianDay), observerLongitude, rightAscension).ToRadians();
            double sinPi = Sin(pi);
            double delta = declination.ToRadians();
            double deltaAlpha = Atan2(-C * sinPi * Sin(H), Cos(delta) - C * sinPi * Cos(H));
            TimeSpan alphaPrime = (rightAscension.ToDegrees() + deltaAlpha.ToDegrees()).ToRightAscension();
            double deltaPrime = Atan2((Sin(delta) - S * sinPi) * Cos(deltaAlpha), Cos(delta) - C * sinPi * Cos(H)).ToDegrees();
            return (alphaPrime, deltaPrime);
        }

        /// <param name="distance">In AU.</param>
        /// <param name="observerElevation">In meters.</param>
        /// <returns>Coordinates corrected for parallax.</returns>
        /// <remarks>p. 279</remarks>
        public static (double longitude, double latitude, double semidiameter)
        ToTopocentric(double longitude, double latitude, double parallax, double semidiameter, double observerLongitude, double observerLatitude, double observerElevation, double julianDay) {
            double lambda = longitude.ToRadians();
            double beta = latitude.ToRadians();
            double s = semidiameter.ToRadians();
            double epsilon = Sky.ObliquityOfEcliptic(julianDay).ToRadians();
            double theta = Sky.LocalSiderealTime(julianDay, observerLongitude).ToRadians();
            double pi = parallax.ToRadians();
            (double S, double C) = Earth.RhoPhiPrimeQuantities(observerLatitude, observerElevation);
            double N = Cos(lambda) * Cos(beta) - C * Sin(pi) * Cos(theta);
            double lambdaPrime = Atan2(Sin(lambda) * Cos(beta) - Sin(pi) * (S * Sin(epsilon) + C * Cos(epsilon) * Sin(theta)), N);
            double betaPrime = Atan(Cos(lambdaPrime) * (Sin(beta) - Sin(pi) * (S * Cos(epsilon) - C * Sin(epsilon) * Sin(theta))) / N);
            double sPrime = Asin(Cos(lambdaPrime) * Cos(betaPrime) * Sin(s) / N);
            return (lambdaPrime.ToDegrees().To360(), betaPrime.ToDegrees(), sPrime.ToDegrees());
        }
    }
}