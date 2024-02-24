using System;
using static System.Math;
using System.Collections.Generic;
using System.Text;

namespace AA.Net {
    public static class Mars  {
        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Mars, julianEphemerisDay);
        }

        public static (TimeSpan rightAscension, double declination) PositionGeocentric(double julianEphemerisDay) {
            return Body.PositionGeocentric(Bodies.Mars, julianEphemerisDay);
        }

        public static double NextPerihelion(double start) {
            return Body.NextApsis(Bodies.Mars, Apsis.Perihelion, start);
        }

        public static double NextAphelion(double start) {
            return Body.NextApsis(Bodies.Mars, Apsis.Aphelion, start);
        }

        public static (double fraction, double positionAngle, double magnitude) Illumination(double julianEphemerisDay) {
            return Body.Illumination(Bodies.Mars, julianEphemerisDay);
        }

        public static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(double julianEphemerisDay, bool J2000 = false) {
            return Body.OrbitalElements(Bodies.Mars, julianEphemerisDay, J2000);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Mars, julianEphemerisDay);
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>
        /// <list type="table">
        /// <item><term>declinationEarth</term>
        /// <description>Areocentric declination of Earth. When positive, the north pole is tilted toward Earth.</description>
        /// </item>
        /// <item><term>declinationSun</term>
        /// <description>Areocentric declination of Sun. When positive, the north pole is illuminated.</description>
        /// </item>
        /// <item><term>positionAngle</term>
        /// <description>Eastward angle of the northern rotation pole from the North Point of the disk.</description>
        /// </item>
        /// <item><term>defectAmount</term>
        /// <description>Greatest defect of illumination in arcseconds.</description>
        /// </item>
        /// <item><term>defectAngle</term>
        /// <description>Position angle of greatest defect of illumination.</description>
        /// </item>
        /// <item><term>longitudeMeridian</term>
        /// <description>Areographic longitude of the central meridian.</description>
        /// </item>
        /// </list>
        /// </returns>
        /// <remarks>pp. 287-290</remarks>
        public static (double declinationEarth, double declinationSun, double positionAngle, double defectAmount, double defectAngle, double longitudeMeridian)
        Disk(double julianEphemerisDay) {
            double T = Time.JulianCentury2000(julianEphemerisDay);
            double lambda0 = (352.9065 + 1.1733 * T).ToRadians();
            double beta0 = (63.2818 - 0.00394 * T).ToRadians();
            var poz = Body.PositionFromLightTime(Bodies.Mars, julianEphemerisDay);
            double De = Asin(-Sin(beta0) * Sin(poz.beta) - Cos(beta0) * Cos(poz.beta) * Cos(lambda0 - poz.lambda)).ToDegrees();
            double N = 49.5581 + 0.77121 * T;
            double lPrime = (poz.l - 0.00697 / poz.r).ToRadians();
            double bPrime = (poz.b - 0.000225 * (Cos((poz.l - N).ToRadians()) / poz.r)).ToRadians();
            double Ds = (-Sin(beta0) * Sin(bPrime) - Cos(beta0) * Cos(bPrime) * Cos(lambda0 - lPrime)).ToDegrees();
            double W = 11.504 + 350.89200025 * (julianEphemerisDay - poz.tau - 2433282.5);
            double epsilon0 = Sky.ObliquityOfEcliptic(julianEphemerisDay, false);
            var equ = Transform.ToEquatorial(lambda0.ToDegrees(), beta0.ToDegrees(), epsilon0);
            double alpha0 = equ.rightAscension.ToRadians();
            double delta0 = equ.declination.ToRadians();
            epsilon0 = epsilon0.ToRadians();
            double u = poz.y * Cos(epsilon0) - poz.z * Sin(epsilon0);
            double v = poz.y * Sin(epsilon0) - poz.z * Cos(epsilon0);
            double alpha = Atan2(u, poz.x);
            double delta = Atan(v / Sqrt(poz.x * poz.x + u * u));
            double zeta = Atan2(Sin(delta0) * Cos(delta) * Cos(alpha0 - alpha) - Sin(delta) * Cos(delta0), Cos(delta) * Sin(alpha0 - alpha));
            double omega = (W - zeta.ToDegrees()).To360();
            var nut = Sky.Nutation(julianEphemerisDay);
            poz.l0 = poz.l0.ToRadians();
            double corrLambda = 0.005693.ToRadians() * Cos(poz.l0 - poz.lambda) / Cos(poz.beta);
            double corrBeta = 0.005693.ToRadians() * Sin(poz.l0 - poz.lambda) * Sin(poz.beta);
            poz.lambda += corrLambda;
            poz.beta += corrBeta;
            nut.longitude = nut.longitude.ToRadians();
            lambda0 += nut.longitude;
            poz.lambda += nut.longitude;
            double epsilon = epsilon0.ToDegrees() + nut.obliquity;
            equ = Transform.ToEquatorial(lambda0.ToDegrees(), beta0.ToDegrees(), epsilon);
            double alpha0Prime = equ.rightAscension.ToRadians();
            double delta0Prime = equ.declination.ToRadians();
            equ = Transform.ToEquatorial(poz.lambda.ToDegrees(), poz.beta.ToDegrees(), epsilon);
            double alphaPrime = equ.rightAscension.ToRadians();
            double deltaPrime = equ.declination.ToRadians();
            double P = Atan2(Cos(delta0Prime) * Sin(alpha0Prime - alphaPrime), Sin(delta0Prime) * Cos(deltaPrime) - Cos(delta0Prime) * Sin(deltaPrime) * Cos(alpha0Prime - alphaPrime)).ToDegrees().To360();
            var sun = Body.Position(Bodies.Sun, julianEphemerisDay);
            var sunEq = Transform.ToEquatorial(sun.longitude, sun.latitude, epsilon);
            double chi = Body.BrightLimbPositionAngle(sunEq.rightAscension.ToRadians(), sunEq.declination.ToRadians(), alphaPrime, deltaPrime);
            double Q = (chi + 180).To360();
            double d = 9.36 / poz.delta;
            double q = (1 - Body.Illumination(Bodies.Mars, julianEphemerisDay).fraction) * d;
            return (De, Ds, P, q, Q, omega);
        }
    }
}
