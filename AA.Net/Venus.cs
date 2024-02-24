using System;
using System.Collections.Generic;

namespace AA.Net {
    public static class Venus {
        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Venus, julianEphemerisDay);
        }

        public static (TimeSpan rightAscension, double declination) PositionGeocentric(double julianEphemerisDay) {
            return Body.PositionGeocentric(Bodies.Venus, julianEphemerisDay);
        }

        public static double NextPerihelion(double start) {
            return Body.NextApsis(Bodies.Venus, Apsis.Perihelion, start);
        }

        public static double NextAphelion(double start) {
            return Body.NextApsis(Bodies.Venus, Apsis.Aphelion, start);
        }

        public static (double fraction, double positionAngle, double magnitude) Illumination(double julianEphemerisDay) {
            return Body.Illumination(Bodies.Venus, julianEphemerisDay);
        }

        public static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(double julianEphemerisDay, bool J2000 = false) {
            return Body.OrbitalElements(Bodies.Venus, julianEphemerisDay, J2000);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Venus, julianEphemerisDay);
        }
    }
}
