using System;
using System.Collections.Generic;
using System.Text;

namespace AA.Net {
    public static class Uranus {
        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Uranus, julianEphemerisDay);
        }

        public static (TimeSpan rightAscension, double declination) PositionGeocentric(double julianEphemerisDay) {
            return Body.PositionGeocentric(Bodies.Uranus, julianEphemerisDay);
        }

        public static double NextPerihelion(double start) {
            return Body.NextApsis(Bodies.Uranus, Apsis.Perihelion, start);
        }

        public static double NextAphelion(double start) {
            return Body.NextApsis(Bodies.Uranus, Apsis.Aphelion, start);
        }

        public static (double fraction, double positionAngle, double magnitude) Illumination(double julianEphemerisDay) {
            return Body.Illumination(Bodies.Uranus, julianEphemerisDay);
        }

        public static (double meanLongitude, double semimajorAxis, double eccentricity, double inclination, double longitudeAscendingNode, double longitudePerihelion)
        OrbitalElements(double julianEphemerisDay, bool J2000 = false) {
            return Body.OrbitalElements(Bodies.Uranus, julianEphemerisDay, J2000);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Uranus, julianEphemerisDay);
        }
    }
}
