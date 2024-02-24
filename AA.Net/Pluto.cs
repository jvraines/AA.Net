using System;
using static System.Math;
using System.Collections.Generic;
using System.Text;

namespace AA.Net {
    public static class Pluto {

        public static (double longitude, double latitude, double distance) Position(double julianEphemerisDay) {
            return Body.Position(Bodies.Pluto, julianEphemerisDay);
        }

        public static (TimeSpan rightAscension, double declination) PositionGeocentric(double julianEphemerisDay) {
            return Body.PositionGeocentric(Bodies.Pluto, julianEphemerisDay);
        }

        public static (double equatorial, double polar) Semidiameter(double julianEphemerisDay) {
            return Body.Semidiameter(Bodies.Pluto, julianEphemerisDay);
        }
    }
}

