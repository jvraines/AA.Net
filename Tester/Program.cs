using System;
using System.Reflection;
using AA.Net;

namespace Tester {
    public class Program {
        public struct Eclipse {
            public DateTime begins;
            public DateTime ends;
            public double gamma;
            public double sigma;
            public bool isAscending;
        }

        static void Main(string[] args) {
            DateTime midnight = DateTime.Parse("13 Feb 2024 21:03:24");
            double jd = midnight.JulianEphemerisDay();
            Console.WriteLine((midnight - TimeSpan.FromMinutes(Sky.EquationOfTime(jd))).ToLocalTime());
        }
    }
}
 