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
            var foo = Time.ToJulian(new DateTime(1582, 10, 15));
            Console.WriteLine(Time.ToGregorian(foo.year, foo.month, foo.day));
        }
    }
}
 