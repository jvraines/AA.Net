
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;

namespace AA.Net {
    public static class Interpolation {

        private static ArgumentOutOfRangeException BadFactor = new ArgumentOutOfRangeException("n", "Interpolation factor must be between -1 and 1.");
        private const double Threshold = 0.000000000000001;

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <param name="n">Interpolation factor.</param>
        /// <returns>Interpolated value.</returns>
        public static double Interpolate(double y1, double y2, double y3, double n) {
            if (n < -1 || n > 1) throw BadFactor;
            double a = y2 - y1;
            double b = y3 - y2;
            double c = b - a;
            return y2 + n / 2 * (a + b + n * c);
        }

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <param name="y4">Input 4 to interpolation.</param>
        /// <param name="y5">Input 5 to interpolation.</param>
        /// <param name="n">Interpolation factor.</param>
        /// <returns>Interpolated value.</returns>
        public static double Interpolate(double y1, double y2, double y3, double y4, double y5, double n) {
            if (n < -1 || n > 1) throw BadFactor;
            double a = y2 - y1;
            double b = y3 - y2;
            double c = y4 - y3;
            double d = y5 - y4;
            double e = b - a;
            double f = c - b;
            double g = d - c;
            double h = f - e;
            double j = g - f;
            double k = j - h;
            double n2 = n * n;
            return y3 + n / 2 * (b + c) + n2 / 2 * f + n * (n2 - 1) / 12 * (h + j) + n2 * (n2 - 1) / 24 * k;
        }

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <returns>Maximum or minimum value and corresponding interpolation factor.</returns>
        public static (double extreme, double n) Extremum(double y1, double y2, double y3) {
            double a = y2 - y1;
            double b = y3 - y2;
            double c = b - a;
            double extreme = y2 - (a + b) * (a + b) / (8 * c);
            double n = -((a + b) / (2 * c));
            return (extreme, n);
        }

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <param name="y4">Input 4 to interpolation.</param>
        /// <param name="y5">Input 5 to interpolation.</param>
        /// <returns>Maximum or minimum value and corresponding interpolation factor.</returns>
        public static (double extreme, double n) Extremum(double y1, double y2, double y3, double y4, double y5) {
            double a = y2 - y1;
            double b = y3 - y2;
            double c = y4 - y3;
            double d = y5 - y4;
            double e = b - a;
            double f = c - b;
            double g = d - c;
            double h = f - e;
            double j = g - f;
            double k = j - h;
            double n = 0, lastN;
            do {
                lastN = n;
                n = (6 * b + 6 * c - h - j + 3 * n * n * (h + j) + 2 * n * n * n * k) / (k - 12 * f);
            }
            while (Math.Abs(n - lastN) > Threshold);
            double extreme = Interpolate(y1, y2, y3, y4, y5, n);
            return (extreme, n);
        }

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <returns>Interpolation factor when value reaches zero.</returns>
        public static double Zero(double y1, double y2, double y3) {
            double a = y2 - y1;
            double b = y3 - y2;
            double c = b - a;
            double n = 0, lastN;
            do {
                lastN = n;
                n -= (2 * y2 + n * (a + b + c * n)) / (a + b + 2 * c * n);
            }
            while (Math.Abs(n - lastN) > Threshold);
            return n;
        }

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <param name="y4">Input 4 to interpolation.</param>
        /// <param name="y5">Input 5 to interpolation.</param>
        /// <returns>Interpolation factor when value reaches zero.</returns>
        public static double Zero(double y1, double y2, double y3, double y4, double y5) {
            double a = y2 - y1;
            double b = y3 - y2;
            double c = y4 - y3;
            double d = y5 - y4;
            double e = b - a;
            double f = c - b;
            double g = d - c;
            double h = f - e;
            double j = g - f;
            double k = j - h;
            double n = 0, lastN;
            do {
                lastN = n;
                n -= Calculate(n);
            }
            while (Math.Abs(n - lastN) > Threshold);
            return n;

            double Calculate(double x) {
                double M = k / 24;
                double N = (h + j) / 12;
                double P = f / 2 - M;
                double Q = (b + c) / 2 - N;
                double x2 = x * x;
                double x3 = x2 * x;
                double x4 = x3 * x;
                return (M * x4 + N * x3 + P * x2 + Q * x + y3) / (4 * M * x3 + 3 * N * x2 + 2 * P * x + Q);
            }
        }

        /// <param name="y1">Input 1 to interpolation.</param>
        /// <param name="y2">Input 2 to interpolation.</param>
        /// <param name="y3">Input 3 to interpolation.</param>
        /// <param name="y4">Input 4 to interpolation.</param>
        /// <returns>Value halfway between y2 and y3.</returns>
        public static double InterpolateHalf(double y1, double y2, double y3, double y4) {
            double da = y2 - y1;
            double db = y3 - y2;
            double dc = y4 - y3;
            if (da != db || da != dc || db != dc) throw new ArgumentException("Y values must be equally spaced.");
            return (9 * (y2 + y3) - y1 - y4) / 16;
        }

        /// <param name="point">Array of known coordinates.</param>
        /// <param name="x">Input value for interpolation.</param>
        /// <returns>Interpolated value.</returns>
        /// <exception cref="ArgumentException"></exception>
        public static double Interpolate((double x, double y)[] point, double x) {
            double n = point.Length;
            if (n < 2) throw new ArgumentException("Array must have at least two elements.");
            if (point.GroupBy(p => p.x).Any(g => g.Count() > 1)) throw new ArgumentException("Array may not contain duplicated X values.");
            double y = 0;
            for (int i = 0; i < n; i++) {
                double p = 1;
                for (int j = 0; j < n; j++) {
                    if (j != i) p *= (x - point[j].x) / (point[i].x - point[j].x);
                }
                y += point[i].y * p;
            }
            return y;
        }        
    }
}
