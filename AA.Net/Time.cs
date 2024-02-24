using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;
using static System.Math;

namespace AA.Net {
    /// <summary>
    /// Methods relating to conversion and manipulation of time quantities.
    /// </summary>
    public static class Time {

        public const double J2000 = 2451545;
        public const double B1950 = 2433282.4235;

        public static List<(DateTime date, double increment)> LeapSecondDates;

        private static string LeapTableErrMessage = "";
        static Time() {
            try {
                XmlSerializer xml = new XmlSerializer(typeof(List<(DateTime, double)>));
                using (FileStream stream = new FileStream("leap_second.xml", FileMode.Open)) {
                    LeapSecondDates = (List<(DateTime, double)>)xml.Deserialize(stream);
                    stream.Close();
                }
            }
            catch(Exception e) {
                LeapTableErrMessage = e.Message;
            }
        }

        /// <param name="time">A UTC date and time.</param>
        /// <param name="t">A positive or negative time span to add.</param>
        /// <returns>A UTC date and time accounting for any leap seconds within the interval.</returns>
        /// <remarks>Leap second data may be edited in the external file leap_second.xml.</remarks>
        public static DateTime AddWithLeapSeconds(this DateTime time, TimeSpan t) {
            if (LeapSecondDates is null) throw new InvalidOperationException($"Leap second table not loaded: {LeapTableErrMessage}");
            DateTime sum = time + t;
            DateTime first = time; 
            DateTime last = sum;
            int sign = -1;
            if (first > last) {
                (last, first) = (first, last);
                sign = 1;
            }
            double adjust = LeapSecondDates.Where(lsd => lsd.date > first && lsd.date <= last).Sum(lsd => lsd.increment);
            return sum.AddSeconds(sign * adjust);
        }

        /// <param name="time">Date and time in the Gregorian calendar.</param>
        /// <returns>Julian Day.</returns>
        public static double JulianDay(this DateTime time) {
            int Y = time.Year;
            int M = time.Month;
            double D = time.Day + time.TimeOfDay.Ticks / (double)TimeSpan.TicksPerDay;
            return JulianDay(Y, M, D);
        }

        private static double DaysInGregorianYear(int year) {
            return IsLeapYear(year) ? 366 : 365;
        }

        /// <param name="year">A fractional year in the Gregorian calendar.</param>
        /// <returns>Julian Day.</returns>
        public static double JulianDay(this double year) {
            int wholeYear = (int)year;
            DateTime t = new DateTime(wholeYear, 1, 1).AddDays(year % 1 * DaysInGregorianYear(wholeYear));
            return t.JulianDay();
        }

        /// <param name="epochYear">A fractional year.</param>
        /// <returns>Julian Day as offset from J2000.0.</returns>
        public static double JulianEpochDay(this double epochYear) {
            return (epochYear - 2000) * 365.25 + J2000;
        }

        /// <param name="year">Year number.</param>
        /// <param name="month">Month of the year 1-12.</param>
        /// <param name="day">Day of the month including fractional clock time.</param>
        /// <returns>Julian Day.</returns>
        /// <remarks>pp. 60-61. This method assumes that dates before 15 October 1582 are in the Julian calendar. Julian Day begins at noon of the same UT day, i.e. 12 hours later.</remarks>
        public static double JulianDay(int year, int month, double day) {
            if (month < 3) {
                year -= 1;
                month += 12;
            }
            double A = Floor(year / 100d);
            double B = year < 1582 || (year == 1582 && (month < 10 || (month == 10 && day < 15))) ? 0 : 2 - A + Floor(A / 4);
            return Floor(365.25 * (year + 4716)) + Floor(30.6001 * (month + 1)) + day + B - 1524.5;
        }

        /// <param name="julianEphemerisDay">Julian Day in dynamical time.</param>
        /// <returns>Julian Day adjusted from Terrestrial Time to Universal Time.</returns>
        public static double JulianUniversalDay(this double julianEphemerisDay) {
            var ymd = julianEphemerisDay.ToDate();
            return julianEphemerisDay - DeltaT(ymd.year, ymd.month) / 86400;
        }

        /// <param name="year">Gregorian year number.</param>
        /// <returns>January 0.0 (JD0).</returns>
        /// <remarks>p. 62</remarks>
        public static double JulianDayZero(int year) {
            int Y = year - 1;
            int A = Y / 100;
            return Math.Floor(365.25 * Y) - A + A / 4 + 1721424.5;
        }

        /// <param name="time">Date and time in the Gregorian calendar.</param>
        /// <returns>Modified Julian Day (MJD).</returns>
        public static double ModifiedJulianDay(this DateTime time) {
            return time.JulianDay() - 2400000.5;
        }

        /// <param name="year">Year number.</param>
        /// <param name="month">Month of the year 1-12.</param>
        /// <param name="day">Day of the month, including fractional clock time.</param>
        /// <returns>Modified Julian Day (MJD).</returns>
        /// <remarks>p. 63</remarks>
        public static double ModifiedJulianDay(int year, int month, double day) {
            return JulianDay(year, month, day) - 2400000.5;
        }

        /// <param name="time">Date and time in the Gregorian calendar.</param>
        /// <returns>Julian Ephemeris Day (JDE).</returns>
        public static double JulianEphemerisDay(this DateTime time) {
            return time.AddSeconds(DeltaT(time.Year, time.Month)).JulianDay(); 
        }

        /// <param name="year">Year number.</param>
        /// <param name="month">Month of the year 1-12.</param>
        /// <param name="day">Day of the month.</param>
        /// <returns>Julian Ephemeris Day (JDE).</returns>
        public static double JulianEphemerisDay(int year, int month, double day) {
            return JulianDay(year, month, day + DeltaT(year, month) / 86400);
        }
        
        /// <param name="julianDay">A Julian Day.</param>
        /// <returns>Julian Ephemeris Day (JDE).</returns>
        /// <remarks>Chapter 10, generally.</remarks>
        public static double JulianEphemerisDay(this double julianDay) {
            var ymd = julianDay.ToDate();
            return julianDay + DeltaT(ymd.year, ymd.month) / 86400;
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <returns>Day of the week where 0-7 represent Sunday-Saturday.</returns>
        /// <remarks>p. 65</remarks>
        public static int DayOfWeek(this double julianDay) {
            return (int)((julianDay + 1.5) % 7);
        }

        private static bool IsLeapYear(int year) {
            return year % 4 == 0 && (year <= 1582 || year % 100 != 0 || year % 400 == 0);
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <returns>Sequential day of the year.</returns>
        /// <remarks>p. 65</remarks>
        public static double DayOfYear(this double julianDay) {
            (int year, int month, double day) = julianDay.ToDate();
            int k = IsLeapYear(year) ? 1 : 2;
            return Floor(275.0 * month / 9) - k * Floor((month + 9.0) / 12) + day - 30;
        }

        /// <remarks>Espenak-Meeus method from https://www.eclipsewise.com/help/deltatpoly2014.html.</remarks>
        /// <param name="year">Year of the (proleptic) Gregorian calendar.</param>
        /// <param name="month">Month of the year 1-12.</param>
        /// <returns>The difference between Terrestrial Time and Universal Time in seconds.</returns>
        public static double DeltaT(int year, int month) {
            double y = year + (month - 0.5) / 12;
            if (y < -1999 || y > 3000) throw new ArgumentOutOfRangeException("Time must be between -1999 and +3000.");
            double u, t;
            if (y < -500) {
                u = (year - 1820) / 100;
                return -20 + 32 * Pow(u, 2);
            }
            if (y < 500) {
                u = y / 100;
                return 10583.6 - 1014.41 * u + 33.78111 * Pow(u, 2) - 5.952053 * Pow(u, 3) - 0.1798452 * Pow(u, 4) + 0.022174192 * Pow(u, 5) + 0.0090316521 * Pow(u, 6);
            }
            if (y < 1600) {
                u = (y - 1000) / 100;
                return 1574.2 - 556.01 * u + 71.23472 * Pow(u, 2) + 0.319781 * Pow(u, 3) - 0.8503463 * Pow(u, 4) - 0.005050998 * Pow(u, 5) + 0.0083572073 * Pow(u, 6);
            }
            if (y < 1700) {
                t = y - 1600;
                return 120 - 0.9808 * t - 0.01532 * Pow(t, 2) + Pow(t, 3) / 7129;
            }
            if (y < 1800) {
                t = y - 1700;
                return 8.83 + 0.1603 * t - 0.0059285 * Pow(t, 2) + 0.00013336 * Pow(t, 3) - Pow(t, 4) / 1174000;
            }
            if (y < 1860) {
                t = y - 1800;
                return 13.72 - 0.332447 * t + 0.0068612 * Pow(t, 2) + 0.0041116 * Pow(t, 3) - 0.00037436 * Pow(t, 4) + 0.0000121272 * Pow(t, 5) - 0.0000001699 * Pow(t, 6) + 0.000000000875 * Pow(t, 7);
            }
            if (y < 1900) {
                t = y - 1860;
                return 7.62 + 0.5737 * t + 0.251754 * Pow(t, 2) + 0.01680668 * Pow(t, 3) - 0.0004473624 * Pow(t, 4) + Pow(t, 5) / 233174;
            }
            if (y < 1920) {
                t = y - 1900;
                return -2.79 + 1.494119 * t - 0.0598939 * Pow(t, 2) + 0.0061966 * Pow(t, 3) - 0.000197 * Pow(t, 4);
            }
            if (y < 1941) {
                t = y - 1920;
                return 21.20 + 0.84493 * t - 0.076100 * Pow(t, 2) + 0.0020936 * Pow(t, 3);
            }
            if (y < 1961) {
                t = y - 1950;
                return 29.07 + 0.407 * t - Pow(t, 2) / 233 + Pow(t, 3) / 2547;
            }
            if (y < 1986) {
                t = y - 1975;
                return 45.45 + 1.067 * t - Pow(t, 2) / 260 - Pow(t, 3) / 718;
            }
            if (y < 2005) {
                t = y - 2000;
                return 63.86 + 0.3345 * t - 0.060374 * Pow(t, 2) + 0.0017275 * Pow(t, 3) + 0.000651814 * Pow(t, 4) + 0.00002373599 * Pow(t, 5);
            }
            if (y < 2015) {
                t = y - 2005;
                return 64.69 + 0.293 * t;
            }
            t = y - 2015;
            return 67.62 + 0.3645 * t + 0.0039755 * Pow(t, 2);
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <returns>Centuries since Epoch 2000.</returns>
        /// <remarks>pp. 87, 143</remarks>
        public static double JulianCentury2000 (double julianDay) {
            return (julianDay - J2000) / 36525;
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <param name="epochYear">A year from which to measure.</param>
        /// <returns>Centuries since the specified epoch.</returns>
        /// <remarks>p. 410</remarks>
        public static double JulianCentury(double julianDay, double epochYear) {
            return (julianDay - epochYear.JulianDay()) / 36525;
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <returns>Year, month, and day (including fractional clock time) in the Julian or Gregorian calendar.</returns>
        /// <remarks>p. 63. This method returns a Julian calendar date for dates before 15 October 1582.</remarks>
        public static (int year, int month, double day)
        ToDate(this double julianDay) {
            if (julianDay < 0) throw new ArgumentOutOfRangeException("JulianDay must be nonnegative.");
            julianDay += 0.5;
            int Z = (int)julianDay;
            double F = julianDay - Z;
            double A;
            if (Z < 2299161) {
                A = Z;
            }
            else {
                double alpha = Floor((Z - 1867216.25) / 36524.25);
                A = Z + 1 + alpha - Floor(alpha / 4d);
            }
            double B = A + 1524;
            double C = Floor((B - 122.1) / 365.25);
            double D = Floor(365.25 * C);
            double E = Floor((B - D) / 30.6001);
            double day = B - D - Floor(30.6001 * E) + F;
            int month = (int)(E < 14 ? E - 1 : E - 13);
            int year = (int)(month > 2 ? C - 4716 : C - 4715);
            return (year, month, day);
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <param name="fromDynamical">Set True to convert from dynamical time.</param>
        /// <returns>Date and time in the Gregorian or Julian calendar.</returns>
        public static DateTime ToDateTime(this double julianDay, bool fromDynamical = false) {
            var ymd = julianDay.ToDate();
            if (ymd.year < 1 || ymd.year > 9999) throw new ArgumentOutOfRangeException("julianDay", "Evaluates to a year before 1 or after 9999.");
            int day = (int)ymd.day;
            double F = ymd.day - day;
            TimeSpan timeOfDay = TimeSpan.FromTicks((long)(F * TimeSpan.TicksPerDay));
            if (fromDynamical) timeOfDay -= TimeSpan.FromSeconds(DeltaT(ymd.year, ymd.month));
            return new DateTime(ymd.year, ymd.month, day, timeOfDay.Hours, timeOfDay.Minutes, timeOfDay.Seconds);
        }

        /// <param name="julianDay">A Julian Day.</param>
        /// <returns>Year and fraction of year.</returns>
        public static double ToFractionalYear(this double julianDay) {
            double n = julianDay.DayOfYear();
            int y = julianDay.ToDate().year;
            return y + n / (IsLeapYear(y) ? 366 : 365);
        }
    }
}