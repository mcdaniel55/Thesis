using System;
using System.Collections.Generic;
using System.Text;

namespace Thesis
{
    static class Interpolation
    {
        /// <summary>
        /// Produces the result of a linear interpolation between the points (a,b) and (c,d) at location x in [a,c]
        /// </summary>
        /// <param name="a">The first coordinate of the first point</param>
        /// <param name="b">The second coordinate of the first point</param>
        /// <param name="c">The first coordinate of the second point</param>
        /// <param name="d">The second coordinate of the second point</param>
        /// <param name="x">The location in [a,c] at which to evaluate the interpolant</param>
        /// <returns></returns>
        public static double Lerp(double a, double b, double c, double d, double x)
        {
            return ((x - a) * d + (c - x) * b) / (c - a);
        }

        public static double[] Linspace(double a, double b, int size)
        {
            double[] values = new double[size];
            double stepSize = (b - a) / (size - 1);
            values[0] = a;
            for (int i = 1; i < values.Length; i++)
            {
                values[i] = a + i * stepSize;
            }
            return values;   
        }
    }
}
