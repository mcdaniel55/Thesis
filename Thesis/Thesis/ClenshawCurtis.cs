using System;
using System.Collections.Generic;
using System.Text;

namespace Thesis.Quadrature
{
    static class ClenshawCurtis
    {
        public static double[] GetEvalPoints(int n)
        {
            double[] output = new double[n + 1];
            double c = Math.PI / n;

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = Math.Cos(i * c);
            }

            return output;
        }
        
        public static double[] GetOddEvalPoints(int n)
        {
            double[] output = new double[(n - 1)/2];
            double c = Math.PI / n;

            for (int i = 0; i < output.Length; i++)
            {
                output[i] = Math.Cos((2 * i + 1) * c);
            }

            return output;
        }

        /// <summary> Computes the weights for an nth-order Clenshaw-Curtis rule </summary>
        /// <param name="n"> The order of the rule </param>
        /// <returns> A double[] of weights </returns>
        /// <remarks> Follows formula (2.4) in Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules (Waldvogel 2005) 
        /// This has been sped up slightly by setting the ith and (n + 1 - i)th entries together, as the weights are symmetric about 0. </remarks>
        public static double[] GetWeights(int n)
        {
            double[] output = new double[n + 1];
            double c = 2 * Math.PI / n;

            int iterations = (output.Length + 1) / 2;

            for (int i = 0; i < iterations; i++)
            {
                output[i] = 0;
                for (int j = 1; j <= n / 2; j++)
                {
                    double term = j == n / 2 ? 1 : 2; // Apply b_j
                    term /= 4 * j * j - 1;
                    term *= Math.Cos(j * i * c); // Cos(2jv_k)
                    output[i] += term;
                }

                output[i] = (1 - output[i]) / n;
                if (i % n != 0) { output[i] *= 2; } // Apply c_k
                // Use symmetry to update the other side of the array (redundant for the middle x=0 element in odd-sized arrays)
                output[output.Length - i - 1] = output[i]; 
            }

            return output;
        }

        public static double Integrate(Func<double, double> f, double intervalStart, double intervalEnd, double[] evalPoints, double[] weights)
        {
            // Compute a linear transformation from the interval to [-1,1]
            double a = (intervalEnd - intervalStart) / 2.0;
            double b = intervalStart + a;
            double xOfz(double z) => a * z + b; // As z ranges from [-1,1], x ranges over [start,end]

            double sum = 0;
            for (int i = 0; i < weights.Length; i++)
            {
                sum += weights[i] * f(xOfz(evalPoints[i]));
            }

            // Multiply by the derivative of x wrt z and return
            return sum * a;
        }

        public static double Integrate(Func<double, double> f, double intervalStart, double intervalEnd, int order)
        {
            double[] evalPoints = GetEvalPoints(order);
            double[] weights = GetWeights(order);
            return Integrate(f, intervalStart, intervalEnd, evalPoints, weights);
        }
    }
}