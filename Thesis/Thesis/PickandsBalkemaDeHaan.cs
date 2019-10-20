using System;
using System.Collections.Generic;

namespace Thesis
{
    static class PickandsBalkemaDeHaan
    {
        // This trichotomy is how it's written in the paper, but these can all be computed by the first case
        private static double ExcessDistributionFunctionDeprecated(double x, double a, double c)
        {
            // When c is in (-epsilon, epsilon), we assume c = 0.
            const double epsilon = 1E-5;

            if (c > epsilon)
            {
                return Math.Pow(1 + c * x / a, -1 / c);
            }
            else if (c < -epsilon)
            {
                if (x < a / c) return 0;
                return Math.Pow(1 - Math.Abs(c) * x / a, 1 / Math.Abs(c));
            }
            // c = 0
            return Math.Exp(-x / a);
        }

        /// <summary> Computes the excess distribution estimate 1 - G(x) for given values of a and c </summary>
        /// <param name="a"></param>
        /// <param name="c"></param>
        public static double ExcessDistributionFunction(double x, double a, double c)
        {
            if (c < 0 && x > -a / c) return 0;
            if (c == 0) return Math.Exp(-x / a);
            return Math.Pow(1 + c * x / a, -1 / c);
        }

        /// <summary> Computes the tail CDF G(x) </summary>
        public static double TailCDF(double x, double a, double c) => 1 - ExcessDistributionFunction(x, a, c);

        /// <summary> Computes the inverse of the conditional excess distribution estimate 1-G(x) </summary>
        public static double ExcessDistributionQuantileFunction(double q, double a, double c)
        {
            return a / c * (Math.Pow(q, -c) - 1);
        }

        /// <summary> Computes the inverse of the tail CDF </summary>
        public static double TailQuantileFunction(double q, double a, double c)
        {
            return a / c * (Math.Pow(1 - q, -c) - 1);
        }

        /// <summary> Estimates the values of c and a based on a choice of M </summary>
        /// <param name="data"> An indexed list or array of observations of the RV, in increasing order </param>
        /// <param name="m"> How many data points back in the list to compute the estimate with. Must be <= count/4 </param>
        public static void EstimateParams(IList<double> data, int m, out double c, out double a)
        {
            const double ln2inv = 1.44269504088896; // 1 / ln(2)
            int n = data.Count; // Will be inlined by the compiler
            c = ln2inv * Math.Log((data[n - m] - data[n - 2 * m]) / (data[n - 2 * m] - data[n - 4 * m]));
            a = (data[n - 2 * m] - data[n - 4 * m]) * c / (Math.Pow(2, c) - 1);
        }

        /// <summary> Computes the sup norm distance d(l) between the estimate G and the empirical distribution F for a given value of l </summary>
        /// <param name="data"> An indexed list or array of observations of the RV, in increasing order </param>
        /// <param name="l"> The value of M to be tested </param>
        /// <param name="a"> The value of a to be used. Pickand's estimate will be used if this is not specified </param>
        /// <param name="c"> The value of c to be used. Pickand's estimate will be used if a is not specified </param>
        /// <returns> The largest deviation between the empirical tail distribution function F and the estimated complement excess distribution function G at any sample point </returns>
        /// <remarks> Pickands' definition of F is inconsistent (as far as I can tell) with the characterization he gives in the proceeding remarks, 
        /// so I used the definition rather than his characterization. </remarks>
        private static double DeviationSupNorm(IList<double> data, int l, double a = -1, double c = 0)
        {
            double largestDeviation = 0;
            if (a < 0) { EstimateParams(data, l, out c, out a); }
            for (int i = 1; i < 4 * l; i++)
            {
                // The largest deviation should occur at one of the step points
                double GHat = TailCDF(data[data.Count - i] - data[data.Count - 4 * l], a, c);
                largestDeviation = Math.Max(largestDeviation, Math.Abs((4.0 * l - i) / (4.0 * l - 1) - GHat));
                //largestDeviation = Math.Max(largestDeviation, Math.Abs((4.0 * l - i - 1) / (4.0 * l - 1) - GHat)); // If you want to include the lower values as well
            }
            return largestDeviation;
        }

        private static double DeviationWeightedSupNorm(IList<double> data, int l, double a = -1, double c = 0)
        {
            double largestDeviation = 0;
            if (a < 0) { EstimateParams(data, l, out c, out a); }
            for (int i = 1; i < 4 * l; i++)
            {
                // The largest deviation should occur at one of the step points
                double GHat = TailCDF(data[data.Count - i] - data[data.Count - 4 * l], a, c);
                double deviation = Math.Abs((4.0 * l - i) / (4.0 * l - 1) - GHat) * Weight(i);
                largestDeviation = Math.Max(largestDeviation, deviation);
            }
            return largestDeviation;

            double Weight(int i) => (data.Count - i + 1) * 1.0 / data.Count;
        }

        private static double DeviationWeightedLeastSquares(IList<double> data, int l, double a = -1, double c = 0)
        {
            double sum = 0;
            double weightSum = 0;
            if (a < 0) { EstimateParams(data, l, out c, out a); }
            for (int i = 1; i < 4 * l; i++)
            {
                // The largest deviation should occur at one of the step points
                double GHat = TailCDF(data[data.Count - i] - data[data.Count - 4 * l], a, c);
                double deviation = (4.0 * l - i) / (4.0 * l - 1) - GHat;
                double weight = Weight(i);
                sum += deviation * deviation * weight;
                weightSum += weight;
                // If you want to include the lower values as well
                /*
                deviation = (4.0 * l - i - 1) / (4.0 * l - 1) - GHat;
                sum += deviation * deviation;
                count += 2;
                */
            }
            return sum / weightSum;

            // Linear weighting that gives full weight to the upper end of the data and half weighting to the middle
            // Basically each point is weighted by what percentile of the data it rests at
            double Weight(int i) => (data.Count - i + 1) * 1.0 / data.Count;
        }

        /// <summary> Pickands' approach to determining a value for M. </summary>
        /// <remarks> Has (n^2 - n) / 32 calls to ComplementG() </remarks>
        public static void ApproximateExcessDistributionParametersPickands(IList<double> data, out double a, out double c, out int m)
        {
            int bestFitM = 1;
            double smallestDeviation = double.PositiveInfinity;
            for (int i = 1; i <= data.Count / 4; i++)
            {
                double deviation = DeviationSupNorm(data, i);
                //Console.WriteLine($"Deviation: {deviation}"); // Testing
                if (deviation < smallestDeviation)
                {
                    smallestDeviation = deviation;
                    bestFitM = i;
                }
            }
            m = bestFitM;
            EstimateParams(data, bestFitM, out c, out a);
        }

        /// <summary> An optimization for finding the best m, a, and c values, which is not designed for performance </summary>
        private static void ApproximateExcessDistributionParametersSlow(IList<double> data, out double a, out double c, out int m)
        {
            // The highest level of optimization here is deciding on M
            int bestFitM = 1;
            double aBest = 0, cBest = 0;

            double bestDev = double.PositiveInfinity;
            for (int i = 1; i < data.Count / 4; i++)
            {
                double al, cl;
                EstimateParams(data, i, out cl, out al);
                double deviation = OptimizeAC(i, al, cl, out al, out cl);
                if (deviation < bestDev)
                {
                    bestDev = deviation;
                    aBest = al;
                    cBest = cl;
                    bestFitM = i;
                }
            }
            a = aBest;
            c = cBest;
            m = bestFitM;

            double OptimizeAC(int l, double aGuess, double cGuess, out double aOpt, out double cOpt, int layers = 20)
            {
                // Choose a metric here
                Func<IList<double>, int, double, double, double> DeviationMetric = DeviationWeightedLeastSquares;
                // Using Pickand's estimates as initial guesses here
                double aCurrent = aGuess;
                double cCurrent = cGuess;
                double smallestDev = DeviationMetric(data, l, aGuess, cGuess);
                double deviation, step;

                bool improving = true;
                while (improving)
                {
                    improving = false;

                    // Optimize a for the current value of c
                    // Since a is a positive constant, and it is likely to have one local optimum, we can use a doubling/halving method to find a good value
                    int layer = 1;
                    while (layer <= layers)
                    {
                        step = aCurrent / Math.Pow(2, layer);
                        if ((deviation = DeviationMetric(data, l, aCurrent - step, cCurrent)) < smallestDev)
                        {
                            smallestDev = deviation;
                            aCurrent -= step;
                            improving = true;
                        }
                        else if ((deviation = DeviationMetric(data, l, aCurrent + step, cCurrent)) < smallestDev)
                        {
                            smallestDev = deviation;
                            aCurrent += step;
                            improving = true;
                        }
                        else
                        {
                            layer++;
                        }
                    }

                    // Optimize c for the current value of a
                    // Give negated and near-zero variations of c a try
                    if ((deviation = DeviationMetric(data, l, aCurrent, -cCurrent)) < smallestDev)
                    {
                        smallestDev = deviation;
                        cCurrent *= -1;
                        improving = true;
                    }
                    if ((deviation = DeviationMetric(data, l, aCurrent, Math.Pow(2, -30))) < smallestDev)
                    {
                        smallestDev = deviation;
                        cCurrent = Math.Pow(2, -30);
                        improving = true;
                    }
                    // Find a local minimum
                    layer = 1;
                    while (layer <= layers)
                    {
                        step = cCurrent / Math.Pow(2, layer);
                        if ((deviation = DeviationMetric(data, l, aCurrent, cCurrent - step)) < smallestDev)
                        {
                            smallestDev = deviation;
                            cCurrent -= step;
                            improving = true;
                        }
                        else if ((deviation = DeviationMetric(data, l, aCurrent, cCurrent + step)) < smallestDev)
                        {
                            smallestDev = deviation;
                            cCurrent += step;
                            improving = true;
                        }
                        else
                        {
                            layer++;
                        }
                    }
                }
                aOpt = aCurrent;
                cOpt = cCurrent;
                //Console.WriteLine($"OptAC(l={l}, a={aOpt}, c={cOpt}: {smallestDev}"); // Testing
                return smallestDev;
            }
        }
    }


    /// <summary>
    /// An approximation of the CDF and QF of the distribution of a dataset, produced by following Pickand's approach in his 1975 paper.
    /// </summary>
    public class PickandsApproximation
    {
        public double a, c; // Parameters, with c corresponding to the gamma or xi parameter of the associated GEV distribution
        public int transitionIndex;
        public double transitionAbscissa;
        List<double> data;

        public PickandsApproximation(IList<double> data)
        {
            if (data.Count < 30) throw new ArgumentException("Insufficient data count for Pickands Balkema De Haan.");
            this.data = new List<double>(data);
            this.data.Sort();
            PickandsBalkemaDeHaan.ApproximateExcessDistributionParametersPickands(this.data, out a, out c, out transitionIndex); // Write m to transitionIndex
            transitionIndex = this.data.Count - 4 * transitionIndex; // Convert from m to the actual transitionIndex
            transitionAbscissa = this.data[transitionIndex]; // m is guaranteed to be > 0
        }
        
        public double CDF(double x)
        {
            if (x <= transitionAbscissa) // ECDF Case
            {
                // Fast search for which elements to interpolate between
                int idx = data.BinarySearch(x);
                if (idx < 0) idx = ~idx;
                else idx++; // Left-continuity here
                return idx * 1.0 / data.Count;
            }
            // Pickands' tail case
            x -= transitionAbscissa;
            double offset = (transitionIndex + 1.0) / data.Count;
            double scale = 1 - offset;
            return scale * PickandsBalkemaDeHaan.TailCDF(x, a, c) + offset;
        }

        public double Quantile(double q)
        {
            double offset = (transitionIndex + 1.0) / data.Count;
            if (q > offset) // Pickands approximation case
            {
                // Renormalize q for the tail QF
                double scale = 1 - offset;
                q = (q - offset) / scale;
                return transitionAbscissa + PickandsBalkemaDeHaan.TailQuantileFunction(q, a, c);
            }
            // ECDF case
            return data[(int)(q * data.Count)];
        }

        public double Sample(Random rand)
        {
            return Quantile(rand.NextDouble());
        }

        public void Samples(double[] array, Random rand)
        {
            for (int i = 0; i < array.Length; i++)
            {
                array[i] = Quantile(rand.NextDouble());
            }
        }
        
        /// <summary> Constructs an alternative ContinuousDistribution version of the approximation with a piecewise-linear ECDF and an upper tail generated using Pickands' algorithm.  </summary>
        /// <param name="data"> An indexed set (array, list, etc.) of observations from a random variable, sorted in increasing order. </param>
        /// <remarks> This is wonderful for testing, but relatively expensive computation and storage-wise. This also uses the right-continuous rather than left-continuous version of the ECDF, though it hardly matters.</remarks>
        public static ContinuousDistribution ApproximatePiecewiseDistributionWithUpperTail(IList<double> data, int resolution = 1000)
        {
            // --- Construct the linear ECDF ---
            // Copy the abscissas from the sample
            List<double> abscissas = new List<double>(data.Count + resolution);
            abscissas.AddRange(data);
            // Evaluate the ECDF
            List<double> cdfVals = new List<double>(data.Count + resolution);
            for (int i = 0; i < data.Count; i++)
            {
                cdfVals.Add(i * 1.0 / data.Count);
            }

            // --- Attach the tail ---
            // Estimate the tail parameters
            int m;
            double a, c;
            PickandsBalkemaDeHaan.ApproximateExcessDistributionParametersPickands(data, out a, out c, out m);

            // Remove the last 4m-1 CDF values so we can replace them with the tail
            abscissas.RemoveRange(abscissas.Count - 4 * m + 1, 4 * m - 1);
            cdfVals.RemoveRange(cdfVals.Count - 4 * m + 1, 4 * m - 1);
            // The last element of the CDF approximation is now Z_4M

            // Generate tail values evenly spaced over the quantiles of the conditional excess distribution function 1-G(x)
            var quantiles = new List<double>(Interpolation.Linspace(0, 1, resolution + 1));
            // Remove the first point, since we already have a point in the CDF at u and the conditional excess will always be 0 there
            quantiles.RemoveAt(0);
            // If the tail is unbounded, replace the quantile at 1 with 1 - 4 * epsilon
            //if (c >= 0) { quantiles.RemoveAt(quantiles.Count - 1); }
            if (c >= 0) { quantiles[quantiles.Count - 1] = 1 - Math.Pow(2, -50); }
            // Replace the proportions with their associated abscissas (eg, the actual quantiles)
            double Z4M = abscissas[abscissas.Count - 1]; // This is where the tail is to be attached
            for (int i = 0; i < quantiles.Count; i++)
            {
                quantiles[i] = PickandsBalkemaDeHaan.TailQuantileFunction(quantiles[i], a, c);
            }
            // Add the CDF values first, then translate and add the abscissas
            double offset = cdfVals[cdfVals.Count - 1]; // Vertical offset
            double scale = 1 - offset; // How much of the full unit probability is left for the tail
            for (int i = 0; i < quantiles.Count; i++)
            {
                cdfVals.Add(scale * PickandsBalkemaDeHaan.TailCDF(quantiles[i], a, c) + offset);
                quantiles[i] += Z4M;
            }
            abscissas.AddRange(quantiles);

            return new ContinuousDistribution(abscissas, cdfVals);
        }






    }
}
