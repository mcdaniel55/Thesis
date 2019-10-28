using System;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

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
        
        // Data must be sorted before using this method
        internal static void ApproximateExcessDistributionParametersBFGS(List<double> data, out double a, out double c, out double u)
        {
            double Fitness(Vector<double> input) // Input is assumed to be (a,c,u)
            { 
                double sum = 0;
                //double weightsum = 0;

                // Compute the index of the next largest element that is at or after u in the data
                int nextLargestIndex = data.BinarySearch(input[2]);
                if (nextLargestIndex < 0) nextLargestIndex = ~nextLargestIndex + 1;
                int knotCount = data.Count - nextLargestIndex;

                //for (int i = nextLargestIndex; i < data.Count; i++)
                for (int i = 0; i < knotCount; i++)
                {
                    // The largest deviation should occur at one of the step points
                    double GHat = TailCDF(data[nextLargestIndex + i] - input[2], input[0], input[1]); // Args: x - u, a, c
                    double residual = i * 1.0 / knotCount - GHat; // Deviation from the top of the step at x_i
                    //double weight = knotCount - i + 1;
                    sum += residual * residual;
                    //sum += Math.Abs(residual) * weight;
                    //weightsum += weight;
                }
                return sum / knotCount; // Consider dividing by n or n^2 here
                //return sum / (weightsum * knotCount);
            }

            // Get Pickands' estimates of a and c for m = n/16 + 1 as starting guesses, consistent with Z4M at the 75th percentile
            //double pickandsEstA, pickandsEstC;
            EstimateParams(data, data.Count / 16 + 1, out double pickandsEstC, out double pickandsEstA);
            double lowerBoundA = 0;
            double upperBoundA = 3 * pickandsEstA + 1;
            double lowerBoundC = Math.Min(3 * pickandsEstC, -3 * pickandsEstC) - 1;
            double upperBoundC = 0; //-lowerBoundC;
            // Initial guess for u is at data[(3/4)n]

            var optimum = FindMinimum.OfFunctionConstrained(Fitness,
                lowerBound: CreateVector.DenseOfArray(new double[] { lowerBoundA, lowerBoundC, data[0] }),
                upperBound: CreateVector.DenseOfArray(new double[] { upperBoundA, upperBoundC, data[data.Count - 3] }),
                initialGuess: CreateVector.DenseOfArray(new double[] { pickandsEstA, pickandsEstC, data[data.Count * 3 / 4] }));

            // Return parameters
            a = optimum[0];
            c = optimum[1];
            u = optimum[2];
        }

        internal static void ApproximateExcessDistributionParametersMoments(IList<double> data, out double a, out double c, out int m)
        {
            a = c = m = 0;
            void Moments(int k, out double moment1, out double moment2)
            {
                double sum1 = 0;
                double sum2 = 0;
                for (int i = 0; i < k; i++)
                {
                    double deviation = Math.Log(data[data.Count - i - 1]) - Math.Log(data[data.Count-1]);
                    sum1 += deviation;
                    sum2 += deviation * deviation;
                }
                moment1 = sum1 / k;
                moment2 = sum2 / k;
            }

            double GammaHat(int k)
            {
                Moments(k, out double m1, out double m2);
                return (m2 - 2 * m1 * m1) / (2 * (m2 - m1 * m1));
                //sigmaHat = data[data.Count - k - 1] * m1 * (1 - gammaHat);
            }
            
            double MSE(int k, double aHat, double gammaHat)
            {
                double sum = 0;
                for (int i = 0; i < k; i++) // Iterating backwards, from Z1 to Zk
                {
                    double x = data[data.Count - i - 1];
                    double ecdfAtX = 1 - i * 1.0 / k;
                    double deviation = TailCDF(x, aHat, gammaHat) - ecdfAtX;
                    sum += deviation * deviation;
                }
                return sum / k;
            }

            // Temp
            List<double> gammas = new List<double>();
            List<double> sigmas = new List<double>();
            // ---

            double bestObservedFitOverK = double.PositiveInfinity;
            for (int k = 4; k < data.Count; k++)
            {
                double gammaHat = GammaHat(k);

                // Find a good aHat scale parameter for this choice of k
                double middleEst = 4 * Math.Sqrt(Statistics.VarianceEstimate(data));
                double bestEst = middleEst;
                double runningMin = MSE(k, bestEst, gammaHat);
                for (int i = 0; i < 20; i++) // 10 iterations is 2-3 correct digits
                {
                    double delta = middleEst * Math.Pow(2, -i);
                    double higherEstFitness = MSE(k, bestEst + delta, gammaHat);
                    double lowerEstFitness = MSE(k, bestEst - delta, gammaHat);
                    double fitMin = Math.Min(runningMin, Math.Min(higherEstFitness, lowerEstFitness));
                    if (higherEstFitness == fitMin)
                    {
                        bestEst = bestEst + delta;
                        runningMin = fitMin;
                    }
                    else if (lowerEstFitness == fitMin)
                    {
                        bestEst = bestEst - delta;
                        runningMin = fitMin;
                    }
                }

                gammas.Add(gammaHat);
                sigmas.Add(bestEst);

                if (runningMin < bestObservedFitOverK)
                {
                    bestObservedFitOverK = runningMin;
                    // Write to a, c, and m
                    m = k;
                    a = bestEst;
                    c = gammaHat;
                }
            }

            // Temp
            double z1, z2;
            Moments(m, out z1, out z2);
            Program.logger.WriteLine($"Gamma hat: {c}");
            Program.logger.WriteLine($"Sigma hat: {data[data.Count - m - 1] * z1 * (1 - c)}");
            Program.logger.WriteLine($"Gammas");
            for (int i = 0; i < gammas.Count; i++)
            {
                Program.logger.WriteLine($"{gammas[i]}");
            }
            Program.logger.WriteLine($"Sigmas");
            for (int i = 0; i < gammas.Count; i++)
            {
                Program.logger.WriteLine($"{sigmas[i]}");
            }
        }

        /// <summary> An optimization for finding the best m, a, and c values, which is not designed for performance </summary>
        internal static void ApproximateExcessDistributionParametersSlow(IList<double> data, out double a, out double c, out int m)
        {
            // The highest level of optimization here is deciding on M
            int bestFitM = 1;
            double aBest = 0, cBest = 0;

            double bestDev = double.PositiveInfinity;
            for (int i = 1; i < data.Count / 4; i++)
            {
                //double al, cl;
                EstimateParams(data, i, out double cl, out double al);
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
        public double transitionProportion;
        public double transitionAbscissa;
        List<double> sortedData;
        Random rand;

        public enum FittingMethod { Pickands_SupNorm, BFGS_MSE, Moments_MSE }

        public PickandsApproximation(IList<double> data, FittingMethod method = FittingMethod.Pickands_SupNorm, Random rand = null)
        {
            if (data.Count < 30) throw new ArgumentException("Insufficient data count for Pickands Balkema De Haan.");
            sortedData = new List<double>(data);
            sortedData.Sort();
            int m; // Transition index

            switch (method)
            {
                case FittingMethod.BFGS_MSE:
                    PickandsBalkemaDeHaan.ApproximateExcessDistributionParametersBFGS(sortedData, out a, out c, out transitionAbscissa);
                    // Compute the index of the closest element that is at or before u in the data
                    int transitionIndex = sortedData.BinarySearch(transitionAbscissa);
                    if (transitionIndex < 0) transitionIndex = ~transitionIndex;
                    transitionProportion = transitionIndex * 1.0 / sortedData.Count;
                    break;

                case FittingMethod.Moments_MSE:
                    PickandsBalkemaDeHaan.ApproximateExcessDistributionParametersMoments(data, out a, out c, out m);
                    transitionProportion = (sortedData.Count - m + 1.0) / data.Count;
                    transitionAbscissa = sortedData[sortedData.Count - m]; // Convert from m to the actual transitionIndex; m is guaranteed to be > 0
                    break;

                default:
                    PickandsBalkemaDeHaan.ApproximateExcessDistributionParametersPickands(sortedData, out a, out c, out m); // Write m to transitionIndex
                    transitionProportion = (sortedData.Count - 4 * m + 1.0) / data.Count;
                    transitionAbscissa = sortedData[sortedData.Count - 4 * m]; // Convert from m to the actual transitionIndex; m is guaranteed to be > 0
                    break;
            }
            
            if (rand == null) rand = Program.rand;
            this.rand = rand;
        }
        
        public double CDF(double x)
        {
            if (x <= transitionAbscissa) // ECDF Case
            {
                // Fast search for which elements to interpolate between
                int idx = sortedData.BinarySearch(x);
                if (idx < 0) idx = ~idx;
                else idx++; // Left-continuity here
                return idx * 1.0 / sortedData.Count;
            }
            // Pickands' tail case
            x -= transitionAbscissa;
            //double offset = (transitionIndex + 1.0) / data.Count;
            double scale = 1 - transitionProportion;
            return scale * PickandsBalkemaDeHaan.TailCDF(x, a, c) + transitionProportion;
        }

        public double Quantile(double q)
        {
            //double offset = (transitionIndex + 1.0) / data.Count;
            if (q > transitionProportion) // Pickands approximation case
            {
                // Renormalize q for the tail QF
                double scale = 1 - transitionProportion;
                q = (q - transitionProportion) / scale;
                return transitionAbscissa + PickandsBalkemaDeHaan.TailQuantileFunction(q, a, c);
            }
            // ECDF case
            return sortedData[(int)(q * sortedData.Count)];
        }

        public double Sample()
        {
            return Quantile(rand.NextDouble());
        }

        public void Samples(double[] array)
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
            PickandsBalkemaDeHaan.ApproximateExcessDistributionParametersPickands(data, out double a, out double c, out int m);

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
