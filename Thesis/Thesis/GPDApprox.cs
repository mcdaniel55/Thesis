using System;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace Thesis
{
    /// <summary>
    /// An approximation of the CDF and QF of the distribution of a dataset, produced by attaching a fitted GPD tail to a linear interpolation of the ECDF
    /// </summary>
    public class GPDApproximation
    {
        public double a, c; // Parameters, with c corresponding to the gamma or xi parameter of the associated GEV distribution
        public double transitionProportion;
        public double transitionAbscissa;
        readonly List<double> sortedData;
        readonly Random rand;
        const double SHAPE_EPSILON = 1E-6;

        public enum FittingMethod { Pickands_SupNorm, BFGS_MSE, Moments_MSE, V4 }

        private GPDApproximation(IList<double> sortedData, double a, double c, double transitionProportion, double transitionAbscissa, Random rand = null)
        {
            this.sortedData = (List<double>)sortedData;
            this.a = a;
            this.c = c;
            this.transitionProportion = transitionProportion;
            this.transitionAbscissa = transitionAbscissa;
            this.rand = rand ?? Program.rand;
        }

        public GPDApproximation(IList<double> data, FittingMethod method = FittingMethod.Pickands_SupNorm, Random rand = null)
        {
            if (data.Count < 30) throw new ArgumentException("Insufficient data count for Pickands Balkema De Haan theorem.");
            sortedData = new List<double>(data);
            sortedData.Sort();
            int m; // Transition index

            switch (method)
            {
                case FittingMethod.BFGS_MSE: // No longer used
                    ApproximateExcessDistributionParametersBFGS(sortedData, out a, out c, out transitionAbscissa);
                    // Compute the index of the closest element that is at or before u in the data
                    int transitionIndex = sortedData.BinarySearch(transitionAbscissa);
                    if (transitionIndex < 0) transitionIndex = ~transitionIndex;
                    transitionProportion = transitionIndex * 1.0 / sortedData.Count; // Shouldn't this be + 1 here?
                    break;

                case FittingMethod.Moments_MSE:
                    ApproximateExcessDistributionParametersMoments(sortedData, out a, out c, out m);
                    transitionProportion = (sortedData.Count - m + 1.0) / sortedData.Count;
                    transitionAbscissa = sortedData[sortedData.Count - m]; // Convert from m to the actual transitionIndex; m is guaranteed to be > 0
                    break;

                case FittingMethod.V4:
                    ApproximateExcessDistributionParametersV4(sortedData, out a, out c, out transitionAbscissa);
                    int transitionIdx = sortedData.BinarySearch(transitionAbscissa);
                    if (transitionIdx < 0) 
                    {
                        transitionIdx = ~transitionIdx; // Now the index of the next-largest element 
                        if (transitionIdx == 0) { transitionProportion = 0; break; }
                        transitionProportion = Interpolation.Lerp(
                            sortedData[transitionIdx - 1], transitionIdx * 1.0 / sortedData.Count, 
                            sortedData[transitionIdx], (transitionIdx + 1) * 1.0 / sortedData.Count, 
                            transitionAbscissa);
                        break;
                    }
                    transitionProportion = (transitionIdx + 1) * 1.0 / sortedData.Count;
                    break;

                default:
                    ApproximateExcessDistributionParametersPickands(sortedData, out a, out c, out m); // Write m to transitionIndex
                    transitionProportion = (sortedData.Count - 4 * m + 1.0) / sortedData.Count;
                    transitionAbscissa = sortedData[sortedData.Count - 4 * m]; // Convert from m to the actual transitionIndex; m is guaranteed to be > 0
                    break;

            }
            
            this.rand = rand ?? Program.rand;
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
            return scale * TailCDF(x, a, c) + transitionProportion;
        }

        public double Quantile(double q)
        {
            //double offset = (transitionIndex + 1.0) / data.Count;
            if (q > transitionProportion) // Pickands approximation case
            {
                // Renormalize q for the tail QF
                double scale = 1 - transitionProportion;
                q = (q - transitionProportion) / scale;
                return transitionAbscissa + TailQuantileFunction(q, a, c);
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

        #region Static Methods
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
            ApproximateExcessDistributionParametersPickands(data, out double a, out double c, out int m);

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
                quantiles[i] = TailQuantileFunction(quantiles[i], a, c);
            }
            // Add the CDF values first, then translate and add the abscissas
            double offset = cdfVals[cdfVals.Count - 1]; // Vertical offset
            double scale = 1 - offset; // How much of the full unit probability is left for the tail
            for (int i = 0; i < quantiles.Count; i++)
            {
                cdfVals.Add(scale * TailCDF(quantiles[i], a, c) + offset);
                quantiles[i] += Z4M;
            }
            abscissas.AddRange(quantiles);

            return new ContinuousDistribution(abscissas, cdfVals);
        }

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
            if (c < -SHAPE_EPSILON && x > -a / c) return 0;
            if (Math.Abs(c) < SHAPE_EPSILON) return Math.Exp(-x / a);
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
            if (c != 0) return a / c * (Math.Pow(1 - q, -c) - 1);
            return -a * Math.Log(q);
        }

        /// <summary> Estimates the values of c and a based on a choice of M </summary>
        /// <param name="data"> An indexed list or array of observations of the RV, in increasing order </param>
        /// <param name="m"> How many data points back in the list to compute the estimate with. Must be less than or equal to count / 4 </param>
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

        private static double DeviationSupNormSmoothed(IList<double> data, int l, double a = -1, double c = 0)
        {
            double largestDeviation = 0;
            if (a < 0) { EstimateParams(data, l, out c, out a); }
            double xOffset = data[data.Count - 4 * l];
            for (int i = 1; i < 4 * l; i++)
            {
                double x = 0.5 * (data[data.Count - i] + data[data.Count - i - 1]) - xOffset;
                double target = (8 * l - 2 * i - 1) * 0.5 / (4 * l - 1); //((2 * (data.Count - i + 1)) - 1) / (2.0 * data.Count);
                double GHat = TailCDF(x, a, c);
                double deviation = Math.Abs(GHat - target);
                largestDeviation = Math.Max(largestDeviation, deviation);
            }
            return largestDeviation;
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
                //double deviation = DeviationSupNorm(data, i);
                double deviation = DeviationSupNormSmoothed(data, i);
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

                /* ECDF version MSE
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
                */
                //return sum / (weightsum * knotCount);

                // Smoothed version MSE
                for (int i = 0; i < knotCount - 1; i++)
                {
                    double GHat = TailCDF(0.5 * (data[nextLargestIndex + i] + data[nextLargestIndex + i + 1]) - input[2], input[0], input[1]);
                    double residual = (2.0 * i + 3) / (2.0 * knotCount) - GHat;
                    sum += residual * residual;
                }
                //return sum / knotCount;
                return (1 + Math.Abs(input[1])) * sum / knotCount; // Weighted so that smaller magnitudes of c are preferred
            }

            // Get Pickands' estimates of a and c for m = n/16 + 1 as starting guesses, consistent with Z4M at the 75th percentile
            //double pickandsEstA, pickandsEstC;
            EstimateParams(data, data.Count / 16 + 1, out double pickandsEstC, out double pickandsEstA);
            double lowerBoundA = 0;
            double upperBoundA = 3 * pickandsEstA + 1;
            double lowerBoundC = Math.Min(3 * pickandsEstC, -3 * pickandsEstC) - 1;
            double upperBoundC = -lowerBoundC;
            // Initial guess for u is at data[(3/4)n]

            var optimum = FindMinimum.OfFunctionConstrained(Fitness,
                lowerBound: CreateVector.DenseOfArray(new double[] { lowerBoundA, lowerBoundC, data[0] }),
                upperBound: CreateVector.DenseOfArray(new double[] { upperBoundA, upperBoundC, data[data.Count - 3] }),
                initialGuess: CreateVector.DenseOfArray(new double[] { pickandsEstA, pickandsEstC /*Math.Min(pickandsEstC, 0)*/, data[data.Count * 3 / 4] }));

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
                    double deviation = Math.Log(data[data.Count - i - 1]) - Math.Log(data[data.Count - 1]);
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
                //return (1 + Math.Abs(gammaHat)) * sum / k; // Weighted version
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
                        bestEst += delta;
                        runningMin = fitMin;
                    }
                    else if (lowerEstFitness == fitMin)
                    {
                        bestEst -= delta;
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
            Moments(m, out double z1, out _);
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

        // Gets the observations above u, and appends u as the first element
        internal static IList<double> GetTailData(IList<double> sortedData, double uVal)
        {
            var list = new List<double>(sortedData);
            list.RemoveAll(x => x <= uVal);
            list.Insert(0, uVal);
            // Testing
            double avgdif = (sortedData[sortedData.Count - 1] - sortedData[sortedData.Count - 3]) / 2;
            list.Add(sortedData[sortedData.Count - 1] + avgdif); // Haven't noticed any differences from this yet
#if DEBUG
            if (list.Count < 2)
            {
                throw new Exception("Illegal state: list must have at least two elements.");
            }
#endif
            return list;
        }

        // Computes the exact first moment of the tail ECDF interpolant
        internal static double TailMean(IList<double> tailData)
        {
            double sum = 0.5 * (tailData[0] + tailData[tailData.Count - 1]); // Half weight to end points
            for (int i = 1; i < tailData.Count - 1; i++) { sum += tailData[i]; }
            return sum / (tailData.Count - 1);
        }

        // Uses the observations that are greater than u to compute the variance of the upper tail under a linear interpolation of the empirical distribution
        // Note: May be less numerically stable than the second moment. This seems not to work as well as E(X^2) - E(X)^2 depending on the location of the data
        internal static double TailVariance(IList<double> tailData, double xBar)
        {
            double Summand(double x1, double x2)
            {
                return (x2 * x2 * x2 / 3.0
                    - xBar * x2 * x2
                    + xBar * xBar * x2
                    - x1 * x1 * x1 / 3.0
                    + xBar * x1 * x1
                    - xBar * xBar * x1)
                    / (x2 - x1);
            }

            double sum = 0;
            for (int i = 0; i < tailData.Count; i++)
            {
                sum += Summand(i, i + 1);
            }

            return sum / tailData.Count; // Divide by n
        }

        // Computes the exact second moment of the tail ECDF interpolant
        internal static double TailSecondMoment(IList<double> tailData) // there must be at least two elements in the upper tail here, and u is naturally < x_n so this is fine
        {
            var x = tailData; // Aliases to shorten the expression
            var n = tailData.Count;
            double sum = x[0] * x[0] + x[n - 1] * x[n - 1] + x[0] * x[1]; // x_0^2 + x_n^2 + x0 * x1 in the labelling on the board
            for (int i = 1; i < n - 1; i++) { sum += x[i] * x[i + 1] + 2 * x[i] * x[i]; }
            return sum / (3 * n - 3);
        }

        // Estimate the shape and scale parameters given a choice of u, using method of moments estimation
        internal static void EstimateParamsMOM(IList<double> tailData, out double scaleParam, out double shapeParam)
        {
            double m1 = TailMean(tailData);
            double m2 = TailSecondMoment(tailData);
            double variance = m2 - m1 * m1;
            double xBarMinusMu = m1 - tailData[0]; // Transition point is tailData[0]; mu is the TP, not the mean here

            //shapeParam = Math.Min(0.5 * (1 - xBarMinusMu * xBarMinusMu / variance), 0); // Clamped to non-positive shapes
            shapeParam = 0.5 * (1 - xBarMinusMu * xBarMinusMu / variance);
            scaleParam = xBarMinusMu * (1 - shapeParam);

#if DEBUG
            if (scaleParam < 0)
            {
                throw new Exception("Scale cannot be negative.");
            }
            //Program.logger.WriteLine($"Variance: {TailVariance(tailData, m1)} m2-m1^2: {variance}");
#endif
        }

        internal static void ApproximateExcessDistributionParametersV4(IList<double> sortedData, out double a, out double c, out double u)
        {
            // The upper tail is defined here by an ECDF interpolating linearly from (u,0) to (x_i, i/n) for data x_1, x_2, ..., x_n all greater than u.
            // This is the model from which we compute the upper tail parameters, using method of moments.
            // This midpoint version works slightly better than the plain ECDF
            double MidpointMSE(IList<double> tailData, double scaleParam, double shapeParam)
            {
                int n = tailData.Count;

                double sum = 0;
                for (int i = 0; i < n - 1; i++)
                {
                    double GHat = TailCDF(0.5 * (tailData[i] + tailData[i + 1]) - tailData[0], scaleParam, shapeParam);
                    double residual = (2.0 * i + 1) / (2.0 * n) - GHat;
                    sum += residual * residual;
                }
                return sum / (n - 1);
            }

            double GetScore(double uval, out double scaleParam, out double shapeParam)
            {
                var tailData = GetTailData(sortedData, uval);
                EstimateParamsMOM(tailData, out double scaleEst, out double shapeEst);
                scaleParam = scaleEst;
                shapeParam = shapeEst;
                double score = MidpointMSE(tailData, scaleEst, shapeEst);
                return score;
            }

            // Try several choices of u evenly spaced over (x_0, x_n-3), and keep the best fit
            var uValues = Interpolation.Linspace(sortedData[0], sortedData[sortedData.Count - 5], sortedData.Count / 4);
            double bestU = 0;
            double bestA = 0;
            double bestC = 0;
            double bestScore = double.PositiveInfinity;
            for (int i = 0; i < uValues.Length; i++)
            {
                double score = GetScore(uValues[i], out double scaleEst, out double shapeEst);
                if (score < bestScore)
                {
                    bestScore = score;
                    bestU = uValues[i];
                    bestA = scaleEst;
                    bestC = shapeEst;
                }
            }
            // --- Refine the best so far by bisection search ---
            double delta = uValues[1] - uValues[0];
            for (int i = 0; i < 10; i++)
            {
                delta *= 0.5;
                double forwardU = Math.Min(bestU + delta, sortedData[sortedData.Count - 3]); // Don't go so high that we don't have data to work with
                double forwardScore = GetScore(forwardU, out double forwardScale, out double forwardShape);
                double backwardScore = GetScore(bestU - delta, out double backwardScale, out double backwardShape);
                if (forwardScore < bestScore)
                {
                    bestScore = forwardScore;
                    bestU = forwardU;
                    bestA = forwardScale;
                    bestC = forwardShape;
                }
                if (backwardScore < bestScore)
                {
                    bestScore = backwardScore;
                    bestU -= delta;
                    bestA = backwardScale;
                    bestC = backwardShape;
                }
            }

            u = bestU;
            a = bestA;
            c = bestC;
        }

        /// <summary> Computes a variety of GPD/ECDF approximates for the provided data, and assigns each one a weight </summary>
        /// <remarks> 
        /// The overall approach is the same as in V4, but this returns all of the generated approximations rather than just the best one. 
        /// The weights are computed according to P(model | data) ~ P(data | model), which by independence is the product of model PDF values over the datapoints.
        /// </remarks>
        internal static void GetTailApproximatesAndWeights(IList<double> sortedData, out List<GPDApproximation> approximations, out List<double> weights)
        {
            approximations = new List<GPDApproximation>(sortedData.Count);
            weights = new List<double>(sortedData.Count);

            // The upper tail is defined here by an ECDF interpolating linearly from (u,0) to (x_i, i/n) for data x_1, x_2, ..., x_n all greater than u.
            // This is the model from which we compute the upper tail parameters, using method of moments.
            // This midpoint version works slightly better than the plain ECDF
            double MidpointMSE(IList<double> tailData, double scaleParam, double shapeParam)
            {
                int n = tailData.Count;

                double sum = 0;
                for (int i = 0; i < n - 1; i++)
                {
                    double GHat = TailCDF(0.5 * (tailData[i] + tailData[i + 1]) - tailData[0], scaleParam, shapeParam);
                    double residual = (2.0 * i + 1) / (2.0 * n) - GHat;
                    sum += residual * residual;
                }
                return sum / (n - 1);
            }

            double GetScore(double uval, out double scaleParam, out double shapeParam)
            {
                var tailData = GetTailData(sortedData, uval);
                EstimateParamsMOM(tailData, out double scaleEst, out double shapeEst);
                scaleParam = scaleEst;
                shapeParam = shapeEst;
                double score = MidpointMSE(tailData, scaleEst, shapeEst);
                return score;
            }

            // Try several choices of u evenly spaced over (x_0, x_n-3), and keep the best fit
            var uValues = Interpolation.Linspace(sortedData[0], sortedData[sortedData.Count - 5], sortedData.Count / 4);
            double bestU = 0;
            double bestA = 0;
            double bestC = 0;
            double bestScore = double.PositiveInfinity;
            for (int i = 0; i < uValues.Length; i++)
            {
                double score = GetScore(uValues[i], out double scaleEst, out double shapeEst);
                if (score < bestScore)
                {
                    bestScore = score;
                    bestU = uValues[i];
                    bestA = scaleEst;
                    bestC = shapeEst;
                }
            }
            // --- Refine the best so far by bisection search ---
            double delta = uValues[1] - uValues[0];
            for (int i = 0; i < 10; i++)
            {
                delta *= 0.5;
                double forwardU = Math.Min(bestU + delta, sortedData[sortedData.Count - 3]); // Don't go so high that we don't have data to work with
                double forwardScore = GetScore(forwardU, out double forwardScale, out double forwardShape);
                double backwardScore = GetScore(bestU - delta, out double backwardScale, out double backwardShape);
                if (forwardScore < bestScore)
                {
                    bestScore = forwardScore;
                    bestU = forwardU;
                    bestA = forwardScale;
                    bestC = forwardShape;
                }
                if (backwardScore < bestScore)
                {
                    bestScore = backwardScore;
                    bestU -= delta;
                    bestA = backwardScale;
                    bestC = backwardShape;
                }
            }

            u = bestU;
            a = bestA;
            c = bestC;




        }

        /// <summary> An optimization for finding the best m, a, and c values, which does not attempt to minimize computation time </summary>
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
        #endregion
    }
}
