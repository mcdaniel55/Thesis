﻿using System;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;
using Thesis.Quadrature;

namespace ThesisQuadratureTests
{
    static class DiscardProbabilityComputations
    {
        #region Pairwise Comparisons
        /// <summary> Takes two normal distributions and returns the probability that a point sample from of the first will be greater than one from the second. That is, the probability of discarding A because of B. </summary>
        /// <param name="A">A normal distribution</param>
        /// <param name="B">Another normal distribution</param>
        /// <returns>P(A>B)</returns>
        public static double PairwiseExact(Normal A, Normal B)
        {
            Normal BMinusA = new Normal(mean: B.Mean - A.Mean,
                stddev: Math.Sqrt(A.Variance + B.Variance));
            return BMinusA.CumulativeDistribution(0);
        }

        /// <summary>
        /// Looks at every relevant pair-wise comparison between normals[idx] and the other elements of normals, and returns the greatest probability for discarding normals[idx]
        /// </summary>
        /// <param name="normals">The set of normal distributions</param>
        /// <param name="idx">The index of the distribution for which we are computing a discard probability</param>
        /// <returns>The highest probability found</returns>
        public static double BestPairwiseDiscardNaive(Normal[] normals, int idx)
        {
            double bestDiscardProbability = 0;
            for (int i = 0; i < normals.Length; i++)
            {
                if (i == idx) continue;
                //Console.WriteLine($"Discarding {idx} vs {i}: {PairwiseExact(normals[idx], normals[i])}"); // Temp
                bestDiscardProbability = Math.Max(bestDiscardProbability, PairwiseExact(normals[idx], normals[i]));
            }
            return bestDiscardProbability;
        }

        #region Deprecated Code
        /// <summary>
        /// Looks at a reduced set of pair-wise comparisons between normals[idx] and the other elements of normals, and returns the greatest probability for discarding normals[idx]
        /// </summary>
        /// <param name="normals">The set of normal distributions</param>
        /// <param name="idx">The index of the distribution for which we are computing a discard probability</param>
        /// <returns>The highest probability found</returns>
        public static double BestPairwiseDiscardEfficientV1(Normal[] normals, int idx)
        {
            double bestDiscardProbability = 0;

            // Find the dominant subset against which we will compare
            List<Normal> dominantSet = new List<Normal>(normals);
            //dominantSet.Sort((n1, n2) => n1.Mean.CompareTo(n2.Mean)); // Does sorting make it faster? No.

            // Standard generic list method                                                                      
            //dominantSet.RemoveAll(n => dominantSet.Exists(other => ((other.Mean < n.Mean && other.Variance <= n.Variance)
            //        || (other.Variance < n.Variance && other.Mean <= n.Mean))));

            // Old-school index-based method, moderate speed-up over the generic list version
            for (int j = 0; j < dominantSet.Count; j++)
            {
                for (int k = 0; k < dominantSet.Count; k++)
                {
                    if (j == k) continue;
                    // If j is dominated by k, remove j
                    if ((dominantSet[k].Mean < dominantSet[j].Mean && dominantSet[k].Variance <= dominantSet[j].Variance)
                        || (dominantSet[k].Variance < dominantSet[j].Variance && dominantSet[k].Mean <= dominantSet[j].Mean))
                    {
                        dominantSet.RemoveAt(j);
                        j--;
                        break;
                    }
                }
            }

            // Index-based iteration is faster than a foreach here
            for (int i = 0; i < dominantSet.Count; i++)
            {
                if (dominantSet[i] == normals[idx]) continue;
                //Console.WriteLine($"Discarding {idx} vs {norm.ToString()}: {PairwiseExact(normals[idx], norm)}"); // Temp
                bestDiscardProbability = Math.Max(bestDiscardProbability, PairwiseExact(normals[idx], dominantSet[i]));
            }

            return bestDiscardProbability;
        }


        /// <summary>
        /// Looks at a reduced set of pair-wise comparisons between normals[idx] and the other elements of normals, and returns the greatest probability for discarding normals[idx]
        /// </summary>
        /// <remarks> This performs significantly faster than the naive algorithm on test sets with 10-100 normal distributions, usually at least doubling in speed </remarks>
        /// <param name="normals">The set of normal distributions</param>
        /// <param name="idx">The index of the distribution for which we are computing a discard probability</param>
        /// <returns>The highest probability found</returns>
        public static double BestPairwiseDiscardEfficientV2(Normal[] normals, int idx)
        {
            double bestDiscardProbability = 0;

            // Find the dominant subset against which we will compare
            Normal[] dominantSet = (Normal[])normals.Clone();

            // Old-school index-based method, moderate speed-up over the generic list version
            for (int j = 0; j < dominantSet.Length; j++)
            {
                for (int k = 0; k < dominantSet.Length; k++)
                {
                    if (dominantSet[k] == null || j == k) continue;
                    // If j is dominated by k, remove j
                    if ((dominantSet[k].Mean < dominantSet[j].Mean && dominantSet[k].Variance <= dominantSet[j].Variance)
                        || (dominantSet[k].Variance < dominantSet[j].Variance && dominantSet[k].Mean <= dominantSet[j].Mean))
                    {
                        dominantSet[j] = null;
                        break;
                    }
                }
            }

            // Index-based iteration is faster than a foreach here
            for (int i = 0; i < dominantSet.Length; i++)
            {
                if (dominantSet[i] == null || dominantSet[i] == normals[idx]) continue;
                //Console.WriteLine($"Discarding {idx} vs {norm.ToString()}: {PairwiseExact(normals[idx], norm)}"); // Temp
                bestDiscardProbability = Math.Max(bestDiscardProbability, PairwiseExact(normals[idx], dominantSet[i]));
            }

            return bestDiscardProbability;
        }

        /// <summary>
        /// Looks at a reduced set of pair-wise comparisons between normals[idx] and the other elements of normals, and returns the greatest probability for discarding normals[idx]
        /// </summary>
        /// <remarks> This performs significantly faster than the naive algorithm on test sets with 10-100 normal distributions, usually at least doubling in speed.
        /// V4 is asymptotically slower than V3, but outperforms in on 10-70</remarks>
        /// <param name="normals">The set of normal distributions</param>
        /// <param name="idx">The index of the distribution for which we are computing a discard probability</param>
        /// <returns>The highest probability found</returns>
        public static double BestPairwiseDiscardEfficientV4(Normal[] normals, int idx)
        {
            double bestDiscardProbability = 0;

            // Old-school index-based method, moderate speed-up over the List.RemoveAll version
            for (int j = 0; j < normals.Length; j++)
            {
                if (j == idx) continue;
                bool dominated = false;
                for (int k = 0; k < normals.Length; k++)
                {
                    if (k == idx || k == j) continue;
                    // If j is dominated by k, ignore j
                    if ((normals[k].Mean < normals[j].Mean && normals[k].Variance <= normals[j].Variance)
                        || (normals[k].Variance < normals[j].Variance && normals[k].Mean <= normals[j].Mean))
                    {
                        dominated = true;
                        break;
                    }
                }
                if (!dominated)
                {
                    bestDiscardProbability = Math.Max(bestDiscardProbability, PairwiseExact(normals[idx], normals[j]));
                }
            }

            return bestDiscardProbability;
        }
        #endregion

        /// <summary>
        /// Looks at a reduced set of pair-wise comparisons between normals[idx] and the other elements of normals, and returns the greatest probability for discarding normals[idx]
        /// </summary>
        /// <remarks> This performs significantly faster than the naive algorithm on test sets with 10-100 normal distributions, usually at least doubling in speed </remarks>
        /// <param name="normals">The set of normal distributions</param>
        /// <param name="idx">The index of the distribution for which we are computing a discard probability</param>
        /// <returns>The highest probability found</returns>
        public static double BestPairwiseDiscard(Normal[] normals, int idx)
        {
            double bestDiscardProbability = 0;

            // Find the dominant subset against which we will compare
            int[] mask = new int[normals.Length];
            // 1 means the normal at that index is dominated, 0 otherwise

            // Ignore the element at idx for comparisons
            mask[idx] = 1;

            // Old-school index-based method, moderate speed-up over the List.RemoveAll version
            for (int j = 0; j < normals.Length; j++)
            {
                if (mask[j] != 0) continue;
                for (int k = 0; k < normals.Length; k++)
                {
                    if (mask[k] != 0 || k == j) continue;
                    // If j is dominated by k, remove j
                    if ((normals[k].Mean < normals[j].Mean && normals[k].Variance <= normals[j].Variance)
                        || (normals[k].Variance < normals[j].Variance && normals[k].Mean <= normals[j].Mean))
                    {
                        mask[j] = 1;
                        break;
                    }
                }
            }

            // Index-based iteration is faster than a foreach here
            for (int i = 0; i < normals.Length; i++)
            {
                if (mask[i] != 0) continue;
                //Console.WriteLine($"Discarding {idx} vs {norm.ToString()}: {PairwiseExact(normals[idx], norm)}"); // Temp
                bestDiscardProbability = Math.Max(bestDiscardProbability, PairwiseExact(normals[idx], normals[i]));
            }

            return bestDiscardProbability;
        }

        /// <summary>
        /// Looks at every pair-wise comparison between elements of normals and returns the greatest probability for discarding each element
        /// </summary>
        /// <param name="normals">The set of normal distributions</param>
        /// <returns>An array describing the highest discard probability found for each element</returns>
        public static double[] BestPairwiseDiscardsNaive(Normal[] normals)
        {
            double[] bestDiscardProbabilities = new double[normals.Length];
            for (int i = 0; i < normals.Length; i++)
            {
                for (int j = 0; j < normals.Length; j++)
                {
                    if (j == i) continue;
                    //Console.WriteLine($"Discarding {idx} vs {i}: {PairwiseExact(normals[idx], normals[i])}"); // Temp
                    bestDiscardProbabilities[i] = Math.Max(bestDiscardProbabilities[i], PairwiseExact(normals[i], normals[j]));
                }
            }
            return bestDiscardProbabilities;
        }

        /// <summary>
        /// Looks at the dominant subset of all pair-wise comparisons between elements of normals and returns the greatest probability for discarding each element
        /// </summary>
        /// <param name="normals">The set of normal distributions</param>
        /// <returns>An array describing the highest discard probability found for each element</returns>
        public static double[] BestPairwiseDiscards(Normal[] normals)
        {
            double[] bestDiscardProbabilities = new double[normals.Length];
            for (int i = 0; i < normals.Length; i++)
            {
                bestDiscardProbabilities[i] = BestPairwiseDiscard(normals, i);
            }
            return bestDiscardProbabilities;
        }

        /// <summary>
        /// Looks at every pair-wise comparison between elements of normals and returns the index of the element with the greatest probability of being discarded
        /// </summary>
        /// <remarks>This is an O(n) operation in the asymptotic case as n -> infty. The n elements are in a cloud, and the optimal and anti-optimal sets
        /// are each half the perimeter of the cloud. If the cloud ~ n, then the perimeter ~ 2sqrt(n), so there are sqrt(n) in each of the optimal and
        /// anti-optimal sets. Each of these will be compared against every element of the other set, so there will be ~ sqrt(n)*sqrt(n) = n operations.
        /// This is the preferred method for finding the next element to discard, as the anti-optimal set may gain elements when a previous element is discarded.</remarks>
        /// <param name="normals">The set of normal distributions</param>
        /// <param name="discardProbability">The probability with which this element can be discarded based on a pairwise comparison</param>
        /// <returns>An array describing the highest discard probability found for each element</returns>
        public static int BestPairwiseDiscard(Normal[] normals, out double discardProbability)
        {
            double bestDiscardProbability = 0;
            int bestIndex = 0;

            // Find the sets of optimal and anti-optimal elements of normals
            List<int> optimalIndices = new List<int>(normals.Length); // elements that are nondominated in having small means and stdevs
            List<int> antiOptimalIndices = new List<int>(normals.Length); // elements that do not dominate any other element by the same criteria

            for (int i = 0; i < normals.Length; i++)
            {
                bool dominates = false;
                bool isDominated = false;

                for (int j = 0; j < normals.Length; j++)
                {
                    if (j == i) continue;
                    if (normals[j].Mean < normals[i].Mean && normals[j].StdDev <= normals[i].StdDev
                        || normals[j].StdDev < normals[i].StdDev && normals[j].Mean <= normals[i].Mean) isDominated = true;
                    if (normals[i].Mean < normals[j].Mean && normals[i].StdDev <= normals[j].StdDev ||
                        normals[i].StdDev < normals[j].StdDev && normals[i].Mean <= normals[j].Mean) dominates = true;
                    if (dominates && isDominated) { break; }
                }
                if (!dominates) { antiOptimalIndices.Add(i); }
                if (!isDominated) { optimalIndices.Add(i); }
            }

            // Compare each of the optimal elements against each of the antioptimal ones
            for (int i = 0; i < antiOptimalIndices.Count; i++)
            {
                for (int j = 0; j < optimalIndices.Count; j++)
                {
                    if (antiOptimalIndices[i] == optimalIndices[j]) continue;
                    double newprob = PairwiseExact(normals[antiOptimalIndices[i]], normals[optimalIndices[j]]);
                    if (newprob > bestDiscardProbability)
                    {
                        bestDiscardProbability = newprob;
                        bestIndex = antiOptimalIndices[i];
                    }
                }
            }

            discardProbability = bestDiscardProbability;
            return bestIndex;
        }

        #endregion

        #region FullComparisons
        /// <summary> Takes a set of normal distributions and computes the probability that each one would produce the minimal value if 
        /// a point sample was taken from each. Integration is performed by Simpson's 3/8 rule. </summary>
        /// <param name="distributions">The array of distributions to compare</param>
        /// <param name="iterations">The number of iterations to use in Simpson's 3/8 Rule. Defaults to 150.</param>
        /// <returns>An array of probabilities that each distribution will produce the minimum value if a point sample was taken from each.
        /// [i] = P(X_i less than min(all X_j))</returns>
        public static double[] Simpsons38RuleWithoutNegation(Normal[] distributions, int iterations = 150)
        {
            double[] complementProbs = new double[distributions.Length];

            for (int i = 0; i < distributions.Length; i++)
            {
                Normal distribution_i = distributions[i];

                Func<double, double> integrand = x =>
                {
                    double product = distribution_i.Density(x);
                    for (int j = 0; j < distributions.Length; j++)
                    {
                        if (j != i) { product *= 1 - distributions[j].CumulativeDistribution(x); }
                    }
                    return product;
                };

                complementProbs[i] = NewtonCotes.Simpsons38Rule(integrand,
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, iterations);
                //Console.WriteLine($"S38R[{i}]: {complementProbs[i]}");
            }
            return complementProbs;
        }

        // The following are methods with alternative integrands derived from negating the distributions' means and finding the probability that each is the max
        // This approach is generally preferable to the other, as it increases numerical stability by eliminating a subtraction for each CDF.
        #region Helper Methods
        private static Normal[] NegateDistributions(Normal[] distributions)
        {
            Normal[] output = new Normal[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                output[i] = new Normal(mean: -1 * distributions[i].Mean, stddev: distributions[i].StdDev);
            }
            return output;
        }
        #endregion

        public static double[] ComplementsTrapezoid(Normal[] distributions, int steps = 150)
        {
            double[] complementProbs = new double[distributions.Length];
            distributions = NegateDistributions(distributions); // This change is local to this method

            for (int i = 0; i < distributions.Length; i++)
            {
                Normal distribution_i = distributions[i];

                Func<double, double> integrand = x =>
                {
                    double product = distribution_i.Density(x);
                    for (int j = 0; j < distributions.Length; j++)
                    {
                        if (j != i) { product *= distributions[j].CumulativeDistribution(x); }
                    }
                    return product;
                };

                complementProbs[i] = NewtonCotes.TrapezoidRule(integrand,
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, steps);
                //Console.WriteLine($"S38[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }

        public static double[] ComplementsSimpsons38(Normal[] distributions, int steps = 150)
        {
            if (steps % 3 != 0) steps += 3 - steps % 3; // Round up to a multiple of 3
            double[] complementProbs = new double[distributions.Length];
            distributions = NegateDistributions(distributions); // This change is local to this method

            for (int i = 0; i < distributions.Length; i++)
            {
                Normal distribution_i = distributions[i];

                Func<double, double> integrand = x =>
                {
                    double product = distribution_i.Density(x);
                    for (int j = 0; j < distributions.Length; j++)
                    {
                        if (j != i) { product *= distributions[j].CumulativeDistribution(x); }
                    }
                    return product;
                };

                complementProbs[i] = NewtonCotes.Simpsons38Rule(integrand,
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, steps);
                //Console.WriteLine($"S38[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }

        // Method with alt form using gauss hermite
        public static double[] ComplementsHermite(Normal[] distributions, int order)
        {
            double[] complementProbs = new double[distributions.Length];
            distributions = NegateDistributions(distributions); // This change is local to this method

            for (int i = 0; i < distributions.Length; i++)
            {
                Normal distribution_i = distributions[i];

                Func<double, double> xOfZ = (double z) => distribution_i.Mean + z * Math.Sqrt(2) * distribution_i.StdDev; // Converts Z to X

                Func<double, double> integrand = z =>
                {
                    double product = 1; // The density is already factored into the method
                    for (int j = 0; j < distributions.Length; j++)
                    {
                        if (j != i) { product *= distributions[j].CumulativeDistribution(xOfZ(z)); }
                    }
                    return product;
                };

                //complementProbs[i] = (1 / Math.Sqrt(Math.PI)) * GaussHermite.Integrate(integrand, order);
                complementProbs[i] = GaussHermite.Integrate(integrand, order);
                //Console.WriteLine($"GH[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        // Method with alt form using gauss legendre
        public static double[] ComplementsLegendre(Normal[] distributions, int order)
        {
            double[] complementProbs = new double[distributions.Length];
            distributions = NegateDistributions(distributions); // This change is local to this method

            for (int i = 0; i < distributions.Length; i++)
            {
                Normal distribution_i = distributions[i];

                Func<double, double> integrand = x =>
                {
                    double product = distribution_i.Density(x);
                    for (int j = 0; j < distributions.Length; j++)
                    {
                        if (j != i) { product *= distributions[j].CumulativeDistribution(x); }
                    }
                    return product;
                };

                complementProbs[i] = GaussLegendre.Integrate(integrand,
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, order);
                //complementProbs[i] = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(integrand,
                //    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, 96);
                //Console.WriteLine($"GL[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }

        public static double[] ComplementsClenshawCurtis(Normal[] distributions, int order)
        {
            // Change integral to alternative form by negating the distributions
            distributions = NegateDistributions(distributions); // This change is local to this method

            // Compute the evaluation points and weights
            double[] evalPoints = ClenshawCurtis.GetEvalPoints(order);
            double[] weights = ClenshawCurtis.GetWeights(order);


            // Compute the interval of integration
            double maxMean = distributions[0].Mean;
            double maxStdev = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                if (distributions[i].Mean > maxMean) { maxMean = distributions[i].Mean; }
                if (distributions[i].StdDev > maxStdev) { maxStdev = distributions[i].StdDev; }
            }
            // 8 standard deviations is just past the threshold beyond which normal PDFs are less than machine epsilon in double precision
            double intervalLowerLimit = maxMean - 8 * maxStdev;
            double intervalUpperLimit = maxMean + 8 * maxStdev;

            // Compute a linear transformation from that interval to [-1,1]
            double a = (intervalUpperLimit - intervalLowerLimit) / 2.0;
            double b = -1 * (2 * intervalLowerLimit / (intervalUpperLimit - intervalLowerLimit) + 1);
            double xOfz(double z) => (z - b) * a; // As z ranges over [-1,1], x will range over [iLL,iUL]

            // Compute the vector of constants
            double[] C = new double[evalPoints.Length];
            double[] X = new double[evalPoints.Length];
            for (int i = 0; i < C.Length; i++)
            {
                C[i] = weights[i];
                X[i] = xOfz(evalPoints[i]);
                for (int j = 0; j < distributions.Length; j++)
                {
                    C[i] *= distributions[j].CumulativeDistribution(X[i]);
                }
            }

            // --- Perform the Integration ---
            double[] complementProbs = new double[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                complementProbs[i] = 0;
                for (int j = 0; j < C.Length; j++)
                {
                    double CDFij = distributions[i].CumulativeDistribution(X[j]);
                    if (CDFij > 0)
                    {
                        complementProbs[i] += distributions[i].Density(X[j]) * C[j] / CDFij;
                    }
                }
                complementProbs[i] *= a; // Multiply by the derivative dx/dz
                Console.WriteLine($"CCAltInv[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }

        /// <summary> Uses nested Clenshaw-Curtis quadrature on the alternative form with an invariant to compute the probability that each element is the minimum for a set of normal distributions to within a user-specified precision </summary>
        /// <param name="distributions"> The set of normal distributions for which you want to compute P(X = min X) </param>
        /// <param name="errorTolerance"> The maximum total error in P(X = min X) over all regions </param>
        /// <param name="maxIterations"> The maximum number of times the quadrature rule will be used, doubling in order each time </param>
        /// <returns></returns>
        public static double[] ComplementsClenshawCurtisAutomatic(Normal[] distributions, double errorTolerance = 10E-14, int maxIterations = 10)
        {
            // Change integral to alternative form by negating the distributions
            distributions = NegateDistributions(distributions); // This change is local to this method

            // Compute the interval of integration
            double maxOfMeanMinus8Stddev = distributions[0].Mean - 8 * distributions[0].StdDev;
            double maxOfMeanPlus8Stddev = distributions[0].Mean + 8 * distributions[0].StdDev;
            for (int i = 0; i < distributions.Length; i++)
            {
                maxOfMeanMinus8Stddev = Math.Max(maxOfMeanMinus8Stddev, distributions[i].Mean - 8 * distributions[i].StdDev);
                maxOfMeanPlus8Stddev = Math.Max(maxOfMeanPlus8Stddev, distributions[i].Mean + 8 * distributions[i].StdDev);
            }
            // 8 standard deviations is just past the threshold beyond which normal PDFs are less than machine epsilon in double precision
            double intervalLowerLimit = maxOfMeanMinus8Stddev;
            double intervalUpperLimit = maxOfMeanPlus8Stddev;

            // Compute a linear transformation from that interval to [-1,1]
            double a = (intervalUpperLimit - intervalLowerLimit) / 2.0;
            //double b = -1 * (2 * intervalLowerLimit / (intervalUpperLimit - intervalLowerLimit) + 1); // TODO: Consider channging this to ILL + a, then z = az + b
            double b = intervalLowerLimit + a;
            //double xOfz(double z) => (z - b) * a; // As z ranges over [-1,1], x will range over [iLL,iUL]
            double xOfz(double z) => z * a + b; // As z ranges over [-1,1], x will range over [iLL,iUL]

            // --- Initialize the Vectors ---
            int order = 32; // Start with a 33-point CC rule
            double errorSum = double.PositiveInfinity;
            double[] errors = new double[distributions.Length];
            double[] complements = new double[distributions.Length];
            double[] weights = ClenshawCurtis.GetWeights(order);
            double[] X = ClenshawCurtis.GetEvalPoints(order); // Eval points in Z
            for (int i = 0; i < X.Length; i++) { X[i] = xOfz(X[i]); } // Convert from Z to X
            double[] C = new double[X.Length]; // The invariant product for each X value, without weights
            bool[] isFinished = new bool[distributions.Length]; // Keeps track of which regions are already at the desired precision
            for (int i = 0; i < C.Length; i++)
            {
                C[i] = 1;
                for (int j = 0; j < distributions.Length; j++)
                { C[i] *= distributions[j].CumulativeDistribution(X[i]); }
            }

            // --- Iterate higher order quadrature rules until desired precision is obtained ---
            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                // We will have three vectors X[], weights[], and C[] instead of two; the weights are now in weights[] instead of C[]
                // Each iteration replaces these vectors with expanded versions. Half + 1 of the entries are the old entries, and the other nearly half are freshly computed.
                // weights[] is the exception: it is completely replaced each time.

                double[] newComplements = new double[distributions.Length];

                // Update discard complement probability vector
                for (int i = 0; i < distributions.Length; i++)
                {
                    // Skip if this element is at the desired accuracy already
                    if (iteration > 1 && isFinished[i]) // errors[i] < errorTolerance / distributions.Length
                    {
                        newComplements[i] = complements[i];
                        continue;
                    }

                    newComplements[i] = 0;
                    for (int j = 0; j < C.Length; j++)
                    {
                        double CDFij = distributions[i].CumulativeDistribution(X[j]);
                        if (CDFij > 0)
                        {
                            newComplements[i] += distributions[i].Density(X[j]) * C[j] * weights[j] / CDFij;
                        }
                    }
                    newComplements[i] *= a; // Multiply by the derivative dx/dz
                }

                // Update the error
                if (iteration > 0)
                {
                    errorSum = 0;
                    for (int i = 0; i < errors.Length; i++)
                    {
                        double newError = Math.Abs(complements[i] - newComplements[i]);
                        // Detect if finished, which requires the error estimate for the ith term to be decreasing and less than its fair share of the total error tolerance
                        if (!isFinished[i]
                            && i > 1
                            && newError < errorTolerance / distributions.Length
                            && newError < errors[i]) { isFinished[i] = true; }
                        errors[i] = newError;
                        errorSum += errors[i];
                    }
                }
                complements = newComplements;

                // Determine if all regions appear to be finished refining their probability estimates
                //bool allRegionsFinished = false;
                //for (int i = 0; i < isFinished.Length; i++) { allRegionsFinished &= isFinished[i]; }

                // Check if all the probabilities add up to one
                double totalProb = 0;
                for (int i = 0; i < complements.Length; i++) { totalProb += complements[i]; }
                bool probsSumToOne = Math.Abs(totalProb - 1.0) < 1E-12;

                // Handle the end of the iteration
                if ((errorSum < errorTolerance && probsSumToOne /*&& allRegionsFinished*/) || iteration == maxIterations - 1)
                {
                    //Console.WriteLine($"Terminating on iteration {iteration} with remaining error estimate {errorSum}");
                    break; // Terminate and return complements
                }

                // Update the vectors for the next iteration
                order *= 2;
                weights = ClenshawCurtis.GetWeights(order);
                // Stretch the old arrays so there are gaps for the new entries
                double[] newX = new double[weights.Length];
                double[] newC = new double[weights.Length];
                for (int i = 0; i < X.Length; i++)
                {
                    newX[2 * i] = X[i];
                    newC[2 * i] = C[i];
                }
                // Add the new entries to X
                double[] entries = ClenshawCurtis.GetOddEvalPoints(order); // New entries in Z
                for (int i = 0; i < entries.Length; i++)
                {
                    int slot = 2 * i + 1;
                    newX[slot] = xOfz(entries[i]); // Convert from Z to X
                    newC[slot] = 1;
                    for (int j = 0; j < distributions.Length; j++)
                    { newC[slot] *= distributions[j].CumulativeDistribution(newX[slot]); }
                }

                X = newX;
                C = newC;
            }

            return complements;
        }


        public static double[] ComplementsSimpson38Automatic(Normal[] distributions, double errorTolerance = 10E-14, int initialSteps = 33, int maxIterations = 10)
        {
            double[] output = ComplementsSimpsons38(distributions, initialSteps);
            for (int i = 1; i < maxIterations; i++)
            {
                double error = 0;
                for (int j = 0; j < output.Length; j++)
                {
                    error += output[j];
                }
                error = Math.Abs(1 - error);
                if (error < errorTolerance) { break; }

                output = ComplementsSimpsons38(distributions, (int)(initialSteps * Math.Pow(2, i)));
            }
            return output;
        }
        #endregion

        // Monte carlo here, don't mind me
        public static double[] MonteCarlo(Normal[] normals, int iterations = 1000000)
        {
            double[] bestCounts = new double[normals.Length];
            for (int i = 0; i < iterations; i++)
            {
                double bestObserved = normals[0].Sample();
                int bestIdx = 0;
                for (int j = 1; j < normals.Length; j++)
                {
                    double observed = normals[j].Sample();
                    if (observed < bestObserved)
                    {
                        bestObserved = observed;
                        bestIdx = j;
                    }
                }
                bestCounts[bestIdx]++;
            }
            // Divide by the number of total observations and convert to discard probs
            for (int i = 0; i < normals.Length; i++)
            {
                bestCounts[i] = 1 - bestCounts[i] / iterations;
            }
            return bestCounts;
        }

    }
}
