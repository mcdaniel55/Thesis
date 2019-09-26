using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Distributions;

namespace ThesisOptNumericalTest
{
    static class NormalComparison
    {
        /// <summary> Takes two normal distributions and returns the probability that a point sample from of the first will be greater than one from the second. </summary>
        /// <param name="A">A normal distribution</param>
        /// <param name="B">Another normal distribution</param>
        /// <returns>P(A>B)</returns>
        public static double ComputeDiscardProbabilityPairwiseExact(Normal A, Normal B)
        {
            Normal BMinusA = new Normal(mean: B.Mean - A.Mean, 
                stddev: Math.Sqrt(A.Variance + B.Variance));
            return BMinusA.CumulativeDistribution(0);
        }

        /// <summary> Takes a set of normal distributions and computes the probability that each one would produce the minimal value if 
        /// a point sample was taken from each. Integration is performed by Simpson's 3/8 rule. </summary>
        /// <param name="distributions">The array of distributions to compare</param>
        /// <param name="iterations">The number of iterations to use in Simpson's 3/8 Rule. Defaults to 150.</param>
        /// <returns>An array of probabilities that each distribution will produce the minimum value if a point sample was taken from each.
        /// [i] = P(X_i less than min(all X_j))</returns>
        public static double[] ComputeDiscardComplementsSimpson(Normal[] distributions, int iterations = 150)
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

                complementProbs[i] = Integration.SimpsonsRule.Integrate(integrand, 
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, iterations);
                //Console.WriteLine($"S38R[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        // Method for direct computation with gauss hermite
        public static double[] ComputeDiscardComplementsGaussHermite(Normal[] distributions)
        {
            double[] complementProbs = new double[distributions.Length];

            for (int i = 0; i < distributions.Length; i++)
            {
                Normal distribution_i = distributions[i];

                Func<double, double> xOfZ = (double z) => distribution_i.Mean + z * Math.Sqrt(2) * distribution_i.StdDev; // Converts Z to X
                
                Func<double, double> integrand = z =>
                {
                    double product = 1; // distribution_i.Density(xOfZ(z)); // The density is already factored into the method
                    for (int j = 0; j < distributions.Length; j++)
                    {
                        if (j != i) { product *= 1 - distributions[j].CumulativeDistribution(xOfZ(z)); }
                    }
                    return product;
                };

                complementProbs[i] = (1 / Math.Sqrt(Math.PI)) * Integration.GaussHermite.Integrate(integrand, 99);
                Console.WriteLine($"GH[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        // Method for direct computation with gauss legendre
        public static double[] ComputeDiscardComplementsGaussLegendre(Normal[] distributions)
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

                complementProbs[i] = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(integrand, 
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, 96);
                Console.WriteLine($"GL[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        // Methods with alternative forms derived from negating the distributions' means and finding the probability that each is the max
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

        public static double[] ComputeDiscardComplementsSimpsonAlt(Normal[] distributions, int iterations = 150)
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

                complementProbs[i] = Integration.SimpsonsRule.Integrate(integrand, 
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, iterations);
                //Console.WriteLine($"S38RAlt[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }

        // Method with alt form using gauss hermite
        public static double[] ComputeDiscardComplementsGaussHermiteAlt(Normal[] distributions)
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

                complementProbs[i] = (1 / Math.Sqrt(Math.PI)) * Integration.GaussHermite.Integrate(integrand, 99);
                Console.WriteLine($"GHA[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        // Method with alt form using gauss legendre
        public static double[] ComputeDiscardComplementsGaussLegendreAlt(Normal[] distributions)
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

                complementProbs[i] = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(integrand, 
                    distribution_i.Mean - 8 * distribution_i.StdDev, distribution_i.Mean + 8 * distribution_i.StdDev, 96);
                Console.WriteLine($"GLA[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        // Methods with alt forms and invariants
        /// <summary> Uses Simpson's 3/8 Rule on the alternative integral form with an invariant to compute the probability that each of the distributions is the minimum </summary>
        /// <param name="distributions"> The set of distributions to consider </param>
        /// <param name="iterations"> The number of steps to use in Simpson's Rule </param>
        /// <returns> An array of probabilities indexed to match the distribution array </returns>
        /// <remarks> This will require up to 3N(1+3T) special function evaluations, where N is the number of distributions and T is the number of iterations. </remarks>
        public static double[] ComputeDiscardComplementsSimpson38AltInvariant(Normal[] distributions, int iterations = 150)
        {
            if (iterations < 3 || iterations % 3 != 0) { throw new ArgumentOutOfRangeException("Error: Iterations must be a positive multiple of three."); }

            distributions = NegateDistributions(distributions); // This change is local to this method

            // --- Compute the interval of integration ---
            // Get min and max mean and stdev for these distributions
            //double minMean = distributions[0].Mean;
            double maxMean = distributions[0].Mean;
            double maxStdev = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                //if (distributions[i].Mean < minMean) { minMean = distributions[i].Mean; }
                if (distributions[i].Mean > maxMean) { maxMean = distributions[i].Mean; }
                if (distributions[i].StdDev > maxStdev) { maxStdev = distributions[i].StdDev; }
            }
            double intervalLowerLimit = maxMean - 8 * maxStdev;
            double intervalUpperLimit = maxMean + 8 * maxStdev;
            double stepSize = (intervalUpperLimit - intervalLowerLimit) / iterations;

            // Compute the vector of constants
            double[] C = new double[1 + 3 * iterations];
            for (int i = 0; i < C.Length; i++)
            {
                C[i] = 1;
                double x = intervalLowerLimit + i * stepSize;
                for (int j = 0; j < distributions.Length; j++)
                {
                    C[i] *= distributions[j].CumulativeDistribution(x);
                }
            }
            // Apply the weights from Simpson's 3/8ths Rule
            for (int i = 1; i < iterations; i += 3)
            {
                C[i] *= 3;
                C[i+1] *= 3;
            }
            for (int i = 3; i < iterations - 1; i += 3)
            {
                C[i] *= 2;
            }

            // --- Perform the Integration ---
            double[] complementProbs = new double[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                complementProbs[i] = 0;
                for (int j = 0; j < C.Length; j++)
                {
                    double x = intervalLowerLimit + j * stepSize;
                    double CDFi = distributions[i].CumulativeDistribution(x);
                    if (CDFi > 0)
                    {
                        complementProbs[i] += distributions[i].Density(x) * C[j] / CDFi;
                    }
                }
                complementProbs[i] *= 3.0 * stepSize / 8.0;
                //Console.WriteLine($"S38RAltInv[{i}]: {complementProbs[i]}");
            }
            return complementProbs;
        }


        /// <summary> Uses Gauss-Legendre quadrature on the alternative integral form with an invariant to compute the probability that each of the distributions is the minimum </summary>
        /// <param name="distributions"> The set of distributions to consider </param>
        /// <param name="evalPoints"> The set of evaluation points </param>
        /// <param name="weights"> The set of weights, in the same order as the eval points </param>
        /// <returns> An array of probabilities indexed to match the distribution array </returns>
        /// <remarks> This will require up to 3NV special function evaluations, where N is the number of distributions and V is the number of evaluation points. </remarks>
        public static double[] ComputeDiscardComplementsGaussLegendreAltInvariant(Normal[] distributions, double[] evalPoints, double[] weights)
        {
            if (evalPoints.Length != weights.Length) { throw new ArgumentException("Error: Evaluation points must have same length as weights."); }
            int fevals = 0; // temp
            distributions = NegateDistributions(distributions); // This change is local to this method

            // Compute the interval of integration
            double minMean = distributions[0].Mean;
            double maxMean = distributions[0].Mean;
            double maxStdev = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                if (distributions[i].Mean < minMean) { minMean = distributions[i].Mean; }
                if (distributions[i].Mean > maxMean) { maxMean = distributions[i].Mean; }
                if (distributions[i].StdDev > maxStdev) { maxStdev = distributions[i].StdDev; }
            }
            double intervalLowerLimit = minMean - 8 * maxStdev;
            double intervalUpperLimit = maxMean + 8 * maxStdev;

            // Compute the transformation to [-1,1]
            double a = (intervalUpperLimit - intervalLowerLimit) / 2.0;
            double b = -1 * (2 * intervalLowerLimit / (intervalUpperLimit - intervalLowerLimit) + 1);
            Func<double, double> xOfz = z => (z - b) * a;

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
                    fevals++;
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
                    fevals++;
                    if (CDFij > 0)
                    {
                        complementProbs[i] += distributions[i].Density(X[j]) * C[j] / CDFij;
                        fevals++;
                    }
                }
                complementProbs[i] *= a; // Multiply by the derivative dx/dz
                Console.WriteLine($"GLAltInv[{i}]: {complementProbs[i]}");
            }
            Console.WriteLine($"Function evaluations: {fevals}");
            return complementProbs;
        }

        /// <summary> Uses Gauss-Hermite quadrature on the alternative integral form with an invariant to compute the probability that each of the distributions is the minimum </summary>
        /// <param name="distributions"> The set of distributions to consider </param>
        /// <param name="evalPoints"> The set of evaluation points </param>
        /// <param name="weights"> The set of weights, in the same order as the eval points </param>
        /// <returns> An array of probabilities indexed to match the distribution array </returns>
        /// <remarks> This will require up to 3NV special function evaluations, where N is the number of distributions and V is the number of evaluation points. </remarks>
        public static double[] ComputeDiscardComplementsGaussHermiteAltInvariant(Normal[] distributions, double[] evalPoints, double[] weights)
        {
            if (evalPoints.Length != weights.Length) { throw new ArgumentException("Error: Evaluation points must have same length as weights."); }
            distributions = NegateDistributions(distributions); // This change is local to this method

            // Compute the interval of integration
            double minMean = distributions[0].Mean;
            double maxMean = distributions[0].Mean;
            double maxStdev = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                if (distributions[i].Mean < minMean) { minMean = distributions[i].Mean; }
                if (distributions[i].Mean > maxMean) { maxMean = distributions[i].Mean; }
                if (distributions[i].StdDev > maxStdev) { maxStdev = distributions[i].StdDev; }
            }
            //double a = (minMean + maxMean) / 2; // Original
            double a = maxMean;
            //double b = (maxMean - minMean) / Math.Sqrt(2); // Original
            //double b = Math.Sqrt(2) * Math.Min((maxMean - minMean) / 2, maxStdev); // Worse
            double b = Math.Sqrt(2) * Math.Min(maxMean, maxStdev); // Better than original
            Normal U = new Normal(a, b); // Original
            //Normal U = new Normal(a, Math.Min((maxMean - minMean) / 2, maxStdev)); // Worse
            //Normal U = new Normal(a, Math.Min(maxMean, maxStdev)); // Better than the original

            // Compute the change of variable function
            Func<double, double> xOfz = z => b * z + a;

            // Compute the vector of constants
            double[] C = new double[evalPoints.Length];
            double[] X = new double[evalPoints.Length];
            for (int i = 0; i < C.Length; i++)
            {
                X[i] = xOfz(evalPoints[i]);
                C[i] = weights[i] / U.Density(X[i]);
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
                complementProbs[i] /= Math.Sqrt(Math.PI);
                Console.WriteLine($"GHAltInv[{i}]: {complementProbs[i]}");
            }

            return complementProbs;
        }


        public static double[] ComputeDiscardComplementsClenshawCurtisAltInvariant(Normal[] distributions, int order)
        {
            // Change integral to alternative form by negating the distributions
            distributions = NegateDistributions(distributions); // This change is local to this method

            // Compute the evaluation points and weights
            double[] evalPoints = Integration.ClenshawCurtis.getEvalPoints(order);
            double[] weights = Integration.ClenshawCurtis.getWeights(order);


            // Compute the interval of integration
            //double minMean = distributions[0].Mean;
            double maxMean = distributions[0].Mean;
            double maxStdev = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                //if (distributions[i].Mean < minMean) { minMean = distributions[i].Mean; }
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


        // Because CC is an easily nested rule, we can also re-iterate and improve the accuracy of the integration by adding more evaluations until we reach a desired level of precision.
        // If we continue to compute P(D_i) for all regions i, we can use the residual of their sum from 1 to estimate the total error in the quadrature.
        // We can get a great deal more accuracy for fewer evaluations if we are picky about which regions to compute the discard probability for and restrict the interval of integration appropriately.
        // As this prevents us from using the residual as an error estimate, we obtain an alternative estimate of the error in each region's P(D_i) by taking the difference of 
        // that value from that of the previous iteration. This also has the advantage of illuminating where in the vector the error is dense.
        // Unfortunately, it also requires us to know ahead of time some (hopefully small) superset of the regions we wish to discard. These may be computed via pairwise considerations, but for now
        // we will simply compute all of the P(D_i) values and use the residual error estimate.
        // Based on all that, we write a method that will use CC to evaluate P(D_i) with a specified level of precision.

        /// <summary> Uses nested Clenshaw-Curtis quadrature on the alternative form with an invariant to compute the probability that each element is the minimum for a set of normal distributions to within a user-specified precision </summary>
        /// <param name="distributions"> The set of normal distributions for which you want to compute P(X = min X) </param>
        /// <param name="errorTolerance"> The maximum total error in P(X = min X) over all regions </param>
        /// <param name="maxIterations"> The maximum number of times the quadrature rule will be used, doubling in order each time </param>
        /// <returns></returns>
        public static double[] ComputeDiscardComplementsClenshawCurtisAltInvariantAutomatic(Normal[] distributions, double errorTolerance = 10E-14, int maxIterations = 10)
        {
            // Change integral to alternative form by negating the distributions
            distributions = NegateDistributions(distributions); // This change is local to this method

            // Compute the interval of integration
            //double minMean = distributions[0].Mean;
            /*
            double maxMean = distributions[0].Mean;
            double maxStdev = 0;
            */
            double maxOfMeanMinus8Stddev = distributions[0].Mean - 8 * distributions[0].StdDev;
            //double minOfMeanPlus8Stddev = distributions[0].Mean + 8 * distributions[0].StdDev;
            double maxOfMeanPlus8Stddev = distributions[0].Mean + 8 * distributions[0].StdDev;
            for (int i = 0; i < distributions.Length; i++)
            {
                //if (distributions[i].Mean < minMean) { minMean = distributions[i].Mean; }
                /*
                if (distributions[i].Mean > maxMean) { maxMean = distributions[i].Mean; }
                if (distributions[i].StdDev > maxStdev) { maxStdev = distributions[i].StdDev; }
                */
                /*
                double mm8sd = distributions[i].Mean - 8 * distributions[i].StdDev;
                if (mm8sd > maxOfMeanMinus8Stddev) { maxOfMeanMinus8Stddev = mm8sd; }
                mm8sd = distributions[i].Mean + 8 * distributions[i].StdDev;
                if (mm8sd < minOfMeanPlus8Stddev) { minOfMeanPlus8Stddev = mm8sd; }
                */
                maxOfMeanMinus8Stddev = Math.Max(maxOfMeanMinus8Stddev, distributions[i].Mean - 8 * distributions[i].StdDev);
                //minOfMeanPlus8Stddev = Math.Min(minOfMeanPlus8Stddev, distributions[i].Mean + 8 * distributions[i].StdDev);
                maxOfMeanPlus8Stddev = Math.Max(maxOfMeanPlus8Stddev, distributions[i].Mean + 8 * distributions[i].StdDev);
            }
            // 8 standard deviations is just past the threshold beyond which normal PDFs are less than machine epsilon in double precision
            //double intervalLowerLimit = maxMean - 8 * maxStdev;
            double intervalLowerLimit = maxOfMeanMinus8Stddev;
            //double intervalUpperLimit = maxMean + 8 * maxStdev;
            double intervalUpperLimit = maxOfMeanPlus8Stddev;

            // Compute a linear transformation from that interval to [-1,1]
            double a = (intervalUpperLimit - intervalLowerLimit) / 2.0;
            double b = -1 * (2 * intervalLowerLimit / (intervalUpperLimit - intervalLowerLimit) + 1);
            double xOfz(double z) => (z - b) * a; // As z ranges over [-1,1], x will range over [iLL,iUL]

            // --- Initialize the Vectors ---
            int order = 32; // Start with a 33-point CC rule
            double errorSum = double.PositiveInfinity;
            double[] errors = new double[distributions.Length];
            double[] complements = new double[distributions.Length];
            double[] weights = Integration.ClenshawCurtis.getWeights(order);
            double[] X = Integration.ClenshawCurtis.getEvalPoints(order); // Eval points in Z
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
                    Console.WriteLine($"Terminating on iteration {iteration} with remaining error estimate {errorSum}");
                    break; // Terminate and return complements
                }

                // Update the vectors for the next iteration
                order *= 2;
                weights = Integration.ClenshawCurtis.getWeights(order);
                // Stretch the old arrays so there are gaps for the new entries
                double[] newX = new double[weights.Length];
                double[] newC = new double[weights.Length];
                for (int i = 0; i < X.Length; i++)
                {
                    newX[2 * i] = X[i];
                    newC[2 * i] = C[i];
                }
                // Add the new entries to X
                double[] entries = Integration.ClenshawCurtis.getOddEvalPoints(order); // New entries in Z
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


        public static double[] ComputeDiscardComplementsSimpson38AltInvariantAutomatic(Normal[] distributions, double errorTolerance = 10E-14, int initialSteps = 33, int maxIterations = 10)
        {
            double[] output = ComputeDiscardComplementsSimpson38AltInvariant(distributions, initialSteps); ;
            for (int i = 1; i < maxIterations; i++)
            {
                double error = 0;
                for (int j = 0; j < output.Length; j++)
                {
                    error += output[j];
                }
                error = Math.Abs(1 - error);
                if (error < errorTolerance) { break; }

                output = ComputeDiscardComplementsSimpson38AltInvariant(distributions, (int)(initialSteps * Math.Pow(2, i)));
            }
            return output;
        }

    }
}