using System;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;
using Thesis;

namespace ThesisQuadratureTests
{
    static class CDFApprox
    {
        #region Old code
        /*
        public static void TestDistribution(IContinuousDistribution dist)
        {*/
        // Mess here, testing an idea
        // Try this with 2nd midpoints, also see what it does on different distributions
        // From (a,b) to (c,d) at x in [a,c]

        /*
        // From (a,c) to (b,d) at x in [a,b]
        double Lerp(double a, double b, double x, double c, double d)
        {
            var deltaXL = x - a;
            var deltaXR = b - x;
            var dist = b - a;
            var y = (deltaXL * d + deltaXR * c) / dist;
            return y;
        }*/
        // This is designed to output in a way that is easy to copy into excel for graphical analysis
        /*
        const int sampleSize = 200;
            const int smoothingPasses = sampleSize / 20;
            const double smoothingCoefficient = 0.4; // Less than 1
            Console.WriteLine("X, CDF, ECDF, ECDFSmooth");
            double[] s = new double[sampleSize];
            double[] cdf = new double[sampleSize];
            double[] ecdf = new double[sampleSize];
            double[] ecdfSmooth = new double[sampleSize];
            dist.Samples(s);
            List<double> sorter = new List<double>(s);
            sorter.Sort();
            s = sorter.ToArray();
            for (int i = 0; i < sampleSize; i++)
            {
                cdf[i] = dist.CumulativeDistribution(s[i]);
                ecdf[i] = (i + 1) * 1.0 / sampleSize;
                ecdfSmooth[i] = ecdf[i];
            }
            for (int i = 0; i < smoothingPasses; i++) // Smooth the ECDF
            {
                ecdfSmooth = Smoothing.SmoothPrototype(s, ecdfSmooth, smoothingCoefficient);
            }
            for (int i = 0; i < sampleSize; i++)
            {
                Console.WriteLine($"{s[i]}, {cdf[i]}, {ecdf[i]}, {ecdfSmooth[i]}");
            }
            
            Console.WriteLine("Mid, Approx1, Approx1Smooth, Midmid, Approx2, Approx2Smooth");
            // Midpoints partition the support into sampleSize many regions, each with about 1/sampleSize probabililty in them
            double[] midpoints = new double[sampleSize - 1];
            double[] midmidpoints = new double[sampleSize - 2]; // Midpoints of midpoints, or meta-midpoints, if you will
            double[] approx1 = new double[sampleSize - 1]; // Approximate CDF using midpoints
            double[] approx2 = new double[sampleSize - 2]; // Approximate CDF using midpoints of midpoints

            for (int i = 0; i < midpoints.Length; i++)
            {
                midpoints[i] = (s[i] + s[i + 1]) / 2.0;
                approx1[i] = (i + 1) * 1.0 / sampleSize;
            }
            for (int i = 0; i < midmidpoints.Length; i++)
            {
                midmidpoints[i] = (midpoints[i] + midpoints[i + 1]) / 2.0;
                approx2[i] = (i + 1) * 1.0 / (sampleSize - 1);
            }

            
            // Apply smoothing
            double[] adiff = Smoothing.SmoothPrototype(midpoints, approx1, smoothingCoefficient);
            double[] adiff2 = Smoothing.SmoothPrototype(midmidpoints, approx2, smoothingCoefficient);
            for (int i = 1; i < smoothingPasses; i++)
            {
                adiff = Smoothing.SmoothPrototype(midpoints, adiff, smoothingCoefficient);
                adiff2 = Smoothing.SmoothPrototype(midmidpoints, adiff2, smoothingCoefficient);
            }

            for (int i = 0; i < midpoints.Length; i++)
            {
                if (i < midmidpoints.Length)
                    Console.WriteLine($"{midpoints[i]}, {approx1[i]}, {adiff[i]}, {midmidpoints[i]}, {approx2[i]}, {adiff2[i]}");
                else
                    Console.WriteLine($"{midpoints[i]}, {approx1[i]}, {adiff[i]}");
            }

            Console.ReadLine();
        }
        */
        #endregion
        

        /*

        /// <summary>
        /// Creates an approximate distribution by constructing the ECDF of a sample
        /// </summary>
        /// <param name="sample"> A list of real values observed from the distribution to be approximated </param>
        /// <param name="rand"> The random number generator to be used by this distribution </param>
        /// <returns></returns>
        public static ContinuousDistribution ECDFFromSample(double[] sample, Random rand = null)
        {
            if (rand == null) rand = Program.rand;
            List<double> sampleSorted = new List<double>(sample);
            sampleSorted.Sort();
            List<double> cumulativeDensities = new List<double>(sample.Length);
            for (int i = 0; i < sample.Length; i++)
            {
                cumulativeDensities.Add((i+1) * 1.0 / sample.Length);
            }

            return new ContinuousDistribution(sampleSorted, cumulativeDensities, rand);
        }
        */
        /*
        /// <summary>
        /// Create an approximate distribution from a sample by placing 1/n over midpoints, smoothing, and extrapolating tails
        /// </summary>
        /// <param name="sample"> The sample values from which to infer a distribution </param>
        /// <param name="smoothingPasses"> The number of smoothing passes to be applied, defaults to sample size / 20 </param>
        /// <param name="smoothingCoefficient"> The strength of each smoothing pass, between 0 and 1 </param>
        /// <param name="extrapolationMode"> What kind of extrapolation, if any, is to be applied to produce the tails of the distribution</param>
        /// <param name="rand"> The random number generator this distribution will use </param>
        /// <returns></returns>
        public static ContinuousDistribution FromSample(double[] sample, int smoothingPasses = -1, double smoothingCoefficient = 0.3, ExtrapolationMode extrapolationMode = ExtrapolationMode.Moments, Random rand = null)
        {
            // Defaults
            if (smoothingPasses < 0) smoothingPasses = sample.Length / 20;
            if (rand == null) rand = Program.rand;

            List<double> sampleSorted = new List<double>(sample);
            sampleSorted.Sort();
            List<double> abscissas = new List<double>(sample.Length - 1);
            List<double> cumulativeDensities = new List<double>(sample.Length - 1);
            for (int i = 1; i < sample.Length; i++)
            {
                abscissas.Add((sampleSorted[i] + sampleSorted[i-1]) / 2); // Midpoints
                cumulativeDensities.Add(i * 1.0 / sample.Length); // Probability estimates before smoothing
            }

            if (smoothingPasses > 0 && smoothingCoefficient > 0) // Apply smoothing if desired
            {
                cumulativeDensities = Smoothing.SmoothIterative(abscissas, cumulativeDensities, smoothingPasses, smoothingCoefficient); 
            }

            switch (extrapolationMode)
            {
                case ExtrapolationMode.Moments:
                    double leftKurt = Statistics.SampleLeftKurtosis(sample); // Eventually change this to an estimator rather than a sample statistic
#if DEBUG
                    Console.WriteLine($"Estimated left-kurtosis: {leftKurt}");
#endif
                    if (leftKurt > 1.8) // If kurtosis is < 1.8, then it's a sheer drop, so there isn't really a tail and we can just extrapolate linearly
                    {
                        Extrapolation.ExtrapolateLeftDistributionTailByMoments(
                        ref abscissas,
                        ref cumulativeDensities,
                        Statistics.Mean(sample),
                        Statistics.LeftVarianceEstimate(sample),
                        leftKurt);
                    }
                    else goto case ExtrapolationMode.Linear;
                    break;

                case ExtrapolationMode.Linear:
                    Extrapolation.ExtrapolateDistributionTailsLinearly(abscissas, cumulativeDensities);
                    break;
                default: break;
            }
            */
            /*
            if (modelLeftTail)
            {
                double leftKurtosis = Statistics.SampleLeftKurtosis(sample); // Combine these into one function? Lkurt does Lvar anyway.
                double leftVariance = Statistics.LeftVarianceEstimate(sample);

                // Get a Tukey Lambda distribution with lambda chosen so that its kurtosis is equal to the left-kurtosis estimate from our data
                TukeyLambda tukey = new TukeyLambda(TukeyLambda.LambdaFromKurtosis(leftKurtosis), 2000);
                List<double> tailValues, tailAbscissas;
                double tailX = tukey.Quantile(1.0 / sample.Length); // Where to slice off the tail from TL
                int cutIndex = tukey.abscissas.BinarySearch(tailX); // Abscissa in TL immediately before tailX

                // If cutIndex < 0, then the cut is happening between two abscissas, and ~xPrev is the larger of their indices
                // Otherwise, the cut is happening at an abscissa, and xPrev is the index of that abscissa
                if (cutIndex < 0) cutIndex = ~cutIndex;
                // Now cutIndex is the next abscissa after the ones we want to keep in any case

                // Copy everything in the tail of TL up to but not including cutIndex
                tailAbscissas = tukey.abscissas.GetRange(0, cutIndex);
                tailValues = tukey.cumulativeDensities.GetRange(0, cutIndex);

                // --- Transform those values to where they need to go for our CDF ---
                // Scale the tail so the left-variance is equal to the estimate
                double scale = Math.Sqrt(leftVariance);
                // Map the ordinate at which we sliced off the tail in TL's support to the first midpoint in our estimate's support space
                double shift = abscissas[0] - tailX * scale; 
                for (int i = 0; i < tailAbscissas.Count; i++) { tailAbscissas[i] = tailAbscissas[i] * scale + shift; }

                // Append that tail to the beginning of our CDF
                tailAbscissas.AddRange(abscissas);
                abscissas = tailAbscissas;
                tailValues.AddRange(cumulativeDensities);
                cumulativeDensities = tailValues;

                // Note: You may want to try different values of lambda and scale to get the resulting distribution's kurtosis and 
                // variance to match the estimate, not just the TL you took the tail from.
            }
            */
            /*
            return new ContinuousDistribution(abscissas, cumulativeDensities, rand);
        }*/
    }
}
