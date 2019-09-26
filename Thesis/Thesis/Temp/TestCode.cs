using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Text;
using ThesisOptNumericalTest.Integration;

namespace ThesisOptNumericalTest
{
    class TestCode
    {
        public static void TestNumericalIntegration(Random rand)
        {
            BogaertGLWrapper.Initialize();

            /*
            double[] nodesAndWeightsRaw = BogaertGLWrapper.GetGLNodesAndWeights(10000000);
           
            // Count how many can be culled
            int count = 0;
            for (int i = 0; i < nodesAndWeightsRaw.Length / 2; i++)
            {
                double x = nodesAndWeightsRaw[2 * i];
                double w = nodesAndWeightsRaw[2 * i + 1];
                if (x < 0 && w * Normal.CDF(0, 1.0 / 8, x) < 10E-18) { count++; }
                else if (x > 0 && w * Normal.PDF(0, 1.0 / 8, x) < 10E-18) { count++; }
            }
            Console.WriteLine($"Could cull {count} out of {nodesAndWeightsRaw.Length / 2} evaluation points.");

            Console.ReadKey();
            */

            // Parameters
            int numberOfDistributions = 20;
            double minMeanFitness = 8;
            double maxMeanFitness = 60;
            double minStDev = 6;
            double maxStDev = 25;

            // Computed values
            double fitnessRange = maxMeanFitness - minMeanFitness;
            double stDevRange = maxStDev - minStDev;

            // Set up the distributions and pick the one with the biggest mean to be i
            Normal[] distributions = new Normal[numberOfDistributions];
            Normal distribution_i = new Normal();

            double minMean = -1 * minMeanFitness, maxMean = -1 * maxMeanFitness; // Starting points for finding min and max in the set (not an error)
            for (int i = 0; i < distributions.Length; i++)
            {
                distributions[i] = new Normal(-1 * (minMeanFitness + fitnessRange * rand.NextDouble()), minStDev + stDevRange * rand.NextDouble());
                if (distributions[i].Mean > maxMean) { maxMean = distributions[i].Mean; }
                if (distributions[i].Mean < minMean) { minMean = distributions[i].Mean; distribution_i = distributions[i]; }
                Console.WriteLine($"Dist {i}: mean {distributions[i].Mean}, stdev {distributions[i].StdDev}");
            }

            Func<double, double> altForm = x =>
            {
                double cdfi = distribution_i.CumulativeDistribution(x);
                if (cdfi == 0 || double.IsNaN(cdfi)) { return 0; }

                double product = distribution_i.Density(x) / cdfi;
                for (int i = 0; i < distributions.Length; i++)
                {
                    product *= distributions[i].CumulativeDistribution(x);
                }
                return product;
            };

            /*
            double correctResult = SimpsonsRule.Integrate(altForm, minMean - 3 * maxStDev, maxMean + 3 * maxStDev, 600);
            Console.WriteLine($"Simp 3/8 (600): 1 - P(D_i) = {correctResult}");

            double gaussLegendreResult = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(altForm, minMean - 3 * maxStDev, maxMean + 3 * maxStDev, 128);
            Console.WriteLine($"Gauss-Legendre 128: 1 - P(D_i) = {gaussLegendreResult}");
            */

            double[] discardProbs = new double[distributions.Length];
            /*
            for (int i = 0; i < distributions.Length; i++)
            {
                distribution_i = distributions[i];
                discardProbs[i] = SimpsonsRule.Integrate(altForm, minMean - 8 * maxStDev, maxMean + 8 * maxStDev, 1500);
                Console.WriteLine($"Simp 3/8 (1500): 1 - P(D_{i}) = {discardProbs[i]}");
            }*/

            List<double> discardProbList = new List<double>(discardProbs);
            discardProbList.Sort();
            double sum = 0;
            for (int i = 0; i < discardProbList.Count; i++)
            {
                //Console.WriteLine($"Sorted: {discardProbList[i]}");
                sum += discardProbList[i];
            }

            // Console.WriteLine($"Sum of probabilities is {sum}");


            /*
            NormalComparison.ComputeDiscardComplementsSimpson(distributions);
            NormalComparison.ComputeDiscardComplementsGaussHermite(distributions);
            NormalComparison.ComputeDiscardComplementsGaussLegendre(distributions);
            NormalComparison.ComputeDiscardComplementsSimpsonAlt(distributions);
            NormalComparison.ComputeDiscardComplementsGaussHermiteAlt(distributions);
            NormalComparison.ComputeDiscardComplementsGaussLegendreAlt(distributions);
            */

            List<double> output;
            sum = 0;
            output = new List<double>(NormalComparison.ComputeDiscardComplementsSimpson38AltInvariant(distributions, 450));
            output.Sort();

            for (int i = 0; i < output.Count; i++)
            {
                sum += output[i];
            }
            Console.WriteLine($"Sum of probabilities is {sum}");

            /*
            Console.WriteLine($"Dist 1 : 1/sqrt3 & 2 : 1/sqrt3");
            distributions = new Normal[] { new Normal(1, 1.0/Math.Sqrt(3)), new Normal(2, 1.0/ Math.Sqrt(3)) };
            Console.WriteLine($"Exact = {NormalComparison.ComputeDiscardProbabilityPairwiseExact(distributions[0], distributions[1])}");
            NormalComparison.ComputeDiscardComplementsGaussLegendreAlt(distributions);
            NormalComparison.ComputeDiscardComplementsSimpson38AltInvariant(distributions, 210);
            Console.WriteLine($"4 Dists 1 : 1 & 2 : 1");
            distributions = new Normal[10];
            for (int i = 0; i < distributions.Length - 1; i++) { distributions[i] = new Normal(1, 1); }
            distributions[distributions.Length - 1] = new Normal(2, 1);
            */
            /*
            NormalComparison.ComputeDiscardComplementsGaussLegendreAlt(distributions);
            output = new List<double>(NormalComparison.ComputeDiscardComplementsSimpson38AltInvariant(distributions, 210));
            output.Sort();
            sum = 0;
            for (int i = 0; i < output.Count; i++)
            {
                sum += output[i];
            }
            Console.WriteLine($"Sum of probabilities is {sum}");
            */


            output = new List<double>(NormalComparison.ComputeDiscardComplementsGaussLegendreAltInvariant(distributions, GaussLegendre.evalPoints75opt, GaussLegendre.weights75opt));
            output.Sort();
            sum = 0;
            for (int i = 0; i < output.Count; i++)
            {
                sum += output[i];
            }
            Console.WriteLine($"Sum of probabilities is {sum}");

            output = new List<double>(NormalComparison.ComputeDiscardComplementsGaussHermiteAltInvariant(distributions, GaussHermite.evaluationPoints70opt, GaussHermite.weights70opt));
            output.Sort();
            sum = 0;
            for (int i = 0; i < output.Count; i++)
            {
                sum += output[i];
            }
            Console.WriteLine($"Sum of probabilities is {sum}");

            output = new List<double>(NormalComparison.ComputeDiscardComplementsClenshawCurtisAltInvariant(distributions, 450));
            output.Sort();
            sum = 0;
            for (int i = 0; i < output.Count; i++)
            {
                sum += output[i];
            }
            Console.WriteLine($"Sum of probabilities is {sum}");


            System.Diagnostics.Stopwatch watch = new System.Diagnostics.Stopwatch();

            watch.Start();
            double[] bigTest = NormalComparison.ComputeDiscardComplementsClenshawCurtisAltInvariantAutomatic(distributions);
            watch.Stop();

            sum = 0;
            for (int i = 0; i < bigTest.Length; i++)
            {
                //Console.WriteLine($"CCAltInvAuto[{i}]: {bigTest[i]}");
                sum += bigTest[i];
            }
            Console.WriteLine($"Sum of probabilities is {sum}");
            Console.WriteLine($"Total Error lower bound: {Math.Abs(sum - 1)}");
            Console.WriteLine($"Time: {watch.Elapsed.TotalMilliseconds}ms");

            discardProbList = new List<double>(bigTest);
            discardProbList.Sort();
            {
                double certainty = 1;
                int idx = 0;
                while (true)
                {
                    double newval = certainty - discardProbList[idx];
                    if (newval < 0.95) break;
                    certainty = newval;
                    idx++;
                }
                Console.WriteLine($"Can discard {idx} distributions with 95% certainty");
            }

            watch.Restart();
            bigTest = NormalComparison.ComputeDiscardComplementsSimpson38AltInvariantAutomatic(distributions);
            watch.Stop();

            output = new List<double>(bigTest);

            output.Sort();
            sum = 0;
            for (int i = 0; i < output.Count; i++)
            {
                sum += output[i];
            }
            Console.WriteLine($"S38 Sum of probabilities is {sum}");
            Console.WriteLine($"S38 Time: {watch.Elapsed.TotalMilliseconds}ms");
        }
    }
}
