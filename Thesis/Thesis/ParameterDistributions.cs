using System;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;

namespace Thesis
{
    public enum GuidingParameter { Mean, Median, LowerMean, OneOverNthQuantile }
    static class ParameterDistributions
    {
        public static Normal MeanCLT(double[] data)
        {
            return new Normal(Statistics.Mean(data), Math.Sqrt(Statistics.VarianceEstimate(data) / data.Length));
        }

        public static Normal QuantileBootstrap(double[] sortedData, double q, Random rand = null)
        {
            if (rand == null) rand = Program.rand;

            // Get a point estimate of the quantile
            double quantileEstimate = Statistics.Quantile(sortedData, q);

            // Bootstrap the variance of the distribution
            double[] observations = new double[250];
            double[] bootstrapSample = new double[sortedData.Length];
            for (int i = 0; i < 250; i++)
            {
                for (int j = 0; j < sortedData.Length; j++)
                {
                    bootstrapSample[j] = sortedData[rand.Next(sortedData.Length)];
                }
                observations[i] = Statistics.Quantile(bootstrapSample, q);
            }
            double varianceEstimate = Statistics.VarianceEstimate(observations);

            return new Normal(quantileEstimate, varianceEstimate);
        }

        public static Normal MedianBootstrapMemoryFriendly(double[] sortedData, double[] bootstrapStorage, Random rand = null)
        {
            if (rand == null) rand = Program.rand;

            // Get a point estimate of the quantile
            double medianEstimate = Statistics.Median(sortedData);

            // Bootstrap the variance of the distribution
            double[] bootstrapSample = new double[sortedData.Length];
            for (int i = 0; i < bootstrapStorage.Length; i++)
            {
                for (int j = 0; j < sortedData.Length; j++)
                {
                    bootstrapSample[j] = sortedData[rand.Next(sortedData.Length)];
                }
                Sorting.Sort(bootstrapSample);
                bootstrapStorage[i] = Statistics.Median(bootstrapSample);
            }
            double varianceEstimate = Statistics.VarianceEstimate(bootstrapStorage);

            return new Normal(medianEstimate, varianceEstimate);
        }

        public static Normal MeanOfLessThanQuantile(double[] data, double q)
        {
            // Sort the data values in increasing order
            List<double> sortedData = new List<double>(data);
            sortedData.Sort();

            // Find the qth quantile and remove everything above that
            double qthQuantile = Statistics.Quantile(sortedData, q);
            sortedData.RemoveAll(x => x > qthQuantile);

            return MeanCLT(sortedData.ToArray());
        }

        public static GEV SampleMinimumMemoryFriendly(double[] data, double[] bootstrapStorage, Random rand = null)
        {
            if (rand == null) rand = Program.rand;

            // Start by computing a tail estimate. The PBDH theorem says this should be GPD shaped. 
            // We are using a small amount of smoothing on the ECDF as well here
            var pickandsApprox = new PickandsApproximation(data, method: PickandsApproximation.FittingMethod.BFGS_MSE);
            // Bootstrap observations of the max under this model
            for (int i = 0; i < bootstrapStorage.Length; i++)
            {
                double max = double.NegativeInfinity;
                for (int j = 0; j < data.Length; j++) // Same number of observations as original sample
                {
                    max = Math.Max(max, pickandsApprox.Sample());
                }
                bootstrapStorage[i] = max;
            }

            // --- Optimize to find the best-fit GEV model for these observations ---

            #region Helper Methods
            double FitnessApproxModel(GEV model)
            {
                double val = 0;
                for (int i = 0; i < bootstrapStorage.Length; i++)
                {
                    val += Math.Pow(model.CumulativeDistribution(bootstrapStorage[i]) - i * 1.0 / bootstrapStorage.Length, 2);
                }
                return val;
            }

            GEV OptimizeBFGS(Func<Vector<double>, double> objectiveFunc, double initialShape, double initialScale, double initialLocation)
            {
                // Formatted by shape, scale, location
                var lowerBounds = CreateVector.DenseOfArray(new double[] { -10, Math.Min(-3 * initialScale, 3 * initialScale), Math.Min(-3 * initialLocation, 3 * initialLocation) });
                var upperBounds = CreateVector.DenseOfArray(new double[] { 10, Math.Max(-3 * initialScale, 3 * initialScale), Math.Max(-3 * initialLocation, 3 * initialLocation) });
                var initialGuess = CreateVector.DenseOfArray(new double[] { initialShape, initialScale, initialLocation });

                var min = FindMinimum.OfFunctionConstrained(objectiveFunc, lowerBounds, upperBounds, initialGuess);
                return new GEV(min[2], min[1], min[0], Program.rand);
            }
            #endregion

            // Initial guesses (based on the case when shape = 0)
            double scaleGuess = Math.Sqrt(6 * Statistics.Variance(bootstrapStorage)) / Math.PI; // 
            double locationGuess = Statistics.Median(bootstrapStorage) + scaleGuess * Math.Log(Math.Log(2));

            double ObjectiveFunction(Vector<double> x) => FitnessApproxModel(new GEV(x[2], x[1], x[0], Program.rand));
            GEV fittedApproxModel = OptimizeBFGS(ObjectiveFunction, pickandsApprox.c, scaleGuess, locationGuess);
            return new GEV()

        }

    }
}
