﻿using System;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;

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
            var pickandsApprox = new PickandsApproximation(data, method: PickandsApproximation.FittingMethod.Pickands_SupNorm);
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
            double FitnessSquaredError(GEV model)
            {
                double val = 0;
                for (int i = 0; i < bootstrapStorage.Length; i++)
                {
                    double deviation = model.CumulativeDistribution(bootstrapStorage[i]) - i * 1.0 / bootstrapStorage.Length;
                    val += deviation * deviation;
                }
                return val;
            }

            GEV OptimizeBFGS(Func<Vector<double>, double> objectiveFunc, double initialShape, double initialScale, double initialLocation)
            {
                // Formatted by shape, scale, location
                var lowerBounds = CreateVector.DenseOfArray(new double[] { -3, Math.Min(-3 * initialScale, 3 * initialScale), Math.Min(-3 * initialLocation, 3 * initialLocation) });
                var upperBounds = CreateVector.DenseOfArray(new double[] { 3, Math.Max(-3 * initialScale, 3 * initialScale), Math.Max(-3 * initialLocation, 3 * initialLocation) });
                var initialGuess = CreateVector.DenseOfArray(new double[] { initialShape, initialScale, initialLocation });

                //var min = FindMinimum.OfFunctionConstrained(objectiveFunc, lowerBounds, upperBounds, initialGuess, 1E-05, 1E-03, 0.1);

                var result = new BfgsBMinimizer(1E-02, 1E-02, 1E-01, 500).FindMinimum(ObjectiveFunction.Value(objectiveFunc), lowerBounds, upperBounds, initialGuess);
                var min = result.MinimizingPoint;

                return new GEV(min[2], min[1], min[0], rand);
            }
            #endregion

            // Initial guesses (based on the case when shape = 0)
            double scaleGuess = Math.Sqrt(6 * Statistics.Variance(bootstrapStorage)) / Math.PI; // 
            double locationGuess = Statistics.Median(bootstrapStorage) + scaleGuess * Math.Log(Math.Log(2));

        double ObjFunction(Vector<double> x) => FitnessSquaredError(new GEV(x[2], x[1], x[0], rand));
        GEV fittedApproxModel = OptimizeBFGS(ObjFunction, Math.Max(-3, Math.Min(pickandsApprox.c, 3)), scaleGuess, locationGuess);
            return new GEV(data[data.Length - 1], fittedApproxModel.scale, fittedApproxModel.shape, rand);
        }

    }
}
