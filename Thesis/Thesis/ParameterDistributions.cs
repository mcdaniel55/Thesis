using System;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;

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
            double qthQuantile = Statistics.Quantile(data, q);
            sortedData.RemoveAll(x => x > qthQuantile);

            return MeanCLT(sortedData.ToArray());
        }

        public static GEV SampleMinimumMemoryFriendly(double[] data, double[] bootstrapStorage, Random rand = null)
        {
            throw new NotImplementedException();
            if (rand == null) rand = Program.rand;
        }

    }
}
