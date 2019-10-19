using System;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;

namespace Thesis
{
    static class ParameterDistributions
    {
        public static Normal MeanCLT(double[] data)
        {
            return new Normal(Statistics.Mean(data), Math.Sqrt(Statistics.VarianceEstimate(data) / data.Length));
        }

        public static Normal QuantileBootstrap(double[] data, double q, Random rand = null)
        {
            if (rand == null) rand = Program.rand;

            // Get a point estimate of the quantile
            double quantileEstimate = Statistics.Quantile(data, q);

            // Bootstrap the variance of the distribution
            double[] observations = new double[250];
            double[] bootstrapSample = new double[data.Length];
            for (int i = 0; i < 250; i++)
            {
                for (int j = 0; j < data.Length; j++)
                {
                    bootstrapSample[j] = data[rand.Next(data.Length)];
                }
                observations[i] = Statistics.Quantile(bootstrapSample, q);
            }
            double varianceEstimate = Statistics.VarianceEstimate(observations);

            return new Normal(quantileEstimate, varianceEstimate);
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

        public static GEV SampleMinimum(double[] data)
        {
            throw new NotImplementedException();
        }

    }
}
