using System;
using System.Collections.Generic;

namespace Thesis
{
    static class Statistics
    {
        public static double Minimum(IList<double> data)
        {
            double min = data[0];
            for (int i = 1; i < data.Count; i++) { min = Math.Min(min, data[i]); }
            return min;
        }

        public static double Mean(IList<double> data)
        {
            double sum = data[0];
            for (int i = 1; i < data.Count; i++) { sum += data[i]; }
            return sum / data.Count;
        }

        public static double MeanAbsoluteDeviation(IList<double> data)
        {
            double sum = 0;
            double mean = Mean(data);
            for (int i = 0; i < data.Count; i++)
            {
                sum += Math.Abs(data[i] - mean);
            }
            return sum / data.Count;
        }

        public static double VarianceEstimate(IList<double> sample)
        {
            double mean = Mean(sample);
            double sum = 0;
            for (int i = 0; i < sample.Count; i++) { sum += Math.Pow(sample[i] - mean, 2); }
            return sum / (sample.Count - 1);
        }

        public static double Variance(IList<double> data)
        {
            double mean = Mean(data);
            double sum = 0;
            for (int i = 0; i < data.Count; i++) { sum += Math.Pow(data[i] - mean, 2); }
            return sum / data.Count;
        }

        public static double Quantile(IList<double> sortedData, double q)
        {
            if (q < 0 || q > 1) { throw new ArgumentOutOfRangeException($"Desired percentile is out of range: {q}"); }
            if (q == 0) { return sortedData[0]; }
            if (q == 1) { return sortedData[sortedData.Count - 1]; }

            double product = (sortedData.Count - 1) * q;
            int idx = (int) product;
            //return data[idx + 1] * (product - idx) + data[idx] * (idx + 1 - product);
            return Interpolation.Lerp(idx, sortedData[idx], idx + 1, sortedData[idx + 1], product);
        }

        public static double MeanOfLowerHalf(IList<double> sortedData)
        {
            double median = Median(sortedData);
            double sum = 0;
            int count = 0;
            for (int i = 0; i < sortedData.Count; i++)
            {
                if (sortedData[i] < median)
                {
                    sum += sortedData[i];
                    count++;
                }
            }
            return sum / count;
        }

        public static double Median(IList<double> sortedData)
        {
            if (sortedData.Count % 2 == 1) return sortedData[sortedData.Count / 2];
            return (sortedData[sortedData.Count / 2 - 1] + sortedData[sortedData.Count / 2]) / 2;
        }
    }
}
