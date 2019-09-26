using System;

namespace Thesis.Quadrature
{
    static class Statistics
    {
        public static double Minimum(double[] data)
        {
            double min = data[0];
            for (int i = 1; i < data.Length; i++) { min = Math.Min(min, data[i]); }
            return min;
        }

        public static double Mean(double[] data)
        {
            double sum = data[0];
            for (int i = 1; i < data.Length; i++) { sum += data[i]; }
            return sum / data.Length;
        }

        public static double VarianceEstimate(double[] sample)
        {
            double mean = Mean(sample);
            double sum = 0;
            for (int i = 0; i < sample.Length; i++) { sum += Math.Pow(sample[i] - mean, 2); }
            return sum / (sample.Length - 1);
        }

        public static double Variance(double[] data)
        {
            double mean = Mean(data);
            double sum = 0;
            for (int i = 0; i < data.Length; i++) { sum += Math.Pow(data[i] - mean, 2); }
            return sum / data.Length;
        }

        public static double Quantile(double[] data, double q)
        {
            if (q < 0 || q > 1) { throw new ArgumentOutOfRangeException($"Desired percentile is out of range: {q}"); }
            if (q == 0) { return data[0]; }
            if (q == 1) { return data[data.Length - 1]; }

            double product = (data.Length - 1) * q;
            int idx = (int) product;
            //return data[idx + 1] * (product - idx) + data[idx] * (idx + 1 - product);
            return Interpolation.Lerp(idx, data[idx], idx + 1, data[idx + 1], product);
        }

        public static double MeanOfLowerHalf(double[] data)
        {
            double median = Median(data);
            double sum = 0;
            int count = 0;
            for (int i = 0; i < data.Length; i++)
            {
                if (data[i] < median)
                {
                    sum += data[i];
                    count++;
                }
            }
            return sum / count;
        }

        public static double Median(double[] data)
        {
            if (data.Length % 2 == 1) return data[data.Length / 2];
            return (data[data.Length / 2 - 1] + data[data.Length / 2]) / 2;
        }
    }
}
