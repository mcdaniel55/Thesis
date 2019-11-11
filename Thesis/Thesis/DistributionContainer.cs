using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Text;

namespace Thesis
{
    /// <summary>
    /// A wrapper that pairs a distribution with its effective upper and lower bounds, as well as potentially applying a transformation
    /// </summary>
    public interface IDistributionWrapper
    {
        double CumulativeDistribution(double x);
        double Density(double x);
        double Sample();
        double GetUpperBound();
        double GetLowerBound();
        IContinuousDistribution GetWrappedDistribution();
    }

    public struct WrappedDistribution : IDistributionWrapper
    {
        private readonly IContinuousDistribution originalDistribution;
        private readonly double upperBound, lowerBound;
        public WrappedDistribution(IContinuousDistribution originalDistribution, double originalLowerBound, double originalUpperBound)
        {
            this.originalDistribution = originalDistribution;
            upperBound = originalUpperBound;
            lowerBound = originalLowerBound;
        }

        public double CumulativeDistribution(double x)
        {
            return originalDistribution.CumulativeDistribution(x);
        }

        public double Density(double x)
        {
            return originalDistribution.Density(x);
        }

        public double Sample()
        {
            return originalDistribution.Sample();
        }

        public static WrappedDistribution[] WrapDistributions(IContinuousDistribution[] distributions, double[] lowerBounds, double[] upperBounds)
        {
            WrappedDistribution[] result = new WrappedDistribution[distributions.Length];

            for (int i = 0; i < distributions.Length; i++)
            {
                result[i] = new WrappedDistribution(distributions[i], lowerBounds[i], upperBounds[i]);
            }

            return result;
        }

        public double GetUpperBound()
        {
            return upperBound;
        }

        public double GetLowerBound()
        {
            return lowerBound;
        }

        public IContinuousDistribution GetWrappedDistribution()
        {
            return originalDistribution;
        }
    }


    public struct NegatedDistribution : IDistributionWrapper
    {
        private readonly IContinuousDistribution originalDistribution;
        private readonly double upperBound, lowerBound;
        public NegatedDistribution(IContinuousDistribution originalDistribution, double originalLowerBound, double originalUpperBound)
        {
            this.originalDistribution = originalDistribution;
            upperBound = -originalLowerBound;
            lowerBound = -originalUpperBound;
        }

        public double CumulativeDistribution(double x)
        {
            return originalDistribution.CumulativeDistribution(-x);
        }

        public double Density(double x)
        {
            return originalDistribution.Density(-x);
        }

        public double Sample()
        {
            return -originalDistribution.Sample();
        }

        public static NegatedDistribution[] NegateDistributions(IContinuousDistribution[] distributions, double[] lowerBounds, double[] upperBounds)
        {
            NegatedDistribution[] result = new NegatedDistribution[distributions.Length];

            for (int i = 0; i < distributions.Length; i++)
            {
                result[i] = new NegatedDistribution(distributions[i], lowerBounds[i], upperBounds[i]);
            }

            return result;
        }

        public static IDistributionWrapper[] NegateNormalDistributions(Normal[] distributions)
        {
            var result = new IDistributionWrapper[distributions.Length];
            double epsilon = Math.Pow(2, -52);
            for (int i = 0; i < distributions.Length; i++)
            {
                result[i] = new NegatedDistribution(distributions[i], distributions[i].InverseCumulativeDistribution(epsilon), distributions[i].InverseCumulativeDistribution(1 - epsilon));
                //-distributions[i].Mean, distributions[i].StdDev, distributions[i].RandomSource);
            }
            return result;
        }

        public double GetUpperBound()
        {
            return upperBound;
        }

        public double GetLowerBound()
        {
            return lowerBound;
        }

        public IContinuousDistribution GetWrappedDistribution()
        {
            return originalDistribution;
        }
    }

    
}
