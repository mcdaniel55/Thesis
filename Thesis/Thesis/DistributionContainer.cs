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
        double Quantile(double q);
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

        public double Quantile(double q) // Currently supports only these two
        {
            if (originalDistribution.GetType() == typeof(Normal))
            {
                return ((Normal)originalDistribution).InverseCumulativeDistribution(q);
            }
            if (originalDistribution.GetType() == typeof(GEV))
            {
                return ((GEV)originalDistribution).Quantile(q);
            }
            else throw new NotImplementedException($"Quantile function not defined for wrapped distribution type: {originalDistribution.GetType()}");
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

        public double Quantile(double q) // Currently supports only these two
        {
            if (originalDistribution.GetType() == typeof(Normal))
            {
                return -((Normal)originalDistribution).InverseCumulativeDistribution(q);
            }
            if (originalDistribution.GetType() == typeof(GEV))
            {
                return -((GEV)originalDistribution).Quantile(q);
            }
            else throw new NotImplementedException($"Quantile function not defined for wrapped distribution type: {originalDistribution.GetType()}");
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

    public struct ParameterDistribution : IDistributionWrapper
    {
        private readonly IContinuousDistribution originalDistribution;
        private readonly double estimate, upperBound, lowerBound;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="errorDistribution"> The probability distribution of the error of estimation thetaHat - theta, with 0 corresponding to no error. </param>
        /// <param name="estimate"> The observed value of thetaHat. </param>
        /// <param name="errorLowerBound"></param>
        /// <param name="errorUpperBound"></param>
        public ParameterDistribution(IContinuousDistribution errorDistribution, double estimate, double errorLowerBound, double errorUpperBound)
        {
            originalDistribution = errorDistribution;
            upperBound = estimate - errorLowerBound;
            lowerBound = estimate - errorUpperBound;
            this.estimate = estimate;
        }

        public double CumulativeDistribution(double x)
        {
            // 1 - F(-x)
            return 1 - originalDistribution.CumulativeDistribution(estimate - x);
        }

        public double Density(double x)
        {
            // f(-x)
            return originalDistribution.Density(estimate - x);
        }

        public double GetLowerBound()
        {
            return lowerBound;
        }

        public double GetUpperBound()
        {
            return upperBound;
        }

        public double GetEstimate()
        {
            return estimate;
        }

        public IContinuousDistribution GetWrappedDistribution()
        {
            return originalDistribution;
        }

        public double Sample()
        {
            return estimate - originalDistribution.Sample();
        }

        public double Quantile(double q) // Currently supports only these two
        {
            if (originalDistribution.GetType() == typeof(Normal))
            {
                return estimate - ((Normal)originalDistribution).InverseCumulativeDistribution(q);
            }
            if (originalDistribution.GetType() == typeof(GEV))
            {
                return estimate - ((GEV)originalDistribution).Quantile(q);
            }
            else throw new NotImplementedException($"Quantile function not defined for wrapped distribution type: {originalDistribution.GetType()}");
        }
    }

}
