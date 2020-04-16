using System;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;
using MathNet.Numerics;

namespace Thesis
{
    public class GEV : IContinuousDistribution
    {
        public readonly double location, scale, shape;
        public Random RandomSource { get; set; }
        private const double SHAPE_EPSILON = 1E-6;

        public GEV(double location, double scale, double shape, Random rand = null)
        {
            RandomSource = rand ?? Program.rand;
            this.location = location;
            this.scale = scale;
            this.shape = shape;
        }

        public double Mode
        {
            get
            {
                if (Math.Abs(shape) < SHAPE_EPSILON) return location;
                return location + scale * (Math.Pow(1 + shape, -shape) - 1) / shape;
            }
        }

        public double Minimum
        {
            get
            {
                if (shape <= SHAPE_EPSILON) return double.NegativeInfinity;
                return location - scale / shape;
            }
        }

        public double Maximum
        {
            get
            {
                if (shape >= -SHAPE_EPSILON) return double.PositiveInfinity;
                return location - scale / shape;
            }
        }

        public double Mean
        {
            get
            {
                if (Math.Abs(shape) < SHAPE_EPSILON) return location + scale * Constants.EulerMascheroni;
                if (shape < 1) return location + scale * (SpecialFunctions.Gamma(1 - shape) - 1) / shape;
                return double.PositiveInfinity;
            }
        }

        public double Variance
        {
            get
            {
                if (Math.Abs(shape) < SHAPE_EPSILON) return scale * scale * Math.PI * Math.PI / 6;
                if (shape >= 0.5) return double.PositiveInfinity;
                double g1 = SpecialFunctions.Gamma(1 - shape);
                return scale * scale * (SpecialFunctions.Gamma(1 - 2 * shape) - g1 * g1) / (shape * shape);
            }
        }

        public double StdDev
        {
            get
            {
                return Math.Sqrt(Variance);
            }
        }

        public double Entropy
        {
            get
            {
                return Math.Log(scale) + Constants.EulerMascheroni * (shape + 1) + 1; // Needs testing; probably Ln, but might be Log10
            }
        }

        public double Skewness
        {
            get
            {
                /*
                if (shape == 0) return 12 * Math.Sqrt(6) * 1.20205690315959428540d / (Math.PI * Math.PI);
                if (shape < 0.3333333333333333d)
                {
                    double g1 = SpecialFunctions.Gamma(1 - shape);
                    double g2 = SpecialFunctions.Gamma(1 - 2 * shape);
                    double g3 = SpecialFunctions.Gamma(1 - 3 * shape);
                    return Math.Sign(shape) * (g3 - 3 * g2 * g1 + 2 * Math.Pow(g1,3)) / Math.Pow(g2 - g1 * g1,1.5);
                }
                throw new InvalidOperationException("The skewness of the GEV distribution is not defined for shape parameters greater than or equal to 1/3");*/
                throw new NotImplementedException();
            }
        }

        public double Median
        {
            get
            {
                if (Math.Abs(shape) < SHAPE_EPSILON) return location - scale * Math.Log(Math.Log(2));
                return location + scale * (Math.Pow(Math.Log(2), -shape) - 1) / shape;
            }
        }

        public double CumulativeDistribution(double x)
        {
            double s = (x - location) / scale;
            if (Math.Abs(shape) < SHAPE_EPSILON) return Math.Exp(-Math.Exp(-s));
            if (shape > 0 && s <= -1.0 / shape) return 0;
            if (shape < 0 && s >= -1.0 / shape) return 1;

            return Math.Exp(-Math.Pow(1 + shape * s, -1.0 / shape));
        }

        public double Density(double x)
        {
            double s = (x - location) / scale;
            if (Math.Abs(shape) < SHAPE_EPSILON) return Math.Exp(-s) * Math.Exp(-Math.Exp(-s)) / scale;
            if (shape >= SHAPE_EPSILON && s <= -1.0 / shape) return 0;
            if (shape <= -SHAPE_EPSILON && s >= -1.0 / shape) return 0;
            return Math.Pow(1 + shape * s, -1.0 / shape - 1) * Math.Exp(-Math.Pow(1 + shape * s, -1.0 / shape)) / scale; // The / scale here is from the chain rule on the transformation S(x)
        }

        public double DensityLn(double x)
        {
            double s = (x - location) / scale;
            if (Math.Abs(shape) < SHAPE_EPSILON) return -Math.Exp(-s) - s - Math.Log(scale);
            if (shape >= SHAPE_EPSILON && s <= -1.0 / shape) return double.NegativeInfinity;
            if (shape <= -SHAPE_EPSILON && s >= -1.0 / shape) return double.NegativeInfinity;
            return (1.0 / 1 - shape) * Math.Log(1 + shape * s) - Math.Pow(1 + shape * s, -1.0 / shape) - Math.Log(scale);
        }

        public double InverseCumulativeDistribution(double q)
        {
            if (Math.Abs(shape) < SHAPE_EPSILON)
            {
                if (q == 0) return double.NegativeInfinity;
                if (q == 1) return double.PositiveInfinity;
                return location - scale * Math.Log(-Math.Log(q));
            }
            if (shape > 0 && q == 1) return double.PositiveInfinity;
            if (shape < 0 && q == 0) return double.NegativeInfinity;
            return location + scale * (Math.Pow(-Math.Log(q), -shape) - 1) / shape;
        }

        public double Sample()
        {
            return InverseCumulativeDistribution(RandomSource.NextDouble());
        }

        public void Samples(double[] values)
        {
            for (int i = 0; i < values.Length; i++)
            {
                values[i] = InverseCumulativeDistribution(RandomSource.NextDouble());
            }
        }

        public IEnumerable<double> Samples()
        {
            throw new NotImplementedException();
        }

        public override string ToString()
        {
            return $"GEV Location {location} Scale {scale} Shape {shape}";
        }
    }
}
