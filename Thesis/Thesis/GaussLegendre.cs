using System;

namespace Thesis.Quadrature
{
    static class GaussLegendre
    {
        public static double Integrate(Func<double, double> f, double intervalStart, double intervalEnd, double[] nodesAndWeights)
        {
            // Get a linear transformation from intervalStart to intervalEnd to [-1,1] 
            double a = (intervalEnd - intervalStart) / 2.0;
            double b = (intervalEnd + intervalStart) / 2.0;
            double xOfz(double z) => a * z + b;

            double sum = 0;
            for (int i = 0; i < nodesAndWeights.Length; i++)
            {
                sum += f(xOfz(nodesAndWeights[2 * i])) * nodesAndWeights[2 * i + 1];
            }

            return sum * a;
        }

        public static double Integrate(Func<double, double> f, double intervalStart, double intervalEnd, int order = 75)
        {
            // Get nodes and weights for the GL rule of appropriate order
            double[] nodeWeightArray = new double[2 * order];
            BogaertGLWrapper.GetGLNodesAndWeightsNonAlloc(nodeWeightArray);

            // Get a linear transformation from intervalStart to intervalEnd to [-1,1] 
            double a = (intervalEnd - intervalStart) / 2.0;
            double b = (intervalEnd + intervalStart) / 2.0;
            double xOfz(double z) => a * z + b;

            double sum = 0;
            for (int i = 0; i < order; i++)
            {
                sum += f(xOfz(nodeWeightArray[2 * i])) * nodeWeightArray[2 * i + 1];
            }

            return sum * a;
        }
    }
}
