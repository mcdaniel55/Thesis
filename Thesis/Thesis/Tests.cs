using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Distributions;

namespace Thesis
{
    static class Tests
    {
        
        // --- Temporary Tests ---

        public static void TestCCQuadrature1()
        {
            double f(double x) => Math.Sqrt(x);

            double outcome = Quadrature.ClenshawCurtis.Integrate(f, 1, 4, 7);
            Console.WriteLine($"Outcome: {outcome}");
        }


        public static void GenerateHardDistributionSet()
        {
            for (int i = 0; i < 30; i++)
            {
                Program.logger.WriteLine("");
                double f() => 90 * Program.rand.NextDouble() + 30; // Means
                double g() => Math.Abs(Normal.Sample(20, 15)) + 0.3; // StdDevs
                for (int j = 0; j < 4; j++)
                {
                    Program.logger.Write($"new Normal({f()}, {g()}), ");
                }
            }
        }


    }
}
