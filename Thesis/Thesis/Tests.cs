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



    }
}
