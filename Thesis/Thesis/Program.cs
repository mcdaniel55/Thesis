using System;
using MathNet.Numerics.Random;
using MathNet.Numerics.Distributions;

namespace Thesis
{
    class Program
    {
        // Global references
        public static Random rand = new Xoshiro256StarStar(); // RNG
        public static Logger logger = new Logger("Output.txt"); // Logging

        static void Main(string[] args)
        {
            //Tests.TestCCQuadrature1();
            //Tables.MakeLegendreTable();
            //Tables.ShowGaussHermiteTable();
            //Tables.MakeCCTable();
            //Tables.ShowTrapTable();
            //Tables.SamplingDistOfMean();
            //Tables.LegendreExampleTable();
            //Tables.ClenshawCurtisExampleTable();
            //Tables.HermiteExampleTable();

            //Tests.GenerateHardDistributionSet();

            //Tables.ComparisonTablePairwise();
            //Tables.ComparisonTableEasySet();
            //Tables.ComparisonTableMediumSet();
            //Tables.ComparisonTableHardSet();

            GEV dist = new GEV(0, 1, 0.5);
            double q = dist.CumulativeDistribution(0.5);
            Console.WriteLine($"q = F(0.5): {q} f(0.5) {dist.Density(0.5)} Quantile(q): {dist.Quantile(q)}");


            Console.WriteLine("Done.");
            Console.ReadLine();
        }
    }
}
