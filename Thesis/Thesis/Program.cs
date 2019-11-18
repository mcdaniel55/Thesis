using System;
using MathNet.Numerics.Random;

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

            //Tables.ComparisonTablePairwise(figure: true);
            //Tables.ComparisonTableEasySet(figure: true);
            //Tables.ComparisonTableMediumSet(figure: true);
            //Tables.ComparisonTableHardSet(figure: true);

            //Tests.TestPickands();

            //Tests.TestGEV();

            //Tests.RunIntroOptimization();
            //Tests.RunWickedCombOptimization();
            Tests.RunEggholderOptimization();

            //Tests.TestNewTailFittingV4();
            //Tests.TestGEVComplementComputations();

            //Tests.TestNegatedParameterDistributionWrapper();

            Console.WriteLine("Done.");
            Console.ReadLine();
        }
    }
}
