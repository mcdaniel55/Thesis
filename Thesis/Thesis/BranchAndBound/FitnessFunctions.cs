using System;
using EMSOptimizationOverhaul;

namespace Thesis.BranchAndBound
{
    public static class FitnessFunctions
    {
        // --- Functions of one variable ---

        // Wicked comb on (-0.15, 0.15)
        /// <summary> Wicked Comb test function.
        /// Interval: [-0.15, 0.15]
        /// </summary>
        /// <param name="difficulty"> The number of peaks and valleys in the interval is directly proportional to this number; dictates most of the difficulty of the problem </param>
        /// <param name="x2Coefficient"> The coefficient of x^2; larger values can reduce the effect of rounding error and make the problem slightly easier </param>
        /// <returns></returns>
        public static double WickedComb(double x, double difficulty, double x2Coefficient)
        {
            return x2Coefficient * x * x - Math.Cos(1.0 / (Math.Abs(x) + 1.0 / (2 * Math.PI * (difficulty + 0.25))));
        }

        // Old test functions from the previous version
        public static double CubeSin(double x) => Math.Sin(x) * (x * x * x - x); // On (-2,2)

        public static double MixedSin(double x) => Math.Sin(5 * x) + 1.0 / 15 * x * x + Math.Cos(3 * x); // On any interval containing 1

        // The function from the beginning of the thesis, on 

        /// <summary>
        /// Example function from chapter 2
        /// Interval: [0, 3.2]
        /// Minima at 0, 2, and 2 * sqrt(2)
        /// </summary>
        /// <returns></returns>
        public static double ExampleFunction(double x) => 1.0 - Math.Cos(0.5 * Math.PI * x * x);

        // --- Functions of two variables ---
        /// <summary>
        /// Interval: -512 to 512 in both x and y
        /// Global minimum at x = 512, y = 404.2319, f(x,y) = -959.6407
        /// </summary>
        public static double EggHolder(double x, double y)
        {
            return -1 * (y + 47) * Math.Sin(Math.Sqrt(Math.Abs(x / 2 + y + 47))) - x * Math.Sin(Math.Sqrt(Math.Abs(x - y - 47)));
        }

        /// <summary>
        /// Interval: -10 to 10 in both x and y
        /// Global minimum at x = 1, y = 1, f(x,y) = 0
        /// </summary>
        public static double LevyN13(double x, double y)
        {
            return Math.Pow(Math.Sin(3 * Math.PI * x), 2)
                + Math.Pow(x - 1, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * y), 2))
                + Math.Pow(y - 1, 2) * (1 + Math.Pow(Math.Sin(2 * Math.PI * y), 2));
        }

        // --- Functions of 20 variables ---
        public static double EMSPlanFitness(int[] fullAmbs, int[] partAmbs)
        {
            // Test the plan
            Simulation sim = new Simulation(
                startTime: new DateTime(2016, 1, 1, 0, 0, 0), // These don't have to change when you use only 2016 or only 2017 calls
                endTime: new DateTime(2018, 1, 1, 0, 0, 0),
                fullAmbsToSpawn: fullAmbs,
                partAmbsToSpawn: partAmbs,
                centralized: true,
                speedMPH: 24f);
            sim.Run();

            // Compute the score
            return sim.MeanResponseTime * sim.Cost;
        }
    }
}
