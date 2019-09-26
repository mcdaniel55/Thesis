using System;
using System.Collections.Generic;
using System.Text;
using ThesisOptNumericalTest.Optimization;

namespace ThesisOptNumericalTest
{
    /// <summary> Represents a two-dimensional region supporting a piecewise continuous function </summary>
    class TwoDimInterval : Region
    {
        public double StartX { get; private set; }
        public double EndX { get; private set; }
        private double rangeX;
        public double StartY { get; private set; }
        public double EndY { get; private set; }
        private double rangeY;

        public TwoDimInterval(double StartX, double EndX, double StartY, double EndY, Random rand) : base(rand)
        {
            this.StartX = StartX;
            this.StartY = StartY;
            this.EndX = EndX;
            this.EndY = EndY;
            rangeX = EndX - StartX;
            rangeY = EndY - StartY;
        }

        public override Region[] Branch()
        {
            // Quartering version
            return new Region[]
            {
                new TwoDimInterval(StartX, StartX + rangeX / 2, StartY, StartY + rangeY / 2, m_rand), // Bottom left
                new TwoDimInterval(StartX + rangeX / 2, EndX, StartY, StartY + rangeY / 2, m_rand), // Bottom Right
                new TwoDimInterval(StartX, StartX + rangeX / 2, StartY + rangeY / 2, EndY, m_rand), // Upper Left
                new TwoDimInterval(StartX + rangeX / 2, EndX, StartY + rangeY / 2, EndY, m_rand), // Upper Right
            };
        }

        protected override double SampleElement()
        {
            double x = StartX + rangeX * m_rand.NextDouble();
            double y = StartY + rangeY * m_rand.NextDouble();
            return LevyN13(x, y);
        }

        // Test functions
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
    }
}
