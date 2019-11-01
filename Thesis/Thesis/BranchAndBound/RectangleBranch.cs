using System;

namespace Thesis.BranchAndBound
{
    public class RectangleBranch : Branch
    {
        /// <summary> The infimum of the x values in the rectangle </summary>
        public readonly double StartX;
        /// <summary> The supremum of the x values in the rectangle </summary>
        public readonly double EndX;
        /// <summary> The range of the x values contained in the rectangle </summary>
        private readonly double RangeX;
        /// <summary> The infimum of the y values in the rectangle </summary>
        public readonly double StartY;
        /// <summary> The supremum of the y values in the rectangle </summary>
        public readonly double EndY;
        /// <summary> The range of the y values contained in the rectangle </summary>
        private readonly double RangeY;

        public RectangleBranch(double StartX, double EndX, double StartY, double EndY, Random rand = null) : base(rand)
        {
            this.StartX = StartX;
            this.EndX = EndX;
            RangeX = EndX - StartX;
            this.StartY = StartY;
            this.EndY = EndY;
            RangeY = EndY - StartY;
        }

        /*
        /// <summary> Creates a new interval with the designated start and end points </summary>
        /// <param name="start"> A finite number less than end </param>
        /// <param name="end"> A finite number greater than start </param>
        /// <param name="rand"> An override for the default random number generator </param>
        public IntervalBranch(double start, double end, Random rand = null) : base(rand)
        {
            Start = start;
            End = end;
            range = end - start;
        }

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
        }*/


        public override Branch[] GetBranches()
        {
            // Quartering version
            return new Branch[]
            {
                new RectangleBranch(StartX, StartX + 0.5 * RangeX, StartY, StartY + 0.5 * RangeY, rand), // Bottom left
                new RectangleBranch(StartX + 0.5 * RangeX, EndX, StartY, StartY + 0.5 * RangeY, rand), // Bottom Right
                new RectangleBranch(StartX, StartX + 0.5 * RangeX, StartY + 0.5 * RangeY, EndY, rand), // Upper Left
                new RectangleBranch(StartX + 0.5 * RangeX, EndX, StartY + 0.5 * RangeY, EndY, rand), // Upper Right
            };
        }

        protected override object RandomElement() => GetRandomElement(); // Handles return type covariance
        public new Tuple<double,double> GetRandomElement()
        {
            double x = StartX + RangeX * rand.NextDouble();
            double y = StartY + RangeY * rand.NextDouble();
            return new Tuple<double, double>(x, y);
        }
    }
}
