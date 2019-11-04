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

        public override string ToString()
        {
            return $"Rectangle X:[{StartX},{EndX}] Y:[{StartY},{EndY}]";
        }
    }
}
