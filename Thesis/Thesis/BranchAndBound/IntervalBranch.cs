using System;

namespace Thesis.BranchAndBound
{
    public class IntervalBranch : Branch
    {
        /// <summary> The infimum of the interval </summary>
        public double Start { get; private set; }
        /// <summary> The supremum of the interval </summary>
        public double End { get; private set; }
        /// <summary> The range of the values contained in the interval </summary>
        private readonly double range;

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

        public override Branch[] GetBranches()
        {
            // Bisecting version
            double midPoint = (End + Start) / 2;
            return new Branch[] { new IntervalBranch(Start, midPoint, rand),
                new IntervalBranch(midPoint, End, rand) };


            // Trisecting version
            /*
            double mid1 = Start + range / 3;
            double mid2 = Start + 2 * range / 3;
            return new Branch[] { new IntervalBranch(Start, mid1, m_rand),
                new IntervalBranch(mid1, mid2, m_rand),
                new IntervalBranch(mid2, End, m_rand) };
            */

            // Quartering version
            /*
            double q1 = start + range / 4;
            double q2 = start + range / 2;
            double q3 = start + 3 * range / 4;
            return new Branch[] { new IntervalBranch(start, q1, m_rand),
                new IntervalBranch(q1, q2, m_rand),
                new IntervalBranch(q2, q3, m_rand),
                new IntervalBranch(q3, end, m_rand) };
            */
        }

        protected override object RandomElement() => GetRandomElement(); // Handles return type covariance
        public new double GetRandomElement() // Hides the base implementation
        {
            return Start + range * rand.NextDouble();
        }

        public override string ToString()
        {
            return $"Interval [{Start},{End}]";
        }
    }
}
