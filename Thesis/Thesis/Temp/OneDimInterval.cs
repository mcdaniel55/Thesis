using System;
using System.Collections.Generic;
using System.Text;
using ThesisOptNumericalTest.Optimization;

namespace ThesisOptNumericalTest
{
    class OneDimInterval : Region
    {
        public double Start { get; private set; }
        public double End { get; private set; }
        private double range;
        public OneDimInterval(double start, double end, Random rand) : base( rand)
        {
            // Trusting the user (me) to provide an end point that is greater than the start value
            this.Start = start;
            this.End = end;
            range = end - start;
        }

        public override Region[] Branch()
        {
            // Bisection version
            /*
            double midPoint = (end + start) / 2;
            return new Region[] { new FiniteInterval(start, midPoint, m_rand),
                new FiniteInterval(midPoint, end, m_rand) };
            */

            
            // Trisection version
            double mid1 = Start + range / 3;
            double mid2 = Start + 2 * range / 3;
            return new Region[] { new OneDimInterval(Start, mid1, m_rand),
                new OneDimInterval(mid1, mid2, m_rand),
                new OneDimInterval(mid2, End, m_rand) };
            

            // Quadrosection version
            /*
            double q1 = start + range / 4;
            double q2 = start + range / 2;
            double q3 = start + 3 * range / 4;
            return new Region[] { new FiniteInterval(start, q1, m_rand),
                new FiniteInterval(q1, q2, m_rand),
                new FiniteInterval(q2, q3, m_rand),
                new FiniteInterval(q3, end, m_rand) };
            */
        }

        protected override double SampleElement()
        {
            double samplePoint = Start + m_rand.NextDouble()*range;
            return MixedSinTest(samplePoint);
        }

        // Test Functions
        private static double CubeSinTest(double x) => Math.Sin(x) * (x*x*x - x); // On (-2,2)
        private static double MixedSinTest(double x) => Math.Sin(5 * x) + 1.0 / 15 * x * x + Math.Cos(3 * x); // On any interval containing 1
    }
}
