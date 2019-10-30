using System;
using MathNet.Numerics.Distributions;

namespace Thesis
{
    /// <summary> Encapsulates a subset of the solution space that can be branched, and on collections of which a bounding operation can be performed </summary>
    abstract class Branch
    {
        #region Properties
        /// <summary> Keeps track of how many times the solution space was branched to produce this subset </summary>
        public int BranchLevel { get; private set; }
        public IContinuousDistribution SamplingDistribution { get; private set; }
        /// <summary> The best (smallest) observation of the fitness function </summary>
        public double MinimumObservedFitness { get; private set; }
        public object BestObservedSolution { get; private set; }
        public readonly Random m_rand;
        #endregion

        #region Methods
        /// <summary> Obtains a single value at random from this region </summary>
        protected abstract double Sample();

        /// <summary> Creates a cover of this subset according to problem-specific structure </summary>
        public abstract Branch[] GetBranches();

        public virtual void Samples(double[] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                array[i] = Sample();
            }
        }

        /// <summary> A constructor that specifies some common initialization logic and 
        /// provides an explicit parameterization of all subclass constructors </summary>
        protected Branch(Random rand)
        {
            m_rand = rand ?? Program.rand;
            MinimumObservedFitness = double.PositiveInfinity;
        }

        #endregion
    }
}
