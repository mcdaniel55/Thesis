using System;
using MathNet.Numerics.Distributions;

namespace Thesis.BranchAndBound
{
    /// <summary> Encapsulates a subset of the solution space that can be branched, and on collections of which a bounding operation can be performed </summary>
    public abstract class Branch
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
        public object GetRandomElement() => RandomElement(); // This should be hidden by a new method in the implementation
        protected abstract object RandomElement(); 

        /// <summary> Creates a cover of this subset according to problem-specific structure </summary>
        public abstract Branch[] GetBranches();

        public virtual void EvaluateRandomElements(object[] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                array[i] = GetRandomElement();
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
