using System;
using System.Collections.Generic;

namespace Thesis.BranchAndBound
{
    /// <summary>
    /// Encapsulates our optimization problem, combining a Branch that spans the solution space with a fitness function that measures elements of that space
    /// </summary>
    public class SIDBranchAndBound
    {
        /// <summary>
        /// A Branch object that contains all the feasible solutions to the problem, and defines the branching function
        /// </summary>
        public Branch SolutionSpace { get; private set; }

        /// <summary>
        /// Defines the fitness function
        /// </summary>
        public Func<object, double> FitnessFunction { get; private set; }



    }
}
