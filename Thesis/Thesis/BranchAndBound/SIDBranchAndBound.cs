using MathNet.Numerics.Random;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace Thesis.BranchAndBound
{
    /// <summary> Encapsulates our optimization problem, combining a Branch that spans the solution space with a fitness function that measures elements of that space. </summary>
    /// <typeparam name="T"> The kind of object that constitutes an element of the search space </typeparam>
    public class SIDBranchAndBound<T>
    {
        /// <summary>
        /// A Branch object that contains all the feasible solutions to the problem, and defines the branching function
        /// </summary>
        public Branch SolutionSpace { get; private set; }

        /// <summary>
        /// Defines the fitness function
        /// </summary>
        public Func<T, double> FitnessFunction { get; private set; }

        public double BestFitnessObserved { get; private set; }
        public object BestElementObserved { get; private set; }

        public SIDBranchAndBound(Branch SolutionSpace, Func<T,double> FitnessFunction)
        {
            this.SolutionSpace = SolutionSpace;
            this.FitnessFunction = FitnessFunction;
        }
        
        public Branch[] BranchAndBound(int sampleSize, int iterations, double confidenceLevel = 0.95, bool cullDuplicates = true)
        {
            // Sanity check the initial sample size
            if (sampleSize < 30) { throw new ArgumentOutOfRangeException("Sample size must be at least 30."); }

            double branchingFactor = 0; // Keeps track of average branches per active region per layer, which may not be constant if there is overlap between regions
            double BestObservedValue = double.PositiveInfinity;
            Branch[] activeBranches = new Branch[1] { SolutionSpace };
            List<Branch> discardedRegions = new List<Branch>(512);

            // Main Loop
            for (int iteration = 0; iteration < iterations; iteration++)
            {
                // --- Branch() Step ---
                var newActiveBranches = new List<Branch>();
                foreach (Branch branch in activeBranches) { newActiveBranches.AddRange(branch?.GetBranches()); }

                Program.logger.WriteLine($"Branched to produce {newActiveBranches.Count} new regions.");
                // Cull duplicates if desired
                if (cullDuplicates)
                {
                    for (int i = 1; i < newActiveBranches.Count; i++)
                    {
                        for (int j = 0; j < i; j++)
                        {
                            if (newActiveBranches[i].Equals(newActiveBranches[j]))
                            {
                                newActiveBranches.RemoveAt(i);
                                i--;
                                break;
                            }
                        }
                    }
                    Program.logger.WriteLine($"Duplicates culled. {newActiveBranches.Count} regions remain.");
                }

                // Update branching factor
                branchingFactor = (branchingFactor * iteration + newActiveBranches.Count / activeBranches.Length) / (iteration + 1);
                Program.logger.WriteLine($"Branching factor revised to {branchingFactor}");
                activeBranches = newActiveBranches.ToArray();

                // Set up random number generators for each task in a deterministic way
                BranchSamplingBatch[] batch = new BranchSamplingBatch(,rand);
                for (int i = 0; i < rands.Length; i++) { rands[i] = new Xoshiro256StarStar(SolutionSpace.rand.Next()); }

                // Sample the regions (now with 1100% more threads!)
                Task[] tasks = new Task[activeBranches.Length];
                for (int i = 0; i < activeBranches.Length; i++)
                {
                    int j = i; // To avoid issues with lambda closure, one of the few quirks of this language
                    tasks[j] = Task.Run(() => { activeBranches[j].GetRandomElement(rands[j]); });
                }
                Task.WaitAll(tasks);
                //foreach (Region reg in activeRegions) { reg.Sample(sampleSize); }

            }

        }


    }
}
