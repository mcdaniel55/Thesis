using MathNet.Numerics.Random;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace Thesis.BranchAndBound
{
    /// <summary> Encapsulates our optimization problem, combining a Branch that spans the solution space with a fitness function that measures elements of that space. </summary>
    /// <typeparam name="T"> The kind of object that constitutes an element of the search space </typeparam>
    public class SIDBranchAndBound<T>
    {
        public enum GuidingParameter { Mean, Median, LowerMean, OneOverNthQuantile }

        #region Properties
        /// <summary>
        /// A Branch object that contains all the feasible solutions to the problem, and defines the branching function
        /// </summary>
        public Branch SolutionSpace { get; set; }

        /// <summary>
        /// Defines the fitness function
        /// </summary>
        public Func<T, double> FitnessFunction { get; set; }

        public double BestFitnessObserved { get; private set; }

        public object BestElementObserved { get; private set; }

        public GuidingParameter GuidingParam { get; set; }
        #endregion

        #region Methods
        public SIDBranchAndBound(Branch SolutionSpace, Func<T, double> FitnessFunction, GuidingParameter GuidingParam)
        {
            this.SolutionSpace = SolutionSpace;
            this.FitnessFunction = FitnessFunction;
            this.GuidingParam = GuidingParam;
        }
        
        public Branch[] BranchAndBound(int sampleSize, int iterations, double confidenceLevel = 0.95, bool cullDuplicates = true)
        {
            throw new NotImplementedException();
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

                // Estimating distributions for the guiding parameter of each branch
                var parameterDistributions = new IContinuousDistribution[activeBranches.Length];

                // Sample the regions (now with more threads!)
                int threadCount = Environment.ProcessorCount; // # of logical processors, including hyperthreading etc.
                var batches = new BranchSamplingBatch<T>[activeBranches.Length];
                // This provides deterministic RNG without race conditions for multithreaded execution at a modest memory cost
                for (int i = 0; i < batches.Length; i++) { batches[i] = new BranchSamplingBatch<T>(activeBranches[i], new Xoshiro256StarStar(Program.rand.Next())); }

                var tasks = new Task[threadCount];
                for (int i = 0; i < tasks.Length; i++) { tasks[i] = null; } // Initialize the task array to null

                // Shared storage for all threads; each chunk of storage holds SampleSize many doubles, and is checked out for temporary use by a BranchSamplingBatch
                var storage = new double[threadCount][];
                for (int i = 0; i < threadCount; i++) { storage[i] = new double[sampleSize]; }

                // Maps thread indices in 'tasks' to arrays in 'storage'. A value of -1 means it's available
                //var storageMap = new int[threadCount]; 
                //for (int i = 0; i < storageMap.Length; i++) { storageMap[i] = -1; }

                int getAvailableTaskIndex()
                {
                    for (int i = 0; i < tasks.Length; i++) { if (tasks[i] == null) { return i; } }
                    return -1;
                }
                
                for (int i = 0; i < activeBranches.Length; i++)
                {
                    int taskIdx = getAvailableTaskIndex();
                    if (taskIdx > -1)
                    {
                        //storageMap[storageIdx] = tasks.Length; // Associate this new task with a storage location
                        tasks[taskIdx] = Task.Run(() => {
                            // Run the sampling
                            batches[i].SampleNonAlloc(storage[taskIdx], FitnessFunction, sampleSize);
                            // Compute the estimating distribution of the guiding parameter
                            
                            
                            
                            
                            
                            // ...
                        });
                    }
                    else
                    {
                        int completedIdx = Task.WaitAny(tasks);
                        if (completedIdx > -1) { tasks[completedIdx] = null; }
                    }
                }
                Task.WaitAll(tasks); // Once this has completed, we have sampled all of the batches.


            }


            



        }
        
        /* These are already defined in ParameterDistributions
        // Estimate the parameter distribution of the mean (based on CLT sampling distribution)
        private static Normal MeanGuidingParameterDistribution(double[] data)
        {
            throw new NotImplementedException();
        }

        // Estimate the parameter distribution of the median (based on bootstrap)
        private static Normal MedianGuidingParameterDistribution(double[] data)
        {
            throw new NotImplementedException();
        }

        // Estimate the parameter distribution of the mean of the lower half of the data (based on CLT sampling distribution)
        private static Normal LowerMeanGuidingParameterDistribution(double[] data)
        {
            throw new NotImplementedException();
        }

        private static GEV OneOverNthQuantileParameterDistribution(double[] data)
        {
            throw new NotImplementedException();
        }
        */

        #endregion
    }
}
