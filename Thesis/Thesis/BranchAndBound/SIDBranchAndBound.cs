using MathNet.Numerics.Random;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics;

namespace Thesis.BranchAndBound
{
    /// <summary> Encapsulates our optimization problem, combining a Branch that spans the solution space with a fitness function that measures elements of that space. </summary>
    /// <typeparam name="T"> The kind of object that constitutes an element of the search space </typeparam>
    public class SIDBranchAndBound<T>
    {
        #region Properties
        private const int BootstrapSampleSize = 250;
        private static readonly double Epsilon16 = Math.Pow(2, -48);
        private static readonly double ComplementEpsilon16 = 1.0 - Math.Pow(2, -48);
        private static readonly double NormalBoundStdDevs = Normal.InvCDF(0, 1, ComplementEpsilon16); // Roughly 8 standard deviations
        private static readonly double PreemptiveDiscardConfidenceThreshold = 0.99999;

        /// <summary>
        /// A Branch object that contains all the feasible solutions to the problem, and defines the branching function
        /// </summary>
        public Branch SolutionSpace { get; set; }

        /// <summary>
        /// Defines the fitness function
        /// </summary>
        public Func<T, double> FitnessFunction { get; set; }

        public double BestFitnessObserved { get; private set; }

        public T BestElementObserved { get; private set; }

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
            // Sanity check the initial sample size
            if (sampleSize < 30) { throw new ArgumentOutOfRangeException("Sample size must be at least 30."); }

            double branchingFactor = 0; // Keeps track of average branches per active region per layer, which may not be constant if there is overlap between regions
            //double BestObservedValue = double.PositiveInfinity;
            Branch[] activeBranches = new Branch[] { SolutionSpace };
            //List<Branch> discardedRegions = new List<Branch>(512); // Tentative

            // === Main Loop ===
            for (int iteration = 0; iteration < iterations; iteration++)
            {
                #region Branch Step
                var newActiveBranches = new List<Branch>();
                foreach (Branch branch in activeBranches) { if (branch != null) newActiveBranches.AddRange(branch.GetBranches()); }

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

                // Storage for the estimating distribution of the guiding parameter of each branch
                var parameterDistributions = new IContinuousDistribution[activeBranches.Length];
                var upperBounds = new double[parameterDistributions.Length];
                var lowerBounds = new double[parameterDistributions.Length];
                #endregion

                #region Sampling
                // --- Sample the regions (now with more threads!) ---
                int threadCount = Environment.ProcessorCount; // # of logical processors, including hyperthreading etc.
                var batches = new BranchSamplingBatch<T>[activeBranches.Length];
                // This provides deterministic RNG without race conditions for multithreaded execution at a modest memory cost
                for (int i = 0; i < batches.Length; i++) { batches[i] = new BranchSamplingBatch<T>(activeBranches[i], new Xoshiro256StarStar(Program.rand.Next())); }

                var tasks = new Task[threadCount];
                for (int i = 0; i < tasks.Length; i++) { tasks[i] = null; } // Initialize the task array to null

                // Reusable storage for each thread; each chunk of storage holds SampleSize many doubles, and is checked out for temporary use by a BranchSamplingBatch
                var fitnessStorage = new double[threadCount][];
                for (int i = 0; i < threadCount; i++) { fitnessStorage[i] = new double[sampleSize]; }
                // Same for bootstrap sampling storage
                var bootstrapStorage = new double[threadCount][];
                for (int i = 0; i < threadCount; i++) { bootstrapStorage[i] = new double[BootstrapSampleSize]; }

                int getAvailableTaskIndex()
                {
                    for (int i = 0; i < tasks.Length; i++) { if (tasks[i] == null || tasks[i].IsCompleted) { return i; } }
                    return -1;
                }
                
                for (int i = 0; i < activeBranches.Length; i++)
                {
                    int taskIdx = getAvailableTaskIndex();
                    if (taskIdx > -1)
                    {
                        // Capture a copy of the current value of i; do not let the compiler use i in the closure of the lambda expression below
                        int idx = i;

                        tasks[taskIdx] = Task.Run(() => {
                            // Run the sampling
                            batches[idx].SampleNonAlloc(fitnessStorage[taskIdx], FitnessFunction, sampleSize);
                            // Compute the estimating distribution of the guiding parameter
                            switch (GuidingParam)
                            {
                                case GuidingParameter.Mean:
                                    {
                                        Normal dist = ParameterDistributions.MeanCLT(fitnessStorage[taskIdx]);
                                        parameterDistributions[idx] = dist;
                                        // Assign upper and lower bounds
                                        upperBounds[idx] = dist.Mean + NormalBoundStdDevs * dist.StdDev;
                                        lowerBounds[idx] = dist.Mean - NormalBoundStdDevs * dist.StdDev;
                                        break;
                                    }
                                case GuidingParameter.Median:
                                    {
                                        // Sort the fitness data in place
                                        Sorting.Sort(fitnessStorage[taskIdx]);
                                        Normal dist = ParameterDistributions.MedianBootstrapMemoryFriendly(fitnessStorage[taskIdx], bootstrapStorage[taskIdx], activeBranches[idx].rand);
                                        parameterDistributions[idx] = dist;
                                        // Assign upper and lower bounds
                                        upperBounds[idx] = dist.Mean + NormalBoundStdDevs * dist.StdDev;
                                        lowerBounds[idx] = dist.Mean - NormalBoundStdDevs * dist.StdDev;
                                        break;
                                    }
                                case GuidingParameter.LowerMean:
                                    {
                                        Normal dist = ParameterDistributions.MeanOfLessThanQuantile(fitnessStorage[taskIdx], 0.5);
                                        parameterDistributions[idx] = dist;
                                        // Assign upper and lower bounds
                                        upperBounds[idx] = dist.Mean + NormalBoundStdDevs * dist.StdDev;
                                        lowerBounds[idx] = dist.Mean - NormalBoundStdDevs * dist.StdDev;
                                        break;
                                    }
                                case GuidingParameter.OneOverNthQuantile:
                                    {
                                        GEV dist = ParameterDistributions.SampleMinimumMemoryFriendly(fitnessStorage[taskIdx], bootstrapStorage[taskIdx], activeBranches[idx].rand);
                                        parameterDistributions[idx] = dist;
                                        // Assign upper and lower bounds
                                        upperBounds[idx] = dist.Quantile(ComplementEpsilon16);
                                        lowerBounds[idx] = dist.Quantile(Epsilon16);
                                        break;
                                    }
                            }
                        });
                    }
                    else
                    {
                        Task.WaitAny(tasks);
                        i--;
                    }
                }
                for (int i = 0; i < tasks.Length; i++) { if (tasks[i] == null) tasks[i] = Task.Run(() => { return; }); } // Filler
                Task.WaitAll(tasks); // Once this has completed, we have sampled all of the branches.
                #endregion

                #region Discarding
                // --- Pre-emptive discarding ---
                bool performedPreEmptive = false;
                
                // Normal Pairwise Discarding
                if (parameterDistributions[0].GetType() == typeof(Normal))
                {
                    // Manual recast to Normal[]
                    Normal[] normals = new Normal[parameterDistributions.Length];
                    for (int i = 0; i < normals.Length; i++) { normals[i] = (Normal)parameterDistributions[i]; }

                    for (int i = 0; i < parameterDistributions.Length - 1; i++)
                    {
                        int discardIndex = DiscardProbabilityComputation.BestPairwiseDiscardNormal(normals, out double discardProb);
                        if (discardProb > PreemptiveDiscardConfidenceThreshold)
                        {
                            parameterDistributions[discardIndex] = null;
                            normals[discardIndex] = null;
                            activeBranches[discardIndex] = null;
                            performedPreEmptive = true;
                            continue;
                        }
                        break;
                    }
                }
                else // Non-normal case
                {
                    /*
                    double maxLowerBound = double.NegativeInfinity;
                    for (int i = 0; i < parameterDistributions.Length; i++) { maxLowerBound = Math.Max(maxLowerBound, lowerBounds[i]); }
                    // Discard any distribution with an upper bound less than max lower bound
                    for (int i = 0; i < parameterDistributions.Length; i++)
                    {
                        if (upperBounds[i] < maxLowerBound)
                        {
                            parameterDistributions[i] = null;
                            performedPreEmptive = true;
                        }
                    }*/

                    // Still in minimizing mode here, as we have not yet negated the distributions
                    double minUpperBound = double.PositiveInfinity;
                    for (int i = 0; i < parameterDistributions.Length; i++) { minUpperBound = Math.Min(minUpperBound, upperBounds[i]); }
                    // Discard any distribution with a lower bound greater than the smallest upper bound
                    for (int i = 0; i < parameterDistributions.Length; i++)
                    {
                        if (lowerBounds[i] > minUpperBound)
                        {
                            parameterDistributions[i] = null;
                            performedPreEmptive = true;
                        }
                    }
                }

                // Compact the arrays if necessary
                if (performedPreEmptive)
                {
                    var compactedUpperBounds = new List<double>(parameterDistributions.Length);
                    var compactedLowerBounds = new List<double>(parameterDistributions.Length);
                    var compactedDistributions = new List<IContinuousDistribution>(parameterDistributions.Length);
                    var compactBranches = new List<Branch>(activeBranches.Length);
                    for (int i = 0; i < activeBranches.Length; i++)
                    {
                        if (activeBranches[i] != null)
                        {
                            compactBranches.Add(activeBranches[i]);
                            compactedDistributions.Add(parameterDistributions[i]);
                            compactedUpperBounds.Add(upperBounds[i]);
                            compactedLowerBounds.Add(lowerBounds[i]);
                        }
                    }
                    activeBranches = compactBranches.ToArray();
                    parameterDistributions = compactedDistributions.ToArray();
                    upperBounds = compactedUpperBounds.ToArray();
                    lowerBounds = compactedLowerBounds.ToArray();
                }
                
                // --- Compute Discard Probabilities ---
                // Negate the distributions
                var negatedDists = DiscardProbabilityComputation.NegateDistributions(parameterDistributions, lowerBounds, upperBounds);

                double[] discardComplements = DiscardProbabilityComputation.ComplementsClenshawCurtisAutomatic(negatedDists);
                // The discard probabilities are in complement form (1 - P(D_i)) in this array

                // --- Discarding ---
                // Note: Nullifying a branch in activeBranches[] will discard it
                // Find the largest discard complement for reference
                double largestComplementDiscard = 0;
                for (int i = 0; i < discardComplements.Length; i++) { largestComplementDiscard = Math.Max(largestComplementDiscard, discardComplements[i]); }

                double confidence = 1.0;
                while (true)
                {
                    // Find the branch with the best probability to discard
                    double smallestComplementDiscard = 1;
                    int smallestCDIndex = 0;
                    for (int i = 0; i < discardComplements.Length; i++)
                    {
                        if (activeBranches[i]!= null && discardComplements[i] < smallestComplementDiscard)
                        {
                            smallestComplementDiscard = discardComplements[i];
                            smallestCDIndex = i;
                        }
                    }
                    // Try to discard the associated branch
                    if (smallestComplementDiscard < 0.5 * largestComplementDiscard // Relative size requirement
                        && confidence - smallestComplementDiscard >= confidenceLevel) // Maintain the confidence level
                    {
                        activeBranches[smallestCDIndex] = null; // Discard the branch
                        confidence -= smallestComplementDiscard; // Account for the cost in confidence
                    }
                    else break; // If we can't discard that one, then we can't do any better, so end the discard loop
                }

                #endregion
            }

            // Compact the result and return
            var compactedBranches = new List<Branch>(activeBranches);
            compactedBranches.RemoveAll(b => b == null); // Remove null
            activeBranches = compactedBranches.ToArray();

            // Assign the best observed
            BestElementObserved = (T)activeBranches[0].BestObservedSolution;
            BestFitnessObserved = activeBranches[0].MinimumObservedFitness;
            for (int i = 1; i < activeBranches.Length; i++)
            {
                if (activeBranches[i].MinimumObservedFitness < BestFitnessObserved)
                {
                    BestFitnessObserved = activeBranches[i].MinimumObservedFitness;
                    BestElementObserved = (T)activeBranches[i].BestObservedSolution;
                }
            }

            return activeBranches;
        }
        
        #endregion
    }
}
