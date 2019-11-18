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
        private static readonly double Epsilon12 = Math.Pow(2, -40);
        private static readonly double ComplementEpsilon12 = 1.0 - Math.Pow(2, -40);
        private static readonly double NormalBoundStdDevs = Normal.InvCDF(0, 1, 1 - Math.Pow(2, -48)); // Roughly 8 standard deviations
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
        
        public Branch[] BranchAndBound(int sampleSize, int iterations, double confidenceLevel = 0.95, bool cullDuplicates = true, bool multiThread = true)
        {
            // Sanity check the initial sample size
            if (sampleSize < 30) { throw new ArgumentOutOfRangeException("Sample size must be at least 30."); }

            BestFitnessObserved = double.PositiveInfinity;
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

                // Print the list to the log
                for (int i = 0; i < activeBranches.Length; i++)
                {
                    Program.logger.WriteLine($"{activeBranches[i].ToString()}");
                }

                // Storage for the estimating distribution of the guiding parameter of each branch
                //var parameterDistributions = new IContinuousDistribution[activeBranches.Length];
                var negatedDistributions = new IDistributionWrapper[activeBranches.Length];
                //var upperBounds = new double[parameterDistributions.Length];
                //var lowerBounds = new double[parameterDistributions.Length];
                #endregion

                #region Sampling
                // --- Sample the regions (now with more threads!) ---
                var batches = new BranchSamplingBatch<T>[activeBranches.Length];
                // This provides deterministic RNG without race conditions for multithreaded execution at a modest memory cost
                for (int i = 0; i < batches.Length; i++) { batches[i] = new BranchSamplingBatch<T>(activeBranches[i], new Xoshiro256StarStar(Program.rand.Next())); }

                if (multiThread)
                {
                    int threadCount = Environment.ProcessorCount; // # of logical processors, including hyperthreading etc.
                    
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
                                batches[idx].SampleNonAlloc(fitnessStorage[taskIdx], FitnessFunction);

                                // Compute the estimating distribution of the guiding parameter
                                switch (GuidingParam)
                                {
                                    case GuidingParameter.Mean:
                                        {
                                            Normal dist = ParameterDistributions.MeanCLT(fitnessStorage[taskIdx]);
                                            negatedDistributions[idx] = new NegatedDistribution(dist, dist.Mean - NormalBoundStdDevs * dist.StdDev, dist.Mean + NormalBoundStdDevs * dist.StdDev);
                                            //parameterDistributions[idx] = dist;
                                            // Assign upper and lower bounds
                                            //upperBounds[idx] = dist.Mean + NormalBoundStdDevs * dist.StdDev;
                                            //lowerBounds[idx] = dist.Mean - NormalBoundStdDevs * dist.StdDev;
                                            break;
                                        }
                                    case GuidingParameter.Median:
                                        {
                                            // Sort the fitness data in place
                                            Sorting.Sort(fitnessStorage[taskIdx]);
                                            Normal dist = ParameterDistributions.MedianBootstrapMemoryFriendly(fitnessStorage[taskIdx], bootstrapStorage[taskIdx], activeBranches[idx].rand);
                                            negatedDistributions[idx] = new NegatedDistribution(dist, dist.Mean - NormalBoundStdDevs * dist.StdDev, dist.Mean + NormalBoundStdDevs * dist.StdDev);
                                            //parameterDistributions[idx] = dist;
                                            // Assign upper and lower bounds
                                            //upperBounds[idx] = dist.Mean + NormalBoundStdDevs * dist.StdDev;
                                            //lowerBounds[idx] = dist.Mean - NormalBoundStdDevs * dist.StdDev;
                                            break;
                                        }
                                    case GuidingParameter.LowerMean:
                                        {
                                            Normal dist = ParameterDistributions.MeanOfLessThanQuantile(fitnessStorage[taskIdx], 0.5);
                                            negatedDistributions[idx] = new NegatedDistribution(dist, dist.Mean - NormalBoundStdDevs * dist.StdDev, dist.Mean + NormalBoundStdDevs * dist.StdDev);
                                            //parameterDistributions[idx] = dist;
                                            // Assign upper and lower bounds
                                            //upperBounds[idx] = dist.Mean + NormalBoundStdDevs * dist.StdDev;
                                            //lowerBounds[idx] = dist.Mean - NormalBoundStdDevs * dist.StdDev;
                                            break;
                                        }
                                    case GuidingParameter.OneOverNthQuantile:
                                        {
                                            // Negate the sample first, so the new upper tail is the old lower tail
                                            /* Deprecated: In-place negation for sorted data. The sorting is handled by Pickands, so it doesn't matter.
                                            double temp;
                                            // Swap elements and negate to preserve ordering
                                            for (int j = 0; j < fitnessStorage[taskIdx].Length / 2; j++)
                                            {
                                                temp = fitnessStorage[taskIdx][j];
                                                fitnessStorage[taskIdx][j] = -fitnessStorage[taskIdx][sampleSize - j - 1];
                                                fitnessStorage[taskIdx][sampleSize - j - 1] = -temp;
                                            }
                                            // Handle middle element if array length is odd
                                            if (sampleSize % 2 == 1) 
                                            {
                                                fitnessStorage[taskIdx][sampleSize / 2] *= -1.0;
                                            }*/
                                            for (int j = 0; j < sampleSize; j++)
                                            {
                                                fitnessStorage[taskIdx][j] *= -1.0;
                                            }

                                            negatedDistributions[idx] = ParameterDistributions.NMinusOneOverNthQuantileViaSampleMinimumParameterDistribution(fitnessStorage[taskIdx], bootstrapStorage[taskIdx], activeBranches[idx].rand);
                                            //GEV dist = ParameterDistributions.SampleMinimumErrorDistMemoryFriendly(fitnessStorage[taskIdx], bootstrapStorage[taskIdx], activeBranches[idx].rand);
                                            //parameterDistributions[idx] = dist;
                                            // Assign upper and lower bounds
                                            //upperBounds[idx] = dist.Quantile(ComplementEpsilon12);
                                            //lowerBounds[idx] = dist.Quantile(Epsilon12);
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
                }
                else // Synchronous case
                {
                    var fitnessStorage = new double[sampleSize];
                    var bootstrapStorage = new double[BootstrapSampleSize];
                    for (int i = 0; i < batches.Length; i++)
                    {
                        batches[i].SampleNonAlloc(fitnessStorage, FitnessFunction);
                        // Assuming 1/nth quantile for testing
                        // Negate the observations to get observations from -X
                        for (int j = 0; j < sampleSize; j++)
                        {
                            fitnessStorage[j] *= -1.0;
                        }

                        Sorting.Sort(fitnessStorage);
                        negatedDistributions[i] = ParameterDistributions.NMinusOneOverNthQuantileViaSampleMinimumParameterDistribution(fitnessStorage, bootstrapStorage, activeBranches[i].rand);
                        //parameterDistributions[i] = dist;
                        // Assign upper and lower bounds
                        //upperBounds[i] = dist.Quantile(ComplementEpsilon12);
                        //lowerBounds[i] = dist.Quantile(Epsilon12);
                    }
                }
                
                #endregion

                // Update the best observation so far (this is not yet implemented at a lower level)
                for (int i = 0; i < batches.Length; i++)
                {
                    if (batches[i].BestObservedFitness < BestFitnessObserved)
                    {
                        BestFitnessObserved = batches[i].BestObservedFitness;
                        BestElementObserved = batches[i].BestObservation;
                    }
                }

                #region Discarding
                // The goal here is to do pre-emptive discarding and put the remaining (negated) sampling distributions in this list
                //var negatedDistributions = new IDistributionWrapper[] { };
                // --- Pre-emptive discarding ---
                //bool performedPreEmptive = false;
                
                // Normal Pairwise Discarding
                if (negatedDistributions[0].GetWrappedDistribution().GetType() == typeof(Normal))
                {
                    // Manual recast to a non-negated Normal[]
                    Normal[] normals = new Normal[negatedDistributions.Length];
                    for (int i = 0; i < normals.Length; i++) { normals[i] = (Normal)(negatedDistributions[i].GetWrappedDistribution()); }

                    for (int i = 0; i < negatedDistributions.Length - 1; i++)
                    {
                        int discardIndex = DiscardProbabilityComputation.BestPairwiseDiscardNormal(normals, out double discardProb);
                        if (discardProb > PreemptiveDiscardConfidenceThreshold)
                        {
                            negatedDistributions[discardIndex] = null;
                            normals[discardIndex] = null;
                            activeBranches[discardIndex] = null;
                            //performedPreEmptive = true;
                            continue;
                        }
                        break;
                    }
                }
                else // Non-normal case
                {
                    // We have negated the distributions at this point, so we need to compare upper bounds against the largest lower bound
                    double maxLowerBound = double.NegativeInfinity;
                    for (int i = 0; i < negatedDistributions.Length; i++) { maxLowerBound = Math.Max(maxLowerBound, negatedDistributions[i].GetLowerBound()); }
                    // Discard any distribution with an upper bound less than max lower bound
                    for (int i = 0; i < negatedDistributions.Length; i++)
                    {
                        if (negatedDistributions[i].GetUpperBound() < maxLowerBound)
                        {
                            negatedDistributions[i] = null;
                            activeBranches[i] = null;
                            //performedPreEmptive = true;
                        }
                    }
                }

                // Compact the arrays of active branches and their corresponding negated distributions, keeping them paired by their indices
                var compactBranches = new List<Branch>(activeBranches.Length);
                var compactNegatedDists = new List<IDistributionWrapper>(activeBranches.Length);
                for (int i = 0; i < activeBranches.Length; i++)
                {
                    if (activeBranches[i] != null)
                    {
                        compactBranches.Add(activeBranches[i]);
                        compactNegatedDists.Add(negatedDistributions[i]);
                    }
                }
                activeBranches = compactBranches.ToArray();
                negatedDistributions = compactNegatedDists.ToArray();

                // --- Compute Discard Probabilities ---
                double[] discardComplements; // The discard probabilities are in complement form (1 - P(D_i)) in this array
                if (GuidingParam == GuidingParameter.OneOverNthQuantile)
                {
                    // If the guiding parameter is the 1/nth quantile, then the integral will become too hard to compute. 
                    // Use a foolproof alternative for now, until we can find a better way to compute it
                    discardComplements = DiscardProbabilityComputation.ComplementsMonteCarloMaximizing(negatedDistributions);
                }
                else
                {
                    discardComplements = DiscardProbabilityComputation.ComplementsClenshawCurtisAutomatic(negatedDistributions, errorTolerance: 1E-8, maxIterations: 10);
                }
                
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

            return activeBranches;
        }
        
        #endregion
    }
}
