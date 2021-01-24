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
        private const int MONTE_CARLO_SAMPLE_SIZE = 250;
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
                Console.WriteLine($"On iteration {iteration + 1} of {iterations}");
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
                var negatedDistributions = new IDistributionWrapper[activeBranches.Length];
                #endregion

                #region Sampling
                // --- Sample the regions (now with more threads!) ---
                var batches = new BranchSamplingBatch<T>[activeBranches.Length];
                // This provides deterministic RNG without race conditions for multithreaded execution at a modest memory cost
                for (int i = 0; i < batches.Length; i++) { batches[i] = new BranchSamplingBatch<T>(activeBranches[i], new Xoshiro256StarStar(Program.rand.Next())); }

                void GetParameterDistribution(int index)
                {
                    // Inefficient allocation, but not a big deal
                    double[] fitnessStorage = new double[sampleSize];
                    double[] bootstrapStorage = new double[MONTE_CARLO_SAMPLE_SIZE];
                    
                    // Run the sampling
                    batches[index].SampleNonAlloc(fitnessStorage, FitnessFunction);

                    // Compute the estimating distribution of the guiding parameter
                    switch (GuidingParam)
                    {
                        case GuidingParameter.Mean:
                            {
                                Normal dist = ParameterDistributions.MeanCLT(fitnessStorage);
                                negatedDistributions[index] = new NegatedDistribution(dist, dist.Mean - NormalBoundStdDevs * dist.StdDev, dist.Mean + NormalBoundStdDevs * dist.StdDev);
                                break;
                            }
                        case GuidingParameter.Median:
                            {
                                // Sort the fitness data in place
                                Sorting.Sort(fitnessStorage);
                                Normal dist = ParameterDistributions.MedianBootstrapMemoryFriendly(fitnessStorage, bootstrapStorage, activeBranches[index].rand);
                                negatedDistributions[index] = new NegatedDistribution(dist, dist.Mean - NormalBoundStdDevs * dist.StdDev, dist.Mean + NormalBoundStdDevs * dist.StdDev);
                                break;
                            }
                        case GuidingParameter.LowerMean:
                            {
                                Normal dist = ParameterDistributions.MeanOfLessThanQuantile(fitnessStorage, 0.1);
                                negatedDistributions[index] = new NegatedDistribution(dist, dist.Mean - NormalBoundStdDevs * dist.StdDev, dist.Mean + NormalBoundStdDevs * dist.StdDev);
                                break;
                            }
                        case GuidingParameter.OneOverNthQuantile:
                            {

                                for (int j = 0; j < sampleSize; j++)
                                {
                                    fitnessStorage[j] *= -1.0;
                                }
                                Sorting.Sort(fitnessStorage);
                                negatedDistributions[index] = ParameterDistributions.OneOverNthQuantileViaSampleMinimumParameterDistribution(fitnessStorage, bootstrapStorage, activeBranches[index].rand);
                                break;
                            }
                    }

                }

                if (multiThread)
                {
                    //int threadCount = Environment.ProcessorCount; // # of logical processors, including hyperthreading etc.
                    int threadCount = 6; // Temp limit
                    Parallel.For(0, batches.Length, new ParallelOptions() { MaxDegreeOfParallelism = threadCount }, (int i) => { GetParameterDistribution(i); } );
                }
                else // Synchronous case
                {
                    var fitnessStorage = new double[sampleSize];
                    var monteCarloStorage = new double[MONTE_CARLO_SAMPLE_SIZE];
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
                        negatedDistributions[i] = ParameterDistributions.OneOverNthQuantileViaSampleMinimumParameterDistribution(fitnessStorage, monteCarloStorage, activeBranches[i].rand);
                    }
#if DEBUG
                    Program.logger.WriteLine("Break SID 172");
#endif
                }
                
                #endregion

                // Update the best observation so far
                for (int i = 0; i < batches.Length; i++)
                {
                    if (batches[i].BestObservedFitness < BestFitnessObserved)
                    {
                        BestFitnessObserved = batches[i].BestObservedFitness;
                        BestElementObserved = batches[i].BestObservation;
                    }
                }

                #region Discarding
                // --- Pre-emptive discarding ---
                
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
                    //discardComplements = DiscardProbabilityComputation.ComplementsMonteCarloMaximizing(negatedDistributions);
                    //discardComplements = DiscardProbabilityComputation.ComplementsQuantileTrapRule(negatedDistributions);
                    discardComplements = DiscardProbabilityComputation.ComplementsClenshawCurtisAutomatic(negatedDistributions);
                }
                else
                {
                    discardComplements = DiscardProbabilityComputation.ComplementsClenshawCurtisAutomatic(negatedDistributions, errorTolerance: 1E-8, maxIterations: 10);
                    //discardComplements = DiscardProbabilityComputation.ComplementsMonteCarloMaximizing(negatedDistributions);
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

                // Print the list to the log
                Program.logger.WriteLine($"Completed iteration {iteration}. Active branches:");
                for (int i = 0; i < activeBranches.Length; i++)
                {
                    if(activeBranches[i] != null) Program.logger.WriteLine($"{activeBranches[i].ToString()}");
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
