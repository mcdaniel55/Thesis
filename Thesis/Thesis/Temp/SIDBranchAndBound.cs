using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using Thesis;

namespace Thesis.Optimization
{
    /// <summary> In the context of this work, an optimization problem has  
    /// * a solution space on which some branching structure can be imposed
    /// * a function to be minimized under which the branches have structure. </summary>
    public class SIDBranchAndBound
    {
        /// <summary> The starting region for the problem, containing all possible solutions. 
        /// Also provides the branching and evaluation methods for the problem. </summary>
        public Region SolutionSpace { get; private set; }

        public SIDBranchAndBound(Region SolutionSpace)
        {
            this.SolutionSpace = SolutionSpace;
        }

        /// <summary> Runs the algorithm in an attempt to narrow down the solution space to a set of regions 
        /// which can then be checked exhaustively. </summary>
        /// <param name="layers"> The number of iterations before terminating </param>
        /// <param name="confidenceLevel"> A number between 0 and 1 describing the confidence level to be used in statistical comparisons </param>
        /// <param name="Conservative"> Whether or not to always keep all regions for which it is not possible to rule out that they contain optimal values with 95% confidence. 
        /// Results can be higher quality with this on, but it may take significantly longer to run </param>
        /// <param name="CullDuplicates"> If true, duplicate regions will be removed after at the end of each layer. Requires the solution space to implement <see cref="IEquatable{T}"/> </param>
        /// <returns> An array of regions that are likely to contain optimal solutions </returns>
        public Region[] RunBranchAndBound(int sampleSize, int layers, double confidenceLevel, bool Conservative = true, bool CullDuplicates = false)
        {
            // Sanity check the initial sample size
            if (sampleSize < 30) { throw new ArgumentOutOfRangeException("Sample size must be at least 30."); }

            double branchingFactor = 0; // Keeps track of average branches per active region per layer, which may not be constant if there is overlap between regions
            double BestObservedValue = double.PositiveInfinity;
            Region[] activeRegions = new Region[1]{ SolutionSpace };
            List<Region> discardedRegions = new List<Region>(512);

            // Determine how many stdevs to use when computing the threshold for keeping otherwise discardable regions because they might contain minima
            double k = Math.Sqrt(1.0 / (1.0 - confidenceLevel) - 1.0); // Using Chebyshev's inequality, with a stardardization for tightening the bound

            // Main loop
            for (int layer = 0; layer < layers; layer++)
            {
                // --- Branch and Sample Active Regions ---
                // Branch the active regions
                var newActiveRegions = new List<Region>();
                foreach (Region reg in activeRegions) { if (reg != null) newActiveRegions.AddRange(reg.Branch()); }

                Program.logger.WriteLine($"Branched to produce {newActiveRegions.Count} new regions.");
                // Cull duplicates if desired
                if (CullDuplicates)
                {
                    for (int i = 1; i < newActiveRegions.Count; i++)
                    {
                        for (int j = 0; j < i; j++)
                        {
                            if (newActiveRegions[i].Equals(newActiveRegions[j]))
                            {
                                newActiveRegions.RemoveAt(i);
                                i--;
                                break;
                            }
                        }
                    }
                    Program.logger.WriteLine($"Duplicates culled. {newActiveRegions.Count} regions remain.");
                }

                // Update branching factor
                branchingFactor = (branchingFactor * layer + newActiveRegions.Count / activeRegions.Length) / (layer + 1);
                Program.logger.WriteLine($"Branching factor revised to {branchingFactor}");
                activeRegions = newActiveRegions.ToArray();

                // Set up random number generators for each task in a deterministic way
                Random[] rands = new Random[activeRegions.Length];
                for (int i = 0; i < rands.Length; i++) { rands[i] = new Xoshiro256StarStar(SolutionSpace.GetRNG().Next()); }

                // Sample the regions (now with 300% more threads!)
                Task[] tasks = new Task[activeRegions.Length];
                for (int i = 0; i < activeRegions.Length; i++) { int j = i; // To avoid issues with lambda closure, one of the few quirks of this language
                    tasks[j] = Task.Run(() => { activeRegions[j].Sample(sampleSize, rands[j]); }); }
                Task.WaitAll(tasks);
                //foreach (Region reg in activeRegions) { reg.Sample(sampleSize); }

                // Compute discard probabilities and number of discardable regions
                Normal[] samplingDistributions = new Normal[activeRegions.Length + discardedRegions.Count];
                for (int i = 0; i < activeRegions.Length; i++) { samplingDistributions[i] = activeRegions[i].SamplingDistribution; }
                for (int i = 0; i < discardedRegions.Count; i++) { samplingDistributions[i + activeRegions.Length] = discardedRegions[i].SamplingDistribution; }
                double[] discardComplements = NormalComparison.ComputeDiscardComplementsClenshawCurtisAltInvariantAutomatic(samplingDistributions);

                // If discard can likely be increased cheaply, increase sample size repeatedly until it can't
                int currentDiscardCount = EnumerateDiscardableRegions(discardComplements);
                int currentSampleSize = sampleSize;
                bool increasing = true;
                while (increasing)
                {
                    int newSampleSize = currentSampleSize + (int) (currentSampleSize * branchingFactor / activeRegions.Length);
                    // Construct a new set of normal distributions estimating the scenario where a larger sample size was used on each
                    var hypotheticals = new Normal[samplingDistributions.Length];
                    for (int i = 0; i < activeRegions.Length; i++) { hypotheticals[i] = activeRegions[i].EstimateDistributionWithDifferentSampleSize(newSampleSize); }
                    for (int i = 0; i < discardedRegions.Count; i++) { hypotheticals[activeRegions.Length + i] = discardedRegions[i].SamplingDistribution; }
                    double[] hypoDiscardProbs = NormalComparison.ComputeDiscardComplementsClenshawCurtisAltInvariantAutomatic(hypotheticals);
                    int discardable = EnumerateDiscardableRegions(hypoDiscardProbs);
                    if (discardable <= currentDiscardCount) { increasing = false; }
                    else {
                        currentSampleSize = newSampleSize;
                        Program.logger.WriteLine($"Revision expected to increase discard count from {currentDiscardCount} to {discardable}");
                        currentDiscardCount = discardable;
                    }
                }
                // If sample size was increased, resample active regions up to that size
                if (currentSampleSize > sampleSize)
                {
                    for (int i = 0; i < activeRegions.Length; i++)
                    {
                        activeRegions[i].Sample(currentSampleSize);
                        samplingDistributions[i] = activeRegions[i].SamplingDistribution;
                    }
                    discardComplements = NormalComparison.ComputeDiscardComplementsClenshawCurtisAltInvariantAutomatic(samplingDistributions);
                    Program.logger.WriteLine($"Layer sample size increased to {currentSampleSize}");
                }

                // Update the best observed value
                foreach (Region reg in activeRegions)
                {
                    BestObservedValue = Math.Min(reg.BestObservation, BestObservedValue);
                }
                Program.logger.WriteLine($"Best value observed in layer: {BestObservedValue}");

                // --- Discard Regions ---
                // Note: The sorting here is inefficient, but conceptually clearer than doing a quicksort and it isn't a performance bottleneck
                double certainty = 1;
                // Handle the discarded regions' contributions first
                for (int i = activeRegions.Length; i < discardComplements.Length; i++)
                {
                    certainty -= discardComplements[i];
                }
                Program.logger.WriteLine($"Uncertainty due to previously discarded regions: {1-certainty}");
                bool discarding = true;
                int discardCountTemp = 0;
                while (discarding)
                {
                    // Find the active region with the largest discard value
                    double min = 1;
                    int minIndex = 0;
                    for (int i = 0; i < activeRegions.Length; i++)
                    {
                        if (discardComplements[i] < min 
                            && activeRegions[i] != null
                            // Skip regions if they have too high of a probability of containing observations that compete with the current best
                            && ((activeRegions[i].SampleMean - k * activeRegions[i].SampleStdDev > BestObservedValue) || !Conservative))
                        {
                            min = discardComplements[i]; minIndex = i;
                        }
                    }
                    // Try to discard it
                    if (certainty - min > confidenceLevel)
                    {
                        certainty -= min;
                        if (min > 1E-18) { discardedRegions.Add(activeRegions[minIndex]); }
                        activeRegions[minIndex] = null;
                        discardCountTemp++;
                    }
                    else { discarding = false; }
                }

                // Discard completely uncompetetive regions with prejudice
                int discarded = 0;
                for (int i = activeRegions.Length; i < discardComplements.Length; i++)
                {
                    if (discardComplements[i] < 1E-13) { discardedRegions.RemoveAt(i - activeRegions.Length - discarded); discarded++; }
                }
                Program.logger.WriteLine($"Forgot {discarded} unimportant discarded regions. {discardedRegions.Count} discarded regions remain.");


                Program.logger.WriteLine($"Layer {layer} complete. Analyzed {activeRegions.Length} active regions, and discarded {discardCountTemp} active regions.");
                Program.logger.WriteLine("--- Active Regions --- ");
                for (int i = 0; i < activeRegions.Length; i++) { if (activeRegions[i] != null)
                    {
                        Program.logger.WriteLine($"Region {i}: {activeRegions[i].ToString()} \n");
                    } }
                Program.logger.WriteLine("--- Discarded Regions Still Affecting Probabilities --- ");
                for (int i = 0; i < discardedRegions.Count; i++)
                {
                    Program.logger.WriteLine($"Region {i}: {discardedRegions[i].ToString()} \n");
                }
            }

            // Shrinkwrap the output
            var output = new List<Region>(activeRegions);
            output.RemoveAll(r => (r == null));
            return output.ToArray();
        }

        private static int EnumerateDiscardableRegions(double[] discardProbabilities)
        {
            var discardProbList = new List<double>(discardProbabilities);
            discardProbList.Sort();

            double certainty = 1;
            for (int count = 0; count < discardProbList.Count; count++)
            {
                double newval = certainty - discardProbList[count];
                if (newval < 0.95) return count;
                certainty = newval;
            }
            throw new Exception("Nonsense Outcome Error: Could discard all regions!");
        }

    }
}
