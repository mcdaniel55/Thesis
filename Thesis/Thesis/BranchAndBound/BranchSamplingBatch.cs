using System;

namespace Thesis.BranchAndBound
{
    /// <summary>
    /// Encapsulates the process of sampling a branch 
    /// </summary>
    /// <typeparam name="T"> The type of object that is an element of the branch </typeparam>
    class BranchSamplingBatch<T>
    {
        public readonly Random rand;
        public Branch Branch { get; private set; }
        public int SampleSize { get; private set; }
        public T BestObservation { get; private set; }
        public double BestObservedFitness { get; private set; }

        public BranchSamplingBatch(Branch Branch, Random rand = null)
        {
            this.Branch = Branch;
            this.rand = rand ?? Program.rand;
        }


        /// <summary> Fills the provided array with fitness values evaluated from elements sampled from the Branch </summary>
        /// <param name="output"></param>
        /// <param name="sampleSize"></param>
        /// <param name="FitnessFunction"></param>
        public void SampleNonAlloc(double[] output, Func<T,double> FitnessFunction, int SampleSize)
        {
            T input;
            Branch.rand = rand; // Point the branch's RNG to the desired one
            for (int i = 0; i < output.Length; i++)
            {
                input = (T)Branch.GetRandomElement();
                output[i] = FitnessFunction(input);
                // Update best observation
                if (output[i] < BestObservedFitness)
                {
                    BestObservedFitness = output[i];
                    BestObservation = input;
                }
            }
            
            this.SampleSize = SampleSize;
        }
    }
}
