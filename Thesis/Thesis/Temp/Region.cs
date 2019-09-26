using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Distributions;

namespace ThesisOptNumericalTest.Optimization
{
    /// <summary> Encapsulates a sampling region. Derived classes must specify how to branch and evaluate elements according to the specific problem </summary>
    public abstract class Region
    {
        #region Properties
        /// <summary> Keeps track of how many times the solution space branched to produce this region </summary>
        public int BranchLevel { get; private set; }
        public int SampleSize { get; private set; }
        public Normal SamplingDistribution { get; private set; }
        public double SampleMean { get; private set; }
        public double SampleStdDev { get; private set; }
        /// <summary> The best (smallest) observation of the fitness function </summary>
        public double BestObservation { get; private set; }

        protected readonly Random m_rand;

        // Used for cheaply updating the mean and variance when the sample size is increased
        private double sumOfSquaredDeviations = 0;
        private double sumOfValues = 0;
        #endregion

        #region Methods
        /// <summary> Obtains a single value at random from this region </summary>
        protected abstract double SampleElement();

        /// <summary> Creates a cover of this region according to problem-specific structure </summary>
        public abstract Region[] Branch();

        /// <summary> A constructor that specifies some common initialization logic and 
        /// provides an explicit parameterization of all subclass constructors </summary>
        protected Region(Random rand)
        {
            m_rand = rand;
            BestObservation = double.PositiveInfinity;
            SampleSize = 0;
        }

        public Random GetRNG()
        {
            return m_rand;
        }

        /// <summary> Takes a sample from the lower half of the region and updates its sampling distribution. Will efficiently increase sample size if re-run with a higher value. </summary>
        /// <param name="totalSize"> The desired total sample size. </param>
        /// <param name="rand"> Override for random number generator for use in threaded sampling </param>
        public void Sample(int totalSize, Random rand = null)
        {
            int additionalSamples = totalSize - SampleSize;
            if (additionalSamples < 1) { return; }

            if (rand == null) rand = m_rand; // Use default rand

            // Figure out how many elements to sample to guarantee the lowest 'additionalSamples' many 
            // elements in the new sample are in the lower half of the population
            int oversampleSize = GetOverSampleSize(additionalSamples);

            var samples = new List<double>(oversampleSize);
            for (int i = 0; i < oversampleSize; i++)
            {
                samples.Add(SampleElement());
            }
            samples.Sort();
            samples.RemoveRange(samples.Count / 2 + 1, samples.Count - samples.Count / 2 - 1);

            // Compute the new mean
            foreach (double val in samples) { sumOfValues += val; }
            SampleMean = sumOfValues / totalSize;

            // Compute the new StdDev, using the old sum of squared deviations to avoid resampling
            foreach (double val in samples) { sumOfSquaredDeviations += Math.Pow(val - SampleMean, 2); }
            SampleStdDev = Math.Sqrt(sumOfSquaredDeviations / (totalSize - 1));

            SamplingDistribution = new Normal(SampleMean, SampleStdDev / Math.Sqrt(totalSize), m_rand);

            // Update the best observed value
            for (int i = 0; i < samples.Count; i++)
            {
                //if (samples[i] < BestObservation) { BestObservation = samples[i]; }
                BestObservation = Math.Min(BestObservation, samples[i]);
            }

            SampleSize = totalSize;
        }

        #region Original Version of Sample()
        /// <summary> Takes a sample from the region and updates its sampling distribution. Will efficiently increase sample size if re-run with a higher value. </summary>
        /// <param name="totalSize"> The desired total sample size. </param>
        /*
        public void Sample(int totalSize)
        {
            int additionalSamples = totalSize - SampleSize;
            if (additionalSamples < 1) { return; }

            double[] samples = new double[additionalSamples];
            for (int i = 0; i < additionalSamples; i++)
            {
                samples[i] = SampleElement();
                sumOfValues += samples[i];
            }

            // Compute the new mean
            SampleMean = sumOfValues / totalSize;

            // Compute the new StdDev, using the old sum of squared deviations to avoid resampling
            for (int i = 0; i < additionalSamples; i++) { sumOfSquaredDeviations += Math.Pow(samples[i] - SampleMean, 2); }
            SampleStdDev = Math.Sqrt(sumOfSquaredDeviations / (totalSize - 1));

            SamplingDistribution = new Normal(SampleMean, SampleStdDev / Math.Sqrt(totalSize), m_rand);

            // Update the best observed value
            for (int i = 0; i < samples.Length; i++)
            {
                //if (samples[i] < BestObservation) { BestObservation = samples[i]; }
                BestObservation = Math.Min(BestObservation, samples[i]);
            }

            SampleSize = totalSize;
        }
        */
        #endregion

        /// <summary> Estimates what the sampling distribution of this region would look like if a different sample size were used </summary>
        /// <param name="newSize"> The new sample size for estimation </param>
        /// <returns> A normal distribution describing the estimated sampling distribution </returns>
        public Normal EstimateDistributionWithDifferentSampleSize(int newSize)
        {
            return new Normal(SamplingDistribution.Mean, Math.Sqrt(sumOfSquaredDeviations / newSize * (newSize - 1)), m_rand);
        }
        #endregion

        /// <summary>
        /// For a given desired sample size N, computes the number of elements that need to be in a sample to ensure there is 95% confidence that 
        /// the lowest N values in the sample are less than the population median
        /// </summary>
        /// <remarks> 
        /// Characterizing an element of the sample as a success if it is less than the median, the number of elements in a sample less than the population
        /// median follows a binomial distribution with p approximately equal to 0.5. Since we always have at least 30 samples, and this p value is well-suited
        /// to such approximations, we can use the normal approximation of the binomial distribution to cheaply compute a one-sided 95% bound. 
        /// </remarks>
        private static int GetOverSampleSize(int targetSize)
        {
            return (int) Math.Ceiling(Math.Pow((1.645 + Math.Sqrt(1.645 * 1.645 + 8 * targetSize)) / 2, 2));
        }
    }
}
