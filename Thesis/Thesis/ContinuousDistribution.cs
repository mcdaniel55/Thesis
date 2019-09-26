using System;
using System.Data;
using System.Collections.Generic;

namespace Thesis.Quadrature
{
    /// <summary>
    /// An approximation of a continuous probability distribution using a piece-wise linear cumulative density function
    /// </summary>
    class ContinuousDistribution
    {
        /// <summary> An array of x-values at which the CDF's value is given explicitly, in increasing order </summary>
        public List<double> abscissas;
        /// <summary> An array of CDF values indexed alike to the array of abscissas </summary>
        public List<double> cumulativeDensities;

        public Random rand = Program.rand;

        protected ContinuousDistribution(){} // Heritable constructor to allow subclasses to easily define their own approximation behaviors

        /// <summary> Create a new distribution with the provided CDF abscissas and values. </summary>
        /// <param name="abscissas"> A list of x-values at which the CDF's value is given explicitly, in increasing order </param>
        /// <param name="cumulativeDensities"> A list of CDF values indexed alike to the array of abscissas </param>
        public ContinuousDistribution(List<double> abscissas, List<double> cumulativeDensities)
        {
            this.abscissas = abscissas;
            this.cumulativeDensities = cumulativeDensities;
            if (abscissas.Count != cumulativeDensities.Count) throw new Exception("CDF: The number of abscissas does not match the number of densities.");
            if (abscissas.Count < 2) throw new Exception($"CDF: Not enough abscissas: {abscissas.Count}");
        }
        
        /// <summary> Create a new distribution with the provided CDF abscissas and values. </summary>
        /// <param name="abscissas"> A list of x-values at which the CDF's value is given explicitly, in increasing order </param>
        /// <param name="cumulativeDensities"> A list of CDF values indexed alike to the array of abscissas </param>
        /// <param name="rand"> Allows the user to specify a random number generator for this distribution to use when sampling </param>
        public ContinuousDistribution(List<double> abscissas, List<double> cumulativeDensities, Random rand) : this(abscissas, cumulativeDensities)
        {
            this.rand = rand;
        }

        /// <summary>
        /// Returns the value of this distribution's CDF at the given location
        /// </summary>
        public double CumulativeDensity(double x)
        {
            // Edge cases
            if (x < abscissas[0]) return 0;
            if (x > abscissas[abscissas.Count - 1]) return 1;
            // Fast search for which elements to interpolate between
            int idx = abscissas.BinarySearch(x);
            // If x is an element of the abscissas, return the corresponding value
            if (idx > -1) return cumulativeDensities[idx];
            // Interpolate
            idx = ~idx; // idx is now the index of the next largest element of abscissas from x
            return Interpolation.Lerp(abscissas[idx - 1], cumulativeDensities[idx - 1], abscissas[idx], cumulativeDensities[idx], x);
        }

        /// <summary>
        /// Returns the value of the quantile function for this distribution at a given quantile
        /// </summary>
        public double Quantile(double q)
        {
            // Edge cases
            if (q < cumulativeDensities[0]) return abscissas[0];
            if (q > cumulativeDensities[cumulativeDensities.Count - 1]) return abscissas[abscissas.Count - 1];
            // Fast search for which elements to interpolate between
            int idx = cumulativeDensities.BinarySearch(q);
            // If q is an element of the cumulative densities list, return the corresponding abscissa
            if (idx > -1) return abscissas[idx];
            // Interpolate
            idx = ~idx; // idx is now the index of the next largest element of cumulative densities from q
            return Interpolation.Lerp(cumulativeDensities[idx - 1], abscissas[idx - 1], cumulativeDensities[idx], abscissas[idx], q);
        }

        /// <summary> Returns the state of the step function describing the density of this distribution at a given point </summary>
        /// <remarks> The left endpoint of a step is inclusive, and the right endpoint is not. This could be improved to approximate the density rather than just 
        /// describing the step function, but evaluating density is already an afterthought as we don't use it. </remarks>
        public double Density(double x)
        {
            if (x < abscissas[0] || x >= abscissas[abscissas.Count - 1]) return 0;
            int idx = abscissas.BinarySearch(x);
            if (idx > -1) return (cumulativeDensities[idx + 1] - cumulativeDensities[idx]) / (abscissas[idx + 1] - abscissas[idx]);
            idx = ~idx;
            return (cumulativeDensities[idx] - cumulativeDensities[idx - 1]) / (abscissas[idx] - abscissas[idx - 1]);
        }

        public double Sample() // This is an elegant property of continuous distributions
        {
            return Quantile(rand.NextDouble());
        }

        public double[] Sample(int size)
        {
            double[] sample = new double[size];
            for (int i = 0; i < sample.Length; i++) { sample[i] = Sample(); }
            return sample;
        }

        public void SampleNonAlloc(double[] arrayToFill)
        {
            for (int i = 0; i < arrayToFill.Length; i++) { arrayToFill[i] = Sample(); }
        }

        public ContinuousDistribution BootstrapSamplingDistribution(Func<double[],double> statistic, int iterations = 250, int smoothingPasses = 0, double smoothingCoefficient = 0.3/*, ExtrapolationMode mode = ExtrapolationMode.Linear*/)
        {
            int sampleSize = abscissas.Count;
            double[] sample = new double[sampleSize];
            double[] observations = new double[iterations];
            for (int i = 0; i < iterations; i++)
            {
                SampleNonAlloc(sample); // Fill the array with fresh sample points
                observations[i] = statistic(sample);
            }
            //return CDFApprox.FromSample(observations, smoothingPasses, smoothingCoefficient, mode, rand);
            return new ContinuousDistribution()
        }

        /// <summary> Writes the values of the provided functions at the provided input values in a comma separated array to the designated textwriter. </summary>
        /// <remarks> Use Console.Out as the first argument to write to the console with this. </remarks>
        public void PrintAtValues(System.IO.TextWriter writer, double[] values, params Func<double, double>[] properties)
        {
            // Header
            writer.Write(nameof(values));
            for (int i = 0; i < properties.Length; i++) { writer.Write($", {properties[i].Method.Name}"); }
            writer.WriteLine();

            // Body
            for (int i = 0; i < values.Length; i++)
            {
                writer.Write(values[i]);
                foreach (Func<double, double> prop in properties) { writer.Write($", {prop(values[i])}"); }
                writer.WriteLine();
            }
        }

        public void PrintAtAbscissas(System.IO.TextWriter writer, params Func<double, double>[] properties)
        {
            // Header 
            writer.Write("Abscissa");
            for (int i = 0; i < properties.Length; i++) { writer.Write($", {properties[i].Method.Name}"); }
            writer.WriteLine();
            // Body
            for (int i = 0; i < abscissas.Count; i++)
            {
                writer.Write(abscissas[i]);
                foreach (Func<double, double> prop in properties) { writer.Write($", {prop(abscissas[i])}"); }
                writer.WriteLine();
            }
        }

        public DataTable TabulateAtValues(string name, double[] inputs, params Func<double,double>[] properties)
        {
            if (name == null) name = this.ToString();
            DataTable table = new DataTable(name);
            table.Columns.Add(nameof(inputs), typeof(double));
            foreach (var function in properties) { table.Columns.Add(function.Method.Name, typeof(double)); }
            for (int i = 0; i < inputs.Length; i++)
            {
                DataRow row = table.NewRow();
                row[nameof(inputs)] = inputs[i];
                foreach (var function in properties) { row[function.Method.Name] = function(inputs[i]); }
            }
            return table;
        }

        public DataTable MakeTableUsingAbscissas(string name, params Func<double, double>[] properties)
        {
            if (name == null) name = this.ToString();
            DataTable table = new DataTable(name);
            table.Columns.Add("Abscissas", typeof(double));
            foreach (var function in properties) { table.Columns.Add(function.Method.Name, typeof(double)); }
            for (int i = 0; i < abscissas.Count; i++)
            {
                DataRow row = table.NewRow();
                row[0] = abscissas[i];
                for (int j = 0; j < properties.Length; j++) { row[j+1] = properties[j](abscissas[i]); }
                table.Rows.Add(row); // This is necessary, as for some reason table.NewRow() only creates a row from the table's schema, and doesn't automatically add it to the table
            }
            return table;
        }
    }
}
