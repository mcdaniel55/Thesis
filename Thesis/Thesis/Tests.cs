using System;
using System.Collections.Generic;
using System.Data;
using MathNet.Numerics.Distributions;

namespace Thesis
{
    static class Tests
    {
        
        // --- Temporary Tests ---

        public static void TestCCQuadrature1()
        {
            double f(double x) => Math.Sqrt(x);

            double outcome = Quadrature.ClenshawCurtis.Integrate(f, 1, 4, 7);
            Console.WriteLine($"Outcome: {outcome}");
        }


        public static void GenerateHardDistributionSet()
        {
            for (int i = 0; i < 30; i++)
            {
                Program.logger.WriteLine("");
                double f() => 90 * Program.rand.NextDouble() + 30; // Means
                double g() => Math.Abs(Normal.Sample(20, 15)) + 0.3; // StdDevs
                for (int j = 0; j < 4; j++)
                {
                    Program.logger.Write($"new Normal({f()}, {g()}), ");
                }
            }
        }


        public static void TestPickands()
        {
            Logger output = new Logger("Pickands Comparison.csv");

            // Generate some data
            ChiSquared chisq = new ChiSquared(4, Program.rand);
            double[] raw = new double[1000];
            chisq.Samples(raw);
            List<double> data = new List<double>(raw);
            data.Sort();

            // Apply the new and old Pickands code and write them to the output for comparison
            ContinuousDistribution oldVersion = PickandsApproximation.ApproximatePiecewiseDistributionWithUpperTail(data, 1000);
            PickandsApproximation newVersion = new PickandsApproximation(data);

            raw = new double[10000];
            // Generate new data from the approximation
            newVersion.Samples(raw, Program.rand);
            var resample = ContinuousDistribution.ECDF(raw, Program.rand);

            // Make the table and write it to a file
            DataTable table = new DataTable("Comparison of Pickands Implementations");
            table.Columns.Add("Abscissas", typeof(double));
            table.Columns.Add("OldCDF", typeof(double));
            table.Columns.Add("NewCDF", typeof(double));
            table.Columns.Add("ResampledECDF", typeof(double));
            table.Columns.Add("TrueCDF", typeof(double));

            for (int i = 0; i < oldVersion.abscissas.Count; i++)
            {
                DataRow row = table.NewRow();
                row["Abscissas"] = oldVersion.abscissas[i];
                row["OldCDF"] = oldVersion.cumulativeDensities[i];
                row["NewCDF"] = newVersion.CDF(oldVersion.abscissas[i]);
                row["ResampledECDF"] = resample.CumulativeDensity(oldVersion.abscissas[i]);
                row["TrueCDF"] = chisq.CumulativeDistribution(oldVersion.abscissas[i]);
                table.Rows.Add(row);
            }

            output.WriteTable(table);

            // Clean up
            output.Dispose();
            table.Dispose();

            /*
            Console.WriteLine($"CDF(1):{newVersion.CDF(1)} QF(that): {newVersion.Quantile(newVersion.CDF(1))}");
            Console.WriteLine($"CDF(2):{newVersion.CDF(2)} QF(that): {newVersion.Quantile(newVersion.CDF(2))}");
            Console.WriteLine($"CDF(3):{newVersion.CDF(3)} QF(that): {newVersion.Quantile(newVersion.CDF(3))}");
            Console.WriteLine($"CDF(4):{newVersion.CDF(4)} QF(that): {newVersion.Quantile(newVersion.CDF(4))}");
            Console.WriteLine($"CDF(11):{newVersion.CDF(11)} QF(that): {newVersion.Quantile(newVersion.CDF(11))}");
            Console.WriteLine($"CDF(13):{newVersion.CDF(13)} QF(that): {newVersion.Quantile(newVersion.CDF(13))}");
            Console.WriteLine($"CDF(15):{newVersion.CDF(15)} QF(that): {newVersion.Quantile(newVersion.CDF(15))}");
            double cdfval = PickandsBalkemaDeHaan.TailCDF(0.4, newVersion.a, newVersion.c);
            Console.WriteLine($"TailCDF(0.4):{cdfval} QF(that): {PickandsBalkemaDeHaan.TailQuantileFunction(cdfval, newVersion.a, newVersion.c)}");
            */
        }

    }
}
