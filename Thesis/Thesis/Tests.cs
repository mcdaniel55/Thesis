using System;
using System.Collections.Generic;
using System.Data;
using MathNet.Numerics.Distributions;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

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

        public static void TestGEV()
        {
            Logger output = new Logger("GEV Test A.csv");
            Logger output2 = new Logger("GEV Test B.csv");
            //var dist = new ChiSquared(4, Program.rand);
            var dist = new Beta(2, 2);
            //var dist = new Exponential(2, Program.rand);
            //var dist = new Gamma(2, 2, Program.rand);
            const int sampleSize = 300;
            
            // Report the distribution 1-1/e quantile
            double upperQuantile = dist.InverseCumulativeDistribution(1 - 1.0 / sampleSize);
            double lowerQuantile = dist.InverseCumulativeDistribution(1.0 / sampleSize);
            output.WriteLine($"1-1/samplesize quantile: {upperQuantile}");
            output.WriteLine($"1/samplesize quantile: {lowerQuantile}");

            // Monte Carlo for the true distribution of the sample maximum
            double[] observations = new double[10000];
            
            for (int i = 0; i < observations.Length; i++)
            {
                double max = double.NegativeInfinity;
                for (int j = 0; j < sampleSize; j++)
                {
                    max = Math.Max(max, dist.Sample());
                }
                observations[i] = max;
            }
            Sorting.Sort(observations);

            ContinuousDistribution MonteCarloDistributionOfTheMaximum = ContinuousDistribution.ECDF(observations, Program.rand);

            // --- Find the best fit GEV distribution for this dataset ---

            #region Old code
            /*
            // Compute location and scale parameter estimates for a given shape parameter Xi using the median and variance
            void EstimateParameters(double shape, double median, double variance, out double location, out double scale)
            {
                if (shape == 0)
                {
                    scale = Math.Sqrt(6 * variance) / Math.PI;
                    location = median + scale * Math.Log(Math.Log(2));
                    return;
                }
                // This scale may or may not work for Xi > 0.5
                scale = Math.Sign(shape) * shape * Math.Sqrt(variance) / Math.Sqrt(SpecialFunctions.Gamma(1 - 2 * shape) - SpecialFunctions.Gamma(1 - shape) * SpecialFunctions.Gamma(1 - shape));
                if (double.IsNaN(scale)) scale = Math.Sqrt(6 * variance) / Math.PI;
                location = median - scale * (Math.Pow(Math.Log(2), -shape) - 1) / shape;
            }*/
            #endregion

            double FitnessExactModel(GEV model)
            {
                double val = 0;
                for (int i = 0; i < observations.Length; i++)
                {
                    val += Math.Pow(model.CumulativeDistribution(observations[i]) - MonteCarloDistributionOfTheMaximum.CumulativeDensity(observations[i]), 2);
                }
                return val;
            }

            #region Old code
            //double medianEst = Statistics.Median(observations);
            //double varianceEst = Statistics.VarianceEstimate(observations);

            /*
            GEV Optimize(double startingval, out double fitness)
            {
                double locationEst;
                double scaleEst;
                double bestScore = double.PositiveInfinity;
                GEV bestSoFar = null;
                bool increasing = false;
                int sinceImproved = 0;
                double shapeEst = startingval; // Neg or pos will stay that way throughout the optimization

                while (true)
                {
                    EstimateParameters(shapeEst, medianEst, varianceEst, out locationEst, out scaleEst);
                    GEV model = new GEV(locationEst, scaleEst, shapeEst, Program.rand);
                    double score = FitnessExactModel(model);
                    if (score < bestScore)
                    {
                        bestScore = score;
                        bestSoFar = model;
                        sinceImproved = 0;
                    }
                    else
                    {
                        increasing ^= true;
                        if (++sinceImproved > 10) break;
                    }
                    if (increasing) shapeEst += 0.3 * startingval;
                    else shapeEst *= 0.5;
                }
                fitness = bestScore;
                return bestSoFar;
            }

            GEV OptimizeV2(double initialGuess, out double fitness)
            {
                double locationEst, scaleEst;
                double bestScore = double.PositiveInfinity;
                GEV bestSoFar = null;
                double shapeEst = initialGuess;
                double bestShapeSoFar = initialGuess;
                // Grow the estimate by doubling until it is no longer improving
                while (true)
                {
                    EstimateParameters(shapeEst, medianEst, varianceEst, out locationEst, out scaleEst);
                    GEV model = new GEV(locationEst, scaleEst, shapeEst, Program.rand);
                    double score = FitnessExactModel(model);
                    if (score < bestScore) // If it improved
                    {
                        bestScore = score;
                        bestSoFar = model;
                        bestShapeSoFar = shapeEst;
                    }
                    else break;
                    shapeEst *= 2;
                }
                double magnitude = bestShapeSoFar;
                for (int i = 0; i < 10; i++) // 10 corresponds to 3 correct digits
                {
                    double delta = magnitude * Math.Pow(2, -(i + 1)); // Half in size for each iteration

                    // Three positions: the current one, one lower by delta, and one higher by delta
                    
                    // Lower Model
                    EstimateParameters(bestShapeSoFar - delta, medianEst, varianceEst, out locationEst, out scaleEst);
                    GEV lowerModel = new GEV(locationEst, scaleEst, bestShapeSoFar - delta, Program.rand);
                    double lowerScore = FitnessExactModel(lowerModel);
                    
                    // Upper Model
                    EstimateParameters(bestShapeSoFar + delta, medianEst, varianceEst, out locationEst, out scaleEst);
                    GEV upperModel = new GEV(locationEst, scaleEst, bestShapeSoFar + delta, Program.rand);
                    double upperScore = FitnessExactModel(upperModel);
                    
                    // Move to the best of the three
                    double bestfitness = Math.Min(bestScore, Math.Min(upperScore, lowerScore));
                    bestScore = bestfitness;
                    if (lowerScore == bestfitness)
                    {
                        bestShapeSoFar = bestShapeSoFar - delta;
                        bestSoFar = lowerModel;
                    }
                    else if (upperScore == bestfitness)
                    {
                        bestShapeSoFar = bestShapeSoFar + delta;
                        bestSoFar = upperModel;
                    }
                }
                fitness = bestScore;
                return bestSoFar;
            }
            */
            #endregion

            GEV OptimizeBFGS(Func<Vector<double>,double> objectiveFunc, double initialShape, double initialScale, double initialLocation)
            {
                // Formatted by shape, scale, location
                var lowerBounds = CreateVector.DenseOfArray(new double[] { -10, Math.Min(-3 * initialScale, 3 * initialScale), Math.Min(-3 * initialLocation, 3 * initialLocation) });
                var upperBounds = CreateVector.DenseOfArray(new double[] { 10, Math.Max(-3 * initialScale, 3 * initialScale), Math.Max(-3 * initialLocation, 3 * initialLocation) });
                var initialGuess = CreateVector.DenseOfArray(new double[] { initialShape, initialScale, initialLocation });

                var min = FindMinimum.OfFunctionConstrained(objectiveFunc, lowerBounds, upperBounds, initialGuess);
                return new GEV(min[2], min[1], min[0], Program.rand);
            }

            #region Old code

            // Optimize for Xi
            /*double fitNeg, fitZero, fitPos;
            GEV bestNeg = OptimizeV2(-1, out fitNeg);
            GEV bestPos = OptimizeV2(1, out fitPos);
            double locZero, scaleZero;
            EstimateParameters(0, medianEst, varianceEst, out locZero, out scaleZero);
            GEV zeroModel = new GEV(locZero, scaleZero, 0, Program.rand);
            fitZero = Fitness(zeroModel);
            // Choose the best model of the three
            double minScore = Math.Min(fitNeg, Math.Min(fitPos, fitZero));
            GEV bestModel = null;
            if (fitNeg == minScore) bestModel = bestNeg;
            if (fitPos == minScore) bestModel = bestPos;
            if (fitZero == minScore) bestModel = zeroModel; // Prefer zero, then pos

            Console.WriteLine($"Best Negative model: shape: {bestNeg.shape} scale: {bestNeg.scale} location: {bestNeg.location} fitness: {fitNeg}");
            Console.WriteLine($"Best Positive model: shape: {bestPos.shape} scale: {bestPos.scale} location: {bestPos.location} fitness: {fitPos}");
            Console.WriteLine($"Zero model: shape: {zeroModel.shape} scale: {zeroModel.scale} location: {zeroModel.location} fitness: {fitZero}");
            */
            #endregion

            double scaleGuess = Math.Sqrt(6 * Statistics.VarianceEstimate(observations)) / Math.PI;
            double locationGuess = Statistics.Median(observations) + scaleGuess * Math.Log(Math.Log(2));
            double shapeGuess = 0.5; // Use Pickands estimator here in the actual model
            Func<Vector<double>, double> objectiveFunction = x => FitnessExactModel(new GEV(x[2], x[1], x[0], Program.rand));
            GEV bestModelMonteCarlo = OptimizeBFGS(objectiveFunction, shapeGuess, scaleGuess, locationGuess);

            output.WriteLine($"MC Exact GEV Model: shape{bestModelMonteCarlo.shape} location{bestModelMonteCarlo.location} scale {bestModelMonteCarlo.scale}");

            double[] sample = new double[sampleSize];
            dist.Samples(sample); // Take a sample from dist
            Sorting.Sort(sample);
            // Report the sample min and max
            output.WriteLine($"Sample maximum: {sample[sample.Length - 1]}");
            //var sorter = new List<double>(sample);
            //sorter.Sort();
            //sample = sorter.ToArray(); 
            var pickandsApprox = new PickandsApproximation(sample, useBFGS: true); // Construct a Pickands tail approx from the sample
            
            // Bootstrap observations of the distribution of the sample maximum from the Pickands model
            double[] approxObservations = new double[observations.Length];
            for (int i = 0; i < approxObservations.Length; i++)
            {
                double max = double.NegativeInfinity;
                for (int j = 0; j < sampleSize; j++)
                {
                    max = Math.Max(max, pickandsApprox.Sample());
                }
                approxObservations[i] = max;
            }

            ContinuousDistribution approxECDF = ContinuousDistribution.ECDF(approxObservations); // ECDF of the bootstrapped observations
            scaleGuess = Math.Sqrt(6 * Statistics.Variance(approxObservations)) / Math.PI;
            locationGuess = Statistics.Median(approxObservations) + scaleGuess * Math.Log(Math.Log(2));
            GEV estimatedGEVUnfitted = new GEV(location: locationGuess, scale: scaleGuess, shape: pickandsApprox.c); // Using the Pickands estimator for shape

            output.WriteLine($"UnfittedGEVModel: shape{estimatedGEVUnfitted.shape} location{estimatedGEVUnfitted.location} scale {estimatedGEVUnfitted.scale}");

            // Fit the model to the data drawn from the Pickands model
            double FitnessApproxModel(GEV model)
            {
                double val = 0;
                for (int i = 0; i < approxObservations.Length; i++)
                {
                    val += Math.Pow(model.CumulativeDistribution(approxObservations[i]) - approxECDF.CumulativeDensity(approxObservations[i]), 2);
                }
                return val;
            }

            objectiveFunction = x => FitnessApproxModel(new GEV(x[2], x[1], x[0], Program.rand));
            GEV fittedApproxModel = OptimizeBFGS(objectiveFunction, pickandsApprox.c, scaleGuess, locationGuess);

            output.WriteLine($"FittedGEVModel: shape{fittedApproxModel.shape} location{fittedApproxModel.location} scale {fittedApproxModel.scale}");

            double[] proportions = Interpolation.Linspace(0.000001, 0.999999, 2000);
            double[] observationQuantiles = Interpolation.Linspace(0.000001, 0.999999, 2000);
            for (int i = 0; i < observationQuantiles.Length; i++) { observationQuantiles[i] = Statistics.Quantile(observations, observationQuantiles[i]); }
            
            output.WriteLine("Abscissas,Monte Carlo ECDF,GEV Fit of MC ECDF,Estimated ECDF,Estimated GEV Unfitted,Estimated GEV Fitted,,ErrDistExactAbscissas,ErrDistExactValues,ErrDistModelAbscissas,ErrDistModelValues");
            for (int i = 0; i < observationQuantiles.Length; i++)
            {
                output.WriteLine($"{observationQuantiles[i]}," +
                    $"{MonteCarloDistributionOfTheMaximum.CumulativeDensity(observationQuantiles[i])}," +
                    $"{bestModelMonteCarlo.CumulativeDistribution(observationQuantiles[i])}," +
                    $"{approxECDF.CumulativeDensity(observationQuantiles[i])}," +
                    $"{estimatedGEVUnfitted.CumulativeDistribution(observationQuantiles[i])}," +
                    $"{fittedApproxModel.CumulativeDistribution(observationQuantiles[i])}," +
                    $"," + // Space
                    $"{observationQuantiles[i] - upperQuantile}," +
                    $"{MonteCarloDistributionOfTheMaximum.CumulativeDensity(observationQuantiles[i])}," +
                    //$"{quantiles[i] - sample[sample.Length - 1]}," +
                    $"{fittedApproxModel.Quantile(proportions[i]) - fittedApproxModel.location}," +
                    $"{proportions[i]}");
            }

            double[] distributionQuantiles = Interpolation.Linspace(0.000001, 0.999999, 2000);
            for (int i = 0; i < distributionQuantiles.Length; i++) { distributionQuantiles[i] = dist.InverseCumulativeDistribution(distributionQuantiles[i]); }
            output2.WriteLine("Abscissas,True CDF,Pickands Estimate");
            for (int i = 0; i < distributionQuantiles.Length; i++)
            {
                output2.WriteLine($"{distributionQuantiles[i]}," +
                    $"{dist.CumulativeDistribution(distributionQuantiles[i])}," +
                    $"{pickandsApprox.CDF(distributionQuantiles[i])}");
            }

            // Clean up
            output.Dispose();
            output2.Dispose();
            //table.Dispose();
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
            newVersion.Samples(raw);
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
