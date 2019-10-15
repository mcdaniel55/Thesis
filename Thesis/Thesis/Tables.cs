using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using Thesis.Quadrature;
using System.Data;

namespace Thesis
{
    static class Tables
    {

        // Print out the Monte Carlo, CLT, and Bootstrap estimates of the sampling distribution the mean of an Exponential(2)
        public static void SamplingDistOfMean()
        {
            Exponential dist = new Exponential(1000, Program.rand);
            const int iterations = 10000;
            const int sampleSize = 30;
            double[] sample1 = new double[sampleSize];
            double[] sample2 = new double[sampleSize];
            Logger logger = new Logger("DistOfMeanComparison.csv");

            // Monte Carlo Samples
            List<double> MCMeans = new List<double>(iterations);
            for (int i = 0; i < iterations; i++)
            {
                dist.Samples(sample1);
                MCMeans.Add(Statistics.Mean(sample1));
            }
            MCMeans.Sort();
            List<double> MCMeansError = new List<double>(iterations);
            for (int i = 0; i < iterations; i++)
            {
                MCMeansError.Add(MCMeans[i] - dist.Mean);
            }

            // CLT Model
            double sampleMean = Statistics.Mean(sample1);
            double sampleVariance = Statistics.VarianceEstimate(sample1);
            Normal CLTModel = new Normal(sampleMean, Math.Sqrt(sampleVariance / sampleSize));
            List<double> CLTNodes = new List<double>(iterations);
            for (int i = 0; i < iterations - 1; i++)
            {
                CLTNodes.Add(CLTModel.InverseCumulativeDistribution((i + 1) * 1.0 / iterations)); // This is already sorted
            }
            CLTNodes.Add(sampleMean + 8 * Math.Sqrt(sampleVariance / sampleSize)); // Model upper bound to map to 1
            List<double> CLTErrors = new List<double>(iterations);
            for (int i = 0; i < iterations; i++)
            {
                CLTErrors.Add(CLTNodes[i] - sampleMean);
            }

            // Bootstrap Model
            List<double> BootstrapMeans = new List<double>(iterations);
            for (int i = 0; i < iterations; i++)
            {
                // Resample from sample1
                for (int j = 0; j < sampleSize; j++)
                {
                    sample2[j] = sample1[Program.rand.Next(sampleSize)];
                }
                BootstrapMeans.Add(Statistics.Mean(sample2));
            }
            BootstrapMeans.Sort();
            List<double> BootstrapErrors = new List<double>(iterations);
            for (int i = 0; i < iterations; i++)
            {
                BootstrapErrors.Add(BootstrapMeans[i] - sampleMean);
            }

            // Output the results to a file
            logger.WriteLine("Monte Carlo,CLT,Bootstrap,MCErr,CLTErr,BSErr,Quantile");
            for (int i = 0; i < iterations; i++)
            {
                logger.WriteLine($"{MCMeans[i]},{CLTNodes[i]},{BootstrapMeans[i]},{MCMeansError[i]},{CLTErrors[i]},{BootstrapErrors[i]},{(i + 1) * 1.0 / iterations}");
            }
            logger.Dispose();
        }

        // Show values for Gauss-Hermite quadrature in the console
        public static void ShowGaussHermiteTable()
        {
            double[] values = GaussHermite.GetNodesAndWeights(7);
            for (int i = 0; i < values.Length; i += 2)
            {
                Console.WriteLine($"Node[{i}]: {values[i]}, Weight[{i}]: {values[i + 1]}");
            }
        }

        // Print out TikZ code for the Legendre Table
        public static void MakeLegendreTable()
        {
            double[] vals = BogaertGLWrapper.GetGLNodesAndWeights(7);
            double[] nodes = new double[7];
            double[] weights = new double[7];
            for (int i = 0; i < 3; i++)
            {
                nodes[i] = -vals[4 * i];
                nodes[6 - i] = vals[4 * i];
                weights[i] = vals[4 * i + 1];
                weights[6 - i] = weights[i];
            }
            nodes[3] = 0;
            weights[3] = vals[13];

            NodeAndWeightTikZTable("Gauss Legendre", nodes, weights);
        }

        // Print out TikZ code for the Clenshaw Curtis Table
        public static void MakeCCTable()
        {
            double[] nodes = ClenshawCurtis.GetEvalPoints(6);
            double[] weights = ClenshawCurtis.GetWeights(6);

            NodeAndWeightTikZTable("Clenshaw Curtis", nodes, weights);
        }

        // Print out nodes and weights for the Trap Rule Table to the console
        public static void ShowTrapTable()
        {
            Console.WriteLine($"1:  {-1.0}   {1.0 / 6}");
            Console.WriteLine($"2:  {-2.0 / 6}   {1.0 / 3}");
            Console.WriteLine($"3:  {-4.0 / 6}   {1.0 / 3}");
            Console.WriteLine($"4:  {0.0}   {1.0 / 3}");
            Console.WriteLine($"5:  {2.0 / 6}   {1.0 / 3}");
            Console.WriteLine($"6:  {4.0 / 6}   {1.0 / 3}");
            Console.WriteLine($"7:  {1.0}   {1.0 / 6}");
        }


        // Print out TikZ code for a quadrature rule table/graph figure
        // xvals has a 0 in the middle, and both xvals and yvals should have 7 values
        public static void NodeAndWeightTikZTable(string quadratureName, double[] xvals, double[] yvals)
        {
            Logger logger = new Logger(quadratureName + ".txt");

            logger.WriteLine(@"\begin{tikzpicture}");
            logger.WriteLine(@"\newcommand{\scale}{3.3}");
            logger.WriteLine(@"\newcommand{\graphColor}{red}");
            logger.WriteLine(@"\newcommand{\dotSize}{3pt}");
            logger.WriteLine(@"\newcommand{\pointAx}{" + xvals[0] + @"}");
            logger.WriteLine(@"\newcommand{\pointBx}{" + xvals[1] + @"}");
            logger.WriteLine(@"\newcommand{\pointCx}{" + xvals[2] + @"}");
            logger.WriteLine(@"\newcommand{\pointDx}{0}");
            logger.WriteLine(@"\newcommand{\pointEx}{" + xvals[4] + @"}");
            logger.WriteLine(@"\newcommand{\pointFx}{" + xvals[5] + @"}");
            logger.WriteLine(@"\newcommand{\pointGx}{" + xvals[6] + @"}");
            logger.WriteLine(@"\newcommand{\pointAy}{" + yvals[0] + @"}");
            logger.WriteLine(@"\newcommand{\pointBy}{" + yvals[1] + @"}");
            logger.WriteLine(@"\newcommand{\pointCy}{" + yvals[2] + @"}");
            logger.WriteLine(@"\newcommand{\pointDy}{" + yvals[3] + @"}");
            logger.WriteLine(@"\newcommand{\pointEy}{" + yvals[4] + @"}");
            logger.WriteLine(@"\newcommand{\pointFy}{" + yvals[5] + @"}");
            logger.WriteLine(@"\newcommand{\pointGy}{" + yvals[6] + @"}");
            logger.WriteLine(@"");
            logger.WriteLine(@"\node (A) at ({\scale*\pointAx},{\scale*\pointAy}) {};");
            logger.WriteLine(@"\node (B) at ({\scale*\pointBx},{\scale*\pointBy}) {};");
            logger.WriteLine(@"\node (C) at ({\scale*\pointCx},{\scale*\pointCy}) {};");
            logger.WriteLine(@"\node (D) at ({\scale*\pointDx},{\scale*\pointDy}) {};");
            logger.WriteLine(@"\node (E) at ({\scale*\pointEx},{\scale*\pointEy}) {};");
            logger.WriteLine(@"\node (F) at ({\scale*\pointFx},{\scale*\pointFy}) {};");
            logger.WriteLine(@"\node (G) at ({\scale*\pointGx},{\scale*\pointGy}) {};");
            logger.WriteLine(@"");
            logger.WriteLine(@"\draw[step={\scale*.25},gray] (-\scale,0) grid (\scale,\scale);");
            logger.WriteLine(@"\draw[line width=1pt,->] ({-1.1*\scale},0) -- ({1.1*\scale},0) node[right] {$x$};");
            logger.WriteLine(@"\draw[line width=1pt,->] (0,0) -- (0,{1.1*\scale}) node[above] {$y$};");
            logger.WriteLine(@"\draw[line width=1pt] ({1*\scale},{.05*\scale}) --++ (0,{-.1*\scale}) node[below] {$1$};");
            logger.WriteLine(@"\draw[line width=1pt] ({-1*\scale},{.05*\scale}) --++ (0,{-.1*\scale}) node[below] {$-1$};");
            logger.WriteLine(@"\draw[line width=1pt] ({.5*\scale},{.05*\scale}) --++ (0,{-.1*\scale}) node[below] {$0.5$};");
            logger.WriteLine(@"\draw[line width=1pt] ({-.5*\scale},{.05*\scale}) --++ (0,{-.1*\scale}) node[below] {$-0.5$};");
            logger.WriteLine(@"\draw[line width=1pt] ({.05*\scale},{1*\scale}) --++ ({-.1*\scale},0) node[left,fill=white] {$1$};");
            logger.WriteLine(@"\draw node[left,fill=white] at ({-.05*\scale},{.5*\scale}) {\rule{10pt}{0pt}};");
            logger.WriteLine(@"\draw[line width=1pt] ({.05*\scale},{.5*\scale}) --++ ({-.1*\scale},0) node[left=-2] {$0.5$};");
            logger.WriteLine(@"\draw[fill=\graphColor] (A) circle (\dotSize);");
            logger.WriteLine(@"\draw[fill=\graphColor] (B) circle (\dotSize);");
            logger.WriteLine(@"\draw[fill=\graphColor] (C) circle (\dotSize);");
            logger.WriteLine(@"\draw[fill=\graphColor] (D) circle (\dotSize);");
            logger.WriteLine(@"\draw[fill=\graphColor] (E) circle (\dotSize);");
            logger.WriteLine(@"\draw[fill=\graphColor] (F) circle (\dotSize);");
            logger.WriteLine(@"\draw[fill=\graphColor] (G) circle (\dotSize);");
            logger.WriteLine(@"");
            logger.WriteLine(@"\draw[\graphColor] (A) --");
            logger.WriteLine(@"(B) --");
            logger.WriteLine(@"(C) --");
            logger.WriteLine(@"(D) --");
            logger.WriteLine(@"(E) --");
            logger.WriteLine(@"(F) --");
            logger.WriteLine(@"(G);");
            logger.WriteLine(@"");
            logger.WriteLine(@"\node[above=-3.5] at (-8,0) {\scalebox{.97}{\begin{tabular}{|S[table-format=2.15]|S[table-format=1.15]|}");
            logger.WriteLine(@"\multicolumn{2}{c}{\rule[-7pt]{0pt}{10pt}\bf Evaluation Points for " + quadratureName + @"} \\");
            logger.WriteLine(@"\hline");
            logger.WriteLine(@"\multicolumn{1}{|c|}{Nodes $x_i$} & \multicolumn{1}{c|}{Weight $w_i$} \\\hline");
            logger.WriteLine(@"\pointAx & \pointAy \\");
            logger.WriteLine(@"\pointBx & \pointBy \\");
            logger.WriteLine(@"\pointCx & \pointCy \\");
            logger.WriteLine(@"0 & \pointDy \\");
            logger.WriteLine(@"\pointEx & \pointEy \\");
            logger.WriteLine(@"\pointFx & \pointFy \\");
            logger.WriteLine(@"\pointGx & \pointGy \\\hline");
            logger.WriteLine(@"\end{tabular}}};");
            logger.WriteLine(@"\end{tikzpicture}");

            logger.Dispose();
        }

        // Print out values for each step of the Gauss Hermite quadrature example
        public static void HermiteExampleTable()
        {
            Logger logger = new Logger("HermiteExampleTable.csv");
            double f(double x) => Math.Sqrt(Math.Abs(x));
            //const int order = 7;
            const int order = 13;

            double[] nodes = new double[order];
            double[] weights = new double[order];
            double[] nodeWeightArray = GaussHermite.GetNodesAndWeights(order);

            for (int i = 0; i <= order / 2; i++)
            {
                nodes[order - i - 1] = -nodeWeightArray[4 * i];
                nodes[i] = nodeWeightArray[4 * i];
                weights[i] = weights[order - i - 1] = nodeWeightArray[4 * i + 1];
            }

            double sum = 0;

            logger.WriteLine("x,w,f(x),wf(x)");
            for (int i = 0; i < order; i++)
            {
                logger.WriteLine($"{nodes[i]},{weights[i]},{f(nodes[i])},{weights[i] * f(nodes[i])}");
                sum += weights[i] * f(nodes[i]);
            }

            logger.WriteLine($"Sum of wf(x): {sum}");
            logger.WriteLine($"GH order 100: {GaussHermite.Integrate(f, 100)}");
            logger.WriteLine($"GH order 500: {GaussHermite.Integrate(f, 500)}");
            logger.WriteLine($"GL order 100000: {GaussLegendre.Integrate(x => Math.Exp(-(x * x)) * f(x), -100, 100, 100000)}");
            logger.Dispose();
        }

        // Print out values for each step of the Gauss Legendre quadrature example
        public static void LegendreExampleTable()
        {
            Logger logger = new Logger("LegendreExampleTable.csv");
            double f(double x) => Math.Sqrt(x);
            const double intervalStart = 1;
            const double intervalEnd = 4;
            const int order = 7;

            // Get nodes and weights for the GL rule of appropriate order
            double[] nodeWeightArray = new double[2 * order];
            BogaertGLWrapper.GetGLNodesAndWeightsNonAlloc(nodeWeightArray);

            double[] nodes = new double[7];
            double[] weights = new double[7];
            for (int i = 0; i < 4; i++)
            {
                nodes[6 - i] = nodeWeightArray[4 * i];
                nodes[i] = -nodeWeightArray[4 * i];
                weights[i] = weights[6 - i] = nodeWeightArray[4 * i + 1];
            }

            double a = (intervalEnd - intervalStart) / 2.0;
            double b = (intervalEnd + intervalStart) / 2.0;
            double xOfz(double z) => a * z + b;

            logger.WriteLine("x,w,h(x),sqrt(h(x)),wsqrt(h(x))");
            for (int i = 0; i < 7; i++)
            {
                logger.WriteLine($"{nodes[i]},{weights[i]},{xOfz(nodes[i])},{f(xOfz(nodes[i]))},{weights[i] * f(xOfz(nodes[i]))}");
            }

            // You can add up the last column and multiply the result by 'a' to get the answer.

            logger.Dispose();
        }

        // Print out values for each step of the Gauss Legendre quadrature example
        public static void ClenshawCurtisExampleTable()
        {
            Logger logger = new Logger("ClenshawCurtisExampleTable.csv");
            double f(double x) => Math.Sqrt(x);
            const double intervalStart = 1;
            const double intervalEnd = 4;
            const int order = 7;

            double[] nodes = ClenshawCurtis.GetEvalPoints(6); // returns a double[7]
            double[] weights = ClenshawCurtis.GetWeights(6); // returns a double[7]

            double a = (intervalEnd - intervalStart) / 2.0;
            double b = (intervalEnd + intervalStart) / 2.0;
            double xOfz(double z) => a * z + b;

            double sum = 0;
            double result = 0;
            for (int i = 0; i < 7; i++)
            {
                sum += weights[i] * f(xOfz(nodes[i]));
            }
            result = sum * a;

            logger.WriteLine("x,w,h(x),sqrt(h(x)),wsqrt(h(x))");
            for (int i = 0; i < 7; i++)
            {
                logger.WriteLine($"{nodes[i]},{weights[i]},{xOfz(nodes[i])},{f(xOfz(nodes[i]))},{weights[i] * f(xOfz(nodes[i]))}");
            }

            // Add up the last column and multiply the result by 'a' to get the answer.
            logger.WriteLine($"Sum of wsqrt(h(x)): {sum}");
            logger.WriteLine($"Result: {result}");

            logger.Dispose();
        }


        // ===== Convergence Tables =====

        // Todo: Have this pick up where it left off, not restart every time
        private static DataTable MonteCarloConvergenceTable(IContinuousDistribution[] dists, int distIndex, int[] iterationList, double exactAnswer)
        {
            DataTable MCTable = new DataTable("Monte Carlo");
            MCTable.Columns.Add("Evaluations", typeof(int));
            MCTable.Columns.Add("Estimate", typeof(double));
            MCTable.Columns.Add("Error", typeof(double));

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48) * exactAnswer;
            int iterations;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                // Figure out what the datapoint should be
                int count = 0;
                for (int j = 0; j < iterations; j++)
                {
                    double maxNj = double.NegativeInfinity;
                    double ni = dists[distIndex].Sample();
                    for (int k = 0; k < dists.Length; k++)
                    {
                        if (k == distIndex) continue;
                        maxNj = Math.Max(maxNj, dists[k].Sample());
                    }
                    if (ni > maxNj) count++;
                }
                double estimate = count * 1.0 / iterations;
                error = Math.Abs(estimate - exactAnswer);

                // Add the datapoint
                DataRow row = MCTable.NewRow();
                row["Evaluations"] = iterations;
                row["Estimate"] = estimate;
                row["Error"] = error;
                MCTable.Rows.Add(row);
            }

            return MCTable;
        }


        private static DataTable TrapRuleConvergenceTable(IContinuousDistribution[] dists, int[] iterationList, double exactAnswer, double lBound, double uBound)
        {
            DataTable TrapTable = new DataTable("Trapezoid Rule");
            TrapTable.Columns.Add("Evaluations", typeof(int));
            TrapTable.Columns.Add("Estimate", typeof(double));
            TrapTable.Columns.Add("Error", typeof(double));

            //double integrand(double x) => N1.Density(x) * N2.CumulativeDistribution(x);
            //double lowerBound = Math.Max(N1.Mean - 8 * N1.StdDev, N2.Mean - 8 * N2.StdDev);
            //double upperBound = Math.Max(N1.Mean + 8 * N1.StdDev, N2.Mean + 8 * N2.StdDev);

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48) * exactAnswer;
            int iterations;

            double Integrand(double x)
            {
                double val = dists[0].Density(x);
                for (int i = 1; i < dists.Length; i++)
                {
                    val *= dists[i].CumulativeDistribution(x);
                }
                return val;
            }

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                // # of eval points will be one greater than step count
                double estimate = NewtonCotes.TrapezoidRule(Integrand, lBound, uBound, iterations - 1);
                error = Math.Abs(estimate - exactAnswer);

                // Add the datapoint
                DataRow row = TrapTable.NewRow();
                row["Evaluations"] = iterations;
                row["Estimate"] = estimate;
                row["Error"] = error;
                TrapTable.Rows.Add(row);
            }

            return TrapTable;
        }

        public static void ComparisonTableEasySet()
        {
            Normal N1 = new Normal(5, 5);
            Normal N2 = new Normal(8, 4);
            var dists = new IContinuousDistribution[] { N1, N2 };
            double QuantileRef(IContinuousDistribution dist, double x) => ((Normal)dist).InverseCumulativeDistribution(x);
            ComparisonTable(dists, QuantileRef);
        }

        // Create a bunch of tables comparing the performances of different quadrature techniques computing 1-P(D_0) on the provided set of distributions
        public static void ComparisonTable(IContinuousDistribution[] dists, Func<IContinuousDistribution, double, double> quantileFunctionRef)
        {
            Logger logger = new Logger("ComparisonTableEasy.csv");

            // Print the distributions in use for record-keeping
            logger.Write($"Distributions:");
            for (int i = 0; i < dists.Length; i++)
            {
                logger.Write($" Dist {i}: {dists[i].GetType().Name}(mean {dists[i].Mean} stdDev {dists[i].StdDev})");
            }
            logger.WriteLine("");

            // Pre-defined numbers of iterations to use
            int[] iterationList = new int[] { 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 46, 52, 58, 64, 70, 76, 82, 88, 94, 100, 112, 124, 136, 148, 160, 184, 208, 232, 256, 280, 328, 376, 424, 472, 520, 616, 712, 808, 904, 1000, 3001, 5002, 8002, 12001, 15001, 18001};

            // === Compute 1 - P(D_0), which is P(N_0 > max N_j) ===
            // Figure out the bounds first
            double lBound = double.NegativeInfinity;
            double uBound = double.NegativeInfinity;
            double epsilonM = Math.Pow(2, -52);
            for (int i = 0; i < dists.Length; i++)
            {
                lBound = Math.Max(lBound, quantileFunctionRef(dists[i], epsilonM));
                uBound = Math.Max(uBound, quantileFunctionRef(dists[i], 1 - epsilonM));
            }

            logger.WriteLine($"Bounds: ({lBound}:{uBound})");

            // We need to figure out what the answer is so we can compute errors
            double exactDiscardComplement = -1;
            // Exact pairwise computation if possible
            if (dists.Length == 2 && dists[0].GetType() == typeof(Normal))
            {
                Normal diff = new Normal(dists[1].Mean - dists[0].Mean, Math.Sqrt(dists[0].Variance + dists[1].Variance));
                exactDiscardComplement = diff.CumulativeDistribution(0);
                logger.WriteLine($"Exact via Difference of Normals: {exactDiscardComplement}");
            }
            else // Otherwise use automatic Clenshaw Curtis to find the exact answer.
            {
                double Integrand(double x)
                {
                    double val = dists[0].Density(x);
                    for (int i = 1; i < dists.Length; i++)
                    {
                        val *= dists[i].CumulativeDistribution(x);
                    }
                    return val;
                }

                // Inefficient Placeholder
                exactDiscardComplement = NewtonCotes.TrapezoidRule(Integrand, lBound, uBound, 10000);

                logger.WriteLine($"Exact via Integral: {exactDiscardComplement}");
            }

            
            // ...


            #region Old Code
            // --- Monte Carlo ---
            /*
            DataTable MCTable = new DataTable("Monte Carlo");
            MCTable.Columns.Add("Evaluations", typeof(int));
            MCTable.Columns.Add("Estimate", typeof(double));
            MCTable.Columns.Add("Error", typeof(double));


            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48) * exactComplementDiscard;
            int iterations;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                // Figure out what the datapoint should be
                int count = 0;
                double n1, n2;
                for (int j = 0; j < iterations; j++)
                {
                    n1 = N1.Sample();
                    n2 = N2.Sample();
                    if (n1 > n2) count++;
                }
                double estimate = count * 1.0 / iterations;
                error = Math.Abs(estimate - exactComplementDiscard);

                // Add the datapoint
                DataRow row = MCTable.NewRow();
                row["Evaluations"] = iterations;
                row["Estimate"] = estimate;
                row["Error"] = error;
                MCTable.Rows.Add(row);
            }*/
            #endregion

            // Generate the tables
            DataTable MCTable = MonteCarloConvergenceTable(dists, 0, iterationList, exactDiscardComplement);
            DataTable TrapTable = TrapRuleConvergenceTable(dists, iterationList, exactDiscardComplement, lBound, uBound);
            

            // --- Simpson's 3/8 Rule ---
            /*
            DataTable SimpTable = new DataTable("Simpsons 3/8 Rule");
            SimpTable.Columns.Add("Evaluations", typeof(int));
            SimpTable.Columns.Add("Estimate", typeof(double));
            SimpTable.Columns.Add("Error", typeof(double));

            error = double.PositiveInfinity;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i] - 1;

                double estimate = NewtonCotes.Simpsons38Rule(integrand, lowerBound, upperBound, iterations);
                error = Math.Abs(estimate - exactComplementDiscard);
                
                // Add the datapoint
                DataRow row = SimpTable.NewRow();
                row["Evaluations"] = iterations + 1; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                SimpTable.Rows.Add(row);
            }

            // --- Gauss Legendre ---
            DataTable GLTable = new DataTable("Gauss Legendre");
            GLTable.Columns.Add("Evaluations", typeof(int));
            GLTable.Columns.Add("Estimate", typeof(double));
            GLTable.Columns.Add("Error", typeof(double));

            error = double.PositiveInfinity;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                double estimate = GaussLegendre.Integrate(integrand, lowerBound, upperBound, iterations);
                error = Math.Abs(estimate - exactComplementDiscard);

                // Add the datapoint
                DataRow row = GLTable.NewRow();
                row["Evaluations"] = iterations; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                GLTable.Rows.Add(row);
            }

            // --- Clenshaw Curtis ---
            DataTable CCTable = new DataTable("Clenshaw Curtis");
            CCTable.Columns.Add("Evaluations", typeof(int));
            CCTable.Columns.Add("Estimate", typeof(double));
            CCTable.Columns.Add("Error", typeof(double));

            error = double.PositiveInfinity;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                double estimate = ClenshawCurtis.Integrate(integrand, lowerBound, upperBound, iterations);
                error = Math.Abs(estimate - exactComplementDiscard);

                // Add the datapoint
                DataRow row = CCTable.NewRow();
                row["Evaluations"] = iterations; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                CCTable.Rows.Add(row);
            }

            // --- Gauss Hermite ---
            DataTable GHTable = new DataTable("Gauss Hermite");
            GHTable.Columns.Add("Evaluations", typeof(int));
            GHTable.Columns.Add("Estimate", typeof(double));
            GHTable.Columns.Add("Error", typeof(double));

            double HermiteIntegrand(double z) => N2.CumulativeDistribution(N1.Mean + Math.Sqrt(2) * N1.StdDev * z);
            double scalar = 1.0 / Math.Sqrt(Math.PI);

            error = double.PositiveInfinity;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                double estimate = scalar * GaussHermite.Integrate(HermiteIntegrand, iterations);
                error = Math.Abs(estimate - exactComplementDiscard);

                // Add the datapoint
                DataRow row = GHTable.NewRow();
                row["Evaluations"] = iterations; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                GHTable.Rows.Add(row);
            }

            logger.WriteLine($"");
            logger.WriteLine($"");
            */

            //logger.WriteTablesSideBySide(MCTable, TrapTable, SimpTable, GLTable, GHTable, CCTable);
            logger.WriteTablesSideBySide(MCTable, TrapTable);
            logger.Dispose();
        }



    }
}
