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

            double[] nodes = ClenshawCurtis.GetEvalPoints(6); // returns a double[7]
            double[] weights = ClenshawCurtis.GetWeights(6); // returns a double[7]

            double a = (intervalEnd - intervalStart) / 2.0;
            double b = (intervalEnd + intervalStart) / 2.0;
            double xOfz(double z) => a * z + b;

            double sum = 0;
            double result;
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

        private static DataTable MonteCarloConvergenceTable(IContinuousDistribution[] dists, int distIndex, int[] iterationList, double exactAnswer)
        {
            DataTable MCTable = new DataTable("Monte Carlo");
            MCTable.Columns.Add("Evaluations", typeof(int));
            MCTable.Columns.Add("Estimate", typeof(double));
            MCTable.Columns.Add("Error", typeof(double));

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48); //* exactAnswer;
            int iterations;
            int itersUsed = 0;
            int count = 0;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                iterations = iterationList[i];

                // Figure out what the datapoint should be
                for (int j = itersUsed; j < iterations; j++)
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

                itersUsed = iterations;
            }

            return MCTable;
        }


        private static DataTable TrapRuleConvergenceTable(IContinuousDistribution[] dists, int[] iterationList, double exactAnswer, double lowerBound, double upperBound)
        {
            DataTable TrapTable = new DataTable("Trapezoid Rule");
            TrapTable.Columns.Add("Evaluations", typeof(int));
            TrapTable.Columns.Add("Estimate", typeof(double));
            TrapTable.Columns.Add("Error", typeof(double));

            //double integrand(double x) => N1.Density(x) * N2.CumulativeDistribution(x);
            //double lowerBound = Math.Max(N1.Mean - 8 * N1.StdDev, N2.Mean - 8 * N2.StdDev);
            //double upperBound = Math.Max(N1.Mean + 8 * N1.StdDev, N2.Mean + 8 * N2.StdDev);

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48); //* exactAnswer;
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
                double estimate = NewtonCotes.TrapezoidRule(Integrand, lowerBound, upperBound, iterations - 1);
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

        private static DataTable Simpsons38ConvergenceTable(IContinuousDistribution[] dists, int[] iterationList, double exactAnswer, double lowerBound, double upperBound)
        {
            DataTable SimpTable = new DataTable("Simpsons 3/8 Rule");
            SimpTable.Columns.Add("Evaluations", typeof(int));
            SimpTable.Columns.Add("Estimate", typeof(double));
            SimpTable.Columns.Add("Error", typeof(double));

            double Integrand(double x)
            {
                double val = dists[0].Density(x);
                for (int i = 1; i < dists.Length; i++)
                {
                    val *= dists[i].CumulativeDistribution(x);
                }
                return val;
            }

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48);// * exactAnswer;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                int iterations = iterationList[i] - 1;

                double estimate = NewtonCotes.Simpsons38Rule(Integrand, lowerBound, upperBound, iterations);
                error = Math.Abs(estimate - exactAnswer);

                // Add the datapoint
                DataRow row = SimpTable.NewRow();
                row["Evaluations"] = iterations + 1; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                SimpTable.Rows.Add(row);
            }

            return SimpTable;
        }

        private static DataTable GaussLegendreConvergenceTable(IContinuousDistribution[] dists, int[] iterationList, double exactAnswer, double lowerBound, double upperBound)
        {
            DataTable GLTable = new DataTable("Gauss Legendre");
            GLTable.Columns.Add("Evaluations", typeof(int));
            GLTable.Columns.Add("Estimate", typeof(double));
            GLTable.Columns.Add("Error", typeof(double));

            double Integrand(double x)
            {
                double val = dists[0].Density(x);
                for (int i = 1; i < dists.Length; i++)
                {
                    val *= dists[i].CumulativeDistribution(x);
                }
                return val;
            }

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48);// * exactAnswer;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                int iterations = iterationList[i];

                double estimate = GaussLegendre.Integrate(Integrand, lowerBound, upperBound, iterations);
                error = Math.Abs(estimate - exactAnswer);

                // Add the datapoint
                DataRow row = GLTable.NewRow();
                row["Evaluations"] = iterations; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                GLTable.Rows.Add(row);
            }

            return GLTable;
        }

        private static DataTable ClenshawCurtisConvergenceTable(IContinuousDistribution[] dists, int[] iterationList, double exactAnswer, double lowerBound, double upperBound)
        {
            DataTable CCTable = new DataTable("Clenshaw Curtis");
            CCTable.Columns.Add("Evaluations", typeof(int));
            CCTable.Columns.Add("Estimate", typeof(double));
            CCTable.Columns.Add("Error", typeof(double));

            double Integrand(double x)
            {
                double val = dists[0].Density(x);
                for (int i = 1; i < dists.Length; i++)
                {
                    val *= dists[i].CumulativeDistribution(x);
                }
                return val;
            }

            double error = double.PositiveInfinity;
            double epsilon = Math.Pow(2, -48);// * exactAnswer;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                int iterations = iterationList[i];

                double estimate = ClenshawCurtis.Integrate(Integrand, lowerBound, upperBound, iterations);
                error = Math.Abs(estimate - exactAnswer);

                // Add the datapoint
                DataRow row = CCTable.NewRow();
                row["Evaluations"] = iterations; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                CCTable.Rows.Add(row);
            }

            return CCTable;
        }

        private static DataTable GaussHermiteConvergenceTable(IContinuousDistribution[] dists, int[] iterationList, double exactAnswer)
        {
            DataTable GHTable = new DataTable("Gauss Hermite");
            GHTable.Columns.Add("Evaluations", typeof(int));
            GHTable.Columns.Add("Estimate", typeof(double));
            GHTable.Columns.Add("Error", typeof(double));

            // This integrand assumes the first error distribution N_0 is normal
            double HermiteIntegrand(double z)  //=> N2.CumulativeDistribution(N1.Mean + Math.Sqrt(2) * N1.StdDev * z);
            {
                double val = 1;
                for (int i = 1; i < dists.Length; i++)
                {
                    val *= dists[i].CumulativeDistribution(dists[0].Mean + Math.Sqrt(2) * dists[0].StdDev * z);
                }
                return val;
            }
            double scalar = 1.0 / Math.Sqrt(Math.PI);

            double epsilon = Math.Pow(2, -48);// * exactAnswer;
            double error = double.PositiveInfinity;

            for (int i = 0; i < iterationList.Length; i++)
            {
                if (error < epsilon) break;
                int iterations = iterationList[i];

                double estimate = scalar * GaussHermite.Integrate(HermiteIntegrand, iterations);
                error = Math.Abs(estimate - exactAnswer);

                // Add the datapoint
                DataRow row = GHTable.NewRow();
                row["Evaluations"] = iterations; // One more than the step count
                row["Estimate"] = estimate;
                row["Error"] = error;
                GHTable.Rows.Add(row);
            }

            return GHTable;
        }

        private static double QuantileRefNormal(IContinuousDistribution dist, double x) => ((Normal)dist).InverseCumulativeDistribution(x);
        
        public static void ComparisonTablePairwise()
        {
            var dists = new IContinuousDistribution[] {
                new Normal(-5, 5, Program.rand),
                new Normal(-8, 4, Program.rand) };
            ComparisonTable("Pairwise", dists, QuantileRefNormal);
        }

        public static void ComparisonTableEasySet()
        {
            var dists = new IContinuousDistribution[] {
            //new Normal(-1.2, 0.8),
            new Normal(-2, 0.6, Program.rand),
            new Normal(-1.2, 0.8, Program.rand),
            new Normal(-2.2, 0.9, Program.rand),
            new Normal(-3.0, 0.8, Program.rand),
            new Normal(-3.1, 0.4, Program.rand),
            new Normal(-4, 0.6, Program.rand),
            new Normal(-4.2, 0.2, Program.rand),
            new Normal(-4.9, 1.2, Program.rand)};
            ComparisonTable("Easy", dists, QuantileRefNormal);
        }

        public static void ComparisonTableMediumSet()
        {
            var dists = new IContinuousDistribution[] {
                new Normal(-60, 3.5, Program.rand),
                new Normal(-47, 0.5, Program.rand),
                new Normal(-50, 2, Program.rand),
                new Normal(-50.4, 0.4, Program.rand),
                new Normal(-54, 6, Program.rand),
                new Normal(-57, 0.6, Program.rand),
                new Normal(-64, 1.3, Program.rand),
                new Normal(-68, 3, Program.rand),
                new Normal(-70, 8, Program.rand),
                new Normal(-82, 0.45, Program.rand),
                new Normal(-90, 0.25, Program.rand),
                new Normal(-51, 3, Program.rand),
                new Normal(-53, 0.3, Program.rand),
                new Normal(-55, 0.35, Program.rand),
                new Normal(-68, 0.5, Program.rand),
                new Normal(-72, 1.15, Program.rand),
                new Normal(-77, 2, Program.rand),
                new Normal(-78, 0.3, Program.rand),
                new Normal(-66, 0.6, Program.rand),
                new Normal(-61, 2.8, Program.rand)};
            ComparisonTable("Medium", dists, QuantileRefNormal);
        }

        public static void ComparisonTableHardSet()
        {
            var dists = new Normal[] {
                new Normal(-42.1414975984606, 38.764397247822), new Normal(-111.152342410918, 8.06294485778047), new Normal(-34.8209126659686, 32.3179648521382), new Normal(-37.3752628301367, 9.94884264425942),
                new Normal(-84.2851671518582, 1.07780947338582), new Normal(-89.0041488068974, 27.7938992044522), new Normal(-32.0960609701481, 52.4931943177281), new Normal(-104.450689227024, 0.660695429971976),
                new Normal(-53.6926366692492, 18.0306845726305), new Normal(-99.1692119858596, 22.0919749704292), new Normal(-102.476676355163, 23.9987238702257), new Normal(-48.5283953581998, 3.99320735489477),
                new Normal(-89.9428836855282, 10.1133383507201), new Normal(-45.508503894306, 33.4127041406401), new Normal(-93.2630289080305, 24.0155809143977), new Normal(-56.8584865116467, 4.74756332965193),
                new Normal(-115.16115293408, 40.966238647701), new Normal(-97.1023270113475, 16.7190823006295), new Normal(-106.577407409437, 47.6207200133792), new Normal(-113.655667470318, 29.8036544163512),
                new Normal(-96.3038251296076, 46.4329531915872), new Normal(-102.061042181121, 14.9695912321809), new Normal(-44.1417794115925, 40.6072679279608), new Normal(-106.764515464249, 20.4390375274701),
                new Normal(-69.673701865441, 5.99234373751478), new Normal(-69.1970309076148, 28.518333057176), new Normal(-89.1144812640981, 4.14076160016157), new Normal(-30.5828645554285, 16.4368661054208),
                new Normal(-89.692327571862, 21.1362066165639), new Normal(-72.9304469295273, 30.5820177932078), new Normal(-77.9145149631647, 13.7869478799317), new Normal(-55.7293215737974, 11.4745798032422),
                new Normal(-98.227204760894, 7.54345466606997), new Normal(-114.798093559778, 25.2036611682361), new Normal(-42.9962215114193, 3.86503182132926), new Normal(-40.0773308982207, 21.1198519192032),
                new Normal(-49.5419941316212, 12.3796196536311), new Normal(-44.8125446051638, 32.8297631596309), new Normal(-36.6018383223455, 29.2112576769018), new Normal(-56.0452768040867, 47.6187587657857),
                new Normal(-34.7192644291755, 37.4987172147137), new Normal(-71.6287119849206, 2.1547458362078), new Normal(-43.781992264399, 0.949923411867221), new Normal(-38.4100036302721, 29.5017662011),
                new Normal(-40.9761644413763, 20.3216294500652), new Normal(-97.4559322806868, 1.01521909012501), new Normal(-100.94439168582, 5.89775222095568), new Normal(-104.796748428794, 16.6183789566177),
                new Normal(-98.6756571639755, 61.5315538526407), new Normal(-47.9352016588942, 21.6483684655025), new Normal(-66.9753787609356, 12.8724286024192), new Normal(-58.6997096195051, 2.84498576046393),
                new Normal(-118.104559452851, 15.6559499994815), new Normal(-60.4931464244742, 12.4650584906304), new Normal(-95.3639071480042, 20.6049492203587), new Normal(-77.8718548360777, 41.3395113151185),
                new Normal(-41.0423814778629, 31.9573273981194), new Normal(-91.1834002591662, 27.6821064023431), new Normal(-33.8129527557981, 25.3967821518823), new Normal(-102.568670680796, 24.6457779897686),
                new Normal(-118.616002364287, 3.25184840843049), new Normal(-40.9472018602293, 3.83653643557231), new Normal(-68.441824710018, 39.1221695762123), new Normal(-75.5522724083559, 33.7467116764872),
                new Normal(-93.69687432838, 23.6280139401113), new Normal(-53.411407457106, 14.1482773526276), new Normal(-101.000890935419, 30.3877787799743), new Normal(-48.0436639430433, 42.6772381314741),
                new Normal(-46.0257224243527, 9.38828057418235), new Normal(-102.150517051881, 1.09124821338376), new Normal(-100.743862930939, 1.39976206339563), new Normal(-71.6518910339223, 32.5808907695109),
                new Normal(-105.624908636368, 27.0022959060466), new Normal(-71.0093159397578, 19.7327782925212), new Normal(-73.9873999135492, 6.18973651781769), new Normal(-96.5654777308266, 61.2957402588965),
                new Normal(-118.49594361343, 29.2127784665106), new Normal(-98.4770858498748, 2.05923701225468), new Normal(-70.9857042154277, 26.6676024592747), new Normal(-46.4071282872186, 4.16192272515231),
                new Normal(-57.6087050140566, 17.0298128526502), new Normal(-67.5500358844254, 8.80211041601733), new Normal(-83.7155035544487, 21.0303912269346), new Normal(-75.2152674337965, 3.15626321805462),
                new Normal(-50.6133144464617, 14.1517886667423), new Normal(-30.9414891728847, 8.95585013658334), new Normal(-74.830577724805, 22.4019676564699), new Normal(-50.1982391983006, 28.2792327465776),
                new Normal(-39.9019854265256, 35.9683387172033), new Normal(-82.0090508001196, 38.4313275120288), new Normal(-75.3148192308055, 0.33940315380585), new Normal(-87.9649018331313, 30.2143281974075),
                new Normal(-31.7966602386788, 40.3910882870193), new Normal(-40.1859876300017, 17.6982023971142), new Normal(-94.8838348416189, 45.9081612906443), new Normal(-91.4406595414439, 25.0566122505966),
                new Normal(-58.5701521288763, 7.30544095261876), new Normal(-81.3562906416854, 24.5449678552224), new Normal(-45.3639427642084, 27.8692027743594), new Normal(-96.5106042991463, 12.0730904905099),
                new Normal(-95.9035982945569, 34.3432425899892), new Normal(-77.8491305270912, 6.08003223435756), new Normal(-40.2607792231843, 29.0742031514528), new Normal(-55.6267061827505, 13.5689329094278),
                new Normal(-57.317210024643, 2.18070397318167), new Normal(-61.4541453958754, 25.219034235059), new Normal(-30.183550974844, 44.6098379834304), new Normal(-115.070982641335, 25.1692148306899),
                new Normal(-101.136961186162, 14.7134801102065), new Normal(-72.9609850691086, 34.6575851080343), new Normal(-44.9818436748492, 37.3315424749977), new Normal(-75.9024579266147, 24.5372021465239),
                new Normal(-33.48792485955, 6.80745151606826), new Normal(-57.5176712689636, 40.6935700741415), new Normal(-86.2632542195183, 38.0864428499461), new Normal(-65.065282256277, 43.8401272532989),
                new Normal(-89.4144461292769, 34.1244135903195), new Normal(-84.182657873283, 25.2365153826144), new Normal(-91.8358709447718, 15.653145434997), new Normal(-58.3501662074713, 12.6953437070671)};
            ComparisonTable("Hard", dists, QuantileRefNormal);

            // Temp
            Program.logger.WriteLine("CC Auto Complements:");
            double[] vals = DiscardProbabilityComputation.ComplementsClenshawCurtisAutomatic(DiscardProbabilityComputation.NegateNormalDistributions(dists));
            double sum = 0;
            for (int i = 0; i < dists.Length; i++)
            {
                Program.logger.WriteLine($"1 - P(D_{i}) = {vals[i]}");
                sum += vals[i];
            }
            Program.logger.WriteLine($"Sum of 1 - P(D_i) = {sum}");
        }


        // Create a bunch of tables comparing the performances of different quadrature techniques computing 1-P(D_0) on the provided set of distributions, assuming they are already negated
        private static void ComparisonTable(string name, IContinuousDistribution[] dists, Func<IContinuousDistribution, double, double> quantileFunctionRef)
        {
            Logger logger = new Logger("ComparisonTable" + name + ".csv");

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
            double lowerBound = double.NegativeInfinity;
            double upperBound = double.NegativeInfinity;
            double epsilonM = Math.Pow(2, -52);
            for (int i = 0; i < dists.Length; i++)
            {
                lowerBound = Math.Max(lowerBound, quantileFunctionRef(dists[i], epsilonM));
                upperBound = Math.Max(upperBound, quantileFunctionRef(dists[i], 1 - epsilonM));
            }

            logger.WriteLine($"Bounds: ({lowerBound}:{upperBound})");

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
                exactDiscardComplement = NewtonCotes.TrapezoidRule(Integrand, lowerBound, upperBound, 5000);

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
            DataTable TrapTable = TrapRuleConvergenceTable(dists, iterationList, exactDiscardComplement, lowerBound, upperBound);
            DataTable SimpTable = Simpsons38ConvergenceTable(dists, iterationList, exactDiscardComplement, lowerBound, upperBound);
            DataTable GLTable = GaussLegendreConvergenceTable(dists, iterationList, exactDiscardComplement, lowerBound, upperBound);
            DataTable CCTable = ClenshawCurtisConvergenceTable(dists, iterationList, exactDiscardComplement, lowerBound, upperBound);
            DataTable GHTable = GaussHermiteConvergenceTable(dists, iterationList, exactDiscardComplement);

            #region Old Code

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
            #endregion

            //logger.WriteTablesSideBySide(MCTable, TrapTable, SimpTable, GLTable, GHTable, CCTable);
            logger.WriteTablesSideBySide(MCTable, TrapTable, SimpTable, GLTable, GHTable, CCTable);

            // Dispose of the tables and logger
            MCTable.Dispose();
            TrapTable.Dispose();
            SimpTable.Dispose();
            GLTable.Dispose();
            GHTable.Dispose();
            CCTable.Dispose();
            logger.Dispose();
        }



    }
}
