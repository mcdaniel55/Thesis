using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Distributions;
using Thesis.Quadrature;

namespace Thesis
{
    static class Tables
    {

        // Print out the Monte Carlo, CLT, and Bootstrap estimates of the sampling distribution of a Gamma(2,2)
        public static void SamplingDistOfMean()
        {
            Gamma gamma = new Gamma(2, 2);

            // Monte Carlo CDF

            
            
            
            
            //Program.logger.WriteTablesSideBySide(Console.Out, );



        }

        // Show values for Gauss-Hermite quadrature in the console
        public static void ShowGaussHermiteTable()
        {
            double[] values = GaussHermite.GetNodesAndWeights(7);
            for (int i = 0; i < values.Length; i+=2)
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
                weights[6-i] = weights[i];
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
        public static void NodeAndWeightTikZTable(string quadratureName, double[] xvals, double[]yvals)
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



    }
}
