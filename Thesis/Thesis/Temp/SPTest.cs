using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Random;
using MathNet.Numerics.Distributions;

namespace ThesisOptNumericalTest
{
    static class SPTest
    {
        public static void P3()
        {
            Xoshiro256StarStar rand = new Xoshiro256StarStar();

            int RunProcess(int time)
            {
                int popsize = 1;
                for (int t = 0; t < time; t++)
                {
                    int newpopSize = popsize;
                    for (int i = 0; i < popsize; i++)
                    {
                        double val = rand.NextDouble();
                        if (val < 0.25) { newpopSize--; }
                        if (val > 0.5) { newpopSize++; }
                        if (val > 0.75) { newpopSize++; }
                    }
                    popsize = newpopSize;
                }
                return popsize;
            }

            int sum = 0;
            int tests = 10000000;
            /*
            for (int i = 0; i < tests; i++)
            {
                sum += RunProcess(5);
            }
            double avgAfter5 = sum *1.0 / tests;

            Console.WriteLine($"E(5) = {avgAfter5}");

            int count = 0;
            tests = 5000;
            int limit = 30;
            for (int i = 0; i < tests; i++)
            {
                if (RunProcess(limit) == 0) { count++; }
            }
            double proportionDead = count * 1.0 / tests;

            Console.WriteLine($"Pi_0 = {proportionDead}");
            */

            tests = 100000000;
            sum = 0;
            double lambda1 = 1.0 / 11;
            double lambda2 = 1.0 / 9;
            double lambda3 = 1.0 / 8;
            for (int i = 0; i < tests; i++)
            {
                double s1, s2, s3;
                s1 = Exponential.Sample(rand, lambda1);
                s2 = Exponential.Sample(rand, lambda2);
                s3 = Exponential.Sample(rand, lambda3);
                if (s1 > 10 && s2 > 10 && s3 > 10) { sum++; }
            }

            double proportionLess = sum * 1.0 / tests;
            Console.WriteLine($"Proportion = {proportionLess}");
            //Console.WriteLine($"L1 / (L1 + L2) = {lambda1 / (lambda1 + lambda2)}");

            Console.ReadLine();
        }

    }
}
