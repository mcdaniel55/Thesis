using System;

namespace Thesis.Quadrature
{
    static class NewtonCotes
    {

        // Number of evaluations is 1 + 3 * steps
        public static double Simpsons38Rule(Func<double, double> f, double start, double end, int steps)
        {
            if (steps % 3 != 0) { throw new ArgumentOutOfRangeException("Step count must be a multiple of 3."); }
            double stepSize = (end - start) / steps;
            steps /= 3;
            double sum = 0;

            // Add the boundary points
            sum += (f(start) + f(end));

            //int counter = 2; // temp
            // Add the middle points
            for (int step = 0; step < steps; step++)
            {
                sum += 3.0 * f(start + (step * 3 + 1) * stepSize);
                sum += 3.0 * f(start + (step * 3 + 2) * stepSize);
                //counter += 2;
            }
            for (int step = 0; step < steps - 1; step++)
            {
                sum += 2.0 * f(start + (step * 3 + 3) * stepSize);
                //counter++;
            }

            sum *= 3.0 * stepSize / 8.0;
            //Console.WriteLine($"Rule of size {steps * 3} used {counter} evals");
            return sum;
        }
        public static double TrapezoidRule(Func<double, double> f, double a, double b, int steps)
        {
            double h = (b - a) / steps;
            double sum = 0;
            double x = a;

            for (int i = 0; i < steps; i++)
            {
                sum += f(x) + f(x += h);
            }

            return h * 0.5 * sum;
        }



    }
}
