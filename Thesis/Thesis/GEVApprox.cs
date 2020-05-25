using System;

namespace Thesis
{
    class GEVApprox
    {
        // Only disagrees with the MLE ordering in a small number of cases that are necessary to avoid gradient issues
        public static double MLEV4(GEV dist, double[] evalPoints)
        {
            double sum = 0;
            const double EPSILON = 1E-6;
            for (int i = 0; i < evalPoints.Length; i++)
            {
                double density = dist.Density(evalPoints[i]);
                sum += Math.Log(EPSILON + density);
            }
            return -sum / evalPoints.Length;
        }

        public static GEV ViaMLE(double[] sortedObservations, double mu, double xi)
        {
            if (xi < -1) return new GEV(mu, (mu - sortedObservations[sortedObservations.Length - 1]) * xi, xi);
            return ViaBinarySearchEstimator(sortedObservations, mu, xi, MLEV4);
        }

        public static GEV ViaBinarySearchEstimator(double[] sortedObservations, double mu, double xi, Func<GEV, double[], double> lossFn, int refinements = 30)
        {
            // Start with the median estimator, since it's decent, reliable, and fast
            double bestSigma = xi > -1E-6 ? (mu - Statistics.Quantile(sortedObservations, 0.5)) / Math.Log(Math.Log(2)) : (Statistics.Quantile(sortedObservations, 0.5) - mu) * xi / (Math.Pow(Math.Log(2), -xi) - 1);
            bestSigma *= 2; // Expand the search range
            double bestLoss = lossFn(new GEV(mu, bestSigma, xi), sortedObservations);
            if (bestSigma > 1) refinements += 3 * (int)Math.Ceiling(Math.Log10(bestSigma));

            double stepSize = 0.5 * bestSigma;
            for (int i = 0; i < refinements; i++)
            {
                // Try forward
                double nextLoss = lossFn(new GEV(mu, bestSigma + stepSize, xi), sortedObservations);
                if (nextLoss < bestLoss)
                {
                    bestLoss = nextLoss;
                    bestSigma += stepSize;
                }
                else // If forward wasn't an improvement, try backward
                {
                    nextLoss = lossFn(new GEV(mu, bestSigma - stepSize, xi), sortedObservations);
                    if (nextLoss < bestLoss)
                    {
                        bestLoss = nextLoss;
                        bestSigma -= stepSize;
                    }
                }
                stepSize *= 0.5; // Narrow the search
            }
            return new GEV(mu, bestSigma, xi);
        }
    }
}
