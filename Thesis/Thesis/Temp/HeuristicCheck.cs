//using Lidgren.Network;
using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.Random;

namespace ThesisOptNumericalTest
{
    static class HeuristicCheck
    {
        public static void Run()
        {

            MersenneTwister rand = new MersenneTwister();
            //NetRandom rand = new XorShiftRandom();
            

            const double xSize = 20;
            const double ySize = 20;
            const int setSize = 30;

            // Create points in the region
            List<Point> points = new List<Point>();
            for (int i = 0; i < setSize; i++)
            {
                points.Add(new Point(rand.NextDouble() * xSize, rand.NextDouble() * ySize));
            }

            #region Brute force check
            // Sort the points into ranks using brute force
            List<List<int>> LLRankSets = new List<List<int>>(20);
            //int[] pointRanks = new int[points.Count];
            List<int> layer = new List<int>();
            List<int> pointIndexList = new List<int>(points.Count);
            for (int i = 0; i < points.Count; i++) { pointIndexList.Add(i); }
            while (pointIndexList.Count > 0)
            {
                layer.Clear();

                // Find solutions that do not dominate other solutions
                for (int i = 0; i < pointIndexList.Count; i++)
                {
                    bool nonDominating = true;

                    for (int j = 0; j < pointIndexList.Count; j++)
                    {
                        if (Point.CompareLL(points[pointIndexList[i]], points[pointIndexList[j]]) < 0)
                        {
                            nonDominating = false; break;
                        }
                    }
                    if (nonDominating)
                    {
                        layer.Add(i);
                    }
                }

                // Remove them from the point index list and put them together in a rank
                List<int> rank = new List<int>(layer.Count);
                LLRankSets.Add(rank);
                for (int i = 0; i < layer.Count; i++)
                {
                    rank.Add(pointIndexList[layer[i]]);
                }
                for (int i = layer.Count - 1; i > -1; i--)
                {
                    pointIndexList.RemoveAt(layer[i]);
                }
            }

            // Sort the resulting ranks by X
            foreach (List<int> list in LLRankSets)
            {
                list.Sort((int a, int b) => points[a].x.CompareTo(points[b].x));
            }

            #endregion

            /*
            // Print the points
            for (int i = 0; i < setSize; i++)
            {
                Point p = points[i];
                Console.WriteLine($"Point {i}: ({p.x}, {p.y}) Rank {pointRanks[i]}");
                if (i < setSize - 1)
                {
                    Console.WriteLine($"Compared to Next: {Point.CompareLL(p,points[i+1])}");
                }
                
            }
            */

            for (int i = 0; i < LLRankSets.Count; i++)
            {
                Console.WriteLine($"Rank {i + 1}:");
                for (int j = 0; j < LLRankSets[i].Count; j++)
                {
                    Point p = points[LLRankSets[i][j]];
                    Console.WriteLine($"Point {j}: ({p.x}, {p.y})");
                }
            }

        }
    }


    class Point
    {
        public double x, y;

        public Point(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        public static Comparison<Point> CompareLL = (Point a, Point b) =>
        {
            if ((a.x < b.x && a.y <= b.y) || (a.x <= b.x && a.y < b.y)) { return -1; }
            if ((a.x > b.x && a.y >= b.y) || (a.x >= b.x && a.y > b.y)) { return 1; }
            else return 0;
        };
    }
}