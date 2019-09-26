using System;

namespace Thesis.Quadrature
{
    static class GaussHermite
    {
        public static double Integrate(Func<double, double> f, double[] nodesAndWeights)
        {
            double sum = 0;
            for (int i = 0; i < nodesAndWeights.Length; i++)
            {
                //sum += f(nodesAndWeights[i]) * nodesAndWeights[++i];
                // Slower, safer version for testing
                double val = f(nodesAndWeights[i]) * nodesAndWeights[++i];
                if (double.IsNaN(val) || double.IsInfinity(val)) continue;
                sum += val;
            }
            return sum;
        }

        public static double Integrate(Func<double, double> f, int order, double center = 0, double scale = 1)
        {
            double[] nodesAndWeights = GetNodesAndWeights(order, center, scale);
            return Integrate(f, nodesAndWeights);
        }

        static Func<double, double> HermitePolynomial0 = (double x) => 1;
        static Func<double, double> HermitePolynomial1 = (double x) => 2 * x;

        // This is an inelegant approach, but this method is only here as a curiosity
        public static Func<double, double> GetHermitePolynomial(int degree)
        {
            if (degree < 0) { throw new ArgumentOutOfRangeException(nameof(degree)); }
            if (degree == 0) { return HermitePolynomial0; }
            if (degree == 1) { return HermitePolynomial1; }
            Func<double, double> HPrev2 = HermitePolynomial0;
            Func<double, double> HPrev1 = HermitePolynomial1;
            Func<double, double> H = (x) => 0;
            for (int i = 2; i <= degree; i++)
            {
                H = (double x) => 2 * x * HPrev1(x) - 2 * i * HPrev2(x);
                HPrev2 = HPrev1;
                HPrev1 = H;
            }

            return H;
        }

        // Below is a C# port of J. Burkardt's C++ port of a Fortran77 routine
        public static double[] GetNodesAndWeights(int n, double centerPoint = 0, double scaleFactor = 1)
        {
            double[] result = new double[2 * n];
            double alpha = 0;
            int i; // Index
            double pi = 3.14159265358979323846264338327950; // Why not use Math.PI?            
            var w = new double[n];
            var x = new double[n];

            // Do the computations
            Function1();

            // Normalize the rule (On second thought, don't!)
            /*
            for (i = 0; i < n; i++)
            {
                w[i] = w[i] * Math.Sqrt(scaleFactor) / Math.Sqrt(pi);
            }
            */

            // Sort from smallest weight to largest to maximize the potential for small 
            // contributions to accumulate into sums that are more noticeable after rounding
            for (i = 0; i < n / 2; i++)
            {
                result[4 * i] = x[i];
                result[4 * i + 1] = w[i];
                result[4 * i + 2] = -x[i];
                result[4 * i + 3] = w[i];
            }
            if (n % 2 == 1)
            {
                result[2 * n - 2] = x[n / 2];
                result[2 * n - 1] = w[n / 2];
            }
            return result;


            // These functions are passed the entire state of this method,
            // and it is more reasonable to write them as internal methods
            // than to pass them copies and references to that state.
            void Function1() // cgqf
            {
                Function2(); // cdgqf
                var mlt = new int[n];
                var ndx = new int[n];
                for (i = 0; i < n; i++)
                {
                    mlt[i] = 1;
                    ndx[i] = i + 1;
                }
                Function6(mlt, ndx);
            }

            void Function2() // cdgqf
            {
                var aj = new double[n];
                var bj = new double[n];
                double zemu = Function3(aj, bj);
                Function4(aj, bj, zemu);
            }

            double Function3(double[] aj, double[] bj) // Class_Matrix
            {
                // #begin facepalm
                double temp = 2.220446049250313E-016;
                double temp2 = 0.5;

                // Gamma compat check
                if (500.0 * temp < Math.Abs(Math.Pow(Gamma(temp2), 2) - pi))
                {
                    throw new Exception("CLASS_MATRIX - Fatal Error: Gamma function does not match machine parameters");
                }

                for (i = 0; i < n; i++) // Unnecessary inititalization
                {
                    aj[i] = 0.0;
                }

                for (i = 1; i <= n; i++)
                {
                    bj[i - 1] = Math.Sqrt((i + alpha * (i % 2)) / 2.0);
                }

                return Gamma((alpha + 1.0) / 2.0); ;
            }

            void Function4(double[] aj, double[] bj, double zemu) // sgqf
            {
                if (zemu <= 0.0)
                {
                    throw new Exception("SGQF - Fatal Error: ZEMU <= 0");
                }

                //  Set up vectors for Function5
                for (i = 0; i < n; i++)
                {
                    x[i] = aj[i];
                }

                w[0] = Math.Sqrt(zemu);
                for (i = 1; i < n; i++) // Unnecessary initialization, could just keep the line above
                {
                    w[i] = 0.0;
                }

                //  Diagonalize the Jacobi matrix.
                Function5(bj);

                // Square the weights, for some reason
                for (i = 0; i < n; i++)
                {
                    w[i] = w[i] * w[i];
                }
            }

            void Function5(double[] bj) // imtqlx
            {
                // d is x
                // e is bj
                // z is w
                double b;
                double c;
                double f;
                double g;
                int ii;
                int itn = 30;
                int j;
                int k;
                int l;
                int m;
                int mml;
                double p;
                double prec;
                double r;
                double s;

                prec = 2.220446049250313E-016;

                if (n == 1)
                {
                    return;
                }

                bj[n - 1] = 0.0;

                for (l = 1; l <= n; l++)
                {
                    j = 0;
                    for (; ; )
                    {
                        for (m = l; m <= n; m++)
                        {
                            if (m == n)
                            {
                                break;
                            }

                            if (Math.Abs(bj[m - 1]) <= prec * (Math.Abs(x[m - 1]) + Math.Abs(x[m])))
                            {
                                break;
                            }
                        }
                        p = x[l - 1];
                        if (m == l)
                        {
                            break;
                        }
                        if (itn <= j)
                        {
                            throw new Exception("IMTQLX - Fatal Error: Iteration limit exceeded");
                        }
                        j = j + 1;
                        g = (x[l] - p) / (2.0 * bj[l - 1]);
                        r = Math.Sqrt(g * g + 1.0);
                        g = x[m - 1] - p + bj[l - 1] / (g + Math.Abs(r) * Sign(g));
                        s = 1.0;
                        c = 1.0;
                        p = 0.0;
                        mml = m - l;

                        for (ii = 1; ii <= mml; ii++)
                        {
                            i = m - ii;
                            f = s * bj[i - 1];
                            b = c * bj[i - 1];

                            if (Math.Abs(g) <= Math.Abs(f))
                            {
                                c = g / f;
                                r = Math.Sqrt(c * c + 1.0);
                                bj[i] = f * r;
                                s = 1.0 / r;
                                c = c * s;
                            }
                            else
                            {
                                s = f / g;
                                r = Math.Sqrt(s * s + 1.0);
                                bj[i] = g * r;
                                c = 1.0 / r;
                                s = s * c;
                            }
                            g = x[i] - p;
                            r = (x[i - 1] - g) * s + 2.0 * c * b;
                            p = s * r;
                            x[i] = g + p;
                            g = c * r - b;
                            f = w[i];
                            w[i] = s * w[i - 1] + c * f;
                            w[i - 1] = c * w[i - 1] - s * f;
                        }
                        x[l - 1] = x[l - 1] - p;
                        bj[l - 1] = g;
                        bj[m - 1] = 0.0;
                    }
                }
                //  Sorting.
                for (ii = 2; ii <= n; ii++)
                {
                    i = ii - 1;
                    k = i;
                    p = x[i - 1];

                    for (j = ii; j <= n; j++)
                    {
                        if (x[j - 1] < p)
                        {
                            k = j;
                            p = x[j - 1];
                        }
                    }

                    if (k != i)
                    {
                        x[k - 1] = x[i - 1];
                        x[i - 1] = p;
                        p = w[i - 1];
                        w[i - 1] = w[k - 1];
                        w[k - 1] = p;
                    }
                }
            }

            void Function6(int[] mlt, int[] ndx) // scqf
            {
                // n, x, mlt, w, n, ndx, w, x, 6, alpha, beta, centerPoint, scaleFactor
                //int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[], double st[], int kind, double alpha, double beta, double a, double b
                double al;
                double be;
                int k;
                int l;
                double p;
                double shft;
                double slp;
                double tmp;

                //parchk(kind, 1, alpha, beta);
                if (scaleFactor <= 0.0)
                {
                    throw new Exception($"SCQF - Fatal Error: B <= 0");
                }

                shft = centerPoint;
                slp = 1.0 / Math.Sqrt(scaleFactor);
                al = alpha;
                be = 0.0;

                p = Math.Pow(slp, al + be + 1.0);

                for (k = 0; k < n; k++)
                {
                    x[k] = shft + slp * x[k];
                    l = Math.Abs(ndx[k]);

                    if (l != 0)
                    {
                        tmp = p;
                        for (i = l - 1; i <= l - 1 + mlt[k] - 1; i++)
                        {
                            w[i] = w[i] * tmp;
                            tmp = tmp * slp;
                        }
                    }
                }
            }
        }

        static double Gamma(double x)
        {
            double[] c = { -1.910444077728E-03,           8.4171387781295E-04,
                           -5.952379913043012E-04,        7.93650793500350248E-04,
                           -2.777777777777681622553E-03,  8.333333333333333331554247E-02,
                            5.7083835261E-03 };
            double eps = 2.22E-16;
            double fact;
            int i;
            int n;
            double[] p = { -1.71618513886549492533811E+00,  2.47656508055759199108314E+01,
                           -3.79804256470945635097577E+02,  6.29331155312818442661052E+02,
                            8.66966202790413211295064E+02, -3.14512729688483675254357E+04,
                           -3.61444134186911729807069E+04,  6.64561438202405440627855E+04 };
            bool parity;
            const double pi = 3.1415926535897932384626434;
            double[] q = { -3.08402300119738975254353E+01,  3.15350626979604161529144E+02,
                           -1.01515636749021914166146E+03, -3.10777167157231109440444E+03,
                            2.25381184209801510330112E+04,  4.75584627752788110767815E+03,
                           -1.34659959864969306392456E+05, -1.15132259675553483497211E+05 };
            double res;
            const double sqrtpi = 0.9189385332046727417803297;
            double sum;
            double value;
            double xbig = 171.624;
            double xden;
            double xinf = 1.79E+308;
            double xminin = 2.23E-308;
            double xnum;
            double y;
            double y1;
            double ysq;
            double z;

            parity = false;
            fact = 1.0;
            n = 0;
            y = x;
            //  Argument is negative.
            if (y <= 0.0)
            {
                y = -x;
                y1 = (int)y;
                res = y - y1;

                if (res != 0.0)
                {
                    if (y1 != (int)(y1 * 0.5) * 2.0)
                    {
                        parity = true;
                    }

                    fact = -pi / Math.Sin(pi * res);
                    y = y + 1.0;
                }
                else
                {
                    res = xinf;
                    value = res;
                    return value;
                }
            }
            //  Argument is positive.
            if (y < eps)
            {
                //  Argument < EPS.
                if (xminin <= y)
                {
                    res = 1.0 / y;
                }
                else
                {
                    res = xinf;
                    value = res;
                    return value;
                }
            }
            else if (y < 12.0)
            {
                y1 = y;
                //  0.0 < argument < 1.0.
                if (y < 1.0)
                {
                    z = y;
                    y = y + 1.0;
                }
                //  1.0 < argument < 12.0.
                //  Reduce argument if necessary.
                else
                {
                    n = (int)(y) - 1;
                    y = y - n;
                    z = y - 1.0;
                }

                //  Evaluate approximation for 1.0 < argument < 2.0.
                xnum = 0.0;
                xden = 1.0;
                for (i = 0; i < 8; i++)
                {
                    xnum = (xnum + p[i]) * z;
                    xden = xden * z + q[i];
                }
                res = xnum / xden + 1.0;
                //
                //  Adjust result for case  0.0 < argument < 1.0.
                //
                if (y1 < y)
                {
                    res = res / y1;
                }
                //
                //  Adjust result for case 2.0 < argument < 12.0.
                //
                else if (y < y1)
                {
                    for (i = 1; i <= n; i++)
                    {
                        res = res * y;
                        y = y + 1.0;
                    }
                }
            }
            else
            {
                //
                //  Evaluate for 12.0 <= argument.
                //
                if (y <= xbig)
                {
                    ysq = y * y;
                    sum = c[6];
                    for (i = 0; i < 6; i++)
                    {
                        sum = sum / ysq + c[i];
                    }
                    sum = sum / y - y + sqrtpi;
                    sum = sum + (y - 0.5) * Math.Log(y);
                    res = Math.Exp(sum);
                }
                else
                {
                    res = xinf;
                    value = res;
                    return value;
                }
            }
            //
            //  Final adjustments and return.
            //
            if (parity)
            {
                res = -res;
            }

            if (fact != 1.0)
            {
                res = fact / res;
            }

            value = res;

            return value;
        }

        // Special sign function
        static double Sign(double x)
        {
            return x < 0 ? -1 : 1;
        }
    }
}
