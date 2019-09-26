using System.Runtime.InteropServices;

namespace Thesis.Quadrature
{
    class BogaertGLWrapper
    {
        /// <summary> Fills an array of length 2 * nodeCount with nodes and weights, with the nodes stored at even indices </summary>
        [DllImport("BogaertFastGL.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void GetNodesAndWeights(double[] result, int nodeCount);

        /// <summary> Loads and initializes the DLL code by running the GL retriever function once on a trivial input </summary>
        public static void Initialize()
        {
            GetNodesAndWeights(new double[2], 1);
        }

        /// <summary> Fills an array of length 2 * order with nodes and weights for Gauss-Legendre quadrature, with the nodes stored at even indices </summary>
        /// <param name="order"> The order of the quadrature rule </param>
        public static double[] GetGLNodesAndWeights(int order)
        {
            double[] result = new double[2 * order];
            GetNodesAndWeights(result, order);
            return result;
        }

        /// <summary> Non-allocating version of the GetGLNodesAndWeights method. Requires a double[] of length 2 * order </summary>
        /// <param name="array"> An array of length 2 * order to hold the nodes and weights </param>
        public static void GetGLNodesAndWeightsNonAlloc(double[] array)
        {
            GetNodesAndWeights(array, array.Length / 2);
        }
    }
}
