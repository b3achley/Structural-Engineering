using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace StructuralCalculator
{
    class Program
    {
        static void Main(string[] args)
        {

            // test 1 w/ cantilever
            BeamSolver test = new BeamSolver(3.5, 1.5, 1600000.00, 5.36, new List<double> { 31.62, 31.62, 252.98 },
                new List<List<object>> {

                    // span 1
                    // pv area
                    new List<object> { 1.896, 1.896, 15.81, 31.62, "UDL", 0 },
                    // nonpv area
                    new List<object> { 1.422, 1.422, 0.00, 15.81, "UDL", 0 },

                    // span 2
                    // pv area
                    new List<object> { 1.896, 1.896, 0.00, 31.62, "UDL", 1 },

                    // span 3
                    // pv area
                    new List<object> { 1.896, 1.896, 0.00, 31.62, "UDL", 2 },
                    // nonpv area
                    new List<object> { 1.422, 1.422, 31.62, 252.98, "UDL", 2 }

                }, true, 50);

            // test 2 w/o cantilever
            //BeamSolver test = new BeamSolver(3.5, 1.5, 1600000.00, 5.36, new List<double> { 63.24, 252.98 },
            //    new List<List<object>> {

            //        // span 1
            //        // pv area
            //        new List<object> { 1.896, 1.896, 15.81, 63.24, "UDL", 0 },
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 0.00, 15.81, "UDL", 0 },

            //        // span 2
            //        // pv area
            //        new List<object> { 1.896, 1.896, 0.00, 31.62, "UDL", 1 },
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 31.62, 252.98, "UDL", 1 }

            //    }, false, 50);

            // test 3 (test 1 but pre-pv scenario)
            //BeamSolver test = new BeamSolver(3.5, 1.5, 1600000.00, 5.36, new List<double> { 31.62, 31.62, 252.98 },
            //    new List<List<object>> {

            //        // span 1
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 0.00, 31.62, "UDL", 0 },

            //        // span 2
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 0.00, 31.62, "UDL", 1 },

            //        // span 3
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 0.00, 252.98, "UDL", 2 }

            //    }, true, 50);

            // test 4 (test 2 but pre-pv scenario)
            //BeamSolver test = new BeamSolver(3.5, 1.5, 1600000.00, 5.36, new List<double> { 63.24, 252.98 },
            //    new List<List<object>> {

            //        // span 1
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 0.00, 63.24, "UDL", 0 },

            //        // span 2
            //        // nonpv area
            //        new List<object> { 1.422, 1.422, 0.00, 252.98, "UDL", 1 }

            //    }, false, 50);





            List<Matrix<double>> d = test.Demand;

            Matrix<double> shear = d[0];
            Matrix<double> moment = d[1];
            Matrix<double> slope = d[2];
            Matrix<double> deflection = d[3];

            foreach (int j in Enumerable.Range(0, 3))
            {
                Console.WriteLine("=========SPAN START=========");
                foreach (int i in Enumerable.Range(0, 51))
                {
                    Console.WriteLine(shear[i, j]);
                }
                Console.WriteLine("=========SPAN END=========");
            }

            //foreach (int j in Enumerable.Range(0, 3))
            //{
            //    Console.WriteLine("=========SPAN START=========");
            //    foreach (int i in Enumerable.Range(0, 51))
            //    {
            //        Console.WriteLine(moment[i, j]);
            //    }
            //    Console.WriteLine("=========SPAN END=========");
            //}

            //foreach (int j in Enumerable.Range(0, 3))
            //{
            //    Console.WriteLine("=========SPAN START=========");
            //    foreach (int i in Enumerable.Range(0, 51))
            //    {
            //        Console.WriteLine(slope[i, j]);
            //    }
            //    Console.WriteLine("=========SPAN END=========");
            //}

            //foreach (int j in Enumerable.Range(0, 3))
            //{
            //    Console.WriteLine("=========SPAN START=========");
            //    foreach (int i in Enumerable.Range(0, 51))
            //    {
            //        Console.WriteLine(deflection[i, j]);
            //    }
            //    Console.WriteLine("=========SPAN END=========");
            //}

            //Matrix<double> bs = test.BeamSegments;

            //foreach (int j in Enumerable.Range(0, 3))
            //{
            //    Console.WriteLine("=========SPAN START=========");
            //    foreach (int i in Enumerable.Range(0, 51))
            //    {
            //        Console.WriteLine(bs[i, j]);
            //    }
            //    Console.WriteLine("=========SPAN END=========");
            //}

            //Vector<double> sr = test.SupportReactions;

            //foreach (double d in sr)
            //{
            //    Console.WriteLine(d);
            //}


















            //// span 1 (eave)

            //    // pv area
            //        // current
            //            UniformDistributedLoad loadPv = new UniformDistributedLoad(1.901, 16.38, 32.77, 31.62);

            //        // should be
            //            UniformDistributedLoad loadPv = new UniformDistributedLoad(1.901, 16.38, 31.62, 31.62);

            //// span 2

            //    // pv area
            //        // current
            //            UniformDistributedLoad loadPv = new UniformDistributedLoad(1.901, 32.77, 65.55, 31.62);

            //        // should be
            //            UniformDistributedLoad loadPv = new UniformDistributedLoad(1.901, 0, 31.62, 31.62);

            //// span 3

            //    // pv area
            //        // current
            //            UniformDistributedLoad loadPv = new UniformDistributedLoad(1.901, 65.55, 98.33, 252.98);


            //        // should be
            //            UniformDistributedLoad loadPv = new UniformDistributedLoad(1.901, 0, 31.62, 252.98);

            //    // nonPV area
            //        // current
            //            UniformDistributedLoad loadNonPv = new UniformDistributedLoad(1.425, 98.33, 327.77, 252.98);

            //        // should be
            //            UniformDistributedLoad loadNonPv = new UniformDistributedLoad(1.425, 31.62, 252.98, 252.98);
        }
    } 
}
