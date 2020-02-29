using Bridge;
using Newtonsoft.Json;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double.MathNet.Numerics.LinearAlgebra;
using System;

namespace MathNetNumerics
{
    public class App
    {
        static VectorBuilder<double> v_builder = BuilderInstance<double>.Vector;
        static MatrixBuilder<double> m_builder = BuilderInstance<double>.Matrix;

        public static Matrix<double> Hat(Vector<double> k)
        {
            Matrix<double> khat = m_builder.Dense(3, 3);
            khat[0, 1] = -k[2];
            khat[0, 2] = k[1];
            khat[1, 0] = k[2];
            khat[1, 2] = -k[0];
            khat[2, 0] = -k[1];
            khat[2, 1] = k[0];
            return khat;
        }

        public static void Main()
        {
            // Write a message to the Console
            Console.WriteLine("Welcome to Bridge.NET");

            Vector<double> k = v_builder.DenseOfArray(new[] { 1.0, 0.0, 0.0 });
            Matrix<double> I = m_builder.DenseIdentity(3);
            Matrix<double> khat = Hat(k);
            Matrix<double> khat2 = khat.Multiply(khat);
            Console.WriteLine(I + khat + khat2);

            // After building (Ctrl + Shift + B) this project, 
            // browse to the /bin/Debug or /bin/Release folder.

            // A new bridge/ folder has been created and
            // contains your projects JavaScript files. 

            // Open the bridge/index.html file in a browser by
            // Right-Click > Open With..., then choose a
            // web browser from the list

            // This application will then run in the browser.
        }
    }
}