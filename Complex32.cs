using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mathnetnumericsbridge
{

    [Serializable]
    [System.Runtime.CompilerServices.TypeForwardedFrom("System.Numerics, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089")]
    public struct Complex32 : IEquatable<Complex32>, IFormattable
    {
        public static readonly Complex32 Zero = new Complex32(0.0f, 0.0f);
        public static readonly Complex32 One = new Complex32(1.0f, 0.0f);
        public static readonly Complex32 ImaginaryOne = new Complex32(0.0f, 1.0f);
        public static readonly Complex32 NaN = new Complex32(float.NaN, float.NaN);
        public static readonly Complex32 Infinity = new Complex32(float.PositiveInfinity, float.PositiveInfinity);

        private const float InverseOfLog10 = (float)0.43429448; // 1 / Log(10)

        // This is the largest x for which (Hypot(x,x) + x) will not overflow. It is used for branching inside Sqrt.
        private static readonly float s_sqrtRescaleThreshold = float.MaxValue /(float) (Math.Sqrt(2.0) + 1.0);

        // This is the largest x for which 2 x^2 will not overflow. It is used for branching inside Asin and Acos.
        private static readonly float s_asinOverflowThreshold = (float)Math.Sqrt(float.MaxValue) / (float)2.0;

        // This value is used inside Asin and Acos.
        private static readonly float s_log2 = (float)Math.Log(2.0);

        // Do not rename, these fields are needed for binary serialization
        private float m_real; // Do not rename (binary serialization)
        private float m_imaginary; // Do not rename (binary serialization)

        public Complex32(float real, float imaginary)
        {
            m_real = real;
            m_imaginary = imaginary;
        }

        public float Real { get { return m_real; } }
        public float Imaginary { get { return m_imaginary; } }

        public float Magnitude { get { return Abs(this); } }
        public float Phase { get { return (float)Math.Atan2(m_imaginary, m_real); } }

        public static Complex32 FromPolarCoordinates(float magnitude, float phase)
        {
            return new Complex32((float)magnitude * (float)Math.Cos(phase), (float)magnitude * (float)Math.Sin(phase));
        }

        public static Complex32 Negate(Complex32 value)
        {
            return -value;
        }

        public static Complex32 Add(Complex32 left, Complex32 right)
        {
            return left + right;
        }

        public static Complex32 Add(Complex32 left, float right)
        {
            return left + right;
        }

        public static Complex32 Add(float left, Complex32 right)
        {
            return left + right;
        }

        public static Complex32 Subtract(Complex32 left, Complex32 right)
        {
            return left - right;
        }

        public static Complex32 Subtract(Complex32 left, float right)
        {
            return left - right;
        }

        public static Complex32 Subtract(float left, Complex32 right)
        {
            return left - right;
        }

        public static Complex32 Multiply(Complex32 left, Complex32 right)
        {
            return left * right;
        }

        public static Complex32 Multiply(Complex32 left, float right)
        {
            return left * right;
        }

        public static Complex32 Multiply(float left, Complex32 right)
        {
            return left * right;
        }

        public static Complex32 Divide(Complex32 dividend, Complex32 divisor)
        {
            return dividend / divisor;
        }

        public static Complex32 Divide(Complex32 dividend, float divisor)
        {
            return dividend / divisor;
        }

        public static Complex32 Divide(float dividend, Complex32 divisor)
        {
            return dividend / divisor;
        }

        public static Complex32 operator -(Complex32 value)  /* Unary negation of a Complex32 number */
        {
            return new Complex32(-value.m_real, -value.m_imaginary);
        }

        public static Complex32 operator +(Complex32 left, Complex32 right)
        {
            return new Complex32(left.m_real + right.m_real, left.m_imaginary + right.m_imaginary);
        }

        public static Complex32 operator +(Complex32 left, float right)
        {
            return new Complex32(left.m_real + right, left.m_imaginary);
        }

        public static Complex32 operator +(float left, Complex32 right)
        {
            return new Complex32(left + right.m_real, right.m_imaginary);
        }

        public static Complex32 operator -(Complex32 left, Complex32 right)
        {
            return new Complex32(left.m_real - right.m_real, left.m_imaginary - right.m_imaginary);
        }

        public static Complex32 operator -(Complex32 left, float right)
        {
            return new Complex32(left.m_real - right, left.m_imaginary);
        }

        public static Complex32 operator -(float left, Complex32 right)
        {
            return new Complex32(left - right.m_real, -right.m_imaginary);
        }

        public static Complex32 operator *(Complex32 left, Complex32 right)
        {
            // Multiplication:  (a + bi)(c + di) = (ac -bd) + (bc + ad)i
            float result_realpart = (left.m_real * right.m_real) - (left.m_imaginary * right.m_imaginary);
            float result_imaginarypart = (left.m_imaginary * right.m_real) + (left.m_real * right.m_imaginary);
            return new Complex32(result_realpart, result_imaginarypart);
        }

        public static Complex32 operator *(Complex32 left, float right)
        {
            if (!float.IsFinite(left.m_real))
            {
                if (!float.IsFinite(left.m_imaginary))
                {
                    return new Complex32(float.NaN, float.NaN);
                }

                return new Complex32(left.m_real * right, float.NaN);
            }

            if (!float.IsFinite(left.m_imaginary))
            {
                return new Complex32(float.NaN, left.m_imaginary * right);
            }

            return new Complex32(left.m_real * right, left.m_imaginary * right);
        }

        public static Complex32 operator *(float left, Complex32 right)
        {
            if (!float.IsFinite(right.m_real))
            {
                if (!float.IsFinite(right.m_imaginary))
                {
                    return new Complex32(float.NaN, float.NaN);
                }

                return new Complex32(left * right.m_real, float.NaN);
            }

            if (!float.IsFinite(right.m_imaginary))
            {
                return new Complex32(float.NaN, left * right.m_imaginary);
            }

            return new Complex32(left * right.m_real, left * right.m_imaginary);
        }

        public static Complex32 operator /(Complex32 left, Complex32 right)
        {
            // Division : Smith's formula.
            float a = left.m_real;
            float b = left.m_imaginary;
            float c = right.m_real;
            float d = right.m_imaginary;

            // Computing c * c + d * d will overflow even in cases where the actual result of the division does not overflow.
            if (Math.Abs(d) < Math.Abs(c))
            {
                float doc = d / c;
                return new Complex32((a + b * doc) / (c + d * doc), (b - a * doc) / (c + d * doc));
            }
            else
            {
                float cod = c / d;
                return new Complex32((b + a * cod) / (d + c * cod), (-a + b * cod) / (d + c * cod));
            }
        }

        public static Complex32 operator /(Complex32 left, float right)
        {
            // IEEE prohibit optimizations which are value changing
            // so we make sure that behaviour for the simplified version exactly match
            // full version.
            if (right == 0)
            {
                return new Complex32(float.NaN, float.NaN);
            }

            if (!float.IsFinite(left.m_real))
            {
                if (!float.IsFinite(left.m_imaginary))
                {
                    return new Complex32(float.NaN, float.NaN);
                }

                return new Complex32(left.m_real / right, float.NaN);
            }

            if (!float.IsFinite(left.m_imaginary))
            {
                return new Complex32(float.NaN, left.m_imaginary / right);
            }

            // Here the actual optimized version of code.
            return new Complex32(left.m_real / right, left.m_imaginary / right);
        }

        public static Complex32 operator /(float left, Complex32 right)
        {
            // Division : Smith's formula.
            float a = left;
            float c = right.m_real;
            float d = right.m_imaginary;

            // Computing c * c + d * d will overflow even in cases where the actual result of the division does not overflow.
            if (Math.Abs(d) < Math.Abs(c))
            {
                float doc = d / c;
                return new Complex32(a / (c + d * doc), (-a * doc) / (c + d * doc));
            }
            else
            {
                float cod = c / d;
                return new Complex32(a * cod / (d + c * cod), -a / (d + c * cod));
            }
        }

        public static float Abs(Complex32 value)
        {
            return Hypot(value.m_real, value.m_imaginary);
        }

        private static float Hypot(float a, float b)
        {
            // Using
            //   sqrt(a^2 + b^2) = |a| * sqrt(1 + (b/a)^2)
            // we can factor out the larger component to dodge overflow even when a * a would overflow.

            a = (float)Math.Abs(a);
            b = (float)Math.Abs(b);

            float small, large;
            if (a < b)
            {
                small = a;
                large = b;
            }
            else
            {
                small = b;
                large = a;
            }

            if (small == 0.0)
            {
                return (large);
            }
            else if (float.IsPositiveInfinity(large) && !float.IsNaN(small))
            {
                // The NaN test is necessary so we don't return +inf when small=NaN and large=+inf.
                // NaN in any other place returns NaN without any special handling.
                return (float.PositiveInfinity);
            }
            else
            {
                float ratio = small / large;
                return (float)(large * Math.Sqrt(1.0 + ratio * ratio));
            }

        }


        private static float Log1P(float x)
        {
            // Compute log(1 + x) without loss of accuracy when x is small.

            // Our only use case so far is for positive values, so this isn't coded to handle negative values.


            float xp1 = 1.0f + x;
            if (xp1 == 1.0)
            {
                return x;
            }
            else if (x < 0.75)
            {
                // This is accurate to within 5 ulp with any floating-point system that uses a guard digit,
                // as proven in Theorem 4 of "What Every Computer Scientist Should Know About Floating-Point
                // Arithmetic" (https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
                return x * (float)Math.Log(xp1) / (float)(xp1 - 1.0);
            }
            else
            {
                return (float)Math.Log(xp1);
            }

        }

        public static Complex32 Conjugate(Complex32 value)
        {
            // Conjugate of a Complex32 number: the conjugate of x+i*y is x-i*y
            return new Complex32(value.m_real, -value.m_imaginary);
        }

        public bool IsZero(Complex32 value)
        {
            // Conjugate of a Complex number: the conjugate of x+i*y is x-i*y
            if (value.Real == 0.0f && value.Imaginary == 0.0f) return true;
            return false;
        }


        public static Complex32 Reciprocal(Complex32 value)
        {
            // Reciprocal of a Complex32 number : the reciprocal of x+i*y is 1/(x+i*y)
            if (value.m_real == 0 && value.m_imaginary == 0)
            {
                return Zero;
            }
            return One / value;
        }

        public static bool operator ==(Complex32 left, Complex32 right)
        {
            return left.m_real == right.m_real && left.m_imaginary == right.m_imaginary;
        }

        public static bool operator !=(Complex32 left, Complex32 right)
        {
            return left.m_real != right.m_real || left.m_imaginary != right.m_imaginary;
        }
        /*
        public override bool Equals(object? obj)
        {
            if (!(obj is Complex32)) return false;
            return Equals((Complex32)obj);
        }
        */
        public bool Equals(Complex32 value)
        {
            return m_real.Equals(value.m_real) && m_imaginary.Equals(value.m_imaginary);
        }

        public override int GetHashCode()
        {
            int n1 = 99999997;
            int realHash = m_real.GetHashCode() % n1;
            int imaginaryHash = m_imaginary.GetHashCode();
            int finalHash = realHash ^ imaginaryHash;
            return finalHash;
        }

        public override string ToString()
        {
            return string.Format( "({0}, {1})", m_real, m_imaginary);
        }

        

        public string ToString(string stringin, IFormatProvider provider)
        {
            return string.Format(provider, "({0}, {1})", m_real, m_imaginary);
        }
        

        public static Complex32 Sin(Complex32 value)
        {
            // We need both sinh and cosh of imaginary part. To avoid multiple calls to Math.Exp with the same value,
            // we compute them both here from a single call to Math.Exp.
            float p = (float)Math.Exp(value.m_imaginary);
            float q = 1.0f / p;
            float sinh = (p - q) * 0.5f;
            float cosh = (p + q) * 0.5f;
            return new Complex32((float)Math.Sin(value.m_real) * cosh, (float)Math.Cos(value.m_real) * sinh);
            // There is a known limitation with this algorithm: inputs that cause sinh and cosh to overflow, but for
            // which sin or cos are small enough that sin * cosh or cos * sinh are still representable, nonetheless
            // produce overflow. For example, Sin((0.01, 711.0)) should produce (~3.0E306, PositiveInfinity), but
            // instead produces (PositiveInfinity, PositiveInfinity).
        }


        public static Complex32 Sinh(Complex32 value)
        {
            // Use sinh(z) = -i sin(iz) to compute via sin(z).
            Complex32 sin = Sin(new Complex32(-value.m_imaginary, value.m_real));
            return new Complex32(sin.m_imaginary, -sin.m_real);
        }
        

        public static Complex32 Asin(Complex32 value)
        {
            float b, bPrime, v;
            Asin_Internal((float)Math.Abs(value.Real), (float)Math.Abs(value.Imaginary), out b, out bPrime, out v);

            float u;
            if (bPrime < 0.0)
            {
                u = (float)Math.Asin(b);
            }
            else
            {
                u = (float)Math.Atan(bPrime);
            }

            if (value.Real < 0.0) u = -u;
            if (value.Imaginary < 0.0) v = -v;

            return new Complex32(u, v);
        }

        public static Complex32 Cos(Complex32 value)
        {
            float p = (float)Math.Exp(value.m_imaginary);
            float q = 1.0f / p;
            float sinh = (p - q) * 0.5f;
            float cosh = (p + q) * 0.5f;
            return new Complex32((float)Math.Cos(value.m_real) * cosh, (float)-Math.Sin(value.m_real) * sinh);
        }


        public static Complex32 Cosh(Complex32 value)
        {
            // Use cosh(z) = cos(iz) to compute via cos(z).
            return Cos(new Complex32(-value.m_imaginary, value.m_real));
        }

        public static Complex32 Acos(Complex32 value)
        {
            float b, bPrime, v;
            Asin_Internal((float)Math.Abs(value.Real), (float)Math.Abs(value.Imaginary), out b, out bPrime, out v);

            float u;
            if (bPrime < 0.0)
            {
                u = (float)Math.Acos(b);
            }
            else
            {
                u = (float)Math.Atan(1.0 / bPrime);
            }

            if (value.Real < 0.0) u = (float)Math.PI - u;
            if (value.Imaginary > 0.0) v = -v;

            return new Complex32(u, v);
        }

        public static Complex32 Tan(Complex32 value)
        {
            // tan z = sin z / cos z, but to avoid unnecessary repeated trig computations, use
            //   tan z = (sin(2x) + i sinh(2y)) / (cos(2x) + cosh(2y))
            // (see Abramowitz & Stegun 4.3.57 or derive by hand), and compute trig functions here.

            // This approach does not work for |y| > ~355, because sinh(2y) and cosh(2y) overflow,
            // even though their ratio does not. In that case, divide through by cosh to get:
            //   tan z = (sin(2x) / cosh(2y) + i \tanh(2y)) / (1 + cos(2x) / cosh(2y))
            // which correctly computes the (tiny) real part and the (normal-sized) imaginary part.

            float x2 = 2.0f * value.m_real;
            float y2 = 2.0f* value.m_imaginary;
            float p = (float)Math.Exp(y2);
            float q = 1.0f / p;
            float cosh = (p + q) * 0.5f;
            if (Math.Abs(value.m_imaginary) <= 4.0)
            {
                float sinh = (p - q) * 0.5f;
                float D = (float)Math.Cos(x2) + cosh;
                return new Complex32((float)Math.Sin(x2) / D, sinh / D);
            }
            else
            {
                float D = 1.0f + (float)Math.Cos(x2) / cosh;
                return new Complex32((float)Math.Sin(x2) / cosh / D, (float)Math.Tanh(y2) / D);
            }
        }


        public static Complex32 Tanh(Complex32 value)
        {
            // Use tanh(z) = -i tan(iz) to compute via tan(z).
            Complex32 tan = Tan(new Complex32(-value.m_imaginary, value.m_real));
            return new Complex32(tan.m_imaginary, -tan.m_real);
        }

        public static Complex32 Atan(Complex32 value)
        {
            Complex32 two = new Complex32(2.0f, 0.0f);
            return (ImaginaryOne / two) * (Log(One - ImaginaryOne * value) - Log(One + ImaginaryOne * value));
        }

        private static void Asin_Internal(float x, float y, out float b, out float bPrime, out float v)
        {

            // This method for the inverse Complex32 sine (and cosine) is described in Hull, Fairgrieve,
            // and Tang, "Implementing the Complex32 Arcsine and Arccosine Functions Using Exception Handling",
            // ACM Transactions on Mathematical Software (1997)
            // (https://www.researchgate.net/profile/Ping_Tang3/publication/220493330_Implementing_the_Complex32_Arcsine_and_Arccosine_Functions_Using_Exception_Handling/links/55b244b208ae9289a085245d.pdf)

            // First, the basics: start with sin(w) = (e^{iw} - e^{-iw}) / (2i) = z. Here z is the input
            // and w is the output. To solve for w, define t = e^{i w} and multiply through by t to
            // get the quadratic equation t^2 - 2 i z t - 1 = 0. The solution is t = i z + sqrt(1 - z^2), so
            //   w = arcsin(z) = - i log( i z + sqrt(1 - z^2) )
            // Decompose z = x + i y, multiply out i z + sqrt(1 - z^2), use log(s) = |s| + i arg(s), and do a
            // bunch of algebra to get the components of w = arcsin(z) = u + i v
            //   u = arcsin(beta)  v = sign(y) log(alpha + sqrt(alpha^2 - 1))
            // where
            //   alpha = (rho + sigma) / 2      beta = (rho - sigma) / 2
            //   rho = sqrt((x + 1)^2 + y^2)    sigma = sqrt((x - 1)^2 + y^2)
            // These formulas appear in DLMF section 4.23. (http://dlmf.nist.gov/4.23), along with the analogous
            //   arccos(w) = arccos(beta) - i sign(y) log(alpha + sqrt(alpha^2 - 1))
            // So alpha and beta together give us arcsin(w) and arccos(w).

            // As written, alpha is not susceptible to cancelation errors, but beta is. To avoid cancelation, note
            //   beta = (rho^2 - sigma^2) / (rho + sigma) / 2 = (2 x) / (rho + sigma) = x / alpha
            // which is not subject to cancelation. Note alpha >= 1 and |beta| <= 1.

            // For alpha ~ 1, the argument of the log is near unity, so we compute (alpha - 1) instead,
            // write the argument as 1 + (alpha - 1) + sqrt((alpha - 1)(alpha + 1)), and use the log1p function
            // to compute the log without loss of accuracy.
            // For beta ~ 1, arccos does not accurately resolve small angles, so we compute the tangent of the angle
            // instead.
            // Hull, Fairgrieve, and Tang derive formulas for (alpha - 1) and beta' = tan(u) that do not suffer
            // from cancelation in these cases.

            // For simplicity, we assume all positive inputs and return all positive outputs. The caller should
            // assign signs appropriate to the desired cut conventions. We return v directly since its magnitude
            // is the same for both arcsin and arccos. Instead of u, we usually return beta and sometimes beta'.
            // If beta' is not computed, it is set to -1; if it is computed, it should be used instead of beta
            // to determine u. Compute u = arcsin(beta) or u = arctan(beta') for arcsin, u = arccos(beta)
            // or arctan(1/beta') for arccos.



            // For x or y large enough to overflow alpha^2, we can simplify our formulas and avoid overflow.
            if ((x > s_asinOverflowThreshold) || (y > s_asinOverflowThreshold))
            {
                b = -1.0f;
                bPrime = x / y;

                float small, big;
                if (x < y)
                {
                    small = x;
                    big = y;
                }
                else
                {
                    small = y;
                    big = x;
                }
                float ratio = small / big;
                v = s_log2 + (float)Math.Log(big) + 0.5f * Log1P(ratio * ratio);
            }
            else
            {
                float r = Hypot((x + 1.0f), y);
                float s = Hypot((x - 1.0f), y);

                float a = (r + s) * 0.5f;
                b = x / a;

                if (b > 0.75)
                {
                    if (x <= 1.0)
                    {
                        float amx = (y * y / (r + (x + 1.0f)) + (s + (1.0f - x))) * 0.5f;
                        bPrime = x / (float)Math.Sqrt((a + x) * amx);
                    }
                    else
                    {
                        // In this case, amx ~ y^2. Since we take the square root of amx, we should
                        // pull y out from under the square root so we don't lose its contribution
                        // when y^2 underflows.
                        float t = (1.0f / (r + (x + 1.0f)) + 1.0f / (s + (x - 1.0f))) * 0.5f;
                        bPrime = x / y / (float)Math.Sqrt((a + x) * t);
                    }
                }
                else
                {
                    bPrime = -1.0f;
                }

                if (a < 1.5)
                {
                    if (x < 1.0)
                    {
                        // This is another case where our expression is proportional to y^2 and
                        // we take its square root, so again we pull out a factor of y from
                        // under the square root.
                        float t = (1.0f / (r + (x + 1.0f)) + 1.0f / (s + (1.0f - x))) * 0.5f;
                        float am1 = y * y * t;
                        v = Log1P(am1 + y * (float)Math.Sqrt(t * (a + 1.0f)));
                    }
                    else
                    {
                        float am1 = (y * y / (r + (x + 1.0f)) + (s + (x - 1.0f))) * 0.5f;
                        v = Log1P(am1 + (float)Math.Sqrt(am1 * (a + 1.0)));
                    }
                }
                else
                {
                    // Because of the test above, we can be sure that a * a will not overflow.
                    v = (float)Math.Log(a + Math.Sqrt((a - 1.0) * (a + 1.0)));
                }
            }
        }

        public static bool IsFinite(Complex32 value) => float.IsFinite(value.m_real) && float.IsFinite(value.m_imaginary);

        public static bool IsInfinity(Complex32 value) => float.IsInfinity(value.m_real) || float.IsInfinity(value.m_imaginary);

        public static bool IsNaN(Complex32 value) => !IsInfinity(value) && !IsFinite(value);

        public static Complex32 Log(Complex32 value)
        {
            return new Complex32((float)Math.Log(Abs(value)), (float)Math.Atan2(value.m_imaginary, value.m_real));
        }

        public static Complex32 Log(Complex32 value, float baseValue)
        {
            return Log(value) / Log(baseValue);
        }

        public static Complex32 Log10(Complex32 value)
        {
            Complex32 tempLog = Log(value);
            return Scale(tempLog, InverseOfLog10);
        }

        public static Complex32 Exp(Complex32 value)
        {
            float expReal = (float)Math.Exp(value.m_real);
            float cosImaginary = expReal * (float)Math.Cos(value.m_imaginary);
            float sinImaginary = expReal * (float)Math.Sin(value.m_imaginary);
            return new Complex32(cosImaginary, sinImaginary);
        }


        public static Complex32 Sqrt(Complex32 value)
        {

            if (value.m_imaginary == 0.0)
            {
                // Handle the trivial case quickly.
                if (value.m_real < 0.0)
                {
                    return new Complex32(0.0f, (float)Math.Sqrt(-value.m_real));
                }
                else
                {
                    return new Complex32((float)Math.Sqrt(value.m_real), 0.0f);
                }
            }
            else
            {

                // One way to compute Sqrt(z) is just to call Pow(z, 0.5), which coverts to polar coordinates
                // (sqrt + atan), halves the phase, and reconverts to cartesian coordinates (cos + sin).
                // Not only is this more expensive than necessary, it also fails to preserve certain expected
                // symmetries, such as that the square root of a pure negative is a pure imaginary, and that the
                // square root of a pure imaginary has exactly equal real and imaginary parts. This all goes
                // back to the fact that Math.PI is not stored with infinite precision, so taking half of Math.PI
                // does not land us on an argument with cosine exactly equal to zero.

                // To find a fast and symmetry-respecting formula for Complex32 square root,
                // note x + i y = \sqrt{a + i b} implies x^2 + 2 i x y - y^2 = a + i b,
                // so x^2 - y^2 = a and 2 x y = b. Cross-substitute and use the quadratic formula to obtain
                //   x = \sqrt{\frac{\sqrt{a^2 + b^2} + a}{2}}  y = \pm \sqrt{\frac{\sqrt{a^2 + b^2} - a}{2}}
                // There is just one complication: depending on the sign on a, either x or y suffers from
                // cancelation when |b| << |a|. We can get aroud this by noting that our formulas imply
                // x^2 y^2 = b^2 / 4, so |x| |y| = |b| / 2. So after computing the one that doesn't suffer
                // from cancelation, we can compute the other with just a division. This is basically just
                // the right way to evaluate the quadratic formula without cancelation.

                // All this reduces our total cost to two sqrts and a few flops, and it respects the desired
                // symmetries. Much better than atan + cos + sin!

                // The signs are a matter of choice of branch cut, which is traditionally taken so x > 0 and sign(y) = sign(b).

                // If the components are too large, Hypot will overflow, even though the subsequent sqrt would
                // make the result representable. To avoid this, we re-scale (by exact powers of 2 for accuracy)
                // when we encounter very large components to avoid intermediate infinities.
                bool rescale = false;
                if ((Math.Abs(value.m_real) >= s_sqrtRescaleThreshold) || (Math.Abs(value.m_imaginary) >= s_sqrtRescaleThreshold))
                {
                    if (float.IsInfinity(value.m_imaginary) && !float.IsNaN(value.m_real))
                    {
                        // We need to handle infinite imaginary parts specially because otherwise
                        // our formulas below produce inf/inf = NaN. The NaN test is necessary
                        // so that we return NaN rather than (+inf,inf) for (NaN,inf).
                        return (new Complex32(float.PositiveInfinity, value.m_imaginary));
                    }
                    else
                    {
                        value.m_real *= 0.25f;
                        value.m_imaginary *= 0.25f;
                        rescale = true;
                    }
                }

                // This is the core of the algorithm. Everything else is special case handling.
                float x, y;
                if (value.m_real >= 0.0f)
                {
                    x = (float)Math.Sqrt((Hypot(value.m_real, value.m_imaginary) + value.m_real) * 0.5);
                    y = value.m_imaginary / (2.0f * x);
                }
                else
                {
                    y = (float)Math.Sqrt((Hypot(value.m_real, value.m_imaginary) - value.m_real) * 0.5);
                    if (value.m_imaginary < 0.0f) y = -y;
                    x = value.m_imaginary / (2.0f * y);
                }

                if (rescale)
                {
                    x *= 2.0f;
                    y *= 2.0f;
                }

                return new Complex32(x, y);

            }

        }

        public static Complex32 Pow(Complex32 value, Complex32 power)
        {
            if (power == Zero)
            {
                return One;
            }

            if (value == Zero)
            {
                return Zero;
            }

            float valueReal = value.m_real;
            float valueImaginary = value.m_imaginary;
            float powerReal = power.m_real;
            float powerImaginary = power.m_imaginary;

            float rho = Abs(value);
            float theta = (float)Math.Atan2(valueImaginary, valueReal);
            float newRho = powerReal * theta + powerImaginary * (float)Math.Log(rho);

            float t = (float)Math.Pow(rho, powerReal) * (float)Math.Pow(Math.E, -powerImaginary * theta);

            return new Complex32(t * (float)Math.Cos(newRho), t * (float)Math.Sin(newRho));
        }

        public static Complex32 Pow(Complex32 value, float power)
        {
            return Pow(value, new Complex32(power, 0));
        }

        private static Complex32 Scale(Complex32 value, float factor)
        {
            float realResult = factor * value.m_real;
            float imaginaryResuilt = factor * value.m_imaginary;
            return new Complex32(realResult, imaginaryResuilt);
        }

        public static implicit operator Complex32(short value)
        {
            return new Complex32(value, 0.0f);
        }

        public static implicit operator Complex32(int value)
        {
            return new Complex32(value, 0.0f);
        }

        public static implicit operator Complex32(long value)
        {
            return new Complex32(value, 0.0f);
        }

        [CLSCompliant(false)]
        public static implicit operator Complex32(ushort value)
        {
            return new Complex32(value, 0.0f);
        }

        [CLSCompliant(false)]
        public static implicit operator Complex32(uint value)
        {
            return new Complex32(value, 0.0f);
        }

        [CLSCompliant(false)]
        public static implicit operator Complex32(ulong value)
        {
            return new Complex32(value, 0.0f);
        }

        [CLSCompliant(false)]
        public static implicit operator Complex32(sbyte value)
        {
            return new Complex32(value, 0.0f);
        }

        public static implicit operator Complex32(byte value)
        {
            return new Complex32(value, 0.0f);
        }

        public static implicit operator Complex32(float value)
        {
            return new Complex32(value, 0.0f);
        }

        


        public static explicit operator Complex32(decimal value)
        {
            return new Complex32((float)value, 0.0f);
        }
    }

}
