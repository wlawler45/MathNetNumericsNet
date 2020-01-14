using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MathNetNumerics
{
    public static class Buffer
    {
        internal static void BlockCopy(byte[] recbuf, int v1, byte[] newbuf, int v2, int v3)
        {
            Array.Copy(recbuf, v1, newbuf, v2, v3);
        }

        internal static int ByteLength(Array a)
        {
            TypeCode t = Type.GetTypeCode(a.GetType().GetElementType());
            switch (t)
            {
                case TypeCode.Double:
                    return 8 * a.Length;
                case TypeCode.Single:
                    return 4 * a.Length;
                case TypeCode.Byte:
                case TypeCode.SByte:
                    return a.Length;
                case TypeCode.UInt16:
                case TypeCode.Int16:
                    return 2 * a.Length;
                case TypeCode.UInt32:
                case TypeCode.Int32:
                    return 4 * a.Length;
                case TypeCode.UInt64:
                case TypeCode.Int64:
                    return 8 * a.Length;
                case TypeCode.Boolean:
                    return a.Length;
                default:
                    throw new ArgumentException("Invalid array type");
            }
        }

        internal static void BlockCopy(Array a, int v, byte[] membuf, int position, int bl)
        {
            var b = new BinaryWriter(new MemoryStream(membuf, position, bl));
            TypeCode t = Type.GetTypeCode(a.GetType().GetElementType());
            switch (t)
            {
                case TypeCode.Double:
                    {
                        var a1 = (double[])a;
                        for (int i = 0; i < bl / 8; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.Single:
                    {
                        var a1 = (float[])a;
                        for (int i = 0; i < bl / 4; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.Byte:
                    {
                        var a1 = (byte[])a;
                        for (int i = 0; i < bl; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.SByte:
                    {
                        var a1 = (sbyte[])a;
                        for (int i = 0; i < bl; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.UInt16:
                    {
                        var a1 = (ushort[])a;
                        for (int i = 0; i < bl / 2; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.Int16:
                    {
                        var a1 = (short[])a;
                        for (int i = 0; i < bl / 2; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.UInt32:
                    {
                        var a1 = (uint[])a;
                        for (int i = 0; i < bl / 4; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.Int32:
                    {
                        var a1 = (int[])a;
                        for (int i = 0; i < bl / 4; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.UInt64:
                    {
                        var a1 = (ulong[])a;
                        for (int i = 0; i < bl / 8; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.Int64:
                    {
                        var a1 = (long[])a;
                        for (int i = 0; i < bl / 8; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                case TypeCode.Boolean:
                    {
                        var a1 = (bool[])a;
                        for (int i = 0; i < bl; i++)
                            b.Write(a1[i + v]);
                        break;
                    }
                default:
                    throw new ArgumentException("Invalid array type");
            }
        }

        internal static void BlockCopy(byte[] membuf, int position, Array a, int v, int bl)
        {
            var b = new BinaryReader(new MemoryStream(membuf, position, bl));
            TypeCode t = Type.GetTypeCode(a.GetType().GetElementType());
            switch (t)
            {
                case TypeCode.Double:
                    {
                        var a1 = (double[])a;
                        for (int i = 0; i < bl / 8; i++)
                            a1[i + v] = b.ReadDouble();
                        break;
                    }
                case TypeCode.Single:
                    {
                        var a1 = (float[])a;
                        for (int i = 0; i < bl / 4; i++)
                            a1[i + v] = b.ReadSingle();
                        break;
                    }
                case TypeCode.Byte:
                    {
                        var a1 = (byte[])a;
                        for (int i = 0; i < bl; i++)
                            a1[i + v] = b.ReadByte();
                        break;
                    }
                case TypeCode.SByte:
                    {
                        var a1 = (sbyte[])a;
                        for (int i = 0; i < bl; i++)
                            a1[i + v] = b.ReadSByte();
                        break;
                    }
                case TypeCode.UInt16:
                    {
                        var a1 = (ushort[])a;
                        for (int i = 0; i < bl / 2; i++)
                            a1[i + v] = b.ReadUInt16();
                        break;
                    }
                case TypeCode.Int16:
                    {
                        var a1 = (short[])a;
                        for (int i = 0; i < bl / 2; i++)
                            a1[i + v] = b.ReadInt16();
                        break;
                    }
                case TypeCode.UInt32:
                    {
                        var a1 = (uint[])a;
                        for (int i = 0; i < bl / 4; i++)
                            a1[i + v] = b.ReadUInt32();
                        break;
                    }
                case TypeCode.Int32:
                    {
                        var a1 = (int[])a;
                        for (int i = 0; i < bl / 4; i++)
                            a1[i + v] = b.ReadInt32();
                        break;
                    }
                case TypeCode.UInt64:
                    {
                        var a1 = (ulong[])a;
                        for (int i = 0; i < bl / 8; i++)
                            a1[i + v] = b.ReadUInt64();
                        break;
                    }
                case TypeCode.Int64:
                    {
                        var a1 = (long[])a;
                        for (int i = 0; i < bl / 8; i++)
                            a1[i + v] = b.ReadInt64();
                        break;
                    }
                case TypeCode.Boolean:
                    {
                        var a1 = (bool[])a;
                        for (int i = 0; i < bl; i++)
                            a1[i + v] = b.ReadBoolean();
                        break;
                    }
                default:
                    throw new ArgumentException("Invalid array type");
            }
        }
    }
}