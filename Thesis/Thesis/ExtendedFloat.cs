using System;

namespace Thesis
{
    struct ExtendedFloat : IComparable<ExtendedFloat>
    {
        double value;
        long exponentOffset;

        public ExtendedFloat(double value, long exponentOffset)
        {
            this.value = value;
            this.exponentOffset = exponentOffset;
        }

        public static unsafe ExtendedFloat operator *(ExtendedFloat a, ExtendedFloat b)
        {
            long cOffset = a.exponentOffset + b.exponentOffset;
            double cVal = a.value * b.value;
            ulong cBits = *(ulong*)&cVal;
            ulong cExp = (cBits >> 52) & 0x7ffL;
            if (cExp != 0x3ffL) // If non-zero exponent
            {
                // Move the exponent into the offset
                cOffset += (long)cExp - 0x3ffL;
                // Clear the exponent bits from c
                cBits &= 0x800fffffffffffffL;
                cVal = *(double*)&cBits;
            }
            return new ExtendedFloat(cVal, cOffset);
        }

        public unsafe double AsDouble()
        {
            if (value == 0) return 0;
            if (exponentOffset < -1023) return 0;
            if (1024 < exponentOffset) return value > 0 ? double.MaxValue : double.MinValue;
            double val = value;
            long asLong = *(long*)&val;
            asLong |= (exponentOffset + 0x3ffL) << 52;
            val = *(double*)&asLong;
            return val;
        }

        public int CompareTo(ExtendedFloat other)
        {
            if (exponentOffset > other.exponentOffset) return value > 0 ? 1 : -1;
            if (exponentOffset < other.exponentOffset) return value > 0 ? -1 : 1;
            return value.CompareTo(other.value);
        }

        public static bool operator <(ExtendedFloat a, ExtendedFloat b)
        {
            return a.CompareTo(b) < 0;
        }

        public static bool operator >(ExtendedFloat a, ExtendedFloat b)
        {
            return a.CompareTo(b) > 0;
        }
    }
}
