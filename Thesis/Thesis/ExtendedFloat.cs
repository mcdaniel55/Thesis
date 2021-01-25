using System;

namespace Thesis
{
    /// <summary>
    /// Provides support for multiplying floating point numbers with very large or small exponents together without clamping to 0 or infinity due to a limited number of bits in the exponent.
    /// </summary>
    struct ExtendedFloat : IComparable<ExtendedFloat>
    {
        double value;
        long exponentOffset;

        public unsafe ExtendedFloat(double value, long exponentOffset = 0)
        {
            this.value = value;
            this.exponentOffset = exponentOffset;
            // Move the exponent part from value into the offset
            ulong vBits = *(ulong*)&value;
            ulong vExp = (vBits >> 52) & 0x7ffL;
            if (vExp != 0x3ffL) // If non-zero exponent
            {
                // Move the exponent into the offset
                this.exponentOffset += (long)vExp - 0x3ffL;
                // Clear the exponent bits from c
                vBits &= 0x800fffffffffffffL;
                // Set the exponent to zero
                vBits |= 0x3ffL << 52;
                this.value = *(double*)&vBits;
            }
        }

        public static unsafe ExtendedFloat operator *(ExtendedFloat a, ExtendedFloat b)
        {
            long cOffset = a.exponentOffset + b.exponentOffset;
            double cVal = a.value * b.value;
            ulong cBits = *(ulong*)&cVal;
            ulong cExp = (cBits >> 52) & 0x7ffL;
            if (cExp != 0x3ffL)
            {
                cOffset += (long)cExp - 0x3ffL;
                cBits &= 0x800fffffffffffffL;
                cBits |= 0x3ffL << 52;
                cVal = *(double*)&cBits;
            }
            return new ExtendedFloat(cVal, cOffset);
        }

        public static explicit operator double(ExtendedFloat val) => val.ToDouble();
        public static implicit operator ExtendedFloat(double val) => new ExtendedFloat(val);

        public unsafe double ToDouble()
        {
            if (value == 0) return 0;
            if (exponentOffset < -1023) return 0;
            if (1024 < exponentOffset) return value > 0 ? double.MaxValue : double.MinValue;
            double val = value;
            ulong asULong = *(ulong*)&val;
            asULong &= 0x800fffffffffffffL;
            asULong |= (ulong)(exponentOffset + 0x3ffL) << 52;
            val = *(double*)&asULong;
            return val;
        }

        public int CompareTo(ExtendedFloat other)
        {
            // Handle zero-value cases
            if (value == 0) return other.value.CompareTo(0);
            if (other.value == 0) return value.CompareTo(0);
            // If exponents are the same, they can be compared as doubles
            if (exponentOffset == other.exponentOffset) return value.CompareTo(other.value);
            // Cases with opposing signs do not require checking the exponents
            if (0 < value && other.value < 0) return 1;
            if (value < 0 && 0 < other.value) return -1;
            // Now they are nonzero and have the same sign, but different exponents
            return exponentOffset > other.exponentOffset ^ value < 0 ? 1 : -1;
        }

        public static bool operator <(ExtendedFloat a, ExtendedFloat b) => a.CompareTo(b) < 0;
        public static bool operator >(ExtendedFloat a, ExtendedFloat b) => a.CompareTo(b) > 0;
        public static bool operator <=(ExtendedFloat a, ExtendedFloat b) => a.CompareTo(b) < 1;
        public static bool operator >=(ExtendedFloat a, ExtendedFloat b) => a.CompareTo(b) > -1;
    }
}
