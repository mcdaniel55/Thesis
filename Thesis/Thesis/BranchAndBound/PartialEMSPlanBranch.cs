using System;
using System.Text;
using EMSOptimizationOverhaul;

namespace Thesis.BranchAndBound
{
    public class PartialEMSPlanBranch : Branch
    {
        public int[] FullAmbs { get; private set; }
        public int[] PartAmbs { get; private set; }
        public int TargetAmbulanceCount { get; private set; }
        public int CurrentAmbulancesCount { get; private set; }

        public PartialEMSPlanBranch(int[] FullAmbs, int[] PartAmbs, int TargetAmbulanceCount, Random rand = null) : base(rand)
        {
            this.FullAmbs = FullAmbs;
            this.PartAmbs = PartAmbs;
            this.TargetAmbulanceCount = TargetAmbulanceCount;
            CurrentAmbulancesCount = 0;
            for (int i = 0; i < FullAmbs.Length; i++) { CurrentAmbulancesCount += FullAmbs[i] + PartAmbs[i]; }
        }

        public override Branch[] GetBranches()
        {
            Branch[] regions = new Branch[2 * FullAmbs.Length];
            for (int i = 0; i < FullAmbs.Length; i++)
            {
                // Add a plan with an extra full time ambulance in each possible way
                int[] newPlanFull = (int[])FullAmbs.Clone();
                newPlanFull[i]++;
                regions[i] = new PartialEMSPlanBranch(newPlanFull, PartAmbs, TargetAmbulanceCount, rand);

                // Do the same for part time
                int[] newPlanPart = (int[])PartAmbs.Clone();
                newPlanPart[i]++;
                regions[i + FullAmbs.Length] = new PartialEMSPlanBranch(FullAmbs, newPlanPart, TargetAmbulanceCount, rand);
            }

            return regions;
        }

        protected override object RandomElement() => GetRandomElement();

        /// <summary>
        /// Gets a random plan filled out from this partial plan
        /// </summary>
        /// <returns> A tuple with the full time amb int[] and the part time amb int[], in that order </returns>
        public new Tuple<int[],int[]> GetRandomElement()
        {
            int[] planToTestFull = (int[])FullAmbs.Clone();
            int[] planToTestPart = (int[])PartAmbs.Clone();
            // Fill out the rest of the plan randomly
            for (int i = 0; i < TargetAmbulanceCount - CurrentAmbulancesCount; i++)
            {
                int idx = rand.Next(FullAmbs.Length);
                if (rand.NextDouble() < 0.5)
                { planToTestFull[idx]++; }
                else { planToTestPart[idx]++; }
            }
            return new Tuple<int[], int[]>(planToTestFull, planToTestPart);
        }

        public void GetRandomElementNonAlloc(Tuple<int[],int[]> existing)
        {
            // Copy in the current plan
            for (int i = 0; i < FullAmbs.Length; i++)
            {
                existing.Item1[i] = FullAmbs[i];
                existing.Item2[i] = PartAmbs[i];
            }
            // Fill out the rest of the plan randomly
            for (int i = 0; i < TargetAmbulanceCount - CurrentAmbulancesCount; i++)
            {
                int idx = rand.Next(FullAmbs.Length);
                if (rand.NextDouble() < 0.5)
                { existing.Item1[idx]++; }
                else { existing.Item2[idx]++; }
            }
        }

        // --- Equality Checking ---

        public override bool Equals(object obj)
        {
            if (obj is PartialEMSPlanBranch) return Equals((PartialEMSPlanBranch)obj);
            return base.Equals(obj);
        }

        // Two partial plans are equal if they agree on how many ambulances have been assigned to each place/time slot
        public bool Equals(PartialEMSPlanBranch other)
        {
            for (int i = 0; i < FullAmbs.Length; i++)
            {
                if (FullAmbs[i] != other.FullAmbs[i]
                    || PartAmbs[i] != other.PartAmbs[i]) { return false; }
            }
            return true;
        }

        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();
            builder.AppendLine();
            for (int i = 0; i < Station.List.Count; i++)
            {
                builder.AppendLine($"{Station.List[i].Name.Truncate(6)} Full {FullAmbs[i]} Part {PartAmbs[i]}");
            }
            return builder.ToString();
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(FullAmbs, PartAmbs);
        }
    }
}
