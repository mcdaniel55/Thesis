using System;
using ThesisOptNumericalTest.Optimization;
using EMSOptimizationOverhaul;
using System.Text;

namespace ThesisOptNumericalTest
{
    class EMSPlanRegion : Region, IEquatable<EMSPlanRegion>
    {
        public int[] CurrentPlanFullAmbs { get; private set; }
        public int[] CurrentPlanPartAmbs { get; private set; }
        public int TargetAmbulanceCount { get; private set; }
        public int AmbulancesAssigned { get; private set; }

        public EMSPlanRegion(int[] CurrentPlanFullAmbs, int[] CurrentPlanPartAmbs, int TargetAmbulanceCount, Random rand) : base(rand)
        {
            this.CurrentPlanFullAmbs = CurrentPlanFullAmbs;
            this.CurrentPlanPartAmbs = CurrentPlanPartAmbs;
            this.TargetAmbulanceCount = TargetAmbulanceCount;
            AmbulancesAssigned = 0;
            for (int i = 0; i < CurrentPlanFullAmbs.Length; i++) { AmbulancesAssigned += CurrentPlanFullAmbs[i]; }
        }

        public override Region[] Branch()
        {
            Region[] regions = new Region[2 * CurrentPlanFullAmbs.Length];
            for (int i = 0; i < CurrentPlanFullAmbs.Length; i++)
            {
                // Add a plan with an extra full time ambulance in each possible way
                int[] newPlanFull = (int[])CurrentPlanFullAmbs.Clone();
                newPlanFull[i]++;
                regions[i] = new EMSPlanRegion(newPlanFull, CurrentPlanPartAmbs, TargetAmbulanceCount, m_rand);

                // Do the same for part time
                int[] newPlanPart = (int[])CurrentPlanPartAmbs.Clone();
                newPlanPart[i]++;
                regions[i + CurrentPlanFullAmbs.Length] = new EMSPlanRegion(CurrentPlanFullAmbs, newPlanPart, TargetAmbulanceCount, m_rand);
            }
            
            return regions;
        }

        protected override double SampleElement()
        {
            for (int attempt = 0; attempt < 5; attempt++)
            {
                try
                {
                    int[] planToTestFull = (int[])CurrentPlanFullAmbs.Clone();
                    int[] planToTestPart = (int[])CurrentPlanPartAmbs.Clone();
                    // Fill out the rest of the plan randomly
                    for (int i = 0; i < TargetAmbulanceCount - AmbulancesAssigned; i++)
                    {
                        int idx = m_rand.Next(CurrentPlanFullAmbs.Length);
                        if (m_rand.NextDouble() < 0.5)
                        { planToTestFull[idx]++; }
                        else { planToTestPart[idx]++; }
                    }

                    // Test the plan
                    Simulation sim = new Simulation(
                        startTime: new DateTime(2016, 1, 1, 0, 0, 0), // These don't have to change when you use only 2016 or only 2017 calls
                        endTime: new DateTime(2018, 1, 1, 0, 0, 0),
                        fullAmbsToSpawn: planToTestFull,
                        partAmbsToSpawn: planToTestPart,
                        centralized: true,
                        speedMPH: 24f);
                    sim.Run();

                    // Compute the score
                    return sim.MeanResponseTime * sim.Cost;
                }
                catch (Exception e) { Logging.Log("SampleElement() ignored an error: " + e.Message + e.StackTrace); }
            }
            throw new Exception("SampleElement: Out of attempts");
        }

        public override bool Equals(object obj)
        {
            if (obj is EMSPlanRegion) return Equals((EMSPlanRegion) obj);
            return base.Equals(obj);
        }

        public bool Equals(EMSPlanRegion other)
        {
            for (int i = 0; i < CurrentPlanFullAmbs.Length; i++)
            {
                if (CurrentPlanFullAmbs[i] != other.CurrentPlanFullAmbs[i] 
                    || CurrentPlanPartAmbs[i] != other.CurrentPlanPartAmbs[i]) { return false; }
            }
            return true;
        }

        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();
            builder.AppendLine();
            for (int i = 0; i < Station.List.Count; i++)
            {
                builder.AppendLine($"{Station.List[i].Name.Truncate(6)} Full {CurrentPlanFullAmbs[i]} Part {CurrentPlanPartAmbs[i]}");
            }
            builder.AppendLine($"Mean: {SampleMean}");
            builder.AppendLine($"StdDev: {SampleStdDev}");
            builder.AppendLine($"Best Observed: {BestObservation}");

            return builder.ToString();
        }
    }
}
