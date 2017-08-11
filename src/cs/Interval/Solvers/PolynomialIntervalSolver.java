package cs.Interval.Solvers;// tambe.usc.DOBSS.gameRepresentation.*;





/**
 * @author 
 * 
 * 
 * A polynomial-time solver for games with interval attacker payoffs
 */
public class PolynomialIntervalSolver {

	public double tolerance = 0.0001;

	private int nTargets;
	private double nResources;
	private double[] coverage;
	private double[] lastGoodCoverage;

	IntervalPayoffs[] payoffs;

	public PolynomialIntervalSolver() {
	}

	//---------- Algorithm 1 -------  ISEGS

	public double[] solve(IntervalSecurityGame isg, int numOrigRes, int[] clustersizes) {
		init(isg);

		double maxPayoff = Double.NEGATIVE_INFINITY;
		double minPayoff = Double.POSITIVE_INFINITY;
		for (int i = 0; i < nTargets; i++) {
			if (payoffs[i].getDefenderUncoveredPayoff() < minPayoff) {
				minPayoff = payoffs[i].getDefenderUncoveredPayoff();
			}
			if (payoffs[i].getDefenderCoveredPayoff() > maxPayoff) {
				maxPayoff = payoffs[i].getDefenderCoveredPayoff();
			}
		}

		//    for(int i = 0; i < nTargets; i++) {
		//      System.out.println("Defender: " + payoffs[i].getDefenderUncoveredPayoff() +
		//                         " Attacker min: " + payoffs[i].getAttackerUncoveredPayoffMin() +
		//                         " Attacker max: " + payoffs[i].getAttackerUncoveredPayoffMax());
		//    }

		while(maxPayoff - minPayoff > tolerance)
		{
			double mid = (maxPayoff + minPayoff) / 2;
			//      System.out.println("MID: " + mid);

			if (determineFeasibility(mid, (double)nResources, numOrigRes, clustersizes)) 
			{
				minPayoff = mid;
				System.arraycopy(coverage, 0, lastGoodCoverage, 0, nTargets);
			} 
			else 
			{
				maxPayoff = mid;
			}
		}

		//    System.out.println("FINAL Coverage: " + Arrays.toString(lastGoodCoverage));
		return lastGoodCoverage;
	}




	//--------- Algorithm 2------------- Feasibility Check

	// returns true if there is a feasible coverage strategy that guarantees the defender payoff,
	// for any realization of the attackers payoffs, and false otherwise


	private boolean determineFeasibility(double defenderPayoff, double maxCoverage,
			int numOrigRes, int[] clustersizes) {
		double[] minCoverageDefenderConstraint = calculateDefenderConstraint(defenderPayoff);

		//    System.out.println("Defender cov constraints: " + Arrays.toString(minCoverageDefenderConstraint));

		for (int target = 0; target < nTargets; target++)
		{
			double requiredCoverage = minCoverageDefenderConstraint[target];
			coverage[target] = requiredCoverage;
			requiredCoverage = minCoverageDefenderConstraint[target]*clustersizes[target];
			
			//      System.out.println("Target " + target + " coverage: " + requiredCoverage);

			// this defines the minimum defender payoff for attacking the chosen target
			//double minAttackerPayoff = (1-requiredCoverage) *payoffs[target].getAttackerUncoveredPayoffMin();
			double minAttackerPayoff = requiredCoverage * payoffs[target].getAttackerCoveredPayoffMin() 
					+
					(1-requiredCoverage) *payoffs[target].getAttackerUncoveredPayoffMin();

			minAttackerPayoff -= 0.0000001; // small epsilon

			//      System.out.println("Min attacker payoff: " + minAttackerPayoff);

			for (int i = 0; i < nTargets; i++) 
			{
				if (i == target) continue;
				//        double coverageNotInAttackSet = Math.max(0, 1 - (minAttackerPayoff / payoffs[i].getAttackerUncoveredPayoffMax()));
				//        double coverageMinConstraint = Math.max(0, 1 - (minAttackerPayoff / payoffs[i].getAttackerUncoveredPayoffMin()));

				double coverageNotInAttackSet = Math.max(0,  (minAttackerPayoff - payoffs[i].getAttackerUncoveredPayoffMax())
						/ 
						(payoffs[i].getAttackerCoveredPayoffMax()-payoffs[i].getAttackerUncoveredPayoffMax()));

				double coverageMinConstraint = Math.max(0,  (minAttackerPayoff - payoffs[i].getAttackerUncoveredPayoffMin())
						/ 
						(payoffs[i].getAttackerCoveredPayoffMin()-payoffs[i].getAttackerUncoveredPayoffMin()));


				double minCoverageInAttackSet = Math.max(coverageMinConstraint, minCoverageDefenderConstraint[i]);
				double min = Math.min(minCoverageInAttackSet, coverageNotInAttackSet);
				
				coverage[i] = min;
				
				requiredCoverage = requiredCoverage +  coverage[i]*clustersizes[i];
				
				 //System.out.println("Updating target "+ i + " coverage "+ coverage[i] + " requirdcov "+ requiredCoverage);

				//        System.out.println("Target: " + i +
				//                           " def constraint: " + minCoverageDefenderConstraint[i] +
				//                           " min constraint: " + coverageMinConstraint +
				//                           " max constraint: " + coverageNotInAttackSet +
				//                           " min: " + min);
			}
			//      System.out.println("Total required: " + requiredCoverage);

			if (requiredCoverage <= numOrigRes) 
				return true;
		}
		return false;
	}

	// returns the minimum coverage necessary for each target to meet the defender payoff constraint
	public double[] calculateDefenderConstraint(double defenderPayoff) {
		double[] tmp = new double[nTargets];
		for (int i = 0; i < nTargets; i++) {
			//tmp[i] = Math.max(0, 1 - (defenderPayoff / payoffs[i].getDefenderUncoveredPayoff()));
			tmp[i] = Math.max(0,((defenderPayoff - payoffs[i].getDefenderUncoveredPayoff() )
					/
					(payoffs[i].getDefenderCoveredPayoff()- payoffs[i].getDefenderUncoveredPayoff())));
		}
		return tmp;
	}

	public void init(IntervalSecurityGame isg) {
		nTargets = isg.getNTargets();
		nResources = isg.getNumDefenders();
		payoffs = isg.getPayoffs();
		coverage = new double[nTargets];
		lastGoodCoverage = new double[nTargets];
	}

}