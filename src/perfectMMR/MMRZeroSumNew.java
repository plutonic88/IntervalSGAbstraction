package perfectMMR;

import java.io.IOException;
import java.util.Vector;

import lpWrapper.LPSolverException;

public class MMRZeroSumNew {
	int nTargets;
	int nRes;
//	double[] def_payoff; //only used for non-zero sum
	
	double[] advPayoffLB;
	double[] advPayoffUB;
	
	// Output
	double[] optDefCov;
	public double[] optAdvPayoff;
	double[] optAdvCov;
	double optRegret;
	
	RelaxedMMRZeroSumNew relaxedMMR;
	MRZeroSumNew mr;
	Vector<double[]> advPayoffs;
	Vector<Double> defUtilities;
	// Just for printing
	public double runtime;
	public MMRZeroSumNew(int nTargets, int nRes, double[] advPayoffLB, double[] advPayoffUB)
	{
		this.nTargets = nTargets;
		this.nRes = nRes;
		this.advPayoffLB = advPayoffLB;
		this.advPayoffUB = advPayoffUB;
		advPayoffs = new Vector<double[]>();
		defUtilities = new Vector<Double>();
		
		
		relaxedMMR = new RelaxedMMRZeroSumNew(advPayoffs, defUtilities, nTargets, nRes);
		mr = new MRZeroSumNew(nTargets, nRes, null, advPayoffLB, advPayoffUB);
	}
	public void loadProblem() throws IOException{
		relaxedMMR.loadProblem();
	}
	public void solve()
	{
		long startTime = System.currentTimeMillis();
		double ub = Double.POSITIVE_INFINITY;
		double lb = 0.0;
		
		// Generate some samples
		for(int i = 0; i < nTargets; i++)
		{
			double[] advPayoff = new double[2 * nTargets];
			advPayoff[i] = advPayoffUB[i]; 
			advPayoff[i + nTargets] = advPayoffUB[i + nTargets]; 
			for(int n = 0; n < nTargets; n++)
				if(n != i)
				{
					advPayoff[n] = advPayoffLB[n]; 
					advPayoff[n + nTargets] = advPayoffLB[n + nTargets]; 
				}
			advPayoffs.add(advPayoff);
			relaxedMMR.update();
		}
		
		// Incremental Payoff Generation
		// Initialize lower bound and upper bound
//		for(int i = 0; i < nTargets; i++)
//		{
//			lb = Math.min(lb, advPayoffLB[nTargets + i]);
//			ub = Math.max(ub, advPayoffUB[i]);
//		}
//		lb *= 2;
//		ub *= 2;
		int numIter = 400;
		int count = 0;
		while(count < numIter && (ub - lb) > 1e-3)
		{
			count++;
			System.out.println(count + "\t" + ub + "\t" + lb);
			// Compute sub-minimaxregret
			try {
				relaxedMMR.solve();
			} catch (LPSolverException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			lb = relaxedMMR.getOptRegret();
			double[] xStar = relaxedMMR.getOptDefCov();
			long endTime = System.currentTimeMillis();
//			for(int i = 0; i < nTargets; i++)
//				System.out.print(xStar[i] + "\t");
//			System.out.println();
			runtime = (endTime - startTime) / 1000.0;
			
			mr.setXStar(xStar);
			mr.solve();
			ub = mr.maxRegret;
			
			if(count == numIter - 1 || ub - lb < 1e-3)
			{
				System.out.println(count + "\t" + ub + "\t" + lb);
				optDefCov = xStar;
				optAdvPayoff = mr.optAdvPayoff;
				optAdvCov = mr.optDefCov;
				optRegret = ub;
				
				System.out.println("Runtime: " + runtime);
			}
			
			advPayoffs.add(mr.optAdvPayoff);
			relaxedMMR.update();
		}
		
	}
	public void end()
	{
		relaxedMMR.end();
		advPayoffs.clear();
		defUtilities.clear();
	}
	public double[] getOptDefCov()
	{
		return optDefCov;
	}
	public double getOptRegret()
	{
		return optRegret;
	}
	public double[] getOptAdvPayoff()
	{
		return optAdvPayoff;
	}
	public double[] getOptAdvCov()
	{
		return optAdvCov;
	}
	
}
