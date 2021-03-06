package perfectMMR;

import java.io.IOException;
import java.util.Vector;

public class MMRZeroSumNew2 {
	int nTargets;
	double nRes;
//	double[] def_payoff; //only used for non-zero sum
	
	double[] advPayoffLB;
	double[] advPayoffUB;
	double[] weight;
	
	// Output
	double[] optDefCov;
	public double[] optAdvPayoff;
	double[] optAdvCov;
	double optRegret;
	
//	RelaxedMMRZeroSumNew relaxedMMR;
	SubMinimaxRegret relaxedMMR;
	MRZeroSumNew mr;
	Vector<double[]> advPayoffs;
	Vector<double[]> defPayoffs;
	Vector<double[]> x;
//	Vector<Double> defUtilities;
	// Just for printing
	public double runtime;
	public MMRZeroSumNew2(int nTargets, double nRes, double[] advPayoffLB, double[] advPayoffUB)
	{
		this.nTargets = nTargets;
		this.nRes = nRes;
		this.advPayoffLB = advPayoffLB;
		this.advPayoffUB = advPayoffUB;
		
		advPayoffs = new Vector<double[]>();
		defPayoffs = new Vector<double[]>();
		x = new Vector<double[]>();
		
		
//		relaxedMMR = new RelaxedMMRZeroSumNew(advPayoffs, defUtilities, nTargets, nRes);
		relaxedMMR = new SubMinimaxRegret(x, advPayoffs, defPayoffs, nTargets, nRes);
		mr = new MRZeroSumNew(nTargets, nRes, null, advPayoffLB, advPayoffUB);
		weight = new double[nTargets];
	}
	public void setWeight(double[] weight)
	{
		for(int i = 0; i < nTargets; i++)
			this.weight[i] = weight[i];
		mr.setWeight(weight);
		relaxedMMR.setWeight(weight);
	}
	public void loadProblem() throws IOException{
//		relaxedMMR.loadProblem();
	}
	public void solve()
	{
		long startTime = System.currentTimeMillis();
		double ub = Double.POSITIVE_INFINITY;
		double lb = 0.0;
		
		// Generate some samples
//		System.out.println("Samples: ");
		for(int i = 0; i < nTargets; i++)
		{
			double[] advPayoff = new double[2 * nTargets];
			double[] defPayoff = new double[2 * nTargets];
			advPayoff[i] = advPayoffUB[i]; 
			advPayoff[i + nTargets] = advPayoffUB[i + nTargets]; 
			defPayoff[i] = -advPayoff[i + nTargets];
			defPayoff[i + nTargets] = -advPayoff[i];
			for(int n = 0; n < nTargets; n++)
				if(n != i)
				{
					advPayoff[n] = advPayoffLB[n]; 
					advPayoff[n + nTargets] = advPayoffLB[n + nTargets]; 
					defPayoff[n] = -advPayoff[n + nTargets];
					defPayoff[n + nTargets] = -advPayoff[n];
				}
			
			double[] defCov = ORIGAMI(nTargets, nRes, advPayoff);
			x.add(defCov);
			advPayoffs.add(advPayoff);
			defPayoffs.add(defPayoff);
			double defEU = Double.POSITIVE_INFINITY;
			for(int n = 0; n < nTargets; n++)
			{
				double temp = defCov[n] * (defPayoff[n] - defPayoff[n + nTargets]) + defPayoff[n + nTargets];
				if(defEU > temp)
					defEU = temp;
			}
//			System.out.println(defEU);
		}
		
		// Incremental Payoff Generation
		// Initialize lower bound and upper bound
		double lbRelaxedMMR = Double.POSITIVE_INFINITY;
		double ubRelaxedMMR = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < nTargets; i++)
		{
			lbRelaxedMMR = Math.min(lbRelaxedMMR, advPayoffLB[nTargets + i]);
			ubRelaxedMMR = Math.max(ubRelaxedMMR, advPayoffUB[i]);
		}
		lbRelaxedMMR *= 2;
		ubRelaxedMMR = 0;
		int numIter = 100;
		int count = 0;
		while(count < numIter && (ub - lb) > 1e-3)
		{
			count++;
//			System.out.println(count + "\t" + ub + "\t" + lb);
			// Compute sub-minimaxregret
			relaxedMMR.binarySearch(lbRelaxedMMR, ubRelaxedMMR);
			ubRelaxedMMR = -relaxedMMR.obj_value;
			lb = relaxedMMR.obj_value;
			double[] xStar = relaxedMMR.def_strategy;
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
				long endTime = System.currentTimeMillis();
				runtime = (endTime - startTime) / 1000.0;
				System.out.println("Runtime: " + runtime);
			}
			double[] advPayoff = new double[2 * nTargets];
			for(int i = 0; i < 2 * nTargets; i++)
				{
					advPayoff[i] = mr.optAdvPayoff[i];
				}
			advPayoffs.add(advPayoff);
			
			double[] defCov = ORIGAMI(nTargets, nRes, mr.optAdvPayoff);
			x.add(defCov);
			double[] defPayoff = new double[2 * nTargets];
			for(int n = 0; n < nTargets; n++)
			{
				defPayoff[n] = -mr.optAdvPayoff[n + nTargets];
				defPayoff[n + nTargets] = -mr.optAdvPayoff[n];
			}
			defPayoffs.add(defPayoff);
		}
		
	}
	public void end()
	{
//		relaxedMMR.end();
		advPayoffs.clear();
		defPayoffs.clear();
		x.clear();
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
	public double getRuntime()
	{
		return runtime;
	}
	public double[] ORIGAMI(int nTargets, double nRes, double[] adv_payoff)
	{
		double[] min_cov = new double[nTargets];
		for(int i = 0; i < nTargets; i++)
			min_cov[i] = 0.0;
		int[] idx = new int[nTargets];
		for(int i = 0; i < nTargets; i++)
			idx[i] = i;
		double[] values = new double[nTargets];
		for(int i = 0; i < nTargets; i++)
			values[i] = adv_payoff[idx[i]];
		QuickSort my_sort = new QuickSort(values, idx, nTargets);
		my_sort.get_sort_id();
		int[] sort_idx = my_sort.sort_id;
		for(int i = 0; i < nTargets / 2; i++)
		{
			int temp = sort_idx[i];
			sort_idx[i] = sort_idx[nTargets - 1 - i];
			sort_idx[nTargets - 1 - i] = temp;
		}
		double r = adv_payoff[sort_idx[0]];
		double res = nRes;
		boolean is_stop = false;
		int final_idx = 1;
		for(int i = 1; i < nTargets; i++)
		{
			final_idx = i;
			double[] dis = new double[i];
			for(int j = 0; j < i; j++)
			{
				dis[j] = (adv_payoff[sort_idx[i]] - r) /
						(adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]);
//				res = res - dis[j];
				res = res - dis[j] * weight[sort_idx[j]];
				min_cov[sort_idx[j]] += dis[j];
				if(res < 0 || min_cov[sort_idx[j]] > 1) 
				{
					for(int k = 0; k <= j; k++)
					{
//						res = res + dis[k];
						res = res + dis[k] * weight[sort_idx[k]];
						min_cov[sort_idx[k]] -= dis[k];
					}
					is_stop = true;
					break;
				}
			}
			if(!is_stop)
				r = adv_payoff[sort_idx[i]];
			else break;
		}
		for(int j = 0; j < nTargets; j++)
			if(min_cov[j] > 1 || min_cov[j] < 0)
				System.out.println("Oops");
		double delta_res = 0.0;
		double sum_res = 0.0;
		if(!is_stop && nTargets > 1) final_idx++;
		for(int j = 0; j < final_idx; j++)
//			sum_res += 1.0 / (adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]);
			sum_res += 1.0 / (adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]) * weight[sort_idx[j]];
		delta_res = res / sum_res;
		
		double delta_cov = Double.NEGATIVE_INFINITY;
		for(int j = 0; j < final_idx; j++)
		{
			double temp = (1 - min_cov[sort_idx[j]]) 
					* (adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]);
			delta_cov = Math.max(delta_cov, temp);
		}
		double delta = Math.max(delta_res, delta_cov);
		for(int j = 0; j < final_idx; j++)
			min_cov[sort_idx[j]] += delta / (adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]);
		for(int j = 0; j < nTargets; j++)
			if(min_cov[j] > 1 || min_cov[j] < 0)
				System.out.println("Oops");
		return min_cov;
	}
}
