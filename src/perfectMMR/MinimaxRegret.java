package perfectMMR;

import java.util.Random;
import java.util.Vector;

public class MinimaxRegret {
	int nTargets;
	double nRes;
//	double[] def_payoff; //only used for non-zero sum
	
	double[] adv_payoff_lb;
	double[] adv_payoff_ub;
	
	// Output
	double[] opt_def_strategy;
	public double[] opt_adv_payoff;
	double[] opt_x;
	double opt_obj;
	
	// Just for printing
	public double runtime;
	public double lb_obj;
	public double ub_obj;
	public double[] def_solution;
	public MinimaxRegret()
	{
		
	}
	public MinimaxRegret(int nTargets, double nRes, double[] adv_payoff_lb, double[] adv_payoff_ub)
	{
		this.nTargets = nTargets;
		this.nRes = nRes;
		this.adv_payoff_lb = adv_payoff_lb;
		this.adv_payoff_ub = adv_payoff_ub;
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
//		for(int i = 0; i < nTargets; i++)
//			System.out.print(adv_payoff[sort_idx[i]] + "\t");
		int final_idx = 1;
		for(int i = 1; i < nTargets; i++)
		{
			final_idx = i;
			double[] dis = new double[i];
			for(int j = 0; j < i; j++)
			{
				dis[j] = (adv_payoff[sort_idx[i]] - r) /
						(adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]);
				res = res - dis[j];
				min_cov[sort_idx[j]] += dis[j];
				if(res < 0 || min_cov[sort_idx[j]] > 1) 
				{
					for(int k = 0; k <= j; k++)
					{
						res = res + dis[k];
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
			sum_res += 1.0 / (adv_payoff[nTargets + sort_idx[j]] - adv_payoff[sort_idx[j]]);
		delta_res = res / sum_res;
		
		double delta_cov = Double.NEGATIVE_INFINITY;
		for(int j = 0; j < final_idx; j++)
		{
//			double x_temp = 1 - min_cov[sort_idx[j]];
//			double u_temp = adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]];
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
	
	
	public void computeMinimaxRegret()
	{
//		long maxBytes = Runtime.getRuntime().maxMemory();
//	       System.out.println("Max memory: " + maxBytes / 1024 / 1024 + "M");
	    long startTime = System.currentTimeMillis();
		double ub = Double.POSITIVE_INFINITY;
		double obj = 0.0;
		Vector<double[]> x = new Vector<double[]>();
		Vector<double[]> adv_payoffs = new Vector<double[]>();
		Vector<double[]> def_payoffs = new Vector<double[]>();
		
		// Initialize
//		int numSamples = 30;
//		int numSamples_add = 3000;
		int numSamples = 0;
		int numSamples_add = 0;
//		MaxRegret g = new MaxRegret();
		for(int i = 0; i < 1; i++)
		{
//			double[] def_strategy = new double[nTargets];
			double[] adv_payoff = new double[2 * nTargets];
			double[] def_payoff = new double[2 * nTargets];
			adv_payoff[i] = adv_payoff_ub[i]; 
			adv_payoff[i + nTargets] = adv_payoff_ub[i + nTargets]; 
//			if(i == 0)
//				adv_payoff[i + nTargets] -= 2;
			def_payoff[i] = -adv_payoff[i + nTargets];
			def_payoff[i + nTargets] = -adv_payoff[i];
			for(int n = 0; n < nTargets; n++)
				if(n != i)
				{
					adv_payoff[n] = adv_payoff_lb[n]; 
					adv_payoff[n + nTargets] = adv_payoff_lb[n + nTargets]; 
					def_payoff[n] = -adv_payoff[n + nTargets];
					def_payoff[n + nTargets] = -adv_payoff[n];
				}
			double[] def_strategy = ORIGAMI(nTargets, nRes, adv_payoff);
//			System.out.println("Payoff " + i);
//			for(int t = 0; t < nTargets; t++)
//				System.out.print(def_strategy[t] + "\t");
//			System.out.println();
			x.add(def_strategy);
			adv_payoffs.add(adv_payoff);
			def_payoffs.add(def_payoff);
		}
//		for(int i = 0; i < nTargets; i++)
//		{
////			double[] def_strategy = new double[nTargets];
//			double[] adv_payoff = new double[2 * nTargets];
//			double[] def_payoff = new double[2 * nTargets];
//			adv_payoff[i] = adv_payoff_ub[i]; 
//			adv_payoff[i + nTargets] = adv_payoff_lb[i + nTargets]; 
////			if(i == 0)
////				adv_payoff[i + nTargets] -= 2;
//			def_payoff[i] = -adv_payoff[i + nTargets];
//			def_payoff[i + nTargets] = -adv_payoff[i];
//			for(int n = 0; n < nTargets; n++)
//				if(n != i)
//				{
//					adv_payoff[n] = adv_payoff_lb[n]; 
//					adv_payoff[n + nTargets] = adv_payoff_lb[n + nTargets]; 
//					def_payoff[n] = -adv_payoff[n + nTargets];
//					def_payoff[n + nTargets] = -adv_payoff[n];
//				}
//			double[] def_strategy = ORIGAMI(nTargets, nRes, adv_payoff);
////			System.out.println("Payoff " + i);
////			for(int t = 0; t < nTargets; t++)
////				System.out.print(def_strategy[t] + "\t");
////			System.out.println();
//			x.add(def_strategy);
//			adv_payoffs.add(adv_payoff);
//			def_payoffs.add(def_payoff);
//		}


		int count = 0;
		Random rand = new Random();
		double lb_old = Double.POSITIVE_INFINITY;
		double ub_old = Double.NEGATIVE_INFINITY;
		// Initialize lower bound and upper bound
		for(int i = 0; i < nTargets; i++)
		{
			lb_old = Math.min(lb_old, adv_payoff_lb[nTargets + i]);
			ub_old = Math.max(ub_old, adv_payoff_ub[i]);
		}
		lb_old *= 2;
		ub_old *= 2;
		int numIter = 10;
//		while (ub - obj > delta && count < 15)
		// Start the iteration process
//		runtime = new double[numIter];
//		lb_obj = new double[numIter];
//		ub_obj = new double[numIter];
		def_solution = new double[nTargets];
		while (count < numIter)
		{
			System.out.println(ub + "\t" + obj);
			// Compute sub-minimaxregret
			SubMinimaxRegret sub_regret = new SubMinimaxRegret(x, adv_payoffs, def_payoffs, nTargets, nRes);
			sub_regret.binarySearch(lb_old, ub_old);
			obj = sub_regret.obj_value;
			ub_old = -obj;
			long endTime = System.currentTimeMillis();
			runtime = endTime - startTime;
			double[] x_star = sub_regret.def_strategy;
			MRZeroSum max_regret = new MRZeroSum(nTargets, nRes, x_star, adv_payoff_lb, adv_payoff_ub);
			// To do: Call the function of computing max_regret
			max_regret.loadProblem();
			max_regret.solve();
			
			ub = max_regret.opt_obj;
			
//			runtime = endTime - startTime;
			lb_obj = obj;
			ub_obj = ub;
			
			if(count == numIter - 1)
			{
				opt_def_strategy = max_regret.opt_def_strategy;
				opt_adv_payoff = max_regret.opt_adv_payoff;
				opt_x = x_star;
				opt_obj = max_regret.opt_obj;
			}
			for(int i = 0; i < nTargets; i++)
				def_solution[i] = x_star[i];
			count++;
		}
		x.clear();
		adv_payoffs.clear();
	}
}
