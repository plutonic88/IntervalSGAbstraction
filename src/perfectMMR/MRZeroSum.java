package perfectMMR;

import lpWrapper.AMIProblem.STATUS_TYPE;
import lpWrapper.LPSolverException;


public class MRZeroSum{
	MRZeroSumLP localSearch;
	int nTargets;
	double nRes;
	double[] x_star;
	double[] adv_payoff_lb;
	double[] adv_payoff_ub;
	
	int numStarts = 100;
	double err = 1e-4;
	
	
	// Output
	public double[] opt_def_strategy;
	public double[] opt_adv_payoff;
	public double opt_obj = Double.NEGATIVE_INFINITY;;
	
	public MRZeroSum(int nTargets, double nRes, double[] x_star, double[] adv_payoff_lb, double[] adv_payoff_ub)
	{
		this.nTargets = nTargets;
		this.nRes = nRes;
		this.x_star = x_star;
		this.adv_payoff_lb = adv_payoff_lb;
		this.adv_payoff_ub = adv_payoff_ub;
		localSearch = new MRZeroSumLP();
		opt_def_strategy = new double[nTargets];
		opt_adv_payoff = new double[2 * nTargets];
	}
	public void end()
	{
		localSearch.end();
	}
	public void loadProblem()
	{
		localSearch.loadProblem();
	}
	/* Function computeMinCov: use Origami 
	 * Input:
	 * 	- resource: number of resources left for other targets
	 * 	- idx: an array of targets which haven't been assigned 
	 * 	- numVar: number of non-assigned targets
	 * Output:
	 * 	- an array of resources to all targets, the not-assigned target will be zero.
	 */
	public double[] computeMinCov(double resource, int[] idx, int numVar)
	{
		double[] min_cov = new double[nTargets];
		for(int i = 0; i < nTargets; i++)
			min_cov[i] = 0.0;
		double[] values = new double[numVar];
		for(int i = 0; i < numVar; i++)
			values[i] = adv_payoff_lb[idx[i]];
		QuickSort my_sort = new QuickSort(values, idx, numVar);
		my_sort.get_sort_id();
		int[] sort_idx = my_sort.sort_id;
		for(int i = 0; i < numVar / 2; i++)
		{
			int temp = sort_idx[i];
			sort_idx[i] = sort_idx[numVar - 1 - i];
			sort_idx[numVar - 1 - i] = temp;
		}
		double r = adv_payoff_lb[sort_idx[0]];
		double res = resource;
		boolean is_stop = false;
		int final_idx = 1;
		for(int i = 1; i < numVar; i++)
		{
			final_idx = i;
			double[] dis = new double[i];
			for(int j = 0; j < i; j++)
			{
				dis[j] = (adv_payoff_lb[sort_idx[i]] - r) /
						(adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]]);
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
				r = adv_payoff_lb[sort_idx[i]];
			else break;
		}
		
		double delta_res = 0.0;
		double sum_res = 0.0;
		if(!is_stop && numVar > 1) final_idx++;
		for(int j = 0; j < final_idx; j++)
			sum_res += 1.0 / (adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]]);
		delta_res = res / sum_res;
		
		double delta_cov = Double.NEGATIVE_INFINITY;
		for(int j = 0; j < final_idx; j++)
		{
			double temp = (1 - min_cov[sort_idx[j]]) 
					* (adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]]);
			delta_cov = Math.max(delta_cov, temp);
		}
		double delta = Math.max(delta_res, delta_cov);
		for(int j = 0; j < final_idx; j++)
			min_cov[sort_idx[j]] += delta / (adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]]);
		
		return min_cov;
	}
	
	public void solve() 
	{
//		double maxValue = Double.NEGATIVE_INFINITY;
		// Iterate all possible attacked targets 
		for(int j = 0; j < nTargets; j++)
		{
//			System.out.println("Target " + j + "................................................................................");
			// Compute cstar;
			double cstar = Double.NEGATIVE_INFINITY;
			for(int t = 0; t < nTargets; t++)
			{
				if(t != j)
				{
					double attEU = x_star[t] * adv_payoff_lb[nTargets + t] + (1 - x_star[t]) * adv_payoff_lb[t];
					cstar = Math.max(cstar, attEU);
				}
			}
			double rLB = adv_payoff_lb[j];
			double rUB = adv_payoff_ub[j];
			double pLB = adv_payoff_lb[j + nTargets];
			double pUB = adv_payoff_ub[j + nTargets];
			// End of compute cstar
			
			for(int iter = 0; iter <= numStarts; iter++)
			{
				double delta = Double.POSITIVE_INFINITY;
				// hill-climbing
				double curOptValue = Double.NEGATIVE_INFINITY;
				double[] reward = new double[nTargets];
				double[] penalty = new double[nTargets];
				double[] cov = new double[nTargets];
				cov[j] = 1.0 / numStarts * iter;
				do
				{
					// Fix x_j------------------------------------------------------------
//					System.out.println("--------------------------------------------------------------------------------------");
					double res = nRes - cov[j];
					int[] idx = new int[nTargets - 1];
					int temp = 0;
					for(int k = 0; k < nTargets; k++)
						if(k != j)
							idx[temp++] = k;
					double[] x = computeMinCov(res, idx, nTargets - 1);
					x[j] = cov[j];
					// Compute c
					double c = Double.POSITIVE_INFINITY;
					for(int t = 0; t < nTargets; t++)
					{
						if(t != j)
						{
							double attEU = x[t] * adv_payoff_lb[nTargets + t] + (1 - x[t]) * adv_payoff_lb[t];
							c = Math.min(c, -attEU);
						}
					}
					// End of compute c
					localSearch.setCoeff(x[j], x_star[j], c, cstar, rLB, rUB, pLB, pUB, j);
					
					try{
						localSearch.solve();
					}catch(LPSolverException e)
					{
						localSearch.writeProb("/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/test");
						break;
					}
					
					if(localSearch.getSolveStatus() != STATUS_TYPE.INFEASIBLE)
					{
//						System.out.println(curOptValue + "\t" + localSearch.getObjValue());
//						double defEU = Double.POSITIVE_INFINITY;
//						for(int t = 0; t < nTargets; t++)
//						{
//							if(t != j)
//							{
//								double attEU = adv_payoff_lb[t] * (1 - x[t]) + adv_payoff_lb[t + nTargets] * x[t];
//								if(defEU > -attEU)
//									defEU = -attEU;
//							}
//							else
//							{
//								double attEU = localSearch.getReward() * (1 - x[t]) + localSearch.getPenalty() * x[t];
//								if(defEU > -attEU)
//									defEU = -attEU;
//							}
//						}
//						if(defEU < localSearch.getDefEU())
//							System.out.println(defEU + "\t" + localSearch.getDefEU() + "\t" + c + "\t" 
//						+ (-localSearch.getReward() * (1 - x[j]) - localSearch.getPenalty() * x[j]));
//						double defEUStar = Double.POSITIVE_INFINITY;
//						for(int t = 0; t < nTargets; t++)
//						{
//							if(t != j)
//							{
//								double attEU = adv_payoff_lb[t] * (1 - x_star[t]) + adv_payoff_lb[t + nTargets] * x_star[t];
//								if(defEUStar > -attEU)
//									defEUStar = -attEU;
//							}
//							else
//							{
//								double attEU = localSearch.getReward() * (1 - x_star[t]) + localSearch.getPenalty() * x_star[t];
//								if(defEUStar > -attEU)
//									defEUStar = -attEU;
//							}
//						}
//						System.out.println(defEU - defEUStar);
						
						curOptValue = localSearch.getObjValue();
						reward[j] = localSearch.getReward();
						penalty[j] = localSearch.getPenalty();
						
					}
					else break;
					
					// Fix R^a_j, P^a_j----------------------------------------------------
					// Use ORIGAMI to compute
					double advRewardLB = adv_payoff_lb[j];
					double advPenaltyLB = adv_payoff_lb[j + nTargets];
					adv_payoff_lb[j] = reward[j];
					adv_payoff_lb[j + nTargets] = penalty[j];
					int[] idx2 = new int[nTargets];
					int temp2 = 0;
					for(int k = 0; k < nTargets; k++)
						idx2[temp2++] = k;
					double[] x2 = computeMinCov(nRes, idx2, nTargets);
					double tempValue = Double.POSITIVE_INFINITY;
					for(int t = 0; t < nTargets; t++)
					{
						double defEU = -(adv_payoff_lb[t] * (1 - x2[t]) + adv_payoff_lb[t + nTargets] * x2[t]);
						if(tempValue > defEU)
							tempValue = defEU;
					}
					tempValue += (x_star[j] * adv_payoff_lb[j + nTargets] + (1 - x_star[j]) * adv_payoff_lb[j]);
					for(int t = 0; t < nTargets; t++)
					{
						reward[t] = adv_payoff_lb[t];
						penalty[t] = adv_payoff_lb[t + nTargets];
						cov[t] = x2[t];
					}
//					System.out.println(curOptValue + "\t" + tempValue);
					// Restore lower bound
					adv_payoff_lb[j] = advRewardLB;
					adv_payoff_lb[j + nTargets] = advPenaltyLB;
					delta = tempValue - curOptValue;
					curOptValue = tempValue;
//					System.out.println("--------------------------------------------------------------------------------------");
				}while(delta > err);
				// end of hill-climbing
				if(opt_obj < curOptValue)
				{
					opt_obj = curOptValue;
					for(int t = 0; t < nTargets; t++)
					{
						opt_adv_payoff[t] = reward[t];
						opt_adv_payoff[t + nTargets] = penalty[t];
						opt_def_strategy[t] = cov[t];
					}
				}
			}
		}
	}
}
