package perfectMMR;

public class Maxmin {
	int nTargets;
	double nRes;
//	double[] def_payoff; //only used for non-zero sum
	
	double[] advPayoffLB;
	double[] advPayoffUB;
	
	double[] defCov;
	public Maxmin(int nTargets, double nRes, double[] advPayoffLB, double[] advPayoffUB)
	{
		this.nTargets = nTargets; 
		this.nRes = nRes;
		this.advPayoffLB = advPayoffLB;
		this.advPayoffUB = advPayoffUB;
	}
	public void solve()
	{
		double[] advPayoff = new double[2 * nTargets];
		for(int i = 0; i < 2 * nTargets; i++)
		{
			advPayoff[i] = (advPayoffLB[i] + advPayoffUB[i]) / 2;
		}
		defCov = ORIGAMI(nTargets, nRes, advPayoff);
//		defCov = ORIGAMI(nTargets, nRes, advPayoffUB);
	}
	public double[] getDefCov()
	{
		return defCov;
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
