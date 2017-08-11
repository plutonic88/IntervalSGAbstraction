package perfectMMR;

import java.util.Vector;

public class SubMinimaxRegret {
	Vector<double[]> x;
	Vector<double[]> adv_payoffs;
	Vector<double[]> def_payoffs;
//	double[] def_payoff;
	int nTargets;
	double nRes;
	double[] weight;
	
	// Output
	public double[] def_strategy;
	public double obj_value;
	
//	public double[] opt_adv_sample;
//	public double[] opt_def_sample;
//	public double value_sample;
	
	public SubMinimaxRegret(Vector<double[]> x, Vector<double[]> adv_payoffs, Vector<double[]> def_payoffs, int nTargets, double nRes)
	{
		this.x = x;
		this.adv_payoffs = adv_payoffs;
		this.def_payoffs = def_payoffs;
		this.nTargets = nTargets;
		this.nRes = nRes;
//		this.def_payoff = def_payoff;
		this.def_strategy = new double[nTargets];
		this.obj_value = Double.NEGATIVE_INFINITY;
		this.weight = new double[nTargets];
	}
	public void setWeight(double[] weight)
	{
		for(int i = 0; i < nTargets; i++)
			this.weight[i] = weight[i];
	}
	public double[] computeLowBound(double v, double[] low_bound_old, double[] adv_payoff, double[] def_payoff)
	{
		boolean yes = false;
		double[] low_bound_new = new double[nTargets];
		double min_cov = Double.POSITIVE_INFINITY;
		double[] temp_cov = new double[nTargets];
		
		for(int i = 0; i < nTargets; i++)
		{
			boolean passed = false;
			double sum_cov = 0.0;
			if (v > def_payoff[i] + 1e-6)
				continue;
			if(v <= def_payoff[nTargets + i] + 1e-6)
				temp_cov[i] = 0;
			else 
				temp_cov[i] = (v - def_payoff[nTargets + i]) / (def_payoff[i] - def_payoff[nTargets  + i]);
			temp_cov[i] = Math.max(temp_cov[i], low_bound_old[i]);
//			sum_cov += temp_cov[i];
			sum_cov += temp_cov[i] * weight[i];
			double key_utility = temp_cov[i] * adv_payoff[nTargets + i] + (1 - temp_cov[i]) * adv_payoff[i];
			for(int j = 0; j < nTargets; j++)
			{
				if(j != i)
				{
					double temp = 0;
					if(key_utility < adv_payoff[nTargets + j])
					{
						passed = false;
						break;
					}
					if(key_utility <= adv_payoff[j])
					{
						temp = (key_utility - adv_payoff[j]) / (adv_payoff[j + nTargets] - adv_payoff[j]);
						temp_cov[j] = Math.max(temp, low_bound_old[j]);
//						sum_cov += temp_cov[j];
						sum_cov += temp_cov[j] * weight[j];
						passed = true;
					}
					else
					{
						temp_cov[j] = 0.0;
						passed = true;
					}
				}
			}
			if(passed && sum_cov < min_cov)
			{
				min_cov = sum_cov;
				for(int j = 0; j < nTargets; j++)
					low_bound_new[j] = temp_cov[j];
				yes = passed;
			}
		}
		if(yes && min_cov <= nRes) 
			return low_bound_new;
		return null;
	}
	public boolean isFeasible(double v, double[] def_strategy, double[] adv_payoff, double[] def_payoff)
	{
		double error = 1e-6;
		boolean is_feasible = false;
		double[] adv_u = new double[nTargets];
		double temp_max = Double.NEGATIVE_INFINITY;
//		double sum_res = 0.0;
//		for(int i = 0; i < nTargets; i++)
//			sum_res += def_strategy[i];
//		if(sum_res - nRes > 1e-4)
//			return is_feasible;
		for(int i = 0; i < nTargets; i++)
		{
			adv_u[i] = def_strategy[i] * adv_payoff[nTargets + i] + (1 - def_strategy[i]) * adv_payoff[i];
			if(adv_u[i] > temp_max)
				temp_max = adv_u[i];
		}
		
		for(int i = 0; i < nTargets; i++)
		{
			double temp = temp_max - adv_u[i];
			if(temp < error)
//			if(temp >= 0.0)
			{
				double temp_u = def_strategy[i] * def_payoff[i] + (1 - def_strategy[i])  * def_payoff[nTargets + i];
				if(v - temp_u < error)
//				if(v - temp_u >= 0)
					return true;
			}
		}
		return is_feasible;
	}
	
	public void binarySearch(double lb_old, double ub_old)
	{
		double error = 1e-3;
//		double lb = Double.POSITIVE_INFINITY;
//		double ub = Double.NEGATIVE_INFINITY;
//		// Initialize lower bound and upper bound
//		for(int i = 0; i < nTargets; i++)
//		{
//			lb = Math.min(lb, def_payoff[nTargets + i]);
//			ub = Math.max(ub, def_payoff[i]);
//		}
//		lb *= 2;
//		ub *= 2;
		double lb = lb_old;
		double ub = ub_old;
		double delta = 1e-6; // tolerance
		
		double[] adv_u = new double[nTargets];
		double[] low_bound_old = new double[nTargets];
		double[] key_def_u = new double[x.size()];
		// Compute the defender's utility of the sample
		for(int k = 0; k < x.size(); k++)
		{
			double[] temp_def = x.get(k);
			double[] temp_adv = adv_payoffs.get(k);
			double[] temp_defPayoff = def_payoffs.get(k);
			double temp_max = Double.NEGATIVE_INFINITY;
			key_def_u[k] = Double.NEGATIVE_INFINITY;
			for(int i = 0; i < nTargets; i++)
			{
				adv_u[i] = temp_def[i] * temp_adv[nTargets + i] + (1 - temp_def[i]) * temp_adv[i];
				if(temp_max < adv_u[i])
					temp_max = adv_u[i];
			}
			for(int i = 0; i < nTargets; i++)
			{
				if(temp_max - adv_u[i] < error)
				{
					double def_u = temp_def[i] * temp_defPayoff[i] + (1 - temp_def[i]) * temp_defPayoff[nTargets + i];
					if(key_def_u[k] < def_u)
						key_def_u[k] = def_u;
				}
			}
		}
		int count = 0;
		double[] v = new double[x.size()];
		while(ub - lb > delta)
		{
			count++;
//			System.out.println(count);
//			System.out.println("Lower bound: " + lb + "\t" + "Upper bound: " + ub);
			double mid = (ub + lb) / 2;
			// Start checking if mid is a feasible defender's utility
			boolean is_stop = false;
			for(int i = 0; i < nTargets; i++)
				low_bound_old[i] = 0.0;
			
//			double[] v = new double[x.size()];
			for(int k = 0; k < x.size(); k++)
				v[k] = mid + key_def_u[k];
//			if(ub - lb <= 2 * delta)
//			{
//				System.out.println("Checking");
//			}
			while(!is_stop)
			{
				for(int k = 0; k < x.size(); k++)
				{
					double[] temp_adv = adv_payoffs.get(k);
					double[] temp_defPayoff = def_payoffs.get(k);
					// Compute the lower bound of def's coverage for every targets
					double[] temp_bound = computeLowBound(v[k], low_bound_old, temp_adv, temp_defPayoff);
					if(temp_bound == null)
					{
						is_stop = true;
						break;
					}
					else
					{
						for(int i = 0; i < nTargets; i++)
							low_bound_old[i] = Math.max(low_bound_old[i], temp_bound[i]);
					}
				}
				double sum_res = 0.0;
				for(int i = 0; i < nTargets; i++)
//					sum_res += low_bound_old[i];
					sum_res += low_bound_old[i] * weight[i];
				if(sum_res - nRes > 1e-6)
					is_stop = true;
				if(is_stop)
				{
					ub = mid;
					break;
				}
//				for(int i = 0; i < nTargets; i++)
//					System.out.print(low_bound_old[i] + "\t");
//				System.out.println("---");
				// Check if we find a feasible solution.
					
				is_stop = true;
				for(int k = 0; k < x.size(); k++)
				{
					boolean is_feasible = isFeasible(v[k], low_bound_old, adv_payoffs.get(k), def_payoffs.get(k));
					if(!is_feasible)
					{
						is_stop = false;
						break;
					}
				}
				if(is_stop)
				{
					lb = mid;
//					def_strategy = low_bound_old;
//					System.out.print("---");
					for(int i = 0; i < nTargets; i++)
					{
						def_strategy[i] = low_bound_old[i];
//						System.out.print(def_strategy[i] + "\t");
					}
//					System.out.println("---");
					obj_value = -lb;
				}
			}
		}
		
		// Testing
//		boolean is_feasible = false;
//		for(int k = 0; k < x.size(); k++)
//		{
//			is_feasible = isFeasible(lb + key_def_u[k], def_strategy, adv_payoffs.get(k));
//			if(!is_feasible)
//				System.out.println("Error");
//		}
//		value_sample = Double.NEGATIVE_INFINITY;
//		for(int k = 0; k < x.size(); k++)
//		{
//			double[] temp_def = def_strategy;
//			double[] temp_adv = adv_payoffs.get(k);
//			double temp_max = Double.NEGATIVE_INFINITY;
//			double key_def = Double.NEGATIVE_INFINITY;
//			for(int i = 0; i < nTargets; i++)
//			{
//				adv_u[i] = temp_def[i] * temp_adv[nTargets + i] + (1 - temp_def[i]) * temp_adv[i];
//				if(temp_max < adv_u[i])
//					temp_max = adv_u[i];
//			}
//			for(int i = 0; i < nTargets; i++)
//			{
//				if(temp_max - adv_u[i] < error)
//				{
//					double def_u = temp_def[i] * def_payoff[i] + (1 - temp_def[i]) * def_payoff[nTargets + i];
//					if(key_def < def_u)
//						key_def = def_u;
//				}
//			}
//			
//			if(value_sample < key_def_u[k] - key_def)
//			{
//				value_sample = key_def_u[k] - key_def;
//				opt_def_sample = x.get(k);
//				opt_adv_sample = temp_adv;
//			}
//		}
//		System.out.println(value_sample + "\t" + lb);
//		is_feasible = false;
//		for(int k = 0; k < x.size(); k++)
//		{
//			is_feasible = isFeasible(-value_sample + key_def_u[k], def_strategy, adv_payoffs.get(k));
//			if(!is_feasible)
//				System.out.println("Error2");
//		}
//		for(int k = 0; k < x.size(); k++)
//		{
//			double[] temp_def = def_strategy;
//			double[] temp_adv = adv_payoffs.get(k);
//			double[] temp = computeLowBound(key_def_u[k] - value_sample, temp_def, temp_adv);
//			if(temp == null)
//				System.out.println("Oops");
//		}
//		System.out.println("Testing");
	}
}

