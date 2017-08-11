package perfectMMR;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
public class MRNonZeroSum {
	// Input
	int nTargets;
	double nRes;
	double[] def_payoff;
	double[] x_star;
	double[] adv_payoff_lb;
	double[] adv_payoff_ub;
	
	// Temp variables
	boolean[][] is_skip_old;
	boolean[][] is_skip;
	
	// Output
	public double[] opt_def_strategy;
	public double[] opt_adv_payoff;
	public double opt_obj;
	
	public MRNonZeroSum()
	{
		
	}
	public MRNonZeroSum(int nTargets, double nRes, double[] def_payoff, double[] x_star, double[] adv_payoff_lb, double[] adv_payoff_ub)
	{
		this.nTargets = nTargets;
		this.nRes = nRes;
		this.def_payoff = def_payoff;
		this.x_star = x_star;
		this.adv_payoff_lb = adv_payoff_lb;
		this.adv_payoff_ub = adv_payoff_ub;
		is_skip = new boolean[nTargets][nTargets];
		is_skip_old = new boolean[nTargets][nTargets];
		for(int i = 0; i < nTargets; i++)
			for(int j = 0; j < nTargets; j++)
			{
				is_skip[i][j] = false;
				is_skip_old[i][j] = false;
			}
//		opt_def_strategy = new double[nTargets];
		opt_adv_payoff = new double[2 * nTargets];
		opt_obj = Double.NEGATIVE_INFINITY;
	}
	
	// idx: index of targets
	// outcome: there will still be nTargets.
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
//			double x_temp = 1 - min_cov[sort_idx[j]];
//			double u_temp = adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]];
			double temp = (1 - min_cov[sort_idx[j]]) 
					* (adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]]);
			delta_cov = Math.max(delta_cov, temp);
		}
		double delta = Math.max(delta_res, delta_cov);
		for(int j = 0; j < final_idx; j++)
			min_cov[sort_idx[j]] += delta / (adv_payoff_lb[nTargets + sort_idx[j]] - adv_payoff_lb[sort_idx[j]]);
		
		return min_cov;
	}
	
	public boolean isFeasible(double obj_value)
	{
		double delta = 1e-4;
		double epsilon = 1e-4;
		boolean is_feasible = false;
//		double[] x_min = new double[nTargets];
		IloCplex cplex;
		IloNumVar Ra_i, Ra_j, Pa_i, Pa_j;
		try {
			cplex = new IloCplex();
			cplex.setOut(null);
		for(int i = 0; i < nTargets; i++) // key targets
		{
			double def_temp_i = x_star[i] * def_payoff[i] + (1 - x_star[i]) * def_payoff[nTargets + i];
			for(int j = 0; j < nTargets; j++) // key targets
			{
				if(!is_skip[i][j]) // if not skip this case
				{
					double res = nRes;
					if(j != i)
					{
						double def_u_star = x_star[j] * def_payoff[j] + (1 - x_star[j]) *  def_payoff[nTargets + j];
						double cov_i = ((obj_value + def_u_star) - def_payoff[i + nTargets]) / (def_payoff[i] - def_payoff[i + nTargets]);
						if(cov_i <= 1.0)
						{
							cov_i = Math.max(0.0, cov_i);
							res = res - cov_i;
							
							double def_temp_j = x_star[j] * def_payoff[j] + (1 - x_star[j]) * def_payoff[nTargets + j];
							
							double max_low_star = Double.NEGATIVE_INFINITY;
							for(int k = 0; k < nTargets; k++)
							{
								if(k != i && k != j)
								{	
									double temp_u_star = x_star[k] * adv_payoff_lb[nTargets + k] + (1 - x_star[k]) * adv_payoff_lb[k];
									double def_temp_k =  x_star[k] * def_payoff[k] + (1 - x_star[k]) * def_payoff[nTargets + k];
									if(def_temp_j >= def_temp_k)
										max_low_star = Math.max(max_low_star, temp_u_star);
									else
										max_low_star = Math.max(max_low_star, temp_u_star + epsilon);
								}
							}
							double max_low_star_test = max_low_star;
							double temp_u_star = x_star[i] * adv_payoff_lb[nTargets + i] + (1 - x_star[i]) * adv_payoff_lb[i];
							if(def_temp_j >= def_temp_i)
								max_low_star_test = Math.max(max_low_star_test, temp_u_star);
							else
								max_low_star_test = Math.max(max_low_star_test, temp_u_star + epsilon);
							temp_u_star = x_star[j] * adv_payoff_ub[nTargets + j] + (1 - x_star[j]) * adv_payoff_ub[j];
							if(temp_u_star < max_low_star_test)
							{
								is_skip[i][j] = true;
								continue;
							}
							
							// do binary search on x_j
							int[] idx = new int[nTargets - 2];
							int temp = 0;
							for(int k = 0; k < nTargets; k++)
								if(k != i && k!= j)
									idx[temp++] = k;
							double lb_x = 0.0;
							double ub_x = Math.min(1.0, res);
							int nSample = 5;
							
							for(int count = 0; count <= nSample; count++)
							{
//								if(count % 1000 == 0)
//									System.out.println(count);
								double mid_x = lb_x + count * (ub_x - lb_x) / nSample;
								
								double res_temp = res - mid_x;
								double[] x_min = computeMinCov(res_temp, idx, nTargets - 2);
								x_min[i] = cov_i;
								x_min[j] = mid_x;	
								
//								double test_cov = 0.0;
//								for(int k = 0; k < nTargets; k++)
//									test_cov += x_min[k];
//								
//								if(test_cov < 3.0 - epsilon)
//									System.out.println(test_cov);
								
								double lb_obj = Double.NEGATIVE_INFINITY;
								double ub_obj = Double.POSITIVE_INFINITY;
								while(ub_obj - lb_obj > delta)
								{	
									double max_low_prime = Double.NEGATIVE_INFINITY;
									for(int k = 0; k < nTargets; k++)
									{
										if(k != i && k != j)
										{
											double temp_u_prime = x_min[k] * adv_payoff_lb[nTargets + k] + (1 - x_min[k]) * adv_payoff_lb[k];
											max_low_prime = Math.max(max_low_prime, temp_u_prime);
										}
									}
									cplex.clearModel();
									Ra_i = cplex.numVar(adv_payoff_lb[i], adv_payoff_ub[i]);
									Ra_i.setName("Ra_i");
									Pa_i = cplex.numVar(adv_payoff_lb[nTargets + i], adv_payoff_ub[nTargets + i]);
									Pa_i.setName("Pa_i");
									Ra_j = cplex.numVar(adv_payoff_lb[j], adv_payoff_ub[j]);
									Ra_j.setName("Ra_j");
									Pa_j = cplex.numVar(adv_payoff_lb[nTargets + j], adv_payoff_ub[nTargets + j]);
									Pa_j.setName("Pa_j");
									
									cplex.addGe(cplex.sum(cplex.prod(x_min[i], Pa_i), cplex.prod(1 - x_min[i], Ra_i)), max_low_prime);
									cplex.addGe(cplex.sum(cplex.prod(x_star[j], Pa_j), cplex.prod(1 - x_star[j], Ra_j)), max_low_star);
									if(def_temp_j >= def_temp_i)
										cplex.addGe(cplex.sum(cplex.prod(x_star[j], Pa_j), cplex.prod(1 - x_star[j], Ra_j)),
												cplex.sum(cplex.prod(x_star[i], Pa_i), cplex.prod(1 - x_star[i], Ra_i)));
									else
										cplex.addGe(cplex.diff(cplex.sum(cplex.prod(x_star[j], Pa_j), cplex.prod(1 - x_star[j], Ra_j)),
												cplex.sum(cplex.prod(x_star[i], Pa_i), cplex.prod(1 - x_star[i], Ra_i))), epsilon);
									
									IloNumExpr f = cplex.diff(cplex.sum(cplex.prod(x_min[i], Pa_i), cplex.prod(1 - x_min[i], Ra_i))
											, cplex.sum(cplex.prod(x_min[j], Pa_j), cplex.prod(1 - x_min[j], Ra_j)));
									cplex.addMaximize(f);
									cplex.solve();
									if(cplex.isPrimalFeasible() && cplex.getObjValue() >= 0)
									{
										// Assign solution
										opt_def_strategy = x_min;
										opt_adv_payoff[i] = cplex.getValue(Ra_i); 
										opt_adv_payoff[nTargets + i] = cplex.getValue(Pa_i); 
										opt_adv_payoff[j] = cplex.getValue(Ra_j); 
										opt_adv_payoff[nTargets + j] = cplex.getValue(Pa_j); 
										for(int k = 0; k < nTargets; k++)
										{
											if(k != i && k != j)
											{
												opt_adv_payoff[k] = adv_payoff_lb[k];
												opt_adv_payoff[nTargets + k] = adv_payoff_lb[nTargets + k];
											}
										}
										// end of assign solution
										cplex.end();
										cplex = null;
										is_feasible = true;
										return is_feasible;
									}
									else if(cplex.isPrimalFeasible())
									{
										lb_obj = cplex.getObjValue();
										double u_prime = cplex.getValue(Ra_i) * (1 - x_min[i]) + cplex.getValue(Pa_i) * x_min[i];
										double temp_res = x_min[i];
										for(int k = 0; k < nTargets; k++)
										{
											if(k != i && k != j)
											{
												x_min[k] = Math.max(0, (u_prime - adv_payoff_lb[k]) / (adv_payoff_lb[nTargets + k] - adv_payoff_lb[k]));
												temp_res += x_min[k];
											}
										}
										x_min[j] = Math.min(1.0, nRes - temp_res);
										ub_obj = u_prime - (x_min[j] * cplex.getValue(Pa_j) + (1 - x_min[j]) * cplex.getValue(Ra_j));
										if(ub_obj >= 0)
										{
											opt_def_strategy = x_min;
											opt_adv_payoff[i] = cplex.getValue(Ra_i); 
											opt_adv_payoff[nTargets + i] = cplex.getValue(Pa_i); 
											opt_adv_payoff[j] = cplex.getValue(Ra_j); 
											opt_adv_payoff[nTargets + j] = cplex.getValue(Pa_j); 
											for(int k = 0; k < nTargets; k++)
											{
												if(k != i && k != j)
												{
													opt_adv_payoff[k] = adv_payoff_lb[k];
													opt_adv_payoff[nTargets + k] = adv_payoff_lb[nTargets + k];
												}
											}
											cplex.end();
											cplex = null;
											is_feasible = true;
											return is_feasible;
										}
									}
									else
									{
										cplex.clearModel();
										break;
									}
									cplex.clearModel();
								}
							} // end doing binary search on x_j
							if(!is_feasible)
								is_skip[i][j] = true;
						}
						else is_skip[i][j] = true;
					}
					else
					{
						int[] idx = new int[nTargets - 1];
						int temp = 0;
						for(int k = 0; k < nTargets; k++)
							if(k != i)
								idx[temp++] = k;
						double def_u_star = x_star[j] * def_payoff[j] + (1 - x_star[j]) *  def_payoff[nTargets + j];
						double cov_i = ((obj_value + def_u_star) - def_payoff[i + nTargets]) / (def_payoff[i] - def_payoff[i + nTargets]);
						if(cov_i <= 1.0)
						{
							cov_i = Math.max(0.0, cov_i);
							res = res - cov_i;
							double[] x_min = computeMinCov(res, idx, nTargets - 1);
							x_min[i] = cov_i;
//							double test_cov = 0.0;
//							for(int k = 0; k < nTargets; k++)
//								test_cov += x_min[k];
//							
//							if(test_cov < 3.0 - epsilon)
//								System.out.println(test_cov);
							double def_temp_j = x_star[j] * def_payoff[j] + (1 - x_star[j]) * def_payoff[nTargets + j];
							double max_low_prime = Double.NEGATIVE_INFINITY;
							double max_low_star = Double.NEGATIVE_INFINITY;
							for(int k = 0; k < nTargets; k++)
							{
								if(k != i)
								{
									double temp_u_prime = x_min[k] * adv_payoff_lb[nTargets + k] + (1 - x_min[k]) * adv_payoff_lb[k];
									max_low_prime = Math.max(max_low_prime, temp_u_prime);
									
									double temp_u_star = x_star[k] * adv_payoff_lb[nTargets + k] + (1 - x_star[k]) * adv_payoff_lb[k];
									double def_temp_k =  x_star[k] * def_payoff[k] + (1 - x_star[k]) * def_payoff[nTargets + k];
									if(def_temp_j >= def_temp_k)
										max_low_star = Math.max(max_low_star, temp_u_star);
									else
										max_low_star = Math.max(max_low_star, temp_u_star + epsilon);
								}
							}
							cplex.clearModel();
							Ra_i = cplex.numVar(adv_payoff_lb[i], adv_payoff_ub[i]);
							Ra_i.setName("Ra_i");
							Pa_i = cplex.numVar(adv_payoff_lb[nTargets + i], adv_payoff_ub[nTargets + i]);
							Pa_i.setName("Pa_i");

							cplex.addGe(cplex.sum(cplex.prod(x_min[i], Pa_i), cplex.prod(1 - x_min[i], Ra_i)), max_low_prime);
				
							cplex.addGe(cplex.sum(cplex.prod(x_star[i], Pa_i), cplex.prod(1 - x_star[i], Ra_i)), max_low_star);
							
							cplex.addMaximize();
							cplex.solve();
							if(cplex.isPrimalFeasible())
							{
								// Assign solution
								opt_def_strategy = x_min;
								opt_adv_payoff[i] = cplex.getValue(Ra_i); 
								opt_adv_payoff[nTargets + i] = cplex.getValue(Pa_i); 
								for(int k = 0; k < nTargets; k++)
								{
									if(k != i)
									{
										opt_adv_payoff[k] = adv_payoff_lb[k];
										opt_adv_payoff[nTargets + k] = adv_payoff_lb[nTargets + k];
									}
								}
								// end of assign solution
								cplex.end();
								is_feasible = true;
								return is_feasible;
							}
							else
							{
								is_skip[i][j] = true;
							}
							cplex.clearModel();
							
						}
						else is_skip[i][j] = true;
					}
				}
			}
		}
		cplex.end();
		cplex = null;
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return is_feasible;
	}
	public void computeMaxRegret()
	{
		double delta = 1e-4;
		double lb_obj = -10, ub_obj = 10; // change later
		// Initialize lower bound and upper bound
		double min_prime = Double.POSITIVE_INFINITY;
		double max_prime = Double.NEGATIVE_INFINITY;
		double min_star = Double.POSITIVE_INFINITY;
		double max_star = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < nTargets; i++)
		{
			min_prime = Math.min(min_prime, def_payoff[nTargets + i]);
			max_prime = Math.max(max_prime, def_payoff[i]);
			double temp = x_star[i] * def_payoff[i] + (1 - x_star[i]) * def_payoff[nTargets + i];
			min_star = Math.min(min_star, temp);
			max_star = Math.max(max_star, temp);
		}
		ub_obj = max_prime - min_star;
		lb_obj = min_prime - max_star;
//		lb_obj = -100;
		// Start doing binary search
		while(ub_obj - lb_obj > delta)
		{
//			System.out.println("Lb: " + lb_obj + "\tUb:" + ub_obj);
			double mid_obj = (ub_obj + lb_obj) / 2;
//			double mid_obj = 7.6;
			if(isFeasible(mid_obj))
			{
				lb_obj = mid_obj;
				for(int i = 0; i < nTargets; i++)
					for(int j = 0; j < nTargets; j++)
						is_skip_old[i][j] = is_skip[i][j]; 
				opt_obj = mid_obj;
			}
			else
			{
				ub_obj = mid_obj;
				for(int i = 0; i < nTargets; i++)
					for(int j = 0; j < nTargets; j++)
						is_skip[i][j] = is_skip_old[i][j]; 
			}
		}
		
	}
}
