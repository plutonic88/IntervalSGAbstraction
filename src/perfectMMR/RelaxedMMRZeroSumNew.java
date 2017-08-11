package perfectMMR;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import lpWrapper.Configuration;
import lpWrapper.MIProblem;

public class RelaxedMMRZeroSumNew extends MIProblem{
	Vector<double[]> advPayoffs; // payoff samples
	Vector<Double> defUtilities; // optimal utilities of defender

	int nTargets;
	int nRes;
	
	HashMap<String, Integer> varMap; 
	HashMap<String, Integer> rowMap;	
	
	public RelaxedMMRZeroSumNew(Vector<double[]> advPayoffs, Vector<Double> defUtilities, int nTargets, int nRes)
	{
//		super();
		this.advPayoffs = advPayoffs;
		this.defUtilities = defUtilities;
		
		this.nTargets = nTargets;
		this.nRes = nRes;
		
		varMap = new HashMap<String, Integer>();
		rowMap = new HashMap<String, Integer>();
	}
	
	@Override
	protected void setProblemType() {
		// TODO Auto-generated method stub
		this.setProblemName("RelaxedMMRZeroSumNew");
		this.setProblemType(PROBLEM_TYPE.LP, OBJECTIVE_TYPE.MIN);
	}

	@Override
	protected void setColBounds() {
		// TODO Auto-generated method stub
		int index = 1;
		// regret
		addAndSetColumn("r", BOUNDS_TYPE.LOWER, 0.0, Configuration.MM, VARIABLE_TYPE.CONTINUOUS, 1.0);
		varMap.put("r", index++);
		
		// defender coverage
		for(int i = 0; i < nTargets; i++)
		{
			addAndSetColumn("x" + i, BOUNDS_TYPE.DOUBLE, 0.0, 1.0, VARIABLE_TYPE.CONTINUOUS, 0.0);
			varMap.put("x" + i, index++);
		}
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
	public void update()
	{
		int curSize = defUtilities.size();
		
		for(int i = curSize; i < advPayoffs.size(); i++)
		{
			double[] advPayoff = advPayoffs.get(i);
			double[] defCov = ORIGAMI(nTargets, nRes, advPayoff);
//			for(int j = 0; j < nTargets; j++)
//				System.out.print(defCov[j] + "\t");
//			System.out.println();
			double defUtility = Double.POSITIVE_INFINITY;
			for(int t = 0; t < nTargets; t++)
			{
				double defEU = -(defCov[t] * (advPayoff[t + nTargets] - advPayoff[t]) + advPayoff[t]);
				if(defUtility > defEU)
					defUtility = defEU;
			}
			defUtilities.add(defUtility);
		}
		int nSamples = advPayoffs.size();
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		for(int k = curSize; k < nSamples; k++)
		{
			double[] advPayoff = advPayoffs.get(k);
			double defUtility = defUtilities.get(k);
			for(int i = 0; i < nTargets; i++)
			{
				addAndSetRow("RC" + k + "_" + i, BOUNDS_TYPE.LOWER, defUtility 
						+ advPayoff[i], Configuration.MM);
				rowMap.put("RC" + k + "_" + i, this.getNumRows());
				ja.add(varMap.get("r"));
				ar.add(1.0);
				ja.add(varMap.get("x" + i));
				ar.add(advPayoff[i] - advPayoff[nTargets + i]);
				this.setMatRow(this.getNumRows(), ja, ar);
				
				ja.clear();
				ar.clear();
			}
		}
	}
	@Override
	protected void setRowBounds() {
		// TODO Auto-generated method stub
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		
		// Resource constraint		
		addAndSetRow("C", BOUNDS_TYPE.UPPER, -Configuration.MM, nRes);
		rowMap.put("C", this.getNumRows());
		for (int i = 0; i < nTargets; i++) {
	
			ja.add(varMap.get("x" + i));
			ar.add(1.0);
		}
		this.setMatRow(this.getNumRows(), ja, ar);
		
		ja.clear();
		ar.clear();
		
		// Regret constraint
		int nSamples = advPayoffs.size();
		for(int k = 0; k < nSamples; k++)
		{
			double[] advPayoff = advPayoffs.get(k);
			double defUtility = defUtilities.get(k);
			for(int i = 0; i < nTargets; i++)
			{
				addAndSetRow("RC" + k + "_" + i, BOUNDS_TYPE.LOWER, defUtility 
						+ advPayoff[i], Configuration.MM);
				rowMap.put("RC" + k + "_" + i, this.getNumRows());
				ja.add(varMap.get("r"));
				ar.add(1.0);
				ja.add(varMap.get("x" + i));
				ar.add(advPayoff[i] - advPayoff[nTargets + i]);
				this.setMatRow(this.getNumRows(), ja, ar);
				
				ja.clear();
				ar.clear();
			}
		}
	}

	@Override
	protected void generateData() {
		// TODO Auto-generated method stub
		
	}
	
	public double getOptRegret()
	{
		return this.getColumnPrimal(varMap.get("r"));
	}
	public double[] getOptDefCov()
	{
		double[] optDefCov = new double[nTargets];
		for(int i = 0; i < nTargets; i++)
		{
			optDefCov[i] = this.getColumnPrimal(varMap.get("x" + i));
		}
		return optDefCov;
	}
	
	public void end()
	{
		super.end();
//		if(!advPayoffs.isEmpty())
//		{
//			advPayoffs.clear();
//			advPayoffs = null;
//		}
//		
//		if(!defUtilities.isEmpty())
//		{
//			defUtilities.clear();
//			defUtilities = null;
//		}
		
		if(!varMap.isEmpty())
		{
			varMap.clear();
			varMap = null;
		}
		
		if(!rowMap.isEmpty())
		{
			rowMap.clear();
			rowMap = null;
		}
	}
}
