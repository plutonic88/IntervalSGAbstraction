package perfectMMR;


public class MRZeroSumNew {
	int nTargets;
	double nRes;
	
	double[] xStar;
	double[] advPayoffLB;
	double[] advPayoffUB;
	double[] weight;
	
	public double maxRegret;
	public double[] optDefCov;
	public double[] optAdvPayoff;
	public MRZeroSumNew(int nTargets, double nRes, double[] xStar, double[] advPayoffLB, double[] advPayoffUB)
	{
		this.xStar = xStar;
		this.advPayoffLB = advPayoffLB;
		this.advPayoffUB = advPayoffUB;
		
		this.nTargets = nTargets;
		this.nRes = nRes;
		
		optDefCov = new double[nTargets];
		optAdvPayoff = new double[2 * nTargets];
		weight = new double[nTargets];
	}
	public void setXStar(double[] xStar)
	{
		this.xStar = xStar;
	}
	public void setWeight(double[] weight)
	{
		for(int i = 0; i < nTargets; i++)
			this.weight[i] = weight[i];
	}
	public void solve()
	{
		maxRegret = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < nTargets; i++)
			solve(i, maxRegret);
	}
	
	public void solve(int targetID, double curRegret)
	{
		double regret = Double.NEGATIVE_INFINITY;
		double[] advPayoff = new double[2 * nTargets];
		
		for(int i = 0; i < nTargets; i++)
		{
			if(i != targetID)
			{
				advPayoff[i] = advPayoffLB[i];
				advPayoff[nTargets + i] = advPayoffLB[nTargets + i];
			}
			else
				advPayoff[i] = advPayoffUB[i];
		}
		// Constant C
		double C = Double.NEGATIVE_INFINITY;
		double advPenaltyLB = 0.0;
		for(int i = 0; i < nTargets; i++)
		{
			if(i != targetID)
			{
				double advEU = xStar[i] * (advPayoffLB[nTargets + i] - advPayoffLB[i])
						+ advPayoffLB[i];
				C = C < advEU ? advEU : C;
			}
		}
		
		if(xStar[targetID] == 0)
		{
			if(C > advPayoffUB[targetID])
				return;
			else advPenaltyLB = advPayoffLB[nTargets + targetID];
		}
		else
		{
			double temp = (C - (1 - xStar[targetID]) * advPayoffUB[targetID]) / xStar[targetID];
			advPenaltyLB = advPayoffLB[nTargets + targetID] < temp 
					? temp : advPayoffLB[nTargets + targetID];
		}
		
		if(advPenaltyLB > advPayoffUB[targetID + nTargets])
			return;
		
		// Narrow down possible attack sets
		int[] sortID = sortID(advPayoff);
		advPayoff[nTargets + targetID] = advPayoffUB[nTargets + targetID];
		int finalIdxUB = computeAttackSet(advPayoff, sortID);
		advPayoff[nTargets + targetID] = advPenaltyLB;
		int finalIdxLB = computeAttackSet(advPayoff, sortID);
		
		// Iterate all possible attack sets
		for(int idx = finalIdxUB; idx <= finalIdxLB; idx++)
		{
//			System.out.println(maxRegret);
			double c = 0;
			double d = nRes;
		
			for(int i = 0; i <= idx; i++)
			{
				if(sortID[i] != targetID)
				{
//					c -= 1 / (advPayoffLB[sortID[i]] - advPayoffLB[nTargets + sortID[i]]);
//					d -= advPayoffLB[sortID[i]] / (advPayoffLB[sortID[i]] - advPayoffLB[nTargets + sortID[i]]);
					c -= 1 / (advPayoffLB[sortID[i]] - advPayoffLB[nTargets + sortID[i]]) * weight[sortID[i]];
					d -= advPayoffLB[sortID[i]] / (advPayoffLB[sortID[i]] - advPayoffLB[nTargets + sortID[i]]) * weight[sortID[i]];
				}
			}
			c /= weight[targetID];
			d /= weight[targetID];
			
			double vLB = Double.NEGATIVE_INFINITY;
			double vUB = Double.POSITIVE_INFINITY;
			for(int i = 0; i <= idx; i++)
			{
				if(sortID[i] != targetID)
				{
					vLB = vLB > -advPayoffLB[sortID[i]] ? vLB : -advPayoffLB[sortID[i]];
					vUB = vUB < -advPayoffLB[sortID[i] + nTargets] ? vUB : -advPayoffLB[sortID[i] + nTargets];
				}
			}
			for(int i = idx + 1; i < nTargets; i++)
			{
				vUB = vUB < -advPayoffLB[sortID[i]] ? vUB : -advPayoffLB[sortID[i]];
			}
			if(vUB < vLB)
				continue;
			// Consider a trivial case when cv + d >= 1
			double tempLB = Math.max((d - 1) / c, advPenaltyLB); // LB of penalty
			tempLB = Math.max(tempLB, -vUB);
			double tempUB = Math.min(-vLB, advPayoffUB[targetID + nTargets]); // UB of penalty
			if(tempLB <= tempUB)
			{
				regret = tempLB * (xStar[targetID] - 1) + (1 - xStar[targetID]) * advPayoffUB[targetID];
				if(regret > curRegret)
				{
//					if(regret > 2.93)
//						System.out.println("Testing");
					for(int i = 0; i <= idx; i++)
					{
						if(sortID[i] != targetID)
						{
							optDefCov[sortID[i]] = (advPayoffLB[sortID[i]] - tempLB) / (advPayoffLB[sortID[i]]
									- advPayoffLB[sortID[i] + nTargets]); 
						}
						else
							optDefCov[targetID] = 1.0;
					}
					for(int i = idx + 1; i < nTargets; i++)
					{
						optDefCov[sortID[i]] = 0.0;
					}
					advPayoff[targetID + nTargets] = tempLB;
					curRegret = regret;
					
					for(int i = 0; i < 2 * nTargets; i++)
						optAdvPayoff[i] = advPayoff[i];
					maxRegret = regret;
				}
			}
			// Consider a trivial case when cv + d = 0
			if(-c * advPayoffUB[targetID] + d == 0 && vLB <= -advPayoffUB[targetID] && -advPayoffUB[targetID] <= vUB)
			{
				regret = -advPayoffUB[targetID] + advPayoffUB[targetID + nTargets] * xStar[targetID] 
						+ (1 - xStar[targetID]) * advPayoffUB[targetID];
				if(regret > curRegret)
				{
//					if(regret > 2.93)
//						System.out.println("Testing");
					for(int i = 0; i <= idx; i++)
					{
						if(sortID[i] != targetID)
						{
							optDefCov[sortID[i]] = (advPayoffLB[sortID[i]] - advPayoffUB[targetID]) / (advPayoffLB[sortID[i]]
									- advPayoffLB[sortID[i] + nTargets]); 
						}
						else
							optDefCov[targetID] = 0.0;
					}
					for(int i = idx + 1; i < nTargets; i++)
					{
						optDefCov[sortID[i]] = 0.0;
					}
					advPayoff[targetID + nTargets] = advPayoffUB[targetID + nTargets];
					curRegret = regret;
					
					for(int i = 0; i < 2 * nTargets; i++)
						optAdvPayoff[i] = advPayoff[i];
					maxRegret = regret;
				}
			}
			// Consider non-trivial case
			if(c < 0)
				vLB = Math.max(vLB, (1 - d) / c);
			else if (d > 1 || d < 0)
				continue;
			double temp1 = c * advPayoffUB[targetID] - 1 - advPenaltyLB * c;
			double temp2 = advPenaltyLB * d + (1 - d) * advPayoffUB[targetID];
			if(temp1 > 0)
				vLB = Math.max(vLB, temp2 / temp1);
			else if(temp1 < 0)
				vUB = Math.min(vUB, temp2 / temp1);
			else if(temp2 > 0)
				continue;
			
			temp1 = c * advPayoffUB[targetID] - 1 - advPayoffUB[targetID + nTargets] * c;
			temp2 = advPayoffUB[targetID + nTargets] * d + (1 - d) * advPayoffUB[targetID];
			if(temp1 < 0)
				vLB = Math.max(vLB, temp2 / temp1);
			else if(temp1 > 0)
				vUB = Math.min(vUB, temp2 / temp1);
			else if(temp2 < 0)
				continue;
			
			if(vLB > vUB)
				continue;
			double a = c * advPayoffUB[targetID] - xStar[targetID];
			double b = (d - xStar[targetID]) * advPayoffUB[targetID];
			double temp = b * c - a * d;
			if(temp <= 0)
			{
				regret = vUB + (a * vUB + b) / (c * vUB + d);
				if(regret > curRegret)
				{
//					if(regret > 2.93)
//						System.out.println("Testing");
					curRegret = regret;
					advPayoff[targetID + nTargets] = (-vUB - (1 - c * vUB - d) * advPayoffUB[targetID]) / (c * vUB + d);
					for(int i = 0; i < 2 * nTargets; i++)
						optAdvPayoff[i] = advPayoff[i];
					maxRegret = regret;
					for(int i = 0; i <= idx; i++)
					{
						if(sortID[i] != targetID)
						{
							optDefCov[sortID[i]] = (advPayoffLB[sortID[i]] + vUB) / (advPayoffLB[sortID[i]]
									- advPayoffLB[sortID[i] + nTargets]); 
						}
						else
							optDefCov[targetID] = c * vUB + d;
					}
					for(int i = idx + 1; i < nTargets; i++)
					{
						optDefCov[sortID[i]] = 0.0;
					}
					
				}
			}
			else
			{
				regret = 0.0;
				double v1 = (Math.sqrt(b * c - a * d) - d) / c;
				double v2 = (-Math.sqrt(b * c - a * d) - d) / c;
				double v = vUB;
				regret = vUB + (a * vUB + b) / (c * vUB + d);
				if(regret < vLB + (a * vLB + b) / (c * vLB + d))
				{
					regret = vLB + (a * vLB + b) / (c * vLB + d);
					v = vLB;
				}
				if(v1 >= vLB && v1 <= vUB && regret < v1 + (a * v1 + b) / (c * v1 + d))
				{
					regret = v1 + (a * v1 + b) / (c * v1 + d);
					v = v1;
				}
				if(v2 >= vLB && v2 <= vUB && regret < v2 + (a * v2 + b) / (c * v2 + d))
				{
					regret = v2 + (a * v2 + b) / (c * v2 + d);
					v = v2;
				}
				if(regret > curRegret)
				{
//					if(regret > 2.93)
//						System.out.println("Testing");
					curRegret = regret;
					advPayoff[targetID + nTargets] = ((c * advPayoffUB[targetID] - 1) * v - 
					                                  (1 - d) * advPayoffUB[targetID]) / (c * v + d);
					for(int i = 0; i < 2 * nTargets; i++)
						optAdvPayoff[i] = advPayoff[i];
					maxRegret = regret;
					for(int i = 0; i <= idx; i++)
					{
						if(sortID[i] != targetID)
						{
							optDefCov[sortID[i]] = (advPayoffLB[sortID[i]] + v) / (advPayoffLB[sortID[i]]
									- advPayoffLB[sortID[i] + nTargets]); 
						}
						else
							optDefCov[targetID] = c * v + d;
					}
					for(int i = idx + 1; i < nTargets; i++)
					{
						optDefCov[sortID[i]] = 0.0;
					}
				}
			}
		}
		
//		return defCov;
	}
	public int[] sortID(double[] advPayoff)
	{
		int[] idx = new int[nTargets];
		for(int i = 0; i < nTargets; i++)
			idx[i] = i;
		double[] values = new double[nTargets];
		for(int i = 0; i < nTargets; i++)
			values[i] = advPayoff[idx[i]];
		QuickSort my_sort = new QuickSort(values, idx, nTargets);
		my_sort.get_sort_id();
		int[] sort_idx = my_sort.sort_id;
		for(int i = 0; i < nTargets / 2; i++)
		{
			int temp = sort_idx[i];
			sort_idx[i] = sort_idx[nTargets - 1 - i];
			sort_idx[nTargets - 1 - i] = temp;
		}
		return sort_idx;
	}
	public int computeAttackSet(double[] advPayoff, int[] sort_idx)
	{
		double[] minCov = new double[nTargets];
		for(int i = 0; i < nTargets; i++)
			minCov[i] = 0.0;

		double r = advPayoff[sort_idx[0]];
		double res = nRes;
		boolean isStop = false;

		int final_idx = 1;
		for(int i = 1; i < nTargets; i++)
		{
			final_idx = i;
			double[] dis = new double[i];
			for(int j = 0; j < i; j++)
			{
				dis[j] = (advPayoff[sort_idx[i]] - r) /
						(advPayoff[nTargets + sort_idx[j]] - advPayoff[sort_idx[j]]);
				if(dis[j] < 0)
					System.out.println("Oops");
//				res = res - dis[j];
				res = res - dis[j] * weight[sort_idx[j]];
				minCov[sort_idx[j]] += dis[j];
				if(res < 0 || minCov[sort_idx[j]] > 1) 
				{
					for(int k = 0; k <= j; k++)
					{
//						res = res + dis[k];
						res = res + dis[k] * weight[sort_idx[k]];
						minCov[sort_idx[k]] -= dis[k];
					}
					isStop = true;
					break;
				}
			}
			if(!isStop)
				r = advPayoff[sort_idx[i]];
			else break;
		}
		for(int j = 0; j < nTargets; j++)
			if(minCov[j] > 1 || minCov[j] < 0)
				System.out.println("Oops");
		
		if(!isStop && nTargets > 1) final_idx++;
		return final_idx - 1;
	}
}
