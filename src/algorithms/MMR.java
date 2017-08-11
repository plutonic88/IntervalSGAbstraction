package algorithms;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;

import lpWrapper.LPSolverException;
import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import models.AttBoundStructure;
import models.PayoffStructure;
import models.SUQRAdversary;
import models.Target;

public class MMR {
	RelaxedMMR relaxMMR;
	MR mr = null;
	List<Target> targetList;
	SUQRAdversary adversary;
	double lb = 0; // update later
	double ub = 10; // update later
	double numRes;
	boolean isZeroSum = true;
	int numSamples = 1;
	double estimatedOpt = 10;
	double thres = 1.0e-2;
	int numPiece = 20;
	int numIter = 5;
	int numRound = 30;
	
	// output
	double[] runtime;
	double[] LB;
	double[] UB;
	double[] trueLB;
	double[][] defStrategy;
	int convergeRound;
	public MMR(List<Target> targetList, SUQRAdversary adversary, double numRes, int numSamples, boolean isZeroSum)
	{
		this.targetList = targetList;
		this.numRes = numRes;
		this.isZeroSum = isZeroSum;
		this.numSamples = numSamples;
		this.adversary = adversary;
		initSample();
		relaxMMR = new RelaxedMMR(targetList, numRes, numPiece, adversary);
		// for max regret only
		double[] advPayoffLB = new double[2 * targetList.size()];
		double[] advPayoffUB = new double[2 * targetList.size()];
		for(Target t : targetList)
		{
			AttBoundStructure attBound = t.getAttBoundStructure();
			advPayoffLB[t.id] = attBound.getAttRewardLB();
			advPayoffLB[t.id + targetList.size()] = attBound.getAttPenaltyLB();
			advPayoffUB[t.id] = attBound.getAttRewardUB();
			advPayoffUB[t.id + targetList.size()] = attBound.getAttPenaltyUB();
			
		}
		mr = new MR(targetList, null, adversary, numRes, numIter);
		
		// Initialize storing
		if(numRound > 0)
			setNumRound(numRound);
		convergeRound = -1;
				
	}
	public void initSample()
	{
		Random rand = new Random();
		if(isZeroSum)
		{
			for(Target t : targetList)
			{
				List<PayoffStructure> payoffList = new ArrayList<PayoffStructure>();
				AttBoundStructure attBoundStructure = t.getAttBoundStructure();
				double attRewardLB = attBoundStructure.getAttRewardLB();
				double attPenaltyLB = attBoundStructure.getAttPenaltyLB();;
				double attRewardUB = attBoundStructure.getAttRewardUB();
				double attPenaltyUB = attBoundStructure.getAttPenaltyUB();
				for(int i = 0; i < numSamples; i++)
				{
//					double pivotReward = rand.nextDouble();
					double pivotReward = 0.5;
//					double attReward = pivotReward < 0.5 ? attRewardLB : attRewardUB;
					double attReward = pivotReward * (attRewardUB - attRewardLB) + attRewardLB;
//					double pivotPenalty = rand.nextDouble();
					double pivotPenalty = 0.5;
//					double attPenalty = pivotPenalty < 0.5 ? attPenaltyLB : attPenaltyUB;
					double attPenalty = pivotPenalty * (attPenaltyUB - attPenaltyLB) + attPenaltyLB;
					PayoffStructure payoffStructure = new PayoffStructure(-attPenalty, -attReward, attReward, attPenalty);
					payoffList.add(payoffStructure);
				}
				if(t.getPayoffList() != null && !t.getPayoffList().isEmpty())
				{
					t.getPayoffList().clear();
				}
				t.setPayoffList(payoffList);
			}
		}
	}
	
	public void setNumRound(int numRound)
	{
		this.numRound = numRound;
		LB = new double[numRound];
		UB = new double[numRound];
		trueLB = new double[numRound];
		runtime = new double[numRound];
		defStrategy = new double[numRound][targetList.size()];
	}
	public void setNumSample(int numSamples)
	{
		this.numSamples = numSamples;
	}
	void setThreshold(double thres)
	{
		this.thres = thres;
	}
	public void loadProblem() throws IOException{
		relaxMMR.loadProblem();
	}
	
	public double initUB()
	{
		double ub = Double.NEGATIVE_INFINITY;
		double maxReward = Double.NEGATIVE_INFINITY;;
		double minPenalty = Double.POSITIVE_INFINITY;;
		
		for(Target t : targetList)
		{
			double tempReward;
			double tempPenalty; 
			if(isZeroSum)
			{
				tempReward = t.getAttBoundStructure().getAttRewardUB();
				tempPenalty = t.getAttBoundStructure().getAttPenaltyLB();
			}
			else
			{
				tempReward = t.getPayoffStructure(0).getDefenderReward();
				tempPenalty = t.getPayoffStructure(0).getDefenderPenalty();
			}
			maxReward = maxReward < tempReward ? tempReward : maxReward;
			minPenalty = minPenalty > tempPenalty ? tempPenalty : minPenalty;
		}
		ub = maxReward - minPenalty;
		return ub;
	}
	public void addConstraint(double[] coeff, double constant)
	{
		// free all current sample (maybe not???)
		// To do: re-initialize all payoff samples if necessary
		initSample();
//		System.out.println("Testing: " +  targetList.get(0).getPayoffList().size());
		// Add constraint
		relaxMMR.addConstraint(coeff, constant);
		mr.addConstraint(coeff, constant);
		
	}
	public void deletePayoffConstraint()
	{
		relaxMMR.deletePayoffConstraint();
	}
	public void solve() throws LPSolverException, MatlabConnectionException, MatlabInvocationException
	{
		double lb = 0;
		double ub = initUB();
		
		relaxMMR.setBounds(lb, ub);
		int iter = 0;
//		while(ub - lb > thres)
		while(iter < numRound)
		{
			
			long start = System.currentTimeMillis();
			System.out.println("Iter: " + iter + "\tLower bound: " + lb + "\tUpper bound: " + ub); 
			relaxMMR.solve();
//			lb = relaxMMR.getCurEstimation();
			lb = relaxMMR.getRegret();
			trueLB[iter] = relaxMMR.getRegret();
			LB[iter] = lb;
			Map<Integer, Double>  tempDefStrategy = relaxMMR.getDefenderCoverageByTargetIds();
			for(int t = 0; t < targetList.size(); t++)
			{
				defStrategy[iter][t] = tempDefStrategy.get(t);
			}
			tempDefStrategy.clear();
			Map<Target, Double> defCov = relaxMMR.getDefenderCoverage();
			double[] xStar = new double[targetList.size()];
			for(Target t : targetList)
				xStar[t.id] = defCov.get(t);
			if(iter == numRound - 1)
			{
				mr.setNumIter(30);
				mr.setDefCov(targetList, defCov);
				mr.solve();
				ub = mr.getMaxRegret();
			}
			else 
			{
				mr.setNumIter(1);
				mr.setDefCov(targetList, defCov);
				mr.solve();
				ub = mr.getMaxRegret();
				for(Target t : targetList)
				{
					PayoffStructure newPayoff = null;
					newPayoff = new PayoffStructure(-mr.getAttPenalty(t), -mr.getAttReward(t)
							, mr.getAttReward(t), mr.getAttPenalty(t));
					t.getPayoffList().add(newPayoff);
				}
				mr.solve();
				for(Target t : targetList)
				{
					PayoffStructure newPayoff = null;
					newPayoff = new PayoffStructure(-mr.getAttPenalty(t), -mr.getAttReward(t)
							, mr.getAttReward(t), mr.getAttPenalty(t));
					t.getPayoffList().add(newPayoff);
				}
				ub = Math.max(ub, mr.getMaxRegret());
			}
			UB[iter] = ub;

			relaxMMR.updateSampledOptList();
			relaxMMR.setBounds(lb, ub);
			if(ub - lb <= thres)
				relaxMMR.setBounds(lb, lb);
			long end = System.currentTimeMillis();
			if(iter == 0)
				runtime[iter] = (end - start) / 1000.0;
			else
				runtime[iter] = runtime[iter - 1] + (end - start) / 1000.0;
			
			if(ub - lb <= thres && convergeRound == -1)
				convergeRound = iter;
			if(iter == numRound - 1)
			{
				System.out.println("Iter: " + iter + "\tLower bound: " + LB[iter] + "\tUpper bound: " + UB[iter]); 
				System.out.println("Runtime: " + runtime[iter]);
			}
			iter++;
		}
	}
	public void saveResults(String folderPath, String runtimeFileName, String regretFileName, String strategyFileName, String convergeFileName)
	{
		try{
			FileOutputStream runtimeFileStream = new FileOutputStream(folderPath + runtimeFileName);
			PrintStream runtimeStream = new PrintStream(runtimeFileStream);
			for(int i = 0; i < numRound; i++)
				runtimeStream.println(runtime[i]);
			
			FileOutputStream regretFileStream = new FileOutputStream(folderPath + regretFileName);
			PrintStream regretStream = new PrintStream(regretFileStream);
			for(int i = 0; i < numRound; i++)
				regretStream.println(UB[i] + "," + LB[i] + "," + trueLB[i]);
			
			FileOutputStream strategyFileStream = new FileOutputStream(folderPath + strategyFileName);
			PrintStream strategyStream = new PrintStream(strategyFileStream);
			for(int i = 0; i < numRound; i++)
			{
				for(int t = 0; t < targetList.size(); t++)
				{
					strategyStream.print(defStrategy[i][t] + ",");
				}
				strategyStream.println();
			}
				
			FileOutputStream convergeFileStream = new FileOutputStream(folderPath + convergeFileName);
			PrintStream convergeStream = new PrintStream(convergeFileStream);
			convergeStream.println(convergeRound);
			
			runtimeStream.close();
			regretStream.close();
			strategyStream.close();
			convergeStream.close();
			
		}catch (FileNotFoundException e1)
		{
			e1.printStackTrace();
		}
	}
	
	public void saveSampleRuntime(String sampleFilePath)
	{
		try{
			FileOutputStream sampleFileStream = new FileOutputStream(sampleFilePath);
			PrintStream sampleStream = new PrintStream(sampleFileStream);
			sampleStream.println(getSampleRuntime());
			sampleStream.close();
			
		}catch (FileNotFoundException e1)
		{
			e1.printStackTrace();
		}
	}
		
	public double getSampleRuntime()
	{
		return relaxMMR.getSampleRuntime();
	}
	public double getRegret()
	{
		return UB[numRound - 1];
	}
	public double getOptCov(int targetID)
	{
		return relaxMMR.getDefenderCoverage(targetID);
	}
	public double[] getOptCov()
	{
		double[] cov = new double[targetList.size()];
		Map<Integer, Double> tempCov = relaxMMR.getDefenderCoverageByTargetIds();
		Iterator it = tempCov.entrySet().iterator();
		while(it.hasNext())
		{
			Entry entry = (Map.Entry)it.next();
			cov[(Integer)entry.getKey()] = (Double)entry.getValue();
		}
		tempCov.clear();
		return cov;
	}
//	public double getOpponentCov(int targetID)
//	{
////		return MRzero.getOptDefCov(targetID);
//	}
	public double getOpponentReward(int targetID)
	{
		return mr.getAttReward(targetID);
	}
	public double getOpponentPenalty(int targetID)
	{
		return mr.getAttPenalty(targetID);
	}
	public double getLB()
	{
		return LB[numRound - 1];
	}
	public void clearPayoffList()
	{
		for(Target t : targetList)
		{
			t.getPayoffList().clear();
		}
	}
	public void end()
	{
		for(Target t : targetList)
		{
			t.getPayoffList().clear();
		}
		relaxMMR.end();
		mr.end();
	}
	public double[] getDefCov() {
		return defStrategy[numRound-1];
	}
}
