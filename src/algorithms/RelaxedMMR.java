package algorithms;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import lpWrapper.LPSolverException;
import models.PayoffStructure;
import models.SUQRAdversary;
import models.Target;

public class RelaxedMMR{
	RelaxedMMR_Piecewise updateBoundsObj;
	PASAQ_SUQR pasaq; // only used for generating samples

	public String modelName = "SUQRModel";
	
	protected double upperBound; 
	protected double lowerBound;

	protected double binarySearchThreshold = 1.0e-4;
	protected int numPiece;
	protected double numResources;
	protected List<Target> targetList;
	public List<Double> sampledOptList; // updated when new payoff is added
	protected SUQRAdversary adversary;
	double runtime = 0.0;
    
	public RelaxedMMR(List<Target> targetList, double numResources, int numPiece, SUQRAdversary adversary){
		this.targetList = targetList;
		this.numPiece = numPiece;
		this.numResources = numResources;
		this.adversary = adversary;
		try{
			pasaq = new PASAQ_SUQR(targetList, numResources, numPiece, adversary);
			initSampledOptList();
			updateBoundsObj = new RelaxedMMR_Piecewise(targetList, sampledOptList, numResources, numPiece, adversary);
		}
		catch(IOException e){
			e.printStackTrace();
		}
		
	}
	public RelaxedMMR(List<Target> targetList, int numResources, int numPiece, SUQRAdversary adversary, boolean random){
		this.targetList = targetList;
		this.numPiece = numPiece;
		this.numResources = numResources;
		this.adversary = adversary;
		try{
			
			pasaq = new PASAQ_SUQR(targetList, numResources, numPiece, adversary);
			if(!random)
				initSampledOptList();
			else
				initSampledOptListRandom();
			updateBoundsObj = new RelaxedMMR_Piecewise(targetList, sampledOptList, numResources, numPiece, adversary);
		}
		catch(IOException e){
			e.printStackTrace();
		}
		
	}
	public void initSampledOptList()
	{
//		System.out.println("Start initialize samples... ");
		long start = System.currentTimeMillis();
		sampledOptList = new ArrayList<Double>();
			
		int numPayoffs = targetList.get(0).getNumPayoff();
		try {
			pasaq.loadProblem();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for(int i = 0; i < numPayoffs; i++)
		{
			// Set up lower bound and upper bound
			double lb = Double.POSITIVE_INFINITY;
			double ub = Double.NEGATIVE_INFINITY;
			for(Target t : targetList)
			{
				double defReward = t.getPayoffList().get(i).getDefenderReward();
				double defPenalty = t.getPayoffList().get(i).getDefenderPenalty();
				if(lb > defPenalty)
					lb = defPenalty;
				if(ub < defReward)
					ub = defReward;
			}
			// End set up lower bound and upper bound
			pasaq.setPayoffIndex(i);
			pasaq.setBounds(lb, ub);
			try {
				pasaq.solve();
			} catch (LPSolverException e) {
				e.printStackTrace();
			}
			sampledOptList.add(pasaq.getDefenderPayoff());
		}
		long end = System.currentTimeMillis();
		runtime = (end - start) / 1000.0;
//		System.out.println("End initialize samples.");
	}
	
	public void initSampledOptListRandom()
	{
//		System.out.println("Start initialize samples... ");
		long start = System.currentTimeMillis();
		sampledOptList = new ArrayList<Double>();
			
		int numPayoffs = targetList.get(0).getNumPayoff();
		Random rand = new Random();
		int nTargets = targetList.size();
		double[] defCov = new double[nTargets];
		double[] attProb = new double[nTargets];
		for(int i = 0; i < numPayoffs; i++)
		{
			double sum = 0.0;
			for(Target t : targetList)
			{
				defCov[t.id] = rand.nextDouble();
				sum += defCov[t.id];
			}
			for(Target t : targetList)
			{
				defCov[t.id] = defCov[t.id] * numResources / sum;
				if(defCov[t.id] > 1.0)
					defCov[t.id] = 1.0;
			}
			
			sum = 0.0;
			for(Target t : targetList)
			{
				double attReward = t.getPayoffList().get(i).getAdversaryReward();
				double attPenalty = t.getPayoffList().get(i).getAdversaryPenalty();
				double attUtility = adversary.w[adversary.DEF_COV] * defCov[t.id] + adversary.w[adversary.ATT_REWARD] * attReward
						+ adversary.w[adversary.ATT_PENALTY] * attPenalty;
				attProb[t.id] = Math.exp(attUtility);
				sum += attProb[t.id];
			}
			
			double defUtility = 0.0;
			for(Target t : targetList)
			{
				double defReward = t.getPayoffList().get(i).getDefenderReward();
				double defPenalty = t.getPayoffList().get(i).getDefenderPenalty();
				double defEU = defReward * defCov[t.id] + defPenalty * (1 - defCov[t.id]);
				defUtility += defEU * attProb[t.id] / sum;
			}
			sampledOptList.add(defUtility);
		}
		long end = System.currentTimeMillis();
		runtime = (end - start) / 1000.0;
//		System.out.println("End initialize samples.");
	}
	public void updateSampledOptList()
	{
		int curSize = sampledOptList.size();
		int numPayoffs = targetList.get(0).getNumPayoff();
		for(int i = curSize; i < numPayoffs; i++)
		{
			// Set up lower bound and upper bound
			double lb = Double.POSITIVE_INFINITY;
			double ub = Double.NEGATIVE_INFINITY;
			for(Target t : targetList)
			{
				double defReward = t.getPayoffList().get(i).getDefenderReward();
				double defPenalty = t.getPayoffList().get(i).getDefenderPenalty();
				if(lb > defPenalty)
					lb = defPenalty;
				if(ub < defReward)
					ub = defReward;
			}
			// End set up lower bound and upper bound
			pasaq.setPayoffIndex(i);
			pasaq.setBounds(lb,ub);
			try {
				pasaq.solve();
			} catch (LPSolverException e) {
				e.printStackTrace();
			}
			sampledOptList.add(pasaq.getDefenderPayoff());
		}
	}
	public List<Double> getSampledOptList()
	{
		return this.sampledOptList;
	}
	
	public void setBinarySearchThreshold(double threshold){
		this.binarySearchThreshold = threshold;		
	}
	
	public double getBinarySearchThreshold(){
		return this.binarySearchThreshold;
	}
	
	public void setBounds(double lowerBound, double upperBound){
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}	
	
	public double getCurEstimation(){
		return updateBoundsObj.estimatedOpt;
	}
	
	public void loadProblem() throws IOException{
		updateBoundsObj.loadProblem();
	}
	public void addConstraint(double[] coeff, double constant)
	{
		// reset all samples
		pasaq.addConstraint(coeff, constant);
		sampledOptList.clear();
		updateSampledOptList();
		updateBoundsObj.addConstraint(coeff, constant);
	}
	public void deletePayoffConstraint()
	{
		updateBoundsObj.deletePayoffConstraint();
	}
	public void solve() throws LPSolverException
	{
		double estimatedOpt = Double.MAX_VALUE;		
		int iter = 0;
		while(upperBound - lowerBound >= binarySearchThreshold)
		{
			iter++;
			System.gc();
//			System.out.println("Iter: " + iter + "\tLower bound: " + lowerBound + "\tUpper bound: " + upperBound); 
			estimatedOpt = 0.5 * (upperBound + lowerBound); 
			updateBoundsObj.setEstimatedOptList(estimatedOpt);
			updateBoundsObj.solve();
			
			// check the minimum value of the binary search objective function.			
			double delta = updateBoundsObj.getLPObjective();
//			double delta = updateBoundsObj.getObjValue();
			if(delta >= 0.0){
				lowerBound = estimatedOpt;
			}
			else{
				upperBound = estimatedOpt;
			}
		}
		
//		if(updateBoundsObj.getObjValue() > 0.0){
		if(updateBoundsObj.getLPObjective() >= 0.0){
			estimatedOpt = upperBound; 
			updateBoundsObj.setEstimatedOptList(estimatedOpt);
			updateBoundsObj.solve();
		}
		
//		System.out.println("Relaxed MMR: " + estimatedOpt + "\t" + getRegret());
//		for(int i = 0; i < sampledOptList.size(); i++)
//			System.out.println("Regret of sample " + i + "\t" + (getDefenderPayoff(i) - sampledOptList.get(i)));
//			
//		System.out.println("Defender optimal strategy: " + getDefenderCoverageByTargetIds());
	}
	
	public double getDefenderPayoff(Map<Integer, Double> defenderCoverage, int payoffIndex){
		Map<Target, Double> mapAttSUQR = updateBoundsObj.getAttackerQR(adversary, defenderCoverage, payoffIndex);
		double defenderEU = 0.0;
		
		for(Target t : targetList){
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			double defCov = 0.0;
			if(defenderCoverage.containsKey(t.id)){
				defCov = defenderCoverage.get(t.id);
			}
			defenderEU += 
				(defCov * payoffs.getDefenderReward() + (1 - defCov) * payoffs.getDefenderPenalty()) * mapAttSUQR.get(t);
		}
		return defenderEU;
	}
	
	public double getDefenderPayoff(int payoffIndex){
		Map<Integer, Double> mapDefCov = updateBoundsObj.getDefenderCoverageByTargetIds();
		Map<Target, Double> mapAttSUQR = updateBoundsObj.getAttackerQR(adversary, payoffIndex);
		double defenderEU = 0.0;
		for(Target t : targetList){
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			defenderEU += 
				(mapDefCov.get(t.id) * payoffs.getDefenderReward() + (1 - mapDefCov.get(t.id)) * payoffs.getDefenderPenalty()) * mapAttSUQR.get(t);
		}
		return defenderEU;
	}
	
	public void writeProb(String fileName){
		this.updateBoundsObj.writeProb(fileName);
	}
		
	public double getDefenderCoverage(int targetId){
		return updateBoundsObj.getDefenderCoverage(targetId) ;
	}

	public double getDefenderCoverage(Target target){
		return updateBoundsObj.getDefenderCoverage(target);		
	}
	
	
	public Map<Integer, Double> getAttackerActionByTargetId(SUQRAdversary type, Map<Integer, Double> defenderCoverage, int payoffIndex){
		Map<Integer, Double> mapAttackerSUQR = new HashMap<Integer, Double>();
		
		double suqr_cum = 0.0;
		double suqr_single = 0.0;
		
		for(Target t : targetList){
			double defCov = 0;
			
			if(defenderCoverage.containsKey(t.id)){
				defCov = defenderCoverage.get(t.id);
			}
			
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			
			double w1 = type.w[SUQRAdversary.DEF_COV] * defCov;
			double w2 = type.w[SUQRAdversary.ATT_REWARD] * payoffs.getAdversaryReward();
			double w3 = type.w[SUQRAdversary.ATT_PENALTY] * payoffs.getAdversaryPenalty();
			
			suqr_single = Math.exp(w1 + w2 + w3);
			mapAttackerSUQR.put(t.id, suqr_single);
			suqr_cum += suqr_single;
		}
		
		for(Target t : targetList){
			suqr_single = mapAttackerSUQR.get(t.id);
			mapAttackerSUQR.remove(t.id);
			mapAttackerSUQR.put(t.id, suqr_single / suqr_cum);
		}
		
		return mapAttackerSUQR;	
	}
	
	public Map<Target, Double> getDefenderCoverage(){
		return updateBoundsObj.getDefenderCoverage();
	}
	
	public Map<Integer, Double> getDefenderCoverageByTargetIds(){
		return updateBoundsObj.getDefenderCoverageByTargetIds();
	}
	
	public Map<Integer, Double> getDefenderPayoffByTargetIds(int payoffIndex){
		return updateBoundsObj.getDefenderPayoffByTargetIds(payoffIndex);
	}
	
	public double getRegret()
	{
		double regret = Double.NEGATIVE_INFINITY;
		int numPayoff = targetList.get(0).getNumPayoff();
		for(int i = 0; i < numPayoff; i++)
		{
			double temp = sampledOptList.get(i) - getDefenderPayoff(i);
			if(regret < temp)
				regret = temp;
		}
		return regret;
	}
	
	public double getSampleRuntime()
	{
		return runtime;
	}
	public void clearOptList()
	{
		sampledOptList.clear();
	}
	public void end(){
		if(sampledOptList != null)
		{
			sampledOptList.clear();
//			sampledOptList = null;
		}
		this.updateBoundsObj.end();
		if(this.pasaq != null)
			this.pasaq.end();
	}
}