package algorithms;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import lpWrapper.LPSolverException;
import models.PayoffStructure;
import models.SUQRAdversary;
import models.Target;

public class PASAQ_SUQR{
	PASAQ_Piecewise updateBoundsObj;

	public String modelName = "SUQRModel";
	
	protected double upperBound; 
	protected double lowerBound;
	protected boolean boundIsSet = false;

	protected double binarySearchThreshold = 0.0001;
	protected int numPiece;
	protected double numResources;
	protected List<Target> targetList;
	int payoffIndex = 0;
	
	protected String lpFolder;
	protected boolean debugPrint = true;
	
	protected SUQRAdversary adversary;
    
	public PASAQ_SUQR(List<Target> targetList, double numResources, int numPiece, SUQRAdversary adversary){
		this.targetList = targetList;
		this.numPiece = numPiece;
		this.numResources = numResources;
		this.adversary = adversary;
		
		try{
			updateBoundsObj = new PASAQ_Piecewise(targetList, numResources, numPiece, adversary);
		}
		catch(IOException e){
			e.printStackTrace();
		}
	}
	public void setPayoffIndex(int payoffIndex)
	{
		this.payoffIndex = payoffIndex;
		updateBoundsObj.setPayoffIndex(payoffIndex);
	}
	public void setBinarySearchThreshold(double threshold){
		this.binarySearchThreshold = threshold;		
	}
	
	public void setLPFolder(String lpFolder) {
		this.lpFolder = lpFolder;
	}
	
	public void setBounds(double lowerBound, double upperBound){
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.boundIsSet = true;		
	}	
	
	public void loadProblem() throws IOException{
		updateBoundsObj.loadProblem();
	}
	
	public void addConstraint(double[] coeff, double constant)
	{
		updateBoundsObj.addConstraint(coeff, constant);
	}
	public void solve() throws LPSolverException{
		double estimatedOpt = Double.NEGATIVE_INFINITY;		
	
		while(upperBound - lowerBound >= binarySearchThreshold){
//			System.gc();
//			System.out.println(lowerBound + "\t" + upperBound);
			estimatedOpt = 0.5 * (upperBound + lowerBound); 
			updateBoundsObj.setEstimatedOpt(estimatedOpt);
			updateBoundsObj.solve();
//			double delta = updateBoundsObj.getObjValue();
			double delta = updateBoundsObj.getLPObjective();
			if(delta < 0.0){
				lowerBound = estimatedOpt;
			}
			else{
				upperBound = estimatedOpt;
			}
		}
//		if(updateBoundsObj.getObjValue() > 0.0){
		if(updateBoundsObj.getLPObjective() >= 0.0){
			estimatedOpt = lowerBound; 
			updateBoundsObj.setEstimatedOpt(estimatedOpt);
			updateBoundsObj.solve();
		}
//		System.out.println("Defender EU: " + estimatedOpt + "\t" + getDefenderPayoff());
//		System.out.println("Defender sampled strategy: " + getDefenderCoverageByTargetIds());
		return;
	}

	public double getDefenderPayoff(Map<Integer, Double> defenderCoverage){
		Map<Target, Double> mapAttSUQR = updateBoundsObj.getAttackerQR(adversary, defenderCoverage);
		
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
	
	public double getDefenderPayoff(){
		Map<Integer, Double> mapDefCov = updateBoundsObj.getDefenderCoverageByTargetIds();
		Map<Target, Double> mapAttSUQR = updateBoundsObj.getAttackerQR(adversary);
		
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
	
	public Map<Target, Double> getAttackerAction(SUQRAdversary type){
		return updateBoundsObj.getAttackerQR(type);		
	}
	
	public Map<Integer, Double> getAttackerActionByTargetId(SUQRAdversary type, Map<Integer, Double> defenderCoverage){
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
	
	public Map<Integer, Double> getDefenderPayoffByTargetIds(){
		return updateBoundsObj.getDefenderPayoffByTargetIds();
	}
	
	public void end(){
		this.updateBoundsObj.end();
	}
}