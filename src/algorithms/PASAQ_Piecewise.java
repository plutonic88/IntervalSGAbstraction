/*
 * @author Matthew Brown
 * 
 */
package algorithms;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import lpWrapper.Configuration;
import lpWrapper.MIProblem;
import models.PayoffStructure;
import models.SUQRAdversary;
import models.Target;

public class PASAQ_Piecewise extends MIProblem {
	SUQRAdversary adversary; // SUQR model
	List<Target> targetList; // list of targets
	HashMap<String, Integer> varMap; // for MILP, column variables
	HashMap<String, Integer> rowMap;	// for MILP, row constraints

	double estimatedOpt = -Configuration.MM;
	double[] ObjCoeff1; // for x * f(x)
	double[] ObjCoeff2; // for f(x)
	int numPiece = 100; // number of segments for piecewise linear approximation
	int payoffIndex = 0; // current sampled payoff
	
	double numResources; // number of defender resources
	
	int nConstraint = 0;
	
	/* Function PASAQ_Piecewise: initialize class variables
	 * Input: 
	 * 	- targetList: list of targets
	 *  - numResources: number of defender resources
	 *  - numPiece: number of segments for piecewise linear approximation
	 *  - adversary: SUQR parameters
	 *  Output: no output
	 *  - varMap & rowMap are initialized
	 *  - cplex are initialized
	 *  - ObjCoeff1 & ObjCoeff2 are computed
	 */
	public PASAQ_Piecewise(List<Target> targetList, double numResources, int numPiece, SUQRAdversary adversary) throws IOException {
		super();
		this.targetList = targetList;
		this.numResources = numResources;
		this.numPiece = numPiece;
		this.adversary = adversary;
		this.varMap = new HashMap<String, Integer>();
		this.rowMap = new HashMap<String, Integer>();
		initObjCoeff();
	}
	
	/* Function initObjCoeff: compute ObjCoeff1 & ObjCoeff2
	 * Input: no input
	 * Output: no output
	 */
	protected void initObjCoeff(){
		ObjCoeff1 = new double[numPiece];
		ObjCoeff2 = new double[numPiece];
		double pieceL = 1.0/numPiece;
		for(int k = 1; k <= numPiece; k++)
		{
			double w1EU1 = this.discountFunction(adversary.w[SUQRAdversary.DEF_COV] * k * pieceL);
			double w1EU2 = this.discountFunction(adversary.w[SUQRAdversary.DEF_COV] * (k - 1) * pieceL);
			ObjCoeff1[k - 1] = w1EU1 * k - w1EU2 * (k - 1);
			ObjCoeff2[k - 1] = (w1EU1 - w1EU2) / pieceL;
		}
	}
	public void setPayoffIndex(int payoffIndex)
	{
		this.payoffIndex = payoffIndex;
	}
	protected void setEstimatedOpt (double estimatedOpt) {
		this.estimatedOpt = estimatedOpt;	
		setObjCoeff(estimatedOpt);
	}
	public void addConstraint(double[] coeff, double constant) // coeff * x <= constant
	{
		nConstraint++;
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		for(Target t : targetList)
		{
			for(int k = 1; k <= numPiece; k++)
			{
				ja.add(varMap.get("x" + t.id + "_" + k));
				ar.add(coeff[t.id]);
			}
		}
		addAndSetRow("Constr" + nConstraint, BOUNDS_TYPE.UPPER, -Configuration.MM, constant);
		rowMap.put("Constr" + nConstraint, this.getNumRows());
		setMatRow(this.getNumRows(), ja, ar);
	}
	public void setObjCoeff(double estimatedOpt) {
		double total = 0.0;
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		ja.add(varMap.get("s"));
		ar.add(1.0);
		for(Target t : targetList)
		{
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			double w2 = adversary.w[SUQRAdversary.ATT_REWARD] * payoffs.getAdversaryReward();
			double w3 = adversary.w[SUQRAdversary.ATT_PENALTY] * payoffs.getAdversaryPenalty();
			double f0 = this.discountFunction(w2 + w3);
			for(int k = 1; k <= numPiece; k++)
			{
				
				double coef = f0 * ((estimatedOpt - payoffs.getDefenderPenalty()) * ObjCoeff2[k - 1]
						- (payoffs.getDefenderReward() - payoffs.getDefenderPenalty()) * ObjCoeff1[k - 1]);
				ja.add(varMap.get("x" + t.id + "_" + k));
				ar.add(-coef);
			}
			total += (estimatedOpt - payoffs.getDefenderPenalty()) * f0;
		}
		if(rowMap.containsKey("payoff_0"))
		{
			setMatRow(rowMap.get("payoff_0"), ja, ar);
			resetRowBound(rowMap.get("payoff_0"), BOUNDS_TYPE.LOWER, total, Configuration.MM);
		}
		else
		{
			addAndSetRow("payoff_0", BOUNDS_TYPE.LOWER, total, Configuration.MM);
			rowMap.put("payoff_0", this.getNumRows());
			setMatRow(this.getNumRows(), ja, ar);
		}
		
	}
	
	// The discounting function can be override by an extended class
	protected double discountFunction(double attEU) {
		// Defaulting discounting function is the exponential function
		return Math.exp(attEU);
	}

	@Override
	protected void setProblemType() {
		this.setProblemName("UpdateBinarySearchBound");
		this.setProblemType(PROBLEM_TYPE.MIP, OBJECTIVE_TYPE.MIN);
	}

	@Override
	protected void setColBounds() {
		// 1, ..., n x[t]
		int index = 1;
		double pieceL = 1.0/numPiece;
		for (Target t : targetList) {

		
			for (int k=1; k<=numPiece; k++){		
				addAndSetColumn("x" + t.id + "_" + k, BOUNDS_TYPE.DOUBLE, 0, pieceL, VARIABLE_TYPE.CONTINUOUS, 0);					
				varMap.put("x" + t.id + "_" + k, index);
				index++;
			}
		}
		
		// integer auxiliary variable for each piece of x_i
		for (Target t : targetList) {
			for (int k=1; k<=numPiece-1; k++){
				addAndSetColumn("l" + t.id + "_" + k, BOUNDS_TYPE.DOUBLE, 0, 1,
						VARIABLE_TYPE.INTEGER, 0.0);
				varMap.put("l" + t.id + "_" + k, index);
				index++;
			}
		}
		
		addAndSetColumn("s", BOUNDS_TYPE.FREE, -Configuration.MM, Configuration.MM, VARIABLE_TYPE.CONTINUOUS, 1.0);
		varMap.put("s", index);
		index++;
		
	}

	protected void setMarginalConstraint() {
		
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		
		// Resource constraint		
		addAndSetRow("C[1]", BOUNDS_TYPE.UPPER, -Configuration.MM, numResources);
		rowMap.put("C1", this.getNumRows());
		for (Target t : targetList) {
			for (int k = 1;k <= numPiece; k++){
				ja.add(varMap.get("x" + t.id + "_" + k));
				ar.add(1.0);
			}
		}
		this.setMatRow(this.getNumRows(), ja, ar);
				
	}

	protected void setPiecewiseLinearConstraints() {
		double M = Configuration.MM;
		double pieceL = 1.0/numPiece;

		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
					
		// valid piecewise break up of x_i and _x_i
		for (Target t : targetList) {	
			
			// l_{ik}*(c_k - c_{k-1}) <= x_{ik}
			for (int k=1; k<=numPiece-1; k++){
				addAndSetRow("C[2]"+ t.id + "_" + k, BOUNDS_TYPE.UPPER, -M, 0);
				
				ja.clear();
				ar.clear();
				
				ja.add(varMap.get("x" + t.id + "_" + k));
				ar.add(-1.0);
				
				ja.add(varMap.get("l" + t.id + "_" + k));
				ar.add(pieceL);
				
				this.setMatRow(this.getNumRows(),ja, ar);
				rowMap.put("C2-" + t.id + "_" + k, this.getNumRows());	
			}
			
						
			// x_{i(k+1)} <= l_{ik}
			for (int k=1; k<=numPiece-1; k++) {
				int kk = k+1; 
				
				addAndSetRow("C[3]"+t.id + "_" + k, BOUNDS_TYPE.UPPER, -M, 0);
				ja.clear();
				ar.clear();
				
				ja.add(varMap.get("x" + t.id + "_" + kk));
				ar.add(1.0);
				
				ja.add(varMap.get("l" + t.id + "_" + k));
				ar.add(-1.0);
				
				this.setMatRow(this.getNumRows(),ja, ar);
				rowMap.put("C3-" + t.id + "_" + k, this.getNumRows());	
			}
		}			
	}
	

		
	@Override
	protected void setRowBounds() {
		this.setMarginalConstraint();
		this.setPiecewiseLinearConstraints();		
	}

	@Override
	protected void generateData() {

	}
	
	public double getObjValue()
	{
		return estimatedOpt - getDefenderPayoff();
	}
	
	public double getDefenderCoverage(int targetId){
		double x_t = 0.0;
		for (int k = 1; k <= numPiece; k++){
			double x_tk = this.getColumnPrimal(varMap.get("x" + targetId + "_" + k));
			x_t = x_t + x_tk;
		}
		return x_t;
	}

	public double getDefenderCoverage(Target target) {
		double x_t = 0.0;
		for (int k = 1; k <= numPiece; k++){
			double x_tk = this.getColumnPrimal(varMap.get("x" + target.id + "_" + k));
			x_t = x_t + x_tk;
		}
		return x_t;		
	}

	public Map<Target, Double> getDefenderCoverage() {
		Map<Target, Double> mapDefCoverage = new HashMap<Target, Double>();
		for (Target t : targetList) {			
			mapDefCoverage.put(t, this.getDefenderCoverage(t.id));
		}
		return mapDefCoverage;
	}
	
	public Map<Integer, Double> getDefenderCoverageByTargetIds() {
		Map<Integer, Double> mapDefCoverage = new HashMap<Integer, Double>();
		for (int tId = 0; tId < targetList.size(); tId++) {			
			mapDefCoverage.put(tId, this.getDefenderCoverage(tId));
		}
		return mapDefCoverage;
	}
	
	public Map<Integer, Double> getDefenderPayoffByTargetIds() {
		Map<Integer, Double> mapDefPayoff = new HashMap<Integer, Double>();
		for (Target t : targetList) {	
			double coverage = this.getDefenderCoverage(t.id);
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			mapDefPayoff.put(t.id, coverage*(payoffs.getDefenderReward() - payoffs.getDefenderPenalty()) + payoffs.getDefenderPenalty());
		}
		return mapDefPayoff;
	}
	public double getDefenderPayoff(){
		Map<Integer, Double> mapDefCov = getDefenderCoverageByTargetIds();
		Map<Target, Double> mapAttSUQR = getAttackerQR(adversary);
		
		double defenderEU = 0.0;
		
		for(Target t : targetList){
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			defenderEU += 
				(mapDefCov.get(t.id) * payoffs.getDefenderReward() + (1 - mapDefCov.get(t.id)) * payoffs.getDefenderPenalty()) * mapAttSUQR.get(t);
		}
		return defenderEU;
	}
	public Map<Target, Double> getAttackerQR(SUQRAdversary adversary){
		Map<Target, Double> mapAttackerSUQR = new HashMap<Target, Double>();
		
		double suqr_cum = 0.0;
		double suqr_single = 0.0;
		
		for(Target t : targetList){
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			double defCov = this.getDefenderCoverage(t);
			
			double w1 = adversary.w[SUQRAdversary.DEF_COV] * defCov;
			double w2 = adversary.w[SUQRAdversary.ATT_REWARD] * payoffs.getAdversaryReward();
			double w3 = adversary.w[SUQRAdversary.ATT_PENALTY] * payoffs.getAdversaryPenalty();
			
			suqr_single = Math.exp(w1 + w2 + w3);
			mapAttackerSUQR.put(t, suqr_single);
			suqr_cum += suqr_single;
		}
		
		for(Target t : targetList){
			suqr_single = mapAttackerSUQR.get(t);
			mapAttackerSUQR.remove(t);
			mapAttackerSUQR.put(t, suqr_single / suqr_cum);
		}
		
		return mapAttackerSUQR;
	}
	
	public Map<Target, Double> getAttackerQR(SUQRAdversary adversary, Map<Integer, Double> defenderCoverage){
		Map<Target, Double> mapAttackerSUQR = new HashMap<Target, Double>();
		
		double suqr_cum = 0.0;
		double suqr_single = 0.0;
		
		for(Target t : targetList){
			PayoffStructure payoffs = t.getPayoffStructure(payoffIndex);
			
			double defCov = 0.0;
			
			if(defenderCoverage.containsKey(t.id)){
				defCov = defenderCoverage.get(t.id);
			}
			
			double w1 = adversary.w[SUQRAdversary.DEF_COV] * defCov;
			double w2 = adversary.w[SUQRAdversary.ATT_REWARD] * payoffs.getAdversaryReward();
			double w3 = adversary.w[SUQRAdversary.ATT_PENALTY] * payoffs.getAdversaryPenalty();
			
			suqr_single = Math.exp(w1 + w2 + w3);
			mapAttackerSUQR.put(t, suqr_single);
			suqr_cum += suqr_single;
		}
		
		for(Target t : targetList){
			suqr_single = mapAttackerSUQR.get(t);
			mapAttackerSUQR.remove(t);
			mapAttackerSUQR.put(t, suqr_single / suqr_cum);
		}
		
		return mapAttackerSUQR;
	}
	
}