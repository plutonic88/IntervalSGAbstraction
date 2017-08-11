package models;

import java.util.ArrayList;
import java.util.List;


/**
 * This class shall track the payoff structure and probability of coverage for a single target.
 * Note that a single Target instance belongs to exactly one GameRound object.
 * @author Ben Ford
 */
public class Target {

	public int id;
	private String targetID;
	private List<PayoffStructure> payoffList; // list of payoff samples
	private PayoffStructure truePayoffStructure;
	private Float coverageProbability;//Valid values: [0,1]
//	private Boolean isCovered;
	AttBoundStructure attBoundStructure; // lower bound and upper bound
	/**
	 * @param targetID The unique ID of the target.
	 */
	public Target(int id, String targetID)
	{
		this.targetID = targetID;
		payoffList = null;
		coverageProbability = null;
//		isCovered = null;
		this.id = id;
//		id = counter;
//		counter++;
	}
	
	/**
	 * @param payoffStructure The associated PayoffStructure object for this target.
	 */
	public Target(int id, List<PayoffStructure> payoffList, AttBoundStructure attBoundStructure)
	{
		this.payoffList = payoffList;
		this.attBoundStructure = attBoundStructure;
		coverageProbability = null;
		this.id = id;
	}
	
	/**
	 * @param payoffStructure The associated PayoffStructure object for this target.
	 * @param coverageProbability The coverage probability for this target.
	 * @throws ValueOutOfBoundsException
	 */
	public Target(int id, ArrayList<PayoffStructure> payoffStructure, float coverageProbability) 
	{
		
		this.payoffList = payoffList;
		this.coverageProbability = coverageProbability;
		
		this.id = id;
//		id = counter;
//		counter++;
	}
	
	public List<PayoffStructure> getPayoffList() {
		return payoffList;
	}
	public PayoffStructure getPayoffStructure(int payoffIndex)
	{
		return payoffList.get(payoffIndex);
	}
	public AttBoundStructure getAttBoundStructure()
	{
		return attBoundStructure;
	}
	public float getCoverageProbability() {
		return coverageProbability;
	}
	
	public void setPayoffList(List<PayoffStructure> payoffList) {
		this.payoffList = payoffList;
	}
	
	public void setCoverageProbability(Float coverageProbability) {
		this.coverageProbability = coverageProbability;
	}
	
	public String getTargetID() {
		return targetID;
	}

	public void setTargetID(String targetID) {
		this.targetID = targetID;
	}
	
	public int getNumPayoff()
	{
		return payoffList.size();
	}
	
	public void setAttBoundStructure(AttBoundStructure newAttBound)
	{
		attBoundStructure.setAttPenaltyLB(newAttBound.getAttPenaltyLB());
		attBoundStructure.setAttPenaltyUB(newAttBound.getAttPenaltyUB());
		attBoundStructure.setAttRewardLB(newAttBound.getAttRewardLB());
		attBoundStructure.setAttRewardUB(newAttBound.getAttRewardUB());
	}
	public void setTruePayoffStructure(PayoffStructure payoffs)
	{
		truePayoffStructure = payoffs;
	}
	
	public void setTrueReward(double trueReward)
	{
		truePayoffStructure.setAdversaryReward(trueReward);
	}
	public PayoffStructure getTruePayoffStructure()
	
	{
		return truePayoffStructure;
	}
	/**
	 * Format: Payoff Structure \n Coverage Probability
	 */
}
