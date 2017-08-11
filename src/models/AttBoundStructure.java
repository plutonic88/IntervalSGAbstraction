package models;

public class AttBoundStructure {
	double attRewardLB;
	double attPenaltyLB;
	double attRewardUB;
	double attPenaltyUB;
	public AttBoundStructure(double attRewardLB, double attPenaltyLB, double attRewardUB, double attPenaltyUB)
	{
		this.attRewardLB = attRewardLB;
		this.attRewardUB = attRewardUB;
		this.attPenaltyLB = attPenaltyLB;
		this.attPenaltyUB = attPenaltyUB;
	}
	
	public AttBoundStructure(AttBoundStructure bound)
	{
		this.attPenaltyLB = bound.getAttPenaltyLB();
		this.attPenaltyUB = bound.getAttPenaltyUB();
		this.attRewardLB = bound.getAttRewardLB();
		this.attRewardUB = bound.getAttRewardUB();
	}
	public double getAttRewardLB()
	{
		return attRewardLB;
	}
	public double getAttRewardUB()
	{
		return attRewardUB;
	}
	public double getAttPenaltyLB()
	{
		return attPenaltyLB;
	}
	public double getAttPenaltyUB()
	{
		return attPenaltyUB;
	}
	
	public void setAttRewardLB(double rewardLB)
	{
		this.attRewardLB = rewardLB;
	}
	public void setAttRewardUB(double rewardUB)
	{
		this.attRewardUB = rewardUB;
	}
	public void setAttPenaltyLB(double penaltyLB)
	{
		this.attPenaltyLB = penaltyLB;
	}
	public void setAttPenaltyUB(double penaltyUB)
	{
		this.attPenaltyUB = penaltyUB;
	}
}
