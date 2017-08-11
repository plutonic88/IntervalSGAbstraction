package cs.Interval.Solvers;

import java.util.Arrays;





/**
 * User: ckiekint
 * Date: Aug 8, 2008
 * Time: 2:29:48 PM
 */
public class IntervalPayoffs {

  private double[] payoffs;
  
  public IntervalPayoffs()
  {
	  payoffs = new double[6]; 
  }

  public IntervalPayoffs(boolean random) {
    payoffs = new double[5];
    if (random) {
      payoffs[0] = 0;
      payoffs[1] = Math.random() * -100;
      payoffs[2] = 0;
      payoffs[4] = Math.random() * 100;
      double intervalSize = Math.random() * 20;
      payoffs[3] = payoffs[4] + intervalSize;
    }
  }
  /*
  public IntervalPayoffs(DistributionalSecurityGame dsg,int target,double parameter) {
	  
	  
		  // Find out mean and stdDev for a target //
		  DistributionSampler attackerUncovered =dsg.getAttackerUncoveredDistribution(target);
		  DistributionSampler attackerCovered=dsg.getAttackerCoveredDistribution(target);
		  //
		 print("in Gaussian");
		 print(dsg.toString());
		 
		 print("mean" + attackerUncovered.getMean());
		 print("std Dev "+attackerUncovered.getStdDev() );
		  payoffs = new double[6];
		    
		      payoffs[0] = dsg.getDefenderCoveredPayoff(target);
		      payoffs[1] = dsg.getDefenderUncoveredPayoff(target);//Math.random() * -100;
		      
		      //--attacker Coverd Payoff
		      
		      //--min
		      payoffs[2] = attackerCovered.getMean() - parameter * attackerCovered.getStdDev();
		      //--max
		      payoffs[3] = attackerCovered.getMean() + parameter * attackerCovered.getStdDev();
		      
		      
		    //--attacker UnCoverd Payoff
		      
		      //min
		      payoffs[4] = attackerUncovered.getMean() - parameter * attackerUncovered.getStdDev();
		      //max
		      payoffs[5] = attackerUncovered.getMean() + parameter * attackerUncovered.getStdDev();
		      //System.out.println("Std Dev : " +attackerUncovered.getStdDev());
	 
	 }*/

  public IntervalPayoffs(double[] payoffs) {
    this.payoffs = payoffs.clone();
  }

  public double[] getPayoffs() {
    return payoffs;
  }

  public void setPayoffs(double[] payoffs) {
    this.payoffs = payoffs.clone();
  }

  public double getDefenderCoveredPayoff() {
    return payoffs[0];
  }

  public double getDefenderUncoveredPayoff() {
    return payoffs[1];
  }

  public double getAttackerCoveredPayoffMin() {
    return payoffs[2];
  }
  
  public double getAttackerCoveredPayoffMax() {
	    return payoffs[3];
	  }
  

  public double getAttackerUncoveredPayoffMin() {
    return payoffs[4];
  }

  public double getAttackerUncoveredPayoffMax() {
    return payoffs[5];
  }

  public void setDefenderCoveredPayoff(double payoff) {
    payoffs[0] = payoff;
  }

  public void setDefenderUncoveredPayoff(double payoff) {
    payoffs[1] = payoff;
  }

  public void setAttackerCoveredPayoffMin(double payoff) {
    payoffs[2] = payoff;
  }
  
  public void setAttackerCoveredPayoffMax(double payoff) {
	payoffs[3] = payoff;
  }

  public void setAttackerUncoveredPayoffMin(double payoff) {
    payoffs[4] = payoff;
  }

  public void setAttackerUncoveredPayoffMax(double payoff) {
    payoffs[5] = payoff;
  }

  public String toString() {
    return Arrays.toString(payoffs);
  }
  
  public void print(String str)
  {
	  System.out.println(str);
  }


}