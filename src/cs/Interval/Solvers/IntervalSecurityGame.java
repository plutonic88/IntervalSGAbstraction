package cs.Interval.Solvers;



/**
 * An SSG with interval payoffs for the attacker rewards
 */
public class IntervalSecurityGame {

  private int nTargets;
  private double nDefenders;
  private IntervalPayoffs[] payoffs;

  public IntervalSecurityGame(int nTargets) {
    this.nTargets = nTargets;
    nDefenders = (int)Math.ceil(nTargets * 0.2);
    payoffs = new IntervalPayoffs[nTargets];
  }

  public int getNTargets() {
    return nTargets;
  }

  public void setNumDefenders(int number) {
    nDefenders = number;
  }
  public void setNumDefenders(double number) {
	    nDefenders = number;
	  }

  public double getNumDefenders() {
    return nDefenders;
  }

  public IntervalPayoffs[] getPayoffs() {
    return payoffs;
  }

  public IntervalPayoffs getPayoffs(int target) {
    return payoffs[target];
  }

  public void setPayoffs(IntervalPayoffs[] payoffs) {
    this.payoffs = payoffs;
  }

  public void setPayoffs(int target, IntervalPayoffs payoffs) {
    this.payoffs[target] = payoffs;
  }

  public static IntervalSecurityGame genRandomGame(int nTargets, int nDefenders) {
	    IntervalSecurityGame ssg = new IntervalSecurityGame(nTargets);
	    ssg.setNumDefenders(nDefenders);

	    for (int target = 0; target < nTargets; target++) 
	    {
	      IntervalPayoffs tmp = new IntervalPayoffs(true);
	      ssg.setPayoffs(target, tmp);
	    }

	    return ssg;
  }
  
  public static IntervalSecurityGame generateAbstractIntervalSecurityGame(int[][][] abstractgame, double ndefenders)
  {
	  int numberoftargets = abstractgame.length;
	  IntervalSecurityGame issg = new IntervalSecurityGame(numberoftargets);
	  issg.setNumDefenders(ndefenders);
	  for (int target = 0; target < numberoftargets; target++) 
	    {
	      IntervalPayoffs tmp = new IntervalPayoffs();
	      tmp.setAttackerCoveredPayoffMin(abstractgame[target][3][0]);
	      tmp.setAttackerCoveredPayoffMax(abstractgame[target][3][1]);
	      
	      tmp.setAttackerUncoveredPayoffMax(abstractgame[target][2][1]);
	      tmp.setAttackerUncoveredPayoffMin(abstractgame[target][2][0]);
	      
	      tmp.setDefenderCoveredPayoff(abstractgame[target][0][0]);
	      tmp.setDefenderUncoveredPayoff(abstractgame[target][1][0]);
	      issg.setPayoffs(target, tmp);
	    }
	  return issg;
	  
	  
	  
	  
  }
}