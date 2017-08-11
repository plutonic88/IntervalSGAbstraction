package src.tambe.usc.DOBSS.gameRepresentation;

/**
 * User: ckiekint
 * Date: Aug 8, 2008
 * Time: 5:58:59 PM
 */

public class SSGSolution {
	private int nTargets;
	  private int nDefenseTypes;
	  private double[][] probabilities;

	  private int attackedTarget;
	  private double defenderPayoff = Double.NEGATIVE_INFINITY;
	  private double attackerPayoff = Double.NEGATIVE_INFINITY;

	  public SSGSolution(int nTargets, int nDefenseTypes) {
	    this.nTargets = nTargets;
	    this.nDefenseTypes = nDefenseTypes;
	    probabilities = new double[nTargets][nDefenseTypes];
	  }

	  public int getAttackedTarget() {
	    return attackedTarget;
	  }

	  public void setAttackedTarget(int attackedTarget) {
	    this.attackedTarget = attackedTarget;
	  }

	  public double getProb(int target, int defenseType) {
	    return probabilities[target][defenseType];
	  }

	  public void setProb(int target, int defenseType, double value) {
		  if(value<0)
		  {
			//  System.out.println("ff");
		  }
	    probabilities[target][defenseType] = value;
	  }

	  public double getDefenderPayoff() {
	    return defenderPayoff;
	  }

	  public void setDefenderPayoff(double defenderPayoff) {
	    this.defenderPayoff = defenderPayoff;
	  }

	  public double getAttackerPayoff() {
	    return attackerPayoff;
	  }

	  public void setAttackerPayoff(double attackerPayoff) {
	    this.attackerPayoff = attackerPayoff;
	  }

	  public double[] getDefenderProbs(int defenseType) {
	   double[] tmp = new double[nTargets];
	   for (int target = 0; target < nTargets; target++) {
	      tmp[target] = probabilities[target][defenseType];
	    }
	    return tmp;
	  }

	  public double[] getTargetProbs(int target) {
	    return probabilities[target];
	  }

	  public double getTotalTargetProb(int target) {
	    double tmp = 0d;
	    for (int defender = 0; defender < nDefenseTypes; defender++) {
	      tmp += probabilities[target][defender];
	    }
	    return tmp;
	  }

	  public double getTotalDefenderProb(int defenseType) {
	    double tmp = 0d;
	    for (int target = 0; target < nTargets; target++) {
	      tmp += probabilities[target][defenseType];
	    }
	    return tmp;
	  }

	  public void computeAttackerStrategy(StructuredSecurityGame ssg) {
	    if (ssg.getNTargets() != nTargets || ssg.getNDefenseTypes() != nDefenseTypes) {
	      System.err.println("SSGsolution::Computing expected payoffs for invalid game!");
	      return;
	    }

	    double TOLERANCE = 0.001d;

	    double bestAttackerPayoff = Double.NEGATIVE_INFINITY;
	    double bestDefenderPayoff = Double.NEGATIVE_INFINITY;
	    int bestTarget = -1;
	    for (int t = 0; t < nTargets; t++) {
	      double prob = 0d;
	      for (int dt = 0; dt < nDefenseTypes; dt++) {
	        prob += probabilities[t][dt];
	      }

	      SSGPayoffs payoffs = ssg.getPayoffs(t);
	      double tmpAttackerPayoff = prob * payoffs.getAttackerCoveredPayoff() + (1-prob) * payoffs.getAttackerUncoveredPayoff();
	      double tmpDefenderPayoff = prob * payoffs.getDefenderCoveredPayoff() + (1-prob) * payoffs.getDefenderUncoveredPayoff();

	      if ( (tmpAttackerPayoff > bestAttackerPayoff + TOLERANCE) ||
	           (tmpAttackerPayoff > bestAttackerPayoff - TOLERANCE && tmpDefenderPayoff > bestDefenderPayoff) ) {
	        bestAttackerPayoff = tmpAttackerPayoff;
	        bestDefenderPayoff = tmpDefenderPayoff;
	        bestTarget = t;
	      }
	    }
	    attackedTarget = bestTarget;
	  }

	  public void computeExpectedPayoffs(StructuredSecurityGame ssg) {
	    if (ssg.getNTargets() != nTargets || ssg.getNDefenseTypes() != nDefenseTypes) {
	      System.err.println("SSGsolution::Computing expected payoffs for invalid game!");
	      return;
	    }

	    computeAttackerStrategy(ssg);

	    double prob = 0d;
	    for (int dt = 0; dt < nDefenseTypes; dt++) {
	      prob += probabilities[attackedTarget][dt];
	    }

	    SSGPayoffs payoffs = ssg.getPayoffs(attackedTarget);
	    defenderPayoff = prob * payoffs.getDefenderCoveredPayoff() + (1-prob) * payoffs.getDefenderUncoveredPayoff();
	    attackerPayoff = prob * payoffs.getAttackerCoveredPayoff() + (1-prob) * payoffs.getAttackerUncoveredPayoff();
	  }

	  public String expectedPayoffString(StructuredSecurityGame ssg) {
	    if (ssg.getNTargets() != nTargets || ssg.getNDefenseTypes() != nDefenseTypes) {
	      System.err.println("SSGsolution::Computing expected payoffs for invalid game!");
	      return null;
	    }

	    StringBuilder sb1 = new StringBuilder();
	    StringBuilder sb2 = new StringBuilder();
	    StringBuilder sb3 = new StringBuilder();

	    sb1.append("Attacker payoffs: ");
	    sb2.append("Defender payoffs: ");
	    sb3.append("Coverage probability: ");

	    for (int t = 0; t < nTargets; t++) {
	      double prob = 0d;
	      for (int dt = 0; dt < nDefenseTypes; dt++) {
	        prob += probabilities[t][dt];
	      }

	      sb3.append(" [").append(t).append("] ").append(prob);

	      SSGPayoffs payoffs = ssg.getPayoffs(t);
	      double tmpAttackerPayoff = prob * payoffs.getAttackerCoveredPayoff() + (1-prob) * payoffs.getAttackerUncoveredPayoff();
	      double tmpDefenderPayoff = prob * payoffs.getDefenderCoveredPayoff() + (1-prob) * payoffs.getDefenderUncoveredPayoff();

	      sb1.append(" [").append(t).append("] ").append(tmpAttackerPayoff);
	      sb2.append(" [").append(t).append("] ").append(tmpDefenderPayoff);
	    }

	    sb3.append("\n").append(sb1.toString()).append("\n");
	    sb3.append("\n").append(sb2.toString()).append("\n");
	    return sb3.toString();
	  }

	  public String toString() {
	    StringBuilder sb = new StringBuilder();
	    sb.append("Attacked target: ").append(attackedTarget+1).append("\n");
	    sb.append("Defender payoff: ").append(defenderPayoff).append("\n");
	    sb.append("Attacker payoff: ").append(attackerPayoff).append("\n");
	    
	    /*for (int target = 0; target < nTargets; target++) 
	    {
	      sb.append("[ ").append(target+1).append("] ");
	      for (int defense = 0; defense < nDefenseTypes; defense++) {
	        sb.append(probabilities[target][defense]).append(" ");
	      }
	      sb.append("   tot: ").append(getTotalTargetProb(target));
	      sb.append("\n");
	    }*/
	    return sb.toString();
	  }


}
