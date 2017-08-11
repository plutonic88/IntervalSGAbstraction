package src.tambe.usc.DOBSS.gameRepresentation;

import java.util.Arrays;

/**
 * User: ckiekint
 * Date: Aug 8, 2008
 * Time: 2:29:48 PM
 */
public class SSGPayoffs {

  private double[] payoffs;

  public SSGPayoffs(boolean random) {
    payoffs = new double[4];
    if (random) {
      for (int i = 0; i < 4; i++) {
        payoffs[i] = Math.random() * 100;
      }
    }
  }

  public SSGPayoffs(double[] payoffs) {
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

  public double getAttackerCoveredPayoff() {
    return payoffs[2];
  }

  public double getAttackerUncoveredPayoff() {
    return payoffs[3];
  }

  public void setDefenderCoveredPayoff(double payoff) {
    payoffs[0] = payoff;
  }

  public void setDefenderUncoveredPayoff(double payoff) {
    payoffs[1] = payoff;
  }

  public void setAttackerCoveredPayoff(double payoff) {
    payoffs[2] = payoff;
  }

  public void setAttackerUncoveredPayoff(double payoff) {
    payoffs[3] = payoff;
  }

  public String toString() {
    return Arrays.toString(payoffs);
  }
}

