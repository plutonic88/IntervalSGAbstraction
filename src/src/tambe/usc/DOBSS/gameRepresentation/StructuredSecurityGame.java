package src.tambe.usc.DOBSS.gameRepresentation;

import java.util.Random;

/**
 * User: ckiekint
 * Date: Aug 8, 2008
 * Time: 2:14:08 PM
 */
public class StructuredSecurityGame {

	private int nTargets;
	private int nDefenseTypes;
	private double[] nDefenders;
	private boolean[][] defenseCapabilities;
	private SSGPayoffs[] payoffs;

	public StructuredSecurityGame(int nTargets, int nDefenseTypes) {
		this.nTargets = nTargets;
		this.nDefenseTypes = nDefenseTypes;
		nDefenders = new double[nDefenseTypes];
		defenseCapabilities = new boolean[nTargets][nDefenseTypes];
		payoffs = new SSGPayoffs[nTargets];
	}

	public int getNTargets() {
		return nTargets;
	}

	public int getNDefenseTypes() {
		return nDefenseTypes;
	}

	public void setNumDefenders(int type, int number) {
		nDefenders[type] = number;
	}

	public double getNumDefenders(int type) {
		return nDefenders[type];
	}

	public double[] getNDefenders() {
		return nDefenders;
	}

	public void setNDefenders(double[] nDefenders) {
		this.nDefenders = nDefenders;
	}

	public boolean[][] getDefenseCapabilities() {
		return defenseCapabilities;
	}

	public boolean getDefenseCapability(int target, int type) {
		return defenseCapabilities[target][type];
	}

	public void setDefenseCapabilities(boolean[][] defenseCapabilities) {
		this.defenseCapabilities = defenseCapabilities;
	}

	public void setDefenseCapability(int target, int type, boolean value) {
		defenseCapabilities[target][type] = value;
	}

	public SSGPayoffs[] getPayoffs() {
		return payoffs;
	}

	public SSGPayoffs getPayoffs(int target) {
		return payoffs[target];
	}

	public void setPayoffs(SSGPayoffs[] payoffs) {
		this.payoffs = payoffs;
	}

	public void setPayoffs(int target, SSGPayoffs payoffs) {
		this.payoffs[target] = payoffs;
	}
	/*  Blocked This Part of the Code */  
	/*public static StructuredSecurityGame genDOBSSGame(SecurityGame sg)
  {

	  StructuredSecurityGame ssg = new StructuredSecurityGame(8, 1);
	  ssg.setNDefenders(new int[] {1});

	    for (int target = 1; target <= 8; target++) {
		      SSGPayoffs tmp = new SSGPayoffs(false);
		      double r1 = sg.getPayoff(new int[] {1,1,target,target}, 1);
		      double r2 = sg.getPayoff(new int[] {1,1,target, target+1}, 1);
		      double r3 = sg.getPayoff(new int[] {1,1,target,target}, 0);
		      double r4 = sg.getPayoff(new int[] {1,1,target, target+1}, 0);
		      tmp.setAttackerCoveredPayoff(r1);
		      tmp.setAttackerUncoveredPayoff(r2);
		      tmp.setDefenderCoveredPayoff(r3);
		      tmp.setDefenderUncoveredPayoff(r4);
		      ssg.setPayoffs(target-1, tmp);
		    }

	    for (int target = 0; target < 8; target++) {
		      ssg.setDefenseCapability(target, 0, true);
		    }

	    return ssg;
  }
	 */
	public static StructuredSecurityGame genRandomGame(int nTargets, int nDefenseTypes, double[] nDefenders) 
	{
		StructuredSecurityGame ssg = new StructuredSecurityGame(nTargets, nDefenseTypes);
		ssg.setNDefenders(nDefenders);

		for (int target = 0; target < nTargets; target++) 
		{
			SSGPayoffs tmp = new SSGPayoffs(false);

			//-------------------------------------------
			//Attacker and Defender's covered payoff is set to zero here
			//-------------------------------------------
			//double r1 = Math.random() * -100;
			double r1=0;
			double r3=0;

			double r2 = Math.random() * -10;
			//double r3 = Math.random() * 100;
			double r4 = Math.random() * 10;

			// System.out.print("R2");
			//System.out.println( r2);

			//---------------------------------------------
			tmp.setAttackerCoveredPayoff(r1);
			tmp.setAttackerUncoveredPayoff(r2);
			tmp.setDefenderCoveredPayoff(r3);
			tmp.setDefenderUncoveredPayoff(r4);
			ssg.setPayoffs(target, tmp);
		}

		Random r = new Random();
		// first, assign at least 1 resource to each schedule (so there are no uncoverable schedules)
		for (int target = 0; target < nTargets; target++) {
			ssg.setDefenseCapability(target, r.nextInt(nDefenseTypes), true);
		}

		// next, randomly add additional resource capabilities up to the density requested
		int cnt = nTargets;
		int max = (int)Math.round((double)(nTargets * nDefenseTypes) * 0.3d);
		while (cnt < max) {
			int tmp1 = r.nextInt(nTargets);
			int tmp2 = r.nextInt(nDefenseTypes);
			if (ssg.getDefenseCapability(tmp1, tmp2)) continue;
			ssg.setDefenseCapability(tmp1, tmp2, true);
			cnt++;
		}

		/*
	    for (int target = 0; target < nTargets; target++) {
	      //ssg.setDefenseCapability(target,0,true);
	      // Set coverage of flights with a random distribution

	      for (int defense = 0; defense < nDefenseTypes; defense++) {
	        ssg.setDefenseCapability(target, defense, Math.random() < 0.3d); 
	      }
	/*
	    	// Set coverage of flights with a random distribution + restriction on the number of offices that
	    	// can cover a single flight
	    	int nCover;
	        for (int defense = 0; defense < nDefenseTypes; defense++) {
	        	nCover = 0;
	        	// Check how many offices can already cover a given target (flight-pair)
	        	for (int i = 0; i < nDefenseTypes; i++) {
	        		if(ssg.getDefenseCapability(target,i))
	        				nCover++;
	        	}
	        	// If less than two offices can cover a target, set the current office to be able to cover it
	        	// with a probability of 40%
	            ssg.setDefenseCapability(target, defense, (Math.random() < 0.4d) & (nCover < 3)); 
	        }

	    }
		 */
		return ssg; 
	}
}

