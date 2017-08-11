package models;

public class SUQRAdversary{
	public static final int DEF_COV = 0;
	public static final int ATT_REWARD = 1;
	public static final int ATT_PENALTY = 2;
	
//	public static int counter = 1;
	
	public double[] w;
	public int id;
	public double prob;
	
	public SUQRAdversary(int id, double w1, double w2, double w3, double prob) {
		
		this.id = id;
		w = new double[3];
		
		w[DEF_COV] = w1;
		w[ATT_REWARD] = w2;
		w[ATT_PENALTY] = w3;
		this.prob = prob;
	}
}
