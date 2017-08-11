package pathCalculator;

public class Step implements Comparable<Step> {
	public int sID;
	public int tID;
	public int tStep;

	public Step(int source, int target, int timeStep) {
		this.sID = source;
		this.tID = target;
		this.tStep = timeStep;
	}
	
	@Override
	public boolean equals(Object obj) {
		assert (obj instanceof Step);
		
		Step s = (Step) obj; 
		
		return this.sID == s.sID && this.tID == s.tID && this.tStep == s.tStep;
	}
	
	@Override
	public String toString() {
		return "<" + tStep + ": "+ sID + " -> " + tID + ">";
	}

	@Override
	public int compareTo(Step o) {
		return (tStep > o.tStep)? 1 : (tStep < o.tStep)? -1 : 0;
	}
	
}
