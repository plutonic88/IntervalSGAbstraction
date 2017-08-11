package pathCalculator;

import java.util.Collection;
import java.util.TreeSet;

public class Path {
	private Collection<Step> stepList;

	public Path(){
		stepList = new TreeSet<Step>();
	}
	
	public void add(Step step) {
		stepList.add(step);
	}

	@Override
	public boolean equals(Object obj) {
		assert (obj instanceof Step);
		
		Path p = (Path) obj;
			
		return stepList.equals(p.stepList);
	}
	
	@Override
	public String toString() {
		return stepList.toString();
	}
}
