package geoInformation;

public class Area {

	public double length;
	public double width;
	public String name;
	
	public int minDensity;
	public int maxDensity;
	
	public Area(String name, double length, double width, int minDensity, int maxDensity) {
		this.name = name;
		
		this.length = length;
		this.width = width;
		
		this.minDensity = minDensity;
		this.maxDensity = maxDensity;
	}

	@Override
	public boolean equals(Object obj) {
		assert(obj instanceof Area);
		
		return this.length == ((Area) obj).length
				&& this.width == ((Area) obj).width
				&& this.minDensity == ((Area) obj).minDensity
				&& this.maxDensity == ((Area) obj).maxDensity;
	}
	
	@Override
	public String toString() {
		return name + ": (" + length + ", " + width +  ", " + minDensity + ", " + maxDensity + ")";
	}
	
	@Override
	public int hashCode() {
		return name.hashCode() + (int) length + (int) width + (int) minDensity + (int) maxDensity;
	}
}
