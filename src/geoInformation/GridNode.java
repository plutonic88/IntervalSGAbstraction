package geoInformation;

public class GridNode {

	public Cell cell;
	
	public GridNode(Cell c) {
		cell = c;
	}
	
	@Override
	public boolean equals(Object obj) {
		assert (obj instanceof GridNode);
			
		GridNode gn = (GridNode) obj;
		
		return this.cell.equals(gn.cell);
	}

	@Override
	public String toString() {
		return cell.toString();
	}

	public boolean isStartingPoint() {
		return cell.start;
	}
}
