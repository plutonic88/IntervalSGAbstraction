package geoInformation;

public class Cell {

	public int id;
	public int x;
	public int y;
	
	// animal density in the specific region
	public double aniDen_lb;
	public double aniDen_ub;

	// true if cell is starting point
	public boolean start;

	public Cell(int cell_id, int coord_i, int coord_j) {
		id = cell_id;
		
		x = coord_i;
		y = coord_j;
	}

	@Override
	public boolean equals(Object obj) {
		assert(obj instanceof Cell);
		
		Cell c1 = (Cell) obj;
		
		return id == c1.id && x == c1.x && y == c1.y;
	}
	
	@Override
	public String toString() {
		return "Cell " + id + ": (" + x + ", " + y + ")";
	}
}
