/*package pathCalculator;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import lpWrapper.Configuration;
import lpWrapper.LPSolverException;

import geoInformation.Cell;
import geoInformation.Grid;
import utilities.JointPatrolSlave;

public class MinNetworkFlow{

	public JointPatrolSlave ps;
	public int cellsToCover;
	private Map<Integer, Cell> targetToCellsMap;

	public MinNetworkFlow(Grid grid, double[] regret, int cellsToCover, boolean returnToSameBase) throws IOException{		
//		Configuration.loadLibrariesCplex(File.separator + "Users"
//										+ File.separator + "dellefav"
//										+ File.separator + "Desktop"
//										+ File.separator + "workspace"
//										+ File.separator + "WildLifeAndUAVs"
//										+ File.separator + "cplexLIB"
//										+ File.separator + "CplexConfig");
				
		this.cellsToCover = cellsToCover;
		this.targetToCellsMap = grid.targetToCellMap;
		
		ps = new JointPatrolSlave(grid, regret, 1, this.cellsToCover, 1, returnToSameBase);		
	}

	public void solve() throws LPSolverException {
//		ps.loadProblem();
		ps.solve();
	}

	public String getBestPath() {
		return ps.getPath().toString();
	}

	public Cell[] getPathCells() {
		int[] cellsIds = ps.getTargets();
		Cell[] output = new Cell[cellsIds.length];
		
		for (int i = 0; i < cellsIds.length; i++) {
			output[i] = targetToCellsMap.get(cellsIds[i]);
		}
		
		return output;
	}
	public void end()
	{
		ps.end();
	}
	
}
*/