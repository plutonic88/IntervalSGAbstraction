/*package geoInformation;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

//import payoffs.AttPayoffsGenerator;
//import payoffs.DefPayoffsGenerator;
//import payoffs.PayoffsGenerator;
//
//import density.DensityGenerator;


public class Grid {
	private Area area;
	
	public double cell_l;
	public double cell_w;
	public int resolution;

	// animal density
	public double[][] aniDen_lb;
	public double[][] aniDen_ub;
	
	// Attacker payoffs
	public double[][] attR_lb;
	public double[][] attR_ub;
	public double[][] attP_lb;
	public double[][] attP_ub;
	
	// Defender payoffs
	public double[][] defR_lb;
	public double[][] defR_ub;
	public double[][] defP_lb;
	public double[][] defP_ub;	

	// graphical tools
	public Collection<Cell> cellList;
	public Collection<GridNode> nodeList;
	public GridNode[][] grid;
	public SimpleWeightedGraph<GridNode, DefaultWeightedEdge> graph;
	public Map<Integer,Cell> targetToCellMap;
	
	//tools
//	DensityGenerator densityGenerator;
//	PayoffsGenerator defPayoffsGenerator, attPayoffsGenerator;

	public boolean graphGenerated;
		
	public Grid(Area area, int grid_resolution){
		this.area = area;
		
		resolution = grid_resolution;
		
		cell_l = area.length / grid_resolution; 
		cell_w = area.width / grid_resolution; 	
		
		// animal density
		aniDen_lb = new double[resolution][resolution];
		aniDen_ub = new double[resolution][resolution];
		
		// Defender payoffs
		defR_lb = new double[resolution][resolution];
		defR_ub = new double[resolution][resolution];
		defP_lb = new double[resolution][resolution];
		defP_ub = new double[resolution][resolution];	

		// Attacker payoffs
		attR_lb = new double[resolution][resolution];
		attR_ub = new double[resolution][resolution];
		attP_lb = new double[resolution][resolution];
		attP_ub = new double[resolution][resolution];
				
		// instantiate graphical tools
		cellList = new ArrayList<Cell>();
		targetToCellMap = new HashMap<Integer, Cell>();
		nodeList = new ArrayList<GridNode>();
		grid = new GridNode[resolution][resolution];
		graph = new SimpleWeightedGraph<GridNode, DefaultWeightedEdge>(GridEdge.class);
	}

//	public void generateAnimalDensity() {		
//		for(int i = 0; i < resolution; i++){
//			for(int j = 0; j < resolution;j++){
//				aniDen_lb[i][j] = densityGenerator.generateLB(area);
//				aniDen_ub[i][j] = densityGenerator.generateUB(area);
//			}
//		}	
//	}
//
//	public void generateDefenderPayoffs() {		
//		for(int i = 0; i < resolution; i++){
//			for(int j = 0; j < resolution; j++){
//				double[] defPayoffs = defPayoffsGenerator.generate(aniDen_lb[i][j], aniDen_ub[i][j]);
//				
//				defR_lb[i][j] = defPayoffs[0];
//				defR_ub[i][j] = defPayoffs[1];
//				defP_lb[i][j] = defPayoffs[2];
//				defP_ub[i][j] = defPayoffs[3];
//			}
//		}
//	}
//
//	public void generateAttackerPayoffs() {
//		for(int i = 0; i < resolution; i++){
//			for(int j = 0; j < resolution; j++){
//				double[] attPayoffs = attPayoffsGenerator.generate(aniDen_lb[i][j], aniDen_ub[i][j]);
//				
//				attR_lb[i][j] = attPayoffs[0];
//				attR_ub[i][j] = attPayoffs[1];
//				attP_lb[i][j] = attPayoffs[2]; 
//				attP_ub[i][j] = attPayoffs[3]; 
//			}
//		}
//	}
//	
//	public void setAnimalGenerator(DensityGenerator generator) {
//		densityGenerator = generator;
//	}

//	public void setPayoffsGenerators(DefPayoffsGenerator defGenerator, AttPayoffsGenerator attGenerator) {
//		attPayoffsGenerator = attGenerator;
//		defPayoffsGenerator = defGenerator;
//	}

	public void instantiateCellList(){
		int cell_id = 1;
		
		for (int i = 0; i < resolution; i++){
			for(int j = 0; j < resolution; j++){
				Cell c = new Cell(cell_id,i,j);
				c.aniDen_lb = aniDen_lb[i][j];
				c.aniDen_ub = aniDen_ub[i][j];
				
				cellList.add(c);
				targetToCellMap.put(cell_id, c);
				
				GridNode gn = new GridNode(c);
				
				nodeList.add(gn);
				grid[i][j] = gn;
				graph.addVertex(gn);
				
				cell_id++;
			}
		}
		
	}
	
	public void readPayoffs(String fileName, char type) {
		try {
			Scanner sc = new Scanner(new FileReader(fileName));
			
			// read payoffs			
			while(sc.hasNext()){
				String[] array = sc.next().split(",");
				
				// we read payoffs for the defender
				int i = Integer.parseInt(array[0]);
				int j = Integer.parseInt(array[1]);

				if(type == 'D'){
					defR_lb[i][j] = Double.parseDouble(array[2]);
					defR_ub[i][j] = Double.parseDouble(array[3]);
					defP_lb[i][j] = Double.parseDouble(array[4]);
					defP_ub[i][j] = Double.parseDouble(array[5]);					
				}

				if(type == 'A'){
					attR_lb[i][j] = Double.parseDouble(array[2]);
					attR_ub[i][j] = Double.parseDouble(array[3]);
					attP_lb[i][j] = Double.parseDouble(array[4]); 
					attP_ub[i][j] = Double.parseDouble(array[5]); 						
				}
			}
			
			sc.close();

		} catch (IOException e) {
			System.err.println("Payoff file not found");
			e.printStackTrace();		
		}
	}

	public void readDensity(String aniDenLoad) {
		try {
			Scanner scAni = new Scanner(new FileReader(aniDenLoad));
			
			// read payoffs			
			while(scAni.hasNext()){
				String[] aniArray = scAni.next().split(",");
				
				// we read payoffs for the defender
				int i = Integer.parseInt(aniArray[0]);
				int j = Integer.parseInt(aniArray[1]);

				aniDen_lb[i][j] = Double.parseDouble(aniArray[2]);
				aniDen_ub[i][j] = Double.parseDouble(aniArray[3]);				
			}
			
			scAni.close();

		} catch (IOException e) {
			System.err.println("Animal Density file not found");
			e.printStackTrace();		
		}		
	}

	public void printAnimalDensity(String aniDenSave) {
		try {
			PrintWriter pwAni = new PrintWriter(new File(aniDenSave));
			
			// write payoffs
			for(int i = 0; i < resolution; i++){
				for(int j = 0; j < resolution; j++){
					pwAni.println(i + "," + j +  "," + aniDen_lb[i][j] + "," + aniDen_ub[i][j]);
				}
			}
			
			pwAni.close();

		} catch (IOException e) {
			System.err.println("Animal Density file not found");
			e.printStackTrace();		
		}		
	}
	
	public void printPayoffs(String fileName, char type) {
		try{
			PrintWriter pw = new PrintWriter(new File(fileName));
			
			// write payoffs
			for(int i = 0; i < resolution; i++){
				for(int j = 0; j < resolution; j++){
					if(type == 'D')
						pw.println(i + "," + j + "," + defR_lb[i][j] + "," + defR_ub[i][j] + "," + defP_lb[i][j] + "," + defP_ub[i][j]);
					
					if(type == 'A')
						pw.println(i + "," + j + "," + attR_lb[i][j] + "," + attR_ub[i][j] + "," + attP_lb[i][j] + "," + attP_ub[i][j]);
				}
			}
			
			pw.close();

		} catch (IOException e) {
			System.err.println("Payoffs files not found");
			e.printStackTrace();		
		}
	}

	public void generateGraphEdges() {		
		for(int i = 0; i < resolution; i++){
			for(int j = 0; j < resolution; j++){
				if(j < resolution - 1){
					GridEdge e1 = (GridEdge) graph.addEdge(grid[i][j], grid[i][j + 1]);
					graph.setEdgeWeight(e1, 1.0);					
				}

				if(i < resolution - 1){
					GridEdge e2 = (GridEdge) graph.addEdge(grid[i][j], grid[i + 1][j]);
					graph.setEdgeWeight(e2, 1.0);
				}
			}
		}
		
		graphGenerated = true; 
	}
	
	public void setStartingPosition(int start_i, int start_j) {
		grid[start_i][start_j].cell.start = true;
	}
	public void resetStartingPosition(int start_i, int start_j)
	{
		grid[start_i][start_j].cell.start = false;
	}

	public double[] readRegret(String fileName, int nbrOfTargets) {
		double[] regret = new double[nbrOfTargets];
		int t = 0;
		
		try {
			Scanner scReg = new Scanner(new FileReader(fileName));
			
			// read payoffs			
			while(scReg.hasNext()){
				String[] regArray = scReg.next().split(",");
				
				for(int i = 0; i < regArray.length; i++){
					regret[t] = Integer.parseInt(regArray[i]);
					
					t++;
				}
			}
			
			scReg.close();

		} catch (IOException e) {
			System.err.println("Regret file not found");
			e.printStackTrace();		
		}				
		
		return regret;
	}	
}*/