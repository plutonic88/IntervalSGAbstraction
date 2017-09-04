package cs.Interval.contraction;

import java.util.ArrayList;
import java.util.HashMap;

public class TargetNode {
	private int targetid;
	private double animaldensity;
	
	public double defenderreward;
	public double defenderpenalty;
	
	public double attackerreward;
	public double attackerpenalty;
	public double distfrombase;
	public boolean isinattackcluster = false;
	
	public double getDistfrombase() {
		return distfrombase;
	}



	public void setDistfrombase(double distfrombase) {
		this.distfrombase = distfrombase;
	}





	double collectedcoinvalue;
	TargetNode maxcoinparent;
	int row;
	int col;
	boolean goal;
	boolean start;
	public TargetNode parent;
	private double distfromstart;
	boolean in;
	boolean out;
	int edgeid;
	boolean usedforconnection = false;
	int tranformedtargetid;
	public double distancecoveredyet;
	


	public boolean isIn() {
		return in;
	}



	public void setIn(boolean in) {
		this.in = in;
	}



	public boolean isOut() {
		return out;
	}



	public void setOut(boolean out) {
		this.out = out;
	}
	
	



	//private boolean isaccessible;
	private ArrayList<TargetNode> neighbors = new ArrayList<TargetNode>();
	private HashMap<TargetNode, Double> distances = new HashMap<TargetNode, Double>();
	private HashMap<TargetNode, ArrayList<TargetNode>> path = new HashMap<TargetNode, ArrayList<TargetNode>>();
	private HashMap<TargetNode, Double> pathutility = new HashMap<TargetNode, Double>();

	
	public  TargetNode(TargetNode node)
	{
		this.row = node.row;
		this.col = node.col;
		this.targetid = node.targetid;
		this.animaldensity = node.animaldensity;
		//this.maxcoinparent = node.maxcoinparent;
		
		this.goal = node.goal;
		this.start = node.start;
		this.parent = new TargetNode();
		this.distancecoveredyet =0;
		//private double distfromstart;
		
		
		
	}

	
	
	public double getDistfromstart() {
		return distfromstart;
	}

	public void setDistfromstart(double distfromstart) {
		this.distfromstart = distfromstart;
	}

	public TargetNode getParent()
	{
		return this.parent;
	}
	
	public void setParent(TargetNode parent)
	{
		this.parent = parent;
	}
	
	
	public static double getAvgDistance(TargetNode node)
	{
		double sum = 0;
		for(TargetNode n: node.neighbors)
		{
			sum += node.getDistance(n);
		}
		sum /= node.neighbors.size();
		return sum;
	}
	
	public void setPathUtility(TargetNode dest, Double utility)
	{
		this.pathutility.put(dest, utility);
	}
	
	public Double getPathUtility(TargetNode dest)
	{
		return this.pathutility.get(dest);
	}
	
	public Double removePathUtility(TargetNode dest)
	{
		return this.pathutility.remove(dest);
	}
	
	
	
	public boolean isGoal() {
		return goal;
	}

	public void setGoal(boolean goal) {
		this.goal = goal;
	}

	public boolean isStart() {
		return start;
	}

	public void setStart(boolean start) {
		this.start = start;
	}
	
	
	public  HashMap<TargetNode, ArrayList<TargetNode>> getPath() {
		return this.path;
	}
	
	public ArrayList<TargetNode> getPath(TargetNode node) {
		return path.get(node);
	}

	public void setPath(TargetNode dest, ArrayList<TargetNode> path) {
		if(this.path.containsKey(dest))
		{
			this.path.remove(dest);
		}
		this.path.put(dest, path);
	}

	public TargetNode()
	{}

	public TargetNode(int targetid, double animaldensity) 
	{
		super();
		this.targetid = targetid;
		this.animaldensity = animaldensity;
		this.collectedcoinvalue = 0;
		this.maxcoinparent = new TargetNode();
		this.start = false;
		this.goal = false;
		

	}
	
	public void setRowCol(int row, int col)
	{
		this.row = row;
		this.col= col;
	}
	
	public int getRow()
	{
		return this.row;
	}
	
	public int getCol()
	{
		return this.col;
	}


	public TargetNode getMaxCoinParent() {
		return maxcoinparent;
	}

	public void setMaxCoinParent(TargetNode parent) {
		this.maxcoinparent = parent;
	}

	public double getCoinvalue() {
		return collectedcoinvalue;
	}

	public void setCoinvalue(double coinvalue) {
		this.collectedcoinvalue = coinvalue;
	}

	public int getTargetid() {
		return targetid;
	}

	public void setTargetid(int targetid) {
		this.targetid = targetid;
	}

	public double getAnimaldensity() {
		return animaldensity;
	}

	public void setAnimaldensity(double animaldensity) {
		this.animaldensity = animaldensity;
	}


	public HashMap<TargetNode, Double> getDistances() {
		return distances;
	}

	public Double getDistance(TargetNode dest) {
		return distances.get(dest);
	}

	public void setDistances(HashMap<TargetNode, Double> distances) {
		this.distances = distances;
	}

	public ArrayList<TargetNode> getNeighbors() {
		return neighbors;
	}

	public void setNeighbors(ArrayList<TargetNode> neighbors)
	{
		for(TargetNode neighbor: neighbors)
		{
			this.neighbors.add(neighbor);
		}

	}
	public void addNeighbor(TargetNode node)
	{
		if(this.neighbors.contains(node))
		{
			this.neighbors.remove(node);
		}
		this.neighbors.add(node);
	}

	public void removeNeighbor(TargetNode neighbor)
	{
		//this.neighbors.remove(this.neighbors.indexOf(neighbor));
		this.neighbors.remove(neighbor);
	}

	public void setDistances(double[] distances)
	{
		int distanceindex = 0;
		for(TargetNode neighbor: this.neighbors)
		{
			this.distances.put(neighbor, distances[distanceindex]);
		}
	}

	public void addDistance(TargetNode node, Double distance)
	{
		if(this.distances.containsKey(node))
		{
			this.distances.remove(node);
		}
		this.distances.put(node, distance);
	}

	/*public void addNodeToPath(TargetNode pathnode)
	{
		this.path.get(this).add(pathnode);
	}*/

	public void removeDistance(TargetNode node)
	{
		this.distances.remove(node);
	}

	public void removePath(TargetNode node)
	{
		this.path.remove(node);
	}



	public double computeDistance(TargetNode src, TargetNode dest) 
	{
		/**
		 * if destination is listed in neighbors
		 */
		
		/*SecurityGameContraction.printNeighbors(src);
		SecurityGameContraction.printDistances(src);
		SecurityGameContraction.printNeighbors(dest);
		SecurityGameContraction.printDistances(dest);
*/
		if(src.getNeighbors().contains(dest))
		{
			/**
			 * get the path
			 *//*
			ArrayList<TargetNode> pathnodes = src.getPath(dest);
			double distsum = 0;
			if(pathnodes.size()>=1)
			{
				distsum += src.getDistances().get(pathnodes.get(0));
				for(int  targetnodeid = 0; targetnodeid< (pathnodes.size()-1); targetnodeid++)
				{
					distsum += pathnodes.get(targetnodeid).getDistances().get(pathnodes.get(targetnodeid+1));
				}
				distsum += pathnodes.get(pathnodes.size()-1).getDistances().get(dest);
			}
			else if(pathnodes.size()==0)
			{
				distsum += src.getDistances().get(dest);
			}*/
			return src.getDistance(dest);
		}
		else
			return Double.MAX_VALUE;


	}



}
