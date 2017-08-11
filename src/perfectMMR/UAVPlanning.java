/*package perfectMMR;

import geoInformation.Area;
import geoInformation.Cell;
import geoInformation.Grid;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;

import lpWrapper.LPSolverException;
import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import models.AttBoundStructure;
import models.PayoffStructure;
import models.SUQRAdversary;
import models.Target;
// Currently, only uncertainty in the animal density (def penalty / att reward)
public class UAVPlanning {
	List<Target> targetList;
	
	//Connectivity information
	Map<Integer, List<Integer>> connectivityMap;
	int pathLength;
	//
	SUQRAdversary adversary;
	boolean isZeroSum = true;
	int numPiece = 5;
	int numSamples = 1;
//	boolean isExactHeuristic = false;
	int numRes;
	int numRound = 2;
	
	int heuristic;
	// Output
	public double[] regret;
	public double[] runtime;
	public double[] trueRegret;
	
	public static final int RANDOM = 0;
	public static final int GREEDY_APPROX = 1;
	public static final int MCNF_APPROX = 2;
	public static final int GREEDY_EXACT = 3;
	public static final int GREEDY_APPROX_SAMPLE = 4;
	public static final int RANDOM_MULTIPLE_APPROX = 5;
	public static final int RANDOM_MULTIPLE_EXACT = 6;
	public static final int GREEDY_UTILITYDIS = 7;
	public static final int GREEDY_UTILITYDIS_PESS = 10;
	public static final int GREEDY_COVDIS = 8;
	public static final int GREEDY_COVDIS_PESS = 12;
	public static final int GREEDY_ATTPROBDIS = 9;
	public static final int GREEDY_ATTPROBDIS_PESS = 11;
	
	public static final int MCNF_UTILITYDIS = 13;
	public static final int MCNF_UTILITYDIS_PESS = 14;
	public static final int MCNF_COVDIS = 15;
	public static final int MCNF_COVDIS_PESS = 16;
	public static final int MCNF_ATTPROBDIS = 17;
	public static final int MCNF_ATTPROBDIS_PESS = 18;
			
	public UAVPlanning(List<Target> targetList, SUQRAdversary adversary, int numRes, int pathLength)
	{
		this.targetList = targetList;
		this.adversary = adversary;
		this.numRes = numRes;
		connectivityMap = loadConnectivityMap(targetList.size());
		this.pathLength = pathLength;
	}
	public void setHeuristic(int heuristic)
	{
		this.heuristic = heuristic;
	}
	public void setNumRound(int numRound)
	{
		this.numRound = numRound;
	}

	
	public int[] selectPathGreedyUtilityDis(int curRound, boolean isOptimistic) throws IOException, LPSolverException, MatlabConnectionException, MatlabInvocationException
	{
		int nTargets = targetList.size();
		int[] pathList = new int[pathLength];
		double[] utilityDis = new double[nTargets];
		for(Target tprime : targetList)
		{
			List<PayoffStructure> payoffList = tprime.getPayoffList();
			if(payoffList != null && !payoffList.isEmpty())
				payoffList.clear();
		}
		double[] adv_lb = new double[2 * nTargets];
		double[] adv_ub = new double[2 * nTargets];
		for(Target t : targetList)
		{
			double curRewardLB = t.getAttBoundStructure().getAttRewardLB();
			double curRewardUB = t.getAttBoundStructure().getAttRewardUB();
			double curPenaltyLB = t.getAttBoundStructure().getAttPenaltyLB();
			double curPenaltyUB = t.getAttBoundStructure().getAttPenaltyUB();
			adv_lb[t.id] = curRewardLB;
			adv_lb[t.id + nTargets] = curPenaltyLB;
			adv_ub[t.id] = curRewardUB;
			adv_ub[t.id + nTargets] = curPenaltyUB;
		}
		MinimaxRegret myTest = new MinimaxRegret(nTargets, numRes, adv_lb, adv_ub);
		myTest.computeMinimaxRegret();
		Map<Target, Double> defCov = new HashMap<Target, Double>();

		for(Target t : targetList)
		{
			defCov.put(t, myTest.def_solution[t.id]);
		}
		MaxRegretZeroSum mr = new MaxRegretZeroSum(targetList, defCov, adversary, numRes, 30);
		mr.solve();
		System.out.println("\nMax Regret:" + mr.getMaxRegret());
		trueRegret[curRound] = mr.getMaxRegret();
		regret[curRound] = mr.getMaxRegret();
		
		
		double[] attProb = new double[nTargets];
		double[] attOpponentProb = new double[nTargets];
	
		double[] optCov = myTest.def_solution;
		double[] opponentCov = new double[nTargets];
		for(Target t : targetList)
		{
			opponentCov[t.id] = mr.getOptDefCov(t.id);
		}

		double maxAttU = Double.NEGATIVE_INFINITY;
		double maxAttOpponentU = Double.NEGATIVE_INFINITY;
		for(Target t : targetList)
		{
			double attUtility = (mr.getAttPenalty(t.id) - mr.getAttReward(t.id)) * optCov[t.id] + mr.getAttReward(t.id);
			attProb[t.id] = attUtility;
			maxAttU = maxAttU < attUtility ? attUtility : maxAttU;
			
			double attOpponentUtility = (mr.getAttPenalty(t.id) - mr.getAttReward(t.id)) * opponentCov[t.id] + mr.getAttReward(t.id);
			attOpponentProb[t.id] = attOpponentUtility;
			maxAttOpponentU = maxAttOpponentU < attOpponentUtility ? attOpponentUtility : maxAttOpponentU;
		}
		boolean isChosen = false;
		boolean isOpponentChosen = false;
		for(Target t : targetList)
		{
			double defUtility = 0.0;
			if(attProb[t.id] == maxAttU && !isChosen)
			{
				defUtility = -maxAttU;
				isChosen = true;
			}
			
			double defOpponentUtility = 0.0;
			if(attOpponentProb[t.id] == maxAttOpponentU && !isOpponentChosen)
			{
				defOpponentUtility = -maxAttOpponentU;
				isOpponentChosen = true;
			}
			if(isOptimistic)
				utilityDis[t.id] = defOpponentUtility - defUtility;
			else
				utilityDis[t.id] = defUtility - defOpponentUtility;
		}
		mr.end();
		
		
		HashSet<Integer> pathSet = new HashSet<Integer>();
		double minRegret = Double.POSITIVE_INFINITY;
		int keyTarget = -1;
		for(int i = 0; i < nTargets; i++)
		{
			if(minRegret > utilityDis[i])
			{
				
				minRegret = utilityDis[i];
				keyTarget = i;
			}
		}
		List<Integer> candidateTargetList = new ArrayList<Integer>();
		for(int i = 0; i < nTargets; i++)
			if(minRegret >= utilityDis[i] - 1e-4)
				candidateTargetList.add(i);
		Random rand = new Random();
		int keyIndex = (int)(rand.nextDouble() * candidateTargetList.size());
		if(keyIndex == candidateTargetList.size())
			keyIndex = candidateTargetList.size() - 1;
		keyTarget = candidateTargetList.get(keyIndex);
		pathList[0] = keyTarget;
		pathSet.add(keyTarget);
		candidateTargetList.clear();
		// End of first target
		int nextTargetID = keyTarget;
		int prevTargetID = -1;
		
		for(int k = 1; k < pathLength; k++)
		{
			System.out.println("Path length: " + k);
			minRegret = Double.POSITIVE_INFINITY;
			keyTarget = -1;
			prevTargetID = nextTargetID;
			List<Integer> neighborList = connectivityMap.get(prevTargetID);
			for(int i = 0; i < neighborList.size(); i++)
			{
				int curTargetID = neighborList.get(i);
				if(!pathSet.contains(curTargetID))
				{
					if(minRegret > utilityDis[curTargetID])
					{
						minRegret = utilityDis[curTargetID];
						keyTarget = curTargetID;
					}
				}				
			}
			nextTargetID = keyTarget;
			pathList[k] = keyTarget;
			pathSet.add(keyTarget);
		}
		pathSet.clear();
		
//		mmr.end();
		return pathList;
	}
	public int[] selectPathMCNFUtilityDis(int pathSize, boolean loops, int curRound, boolean isOptimistic) 
			throws IOException, LPSolverException, MatlabConnectionException, MatlabInvocationException
	{
		int[] pathList = new int[pathSize];
		int nTargets = targetList.size();
		double[] utilityDis = new double[nTargets];
		for(Target tprime : targetList)
		{
			List<PayoffStructure> payoffList = tprime.getPayoffList();
			if(payoffList != null && !payoffList.isEmpty())
				payoffList.clear();
		}
		double[] adv_lb = new double[2 * nTargets];
		double[] adv_ub = new double[2 * nTargets];
		for(Target t : targetList)
		{
			double curRewardLB = t.getAttBoundStructure().getAttRewardLB();
			double curRewardUB = t.getAttBoundStructure().getAttRewardUB();
			double curPenaltyLB = t.getAttBoundStructure().getAttPenaltyLB();
			double curPenaltyUB = t.getAttBoundStructure().getAttPenaltyUB();
			adv_lb[t.id] = curRewardLB;
			adv_lb[t.id + nTargets] = curPenaltyLB;
			adv_ub[t.id] = curRewardUB;
			adv_ub[t.id + nTargets] = curPenaltyUB;
		}
		MinimaxRegret myTest = new MinimaxRegret(nTargets, numRes, adv_lb, adv_ub);
		myTest.computeMinimaxRegret();
		Map<Target, Double> defCov = new HashMap<Target, Double>();

		for(Target t : targetList)
		{
			defCov.put(t, myTest.def_solution[t.id]);
		}
		MaxRegretZeroSum mr = new MaxRegretZeroSum(targetList, defCov, adversary, numRes, 30);
		mr.solve();
		System.out.println("\nMax Regret:" + mr.getMaxRegret());
		trueRegret[curRound] = mr.getMaxRegret();
		regret[curRound] = mr.getMaxRegret();
		
		
		double[] attProb = new double[nTargets];
		double[] attOpponentProb = new double[nTargets];
	
		double[] optCov = myTest.def_solution;
		double[] opponentCov = new double[nTargets];
		for(Target t : targetList)
		{
			opponentCov[t.id] = mr.getOptDefCov(t.id);
		}

		double maxAttU = Double.NEGATIVE_INFINITY;
		double maxAttOpponentU = Double.NEGATIVE_INFINITY;
		for(Target t : targetList)
		{
			double attUtility = (mr.getAttPenalty(t.id) - mr.getAttReward(t.id)) * optCov[t.id] + mr.getAttReward(t.id);
			attProb[t.id] = attUtility;
			maxAttU = maxAttU < attUtility ? attUtility : maxAttU;
			
			double attOpponentUtility = (mr.getAttPenalty(t.id) - mr.getAttReward(t.id)) * opponentCov[t.id] + mr.getAttReward(t.id);
			attOpponentProb[t.id] = attOpponentUtility;
			maxAttOpponentU = maxAttOpponentU < attOpponentUtility ? attOpponentUtility : maxAttOpponentU;
		}
		boolean isChosen = false;
		boolean isOpponentChosen = false;
		for(Target t : targetList)
		{
			double defUtility = 0.0;
			if(attProb[t.id] == maxAttU && !isChosen)
			{
				defUtility = -maxAttU;
				isChosen = true;
			}
			
			double defOpponentUtility = 0.0;
			if(attOpponentProb[t.id] == maxAttOpponentU && !isOpponentChosen)
			{
				defOpponentUtility = -maxAttOpponentU;
				isOpponentChosen = true;
			}
			if(isOptimistic)
				utilityDis[t.id] = defOpponentUtility - defUtility;
			else
				utilityDis[t.id] = defUtility - defOpponentUtility;
		}
		mr.end();
		
		// we fix the starting position
		int start_j = 0;
		int start_i = 0;
		double minRegret = Double.POSITIVE_INFINITY;
		// we calculate the resolution
		int resolution = (int) Math.sqrt(targetList.size());
		Grid grid = new Grid(new Area("Test", 1000 * resolution, 1000 * resolution, 0, 30), resolution);
		
		grid.instantiateCellList();
		
		for(start_j = 0; start_j < resolution; start_j++)
		{
			for(start_i = 0; start_i < resolution; start_i++)
			{
				grid.setStartingPosition(start_i, start_j);
			}
		}
		
		double avgRegret = 0.0;
		
		MinNetworkFlow flow = new MinNetworkFlow(grid, utilityDis, pathSize, loops);		
		flow.solve();
		
		Cell[] cells = flow.getPathCells();
		int ids[] = new int[cells.length];
		for(int i = 0; i < ids.length; i++) 
			avgRegret += utilityDis[cells[i].id - 1];
		avgRegret /= pathSize;
		if(minRegret > avgRegret)
		{
			for(int i = 0; i < ids.length; i++)
				pathList[i] = cells[i].id - 1;
			minRegret = avgRegret;
		}
		flow.end();
//		
//		mmr.end();
		return pathList;
	}
	

	public int[] selectPathRandom()
	{
		int[] pathList = new int[pathLength];
		int numTargets = targetList.size();
		int prevTargetID = -1;
		int nextTargetID = -1;
		Random rand = new Random();
		nextTargetID = rand.nextInt(numTargets - 1);
		pathList[0] = nextTargetID;
		HashSet<Integer> pathSet = new HashSet<Integer>();
		pathSet.add(nextTargetID);
		for(int k = 1; k < pathLength; k++)
		{
			prevTargetID = nextTargetID;
			List<Integer> neighborList = connectivityMap.get(prevTargetID);
			int index = (int)(rand.nextDouble() * neighborList.size());
			if(index == neighborList.size())
				index = neighborList.size() - 1;
			nextTargetID = neighborList.get(index);
			
			while(pathSet.contains(nextTargetID))
			{
				index = (int)(rand.nextDouble() * neighborList.size());
				if(index == neighborList.size())
					index = neighborList.size() - 1;
				nextTargetID = neighborList.get(index);
			}
			pathList[k] = nextTargetID;
			pathSet.add(nextTargetID);
		}
		return pathList;
	}
	public double computeTrueRegret() throws IOException, LPSolverException, MatlabConnectionException, MatlabInvocationException
	{
		double regret = 0;
		int nTargets = targetList.size();
		double[] adv_lb = new double[2 * nTargets];
		double[] adv_ub = new double[2 * nTargets];
		for(Target t : targetList)
		{
			double curRewardLB = t.getAttBoundStructure().getAttRewardLB();
			double curRewardUB = t.getAttBoundStructure().getAttRewardUB();
			double curPenaltyLB = t.getAttBoundStructure().getAttPenaltyLB();
			double curPenaltyUB = t.getAttBoundStructure().getAttPenaltyUB();
			adv_lb[t.id] = curRewardLB;
			adv_lb[t.id + nTargets] = curPenaltyLB;
			adv_ub[t.id] = curRewardUB;
			adv_ub[t.id + nTargets] = curPenaltyUB;
		}
		MinimaxRegret myTest = new MinimaxRegret(nTargets, numRes, adv_lb, adv_ub);
		myTest.computeMinimaxRegret();
		Map<Target, Double> defCov = new HashMap<Target, Double>();

		for(Target t : targetList)
		{
			defCov.put(t, myTest.def_solution[t.id]);
		}
		MaxRegretZeroSum mr = new MaxRegretZeroSum(targetList, defCov, adversary, numRes, 30);
		mr.solve();
		regret = mr.getMaxRegret();
		mr.end();
		return regret;
	}
	public void PESimple(int numRound, int heuristic) throws IOException, LPSolverException, MatlabConnectionException, MatlabInvocationException
	{
		regret = new double[numRound];
		runtime = new double[numRound];
		trueRegret = new double[numRound];
		long start = 0;
		long end = 0;
		for(int i = 0; i < numRound; i++)
		{
			int[] pathList = null;
			start = System.currentTimeMillis();
			if(heuristic == RANDOM)
			{
				regret[i] = trueRegret[i] = computeTrueRegret();
				pathList = selectPathRandom();
			}
			if(heuristic == GREEDY_UTILITYDIS)
			{
				boolean isOptimistic = true;
				pathList = selectPathGreedyUtilityDis(i, isOptimistic);
			}
			if(heuristic == GREEDY_UTILITYDIS_PESS)
			{
				boolean isOptimistic = false;
				pathList = selectPathGreedyUtilityDis(i, isOptimistic);
			}
			
			if(heuristic == MCNF_UTILITYDIS)
			{
				boolean isOptimistic = true;
				pathList = selectPathMCNFUtilityDis(pathLength, false, i, isOptimistic);
			}
			if(heuristic == MCNF_UTILITYDIS_PESS)
			{
				boolean isOptimistic = false;
				pathList = selectPathMCNFUtilityDis(pathLength, false, i, isOptimistic);
			}
			
			// Do simulation
			for(int k = 0; k < pathList.length; k++)
			{
				int targetID = pathList[k];
				Target t = targetList.get(targetID);
				double curRewardLB = t.getAttBoundStructure().getAttRewardLB();
				double curRewardUB = t.getAttBoundStructure().getAttRewardUB();
				double curPenaltyLB = t.getAttBoundStructure().getAttPenaltyLB();
				
				double avgB = (curRewardLB + curRewardUB) / 2;
				PayoffStructure truePayoff = t.getTruePayoffStructure();
				if(truePayoff.getAdversaryReward() < avgB)
					t.setAttBoundStructure(new AttBoundStructure(curRewardLB, curPenaltyLB, avgB, curPenaltyLB));
				else
					t.setAttBoundStructure(new AttBoundStructure(avgB, curPenaltyLB, curRewardUB, curPenaltyLB));
				
			}
			
			end = System.currentTimeMillis();
			if(i > 0)
				runtime[i] = runtime[i - 1] + (end - start) / 1000.0;
			else runtime[i] = (end - start) / 1000.0;
			System.out.print("Round " + i + ": " + regret[i] + "\t");
			for(int k = 0; k < pathList.length; k++)
			{
				System.out.print(pathList[k] + "\t");
			}
			System.out.println();
		}
	}
	
	
	public void testPE() throws IOException, LPSolverException, MatlabConnectionException, MatlabInvocationException
	{
		PESimple(numRound, heuristic);
	}
	
	public Map<Integer, List<Integer>> createConnectivityMap(int numTargets)
	{
		Map<Integer, List<Integer>> areaMap = new HashMap<Integer, List<Integer>>();
		int size = (int) Math.sqrt(numTargets);
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < size; j++)
			{
				int targetID = i * size + j;
				List<Integer> neighborList = new ArrayList<Integer>();
				if(i > 0)
				{
					int neighborID = (i - 1) * size + j;
					neighborList.add(neighborID);
				}
				if(i < size - 1)
				{
					int neighborID = (i + 1) * size + j;
					neighborList.add(neighborID);
				}
				if(j > 0)
				{
					int neighborID = i * size + j - 1;
					neighborList.add(neighborID);
				}
				if(j < size - 1)
				{
					int neighborID = i * size + j + 1;
					neighborList.add(neighborID);
				}
				areaMap.put(targetID, neighborList);
			}
		}
		// Store map
		String filePath = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/UAVPLANNING/" + numTargets + "T/areaMap.csv";
		
			FileOutputStream mapFileStream = null;
			try {
				mapFileStream = new FileOutputStream(filePath);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			PrintStream mapStream = new PrintStream(mapFileStream);
			for(int i = 0; i < numTargets; i++)
			{
				List<Integer> neighborList = areaMap.get(i);
				mapStream.print(i + ",");
				for(int j = 0; j < neighborList.size(); j++)
				{
					mapStream.print(neighborList.get(j) + ",");
				}
				mapStream.println();
			} 
			mapStream.close();
			
		return areaMap;
	}
	
	public Map<Integer, List<Integer>> loadConnectivityMap(int numTargets)
	{
		Map<Integer, List<Integer>> areaMap = new HashMap<Integer, List<Integer>>();
		FileInputStream in = null;
		FileReader fin = null;
		Scanner src = null;
		
		String filePath = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/UAVPLANNING/" + numTargets + "T/areaMap.csv";
		try {
			fin = new FileReader(filePath);
			src = new Scanner(fin);
			for(int i = 0; i < numTargets; i++)
			{
				String[] data = src.nextLine().split(",");
				int targetID = Integer.parseInt(data[0]);
				List<Integer> neighborList = new ArrayList<Integer>();
				for(int j = 1; j < data.length; j++)
				{
					int neighborID = Integer.parseInt(data[j]);
					neighborList.add(neighborID);
				}
				areaMap.put(targetID, neighborList);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.out.println("Couldn't open file for reading.");
			e.printStackTrace();
		}
		src.close();
		return areaMap;
	}
	public void end()
	{
		for(int i = 0; i < targetList.size(); i++)
		{
			connectivityMap.get(i).clear();
		}
		connectivityMap.clear();
		for(Target t : targetList)
		{
//			if(!t.getPayoffList().isEmpty())
//				t.getPayoffList().clear();
		}
		targetList.clear();
	}
}
*/