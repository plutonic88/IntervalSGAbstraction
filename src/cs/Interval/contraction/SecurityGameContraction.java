package cs.Interval.contraction;





import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import com.cs.kruskal.KruskalAlgorithm;

import cs.Interval.Abstraction.SecurityGameAbstraction;
import cs.Interval.ILP.MIPSolver3;
import cs.Interval.ILP.MIPSolver4;
import cs.com.allpair.AllPairShortestPath;
import cs.com.realworld.ReadData;
import groupingtargets.GroupingTargets;
import groupingtargets.SuperTarget;

public class SecurityGameContraction 
{
	private static final int INFINITY = 9999999;
	//public HashMap<Integer, Integer> unaccessibltargets = new HashMap<Integer, Integer>();
	public static ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	public static ArrayList<TargetNode> spanningtree = new ArrayList<TargetNode>();
	public static ArrayList<TargetNode> origtargets = new ArrayList<TargetNode>();
	public static ArrayList<TargetNode> duplicatetargets = new ArrayList<TargetNode>();
	public static TargetNode graph = new TargetNode();
	public static Random rand1 = new Random(100);
	//public static ArrayList<Integer> unimportanttargets = new ArrayList<Integer>();
	public static double instantcontractionshortestdist = 0;

	public SecurityGameContraction(int nrow, int ncol, double[][] u, double[][] e, int[][] gamedata)
	{
		super();

		buildGraph(nrow, ncol, u, e, gamedata);
	}

	public SecurityGameContraction(int numRow, int numCol, int[][] gamedata) 
	{
		super();

		//buildGraph(numRow, numCol, gamedata);

	}
	
	
	public SecurityGameContraction(int nTargets, int[][] gamedata) 
	{
		super();

		//buildGraph(numRow, numCol, gamedata);

	}
	
	public SecurityGameContraction() 
	{
		super();

		//buildGraph(numRow, numCol, gamedata);

	}


	public static ArrayList<TargetNode> makeSpanningTree(ArrayList<TargetNode> targets, int nTargets, ArrayList<TargetNode> domindatednodes)
	{
		ArrayList<TargetNode> spanningtree = new ArrayList<TargetNode>();
		/**
		 * make the adjacency matrix
		 */
		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		int icount =1;

		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}



		makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);


		System.out.println("adjacency matrix is ");
		for (int i = 1; i <= nTargets; i++)
			System.out.print("\t" + i);
		System.out.println();
		for (int source = 1; source <= nTargets; source++)
		{
			System.out.print(source + "\t");
			for (int destination = 1; destination <= nTargets; destination++)
			{
				System.out.print(adjacencymatrix[source][destination] + "\t");
			}
			System.out.println();
		}

		KruskalAlgorithm kruskalAlgorithm = new KruskalAlgorithm(nTargets);
		int[][] sanningmatrix= kruskalAlgorithm.kruskalAlgorithm(adjacencymatrix);


		/**
		 * convert from matrix to set of targets
		 */

		convertMatrixToMST(sanningmatrix,spanningtree, targets, domindatednodes, map, mapback);





		return spanningtree;
	}



	private static void convertMatrixToMST(int[][] adjacencymatrix,
			ArrayList<TargetNode> spanningtree, ArrayList<TargetNode> targets, ArrayList<TargetNode> dominatedtargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) 
	{
		int nTargets = targets.size();


		for(int i=0; i<nTargets; i++)
		{
			TargetNode s = targets.get(i);
			TargetNode tmp= new TargetNode(s.getTargetid(), s.getAnimaldensity());
			tmp.attackerpenalty = s.attackerpenalty;
			tmp.attackerreward = s.attackerreward;
			tmp.defenderpenalty= s.defenderpenalty;
			tmp.defenderreward = s.defenderreward;
			spanningtree.add(tmp);
		}

		/**
		 * add the dominated targets too
		 */

		for(int i=0; i<dominatedtargets.size(); i++)
		{
			TargetNode s = dominatedtargets.get(i);
			TargetNode tmp= new TargetNode(s.getTargetid(), s.getAnimaldensity());
			tmp.attackerpenalty = s.attackerpenalty;
			tmp.attackerreward = s.attackerreward;
			tmp.defenderpenalty= s.defenderpenalty;
			tmp.defenderreward = s.defenderreward;
			spanningtree.add(tmp);
		}







		for (int source = 1; source <= (nTargets); source++)
		{

			for (int destination = 1; destination <= nTargets; destination++)
			{
				if(adjacencymatrix[source][destination]>0 && adjacencymatrix[source][destination]!=999)
				{
					/**
					 * create 
					 */

					System.out.println("src "+source+", destination "+destination + "="+adjacencymatrix[source][destination]);


					int srcid = mapback.get(source);
					int destid = mapback.get(destination);
					System.out.println("orig src "+srcid+", orig destination "+destid);

					TargetNode src = getTargetNode(srcid, spanningtree); // one less in the spanning tree
					TargetNode dest = getTargetNode(destid, spanningtree); // one less in the spanning tree
					src.addNeighbor(dest);
					src.addDistance(dest, (double)adjacencymatrix[source][destination]);

					TargetNode origsrc = getTargetNode(srcid, targets); // one less in the spanning tree
					TargetNode origdest = getTargetNode(destid, targets); // one less in the spanning tree

					ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();

					System.out.println("PATH N0DES ");
					for(TargetNode pnode:	origsrc.getPath(origdest))
					{
						TargetNode spanningpnde = getTargetNode(pnode.getTargetid(), spanningtree);
						System.out.print(spanningpnde.getTargetid() +" ");
						pathnodes.add(spanningpnde);
					}
					src.setPath(dest, pathnodes);
					src.setPathUtility(dest, origsrc.getPathUtility(origdest));
					System.out.println("PATH utility "+ origsrc.getPathUtility(origdest));

				}

			}

		}


	}


	
	private static void purifyAPSPMatrixWithNeighborAndZero(int[][] apsmat,
			ArrayList<TargetNode> targets, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		/*int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}*/


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();




				if (source == destination)
				{
					apsmat[source][destination] = 0;  // reason ache
					continue;
				}
				if (apsmat[source][destination] == 0)
				{
					apsmat[source][destination] = INFINITY;
				}
			}
		}



		//int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(int j=1; j<= nTargets; j++)
			{
				int target = mapback.get(j);
				TargetNode tmp = getTargetNode(target, targets);
				if(!n.getNeighbors().contains(tmp))
				{

					if(n.getTargetid() != target)
					{

						/*System.out.print("["+n.getTargetid()+"]["+target+"] ---> ");
						System.out.println("["+map.get(n.getTargetid())+"]["+map.get(target)+"]=INF");*/
						apsmat[map.get(n.getTargetid())][map.get(target)]=  INFINITY;
					}

				}


			}
			//i++;
		}






	}
	

	private static void purifyAPSPSuperMatrixWithNeighborAndZero(int[][] apsmat,
			HashMap<Integer, TargetNode> targetmaps, int nTargets, 
			HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback, HashMap<Integer, SuperTarget> sts) {


		/*int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}*/


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();




				if (source == destination)
				{
					apsmat[source][destination] = 0;  // reason ache
					continue;
				}
				if (apsmat[source][destination] == 0)
				{
					apsmat[source][destination] = INFINITY;
				}
			}
		}



		//int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(int j=1; j<= nTargets; j++)
			{
				int target = mapback.get(j);
				TargetNode tmp = getTargetNode(target, targets);
				if(!n.getNeighbors().contains(tmp))
				{

					if(n.getTargetid() != target)
					{

						/*System.out.print("["+n.getTargetid()+"]["+target+"] ---> ");
						System.out.println("["+map.get(n.getTargetid())+"]["+map.get(target)+"]=INF");*/
						apsmat[map.get(n.getTargetid())][map.get(target)]=  INFINITY;
					}

				}


			}
			//i++;
		}






	}
	
	
	private static void purifyAPSPSTMatrixZero(int[][] adjacencymatrix,
			HashMap<Integer, SuperTarget> sts, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		/*int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}*/


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();




				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
				else if(adjacencymatrix[source][destination] < 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}



		/*//int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(int j=1; j<= nTargets; j++)
			{
				int target = mapback.get(j);
				TargetNode tmp = getTargetNode(target, targets);
				if(!n.getNeighbors().contains(tmp))
				{

					if(n.getTargetid() != target)
					{

						//System.out.print("["+n.getTargetid()+"]["+target+"] ---> ");
						//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(target)+"]=1");
						adjacencymatrix[map.get(n.getTargetid())][map.get(target)]=  INFINITY;
					}

				}


			}
			//i++;
		}
*/





	}
	
	
	
	public static void purifyAPSPMatrixZero(int[][] adjacencymatrix,
			ArrayList<TargetNode> targets, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		/*int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}*/


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();




				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
				else if(adjacencymatrix[source][destination] < 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}



		/*//int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(int j=1; j<= nTargets; j++)
			{
				int target = mapback.get(j);
				TargetNode tmp = getTargetNode(target, targets);
				if(!n.getNeighbors().contains(tmp))
				{

					if(n.getTargetid() != target)
					{

						//System.out.print("["+n.getTargetid()+"]["+target+"] ---> ");
						//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(target)+"]=1");
						adjacencymatrix[map.get(n.getTargetid())][map.get(target)]=  INFINITY;
					}

				}


			}
			//i++;
		}
*/





	}
	
	
	public static void purifyAPSPMatrixZeroGT(int[][] adjacencymatrix,
			ArrayList<TargetNode> targets, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		/*int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}*/


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();




				if (source == destination)
				{
					adjacencymatrix[source][destination] = 0;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
				else if(adjacencymatrix[source][destination] < 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}



		/*//int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(int j=1; j<= nTargets; j++)
			{
				int target = mapback.get(j);
				TargetNode tmp = getTargetNode(target, targets);
				if(!n.getNeighbors().contains(tmp))
				{

					if(n.getTargetid() != target)
					{

						//System.out.print("["+n.getTargetid()+"]["+target+"] ---> ");
						//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(target)+"]=1");
						adjacencymatrix[map.get(n.getTargetid())][map.get(target)]=  INFINITY;
					}

				}


			}
			//i++;
		}
*/





	}



	
	

	private static void makeAdjacencyMatrixST(int[][] adjacencymatrix,
			HashMap<Integer, SuperTarget> supertargets, int nSTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		int i=1; 
		for(SuperTarget n: supertargets.values())
		{
			//int j=1;
			for(SuperTarget nei: n.neighbors.values())
			{
				//System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.stid)][map.get(nei.stid)]= (int)GroupingTargets.minimumDist(n, nei);
			}
			i++;
		}


		for (int source = 1; source <=nSTargets; source++)
		{
			for (int destination = 1; destination <= nSTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();
				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}





	}


	public static void makeAdjacencyMatrix(int[][] adjacencymatrix,
			ArrayList<TargetNode> targets, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		int i=1; 
		for(TargetNode n: targets)
		{
			//int j=1;
			for(TargetNode nei: n.getNeighbors())
			{
				//System.out.print("["+n.getTargetid()+"]["+nei.getTargetid()+"] ---> ");
				//System.out.println("["+map.get(n.getTargetid())+"]["+map.get(nei.getTargetid())+"]=1");



				adjacencymatrix[map.get(n.getTargetid())][map.get(nei.getTargetid())]=  n.getDistance(nei).intValue();
			}
			i++;
		}


		for (int source = 1; source <=nTargets; source++)
		{
			for (int destination = 1; destination <= nTargets; destination++)
			{
				//adjacencymatrix[source][destination] = scan.nextInt();
				if (source == destination)
				{
					adjacencymatrix[source][destination] = INFINITY;
					continue;
				}
				if (adjacencymatrix[source][destination] == 0)
				{
					adjacencymatrix[source][destination] = INFINITY;
				}
			}
		}





	}



	public static void instantContractionWithExtreamPruningTest(double[][] density, int ITER, int nrow, int ncol, int[] percentages, int[] thresholds, double dmax)
	{



		int nRes = 2;
		//double perc = .20;
		int nTargets = nrow*ncol;
		//double thresholds[]={2,3,4,5};
		//int ITER=4;

		//double percentages[] = {10};
		double[] result = new double[percentages.length];
		int rindex=0;
		double defexp=0;

		//int lstart=1,  lend=4, hstart=8, hend=10;
		//double[][] density=generateRandomDensity( percentages[0], ITER, lstart, lend,  hstart, hend, nTargets);


		for(double perc: percentages)
		{

			for(double threshold: thresholds)
			{
				double sumsol = 0;
				//long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				int targetsize=0;


				for(int iter=0; iter<ITER; iter++)
				{

					targets.clear();

					int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

					//makeZeroSum(gamedata,nTargets);
					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);

					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
					//printEdges(targets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					selectDominatedTargets(targets, domindatednodes, threshold);
					//printNodesWithNeighborsAndPath(domindatednodes, targets);

					preComputeShortestPaths(0, targets, domindatednodes);
					Date start = new Date();
					long l1 = start.getTime();



					instantContractionWithExtreamPrunning(domindatednodes, targets, dmax);


					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;

					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					targetsize+= targets.size();
					int[][] p = new int[targets.size()][];
					try 
					{
						start = new Date();
						l1 = start.getTime();
						ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
						ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
						/**
						 * map has present id
						 * mapback gives the original ids
						 */
						HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
						HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
						makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
						System.out.println("Total path with duplicates "+pathseq.size());
						//printPaths(pathseq);
						pathseq = removeDuplicatePathSimple(pathseq);
						//printPaths(pathseq);
						System.out.println("Total path without duplicates "+pathseq.size());

						//int k = 2;
						Integer[] input = new Integer[pathseq.size()];
						int[] branch = new int[nRes];//{0,0};//new char[k];


						for(int i=0; i<input.length; i++)
						{
							input[i] = i;
						}
						HashSet jSet=new HashSet();
						if(pathseq.size()==0)
						{
							//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
							//choose the worst payoff for defender

							Double mAxpayoff = Double.MIN_VALUE;
							Double defpayoff = 0.0;
							for(int i=0; i<domindatednodes.size(); i++)
							{
								targets.add(domindatednodes.get(i));
							}
							for(TargetNode x: targets)
							{
								if(x.attackerreward>mAxpayoff)
								{
									mAxpayoff= x.attackerreward;
									defpayoff = x.defenderpenalty;
								}
							}

							sumsol += defpayoff;
							System.out.println("Defender expected payoff "+ defpayoff);
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;


						}
						else
						{
							//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
							if(pathseq.size()<nRes)
							{

								branch = new int[pathseq.size()];
								jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
							}
							else
							{
								jSet=combine(input, nRes, 0, branch, 0, jSet);
							}

							List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
							/**
							 * columns will be combination of paths for each resources. 
							 */
							/**
							 * pmat, where columns will be combination of paths. 
							 * rows are targets. 
							 * each entry will say whether the target is in the joint schedule
							 */
							//jSet.

							//printJointSchedule(jset);

							//printNodesAsNeighbors(dominatednodes);

							p = makePmat(pathseq, jset, mapback, targets);
							//printPathMat(p);

							/**
							 * remove duplicates from p
							 */
							//removeDuplicatesFromP(p);
							//System.out.println();
							//printPathMat(p);

							//System.out.println("Number of targets after contraction "+ targets.size());
							//System.out.println("mip in... ");
							System.out.println("Iter "+iter+", mip in... ");
							double[] probdistribution = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;
							if(probdistribution.equals(null))
							{
								throw new Exception("Prob null...");
							}
							System.out.println("Iter "+iter+", mip out... ");
							//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


							start = new Date();
							l1 = start.getTime();
							int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
							//removeDuplicatesFromP(origpmat);
							//printPathMat(origpmat);
							//System.out.println("\n after mapping back");
							//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);

							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							if(threshold>0)
							{
								revmaptime += diff;
							}


							/*for(int i=0; i<coverage.length; i++)
	{
		for(int j=0; j<coverage[i].length; j++)
		{
			if(coverage[i][j]>0)
			{
				System.out.println("selected path : " + j);

				for(int k=0; k<origpmat.length; k++)
				{
					System.out.print(origpmat[k][j] + " ");
				}
				System.out.println();

			}
		}
	}*/			

							int maxtargetforattacker = findAttackTarget(origpmat, probdistribution, gamedata);

							double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);
							defexp = defexpectedpayoff;
							System.out.println("Attacked target is "+ maxtargetforattacker);
							System.out.println("Defender expected payoff "+ defexpectedpayoff);


							Logger.logit("instant contraction :  \n");
							for(TargetNode n: domindatednodes)
							{
								Logger.logit(" dominated node "+ n.getTargetid()+"\n");
							}
							Logger.logit("expected payoff "+ defexpectedpayoff+ "\n" + "threshold "+ threshold+"\n");



							sumsol += defexp;
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defexpectedpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/

						}

					} 
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

					//	writeInFile(iter,(int) threshold, targetsize/(iter+1), defexp, contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));

					//writeInFile(iter,(int) threshold, targetsize/(iter+1), sumsol/(iter+1), contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));
					//writeInFile(iter,(int) threshold, targetsize/(ITER), sumsol/(ITER), contractiontime/(ITER), solvingtime/(ITER), revmaptime/(ITER));


				}
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				//String x = df.format(revtime);
				writeInFile(targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, (long)revtime);
			}
		}


		//System.out.println("Final Targets size "+ targets.size());


	}


	
	


	public static void instantContractionTestAPSP(double[][] density, int ITER, int nrow, int ncol, 
			int[] percentages, int[] thresholds, double dmax, int nRes)
	{



		//int nRes = 2;
		//double perc = .20;
		int nTargets = nrow*ncol;
		//double thresholds[]={2,3,4,5};
		//int ITER=4;

		//double percentages[] = {10};
		double[] result = new double[percentages.length];
		int rindex=0;
		double defexp=0;

		//int lstart=1,  lend=4, hstart=8, hend=10;
		//double[][] density=generateRandomDensity( percentages[0], ITER, lstart, lend,  hstart, hend, nTargets);


		for(double perc: percentages)
		{

			for(double threshold: thresholds)
			{
				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				int targetsize=0;


				for(int iter=0; iter<ITER; iter++)
				{

					targets.clear();

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

					//makeZeroSum(gamedata,nTargets);
					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);

					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
					//printEdges(targets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					selectDominatedTargets(targets, domindatednodes, threshold);
					//printNodesWithNeighborsAndPath(domindatednodes, targets);

					Date start = new Date();
					long l1 = start.getTime();


					Date tstart = new Date();
					long tl1 = tstart.getTime();


					//instantContraction(domindatednodes, targets, dmax);
					instantContractionWithAPSP(domindatednodes, targets, dmax);


					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;

					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					targetsize+= targets.size();
					int[][] p = new int[targets.size()][];
					try 
					{
						start = new Date();
						l1 = start.getTime();
						ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
						ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
						/**
						 * map has present id
						 * mapback gives the original ids
						 */
						HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
						HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
						makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
						System.out.println("Total path with duplicates "+pathseq.size());
						//printPaths(pathseq);
						pathseq = removeDuplicatePathSimple(pathseq);
						//printPaths(pathseq);
						System.out.println("Total path without duplicates "+pathseq.size());

						//int k = 2;
						Integer[] input = new Integer[pathseq.size()];
						int[] branch = new int[nRes];//{0,0};//new char[k];


						for(int i=0; i<input.length; i++)
						{
							input[i] = i;
						}
						HashSet jSet=new HashSet();
						if(pathseq.size()==0)
						{
							//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
							//choose the worst payoff for defender

							Double mAxpayoff = Double.MIN_VALUE;
							Double defpayoff = 0.0;
							for(int i=0; i<domindatednodes.size(); i++)
							{
								targets.add(domindatednodes.get(i));
							}
							for(TargetNode x: targets)
							{
								if(x.attackerreward>mAxpayoff)
								{
									mAxpayoff= x.attackerreward;
									defpayoff = x.defenderpenalty;
								}
							}

							sumsol += defpayoff;
							System.out.println("Defender expected payoff "+ defpayoff);
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;


						}
						else
						{
							//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
							if(pathseq.size()<nRes)
							{

								branch = new int[pathseq.size()];
								jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
							}
							else
							{
								jSet=combine(input, nRes, 0, branch, 0, jSet);
							}

							List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
							/**
							 * columns will be combination of paths for each resources. 
							 */
							/**
							 * pmat, where columns will be combination of paths. 
							 * rows are targets. 
							 * each entry will say whether the target is in the joint schedule
							 */
							//jSet.

							//printJointSchedule(jset);

							//printNodesAsNeighbors(dominatednodes);

							p = makePmat(pathseq, jset, mapback, targets);
							//printPathMat(p);

							/**
							 * remove duplicates from p
							 */
							//removeDuplicatesFromP(p);
							//System.out.println();
							//printPathMat(p);

							//System.out.println("Number of targets after contraction "+ targets.size());
							//System.out.println("mip in... ");
							System.out.println("Iter "+iter+", mip in... ");
							double[] probdistribution = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;
							if(probdistribution.equals(null))
							{
								throw new Exception("Prob null...");
							}
							System.out.println("Iter "+iter+", mip out... ");
							//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


							start = new Date();
							l1 = start.getTime();
							int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
							//removeDuplicatesFromP(origpmat);
							//printPathMat(origpmat);
							//System.out.println("\n after mapping back");
							//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);


							Date tstop = new Date();
							long tl2 = tstop.getTime();
							long tdiff = tl2 - tl1;

							sumtime += tdiff;




							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							if(threshold>0)
							{
								revmaptime += diff;
							}



							int maxtargetforattacker = findAttackTarget(origpmat, probdistribution, gamedata);

							double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);
							defexp = defexpectedpayoff;
							System.out.println("Attacked target is "+ maxtargetforattacker);
							System.out.println("Defender expected payoff "+ defexpectedpayoff);


							Logger.logit("instant contraction :  \n");
							for(TargetNode n: domindatednodes)
							{
								Logger.logit(" dominated node "+ n.getTargetid()+"\n");
							}
							Logger.logit("expected payoff "+ defexpectedpayoff+ "\n" + "threshold "+ threshold+"\n");



							sumsol += defexp;
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defexpectedpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/

						}

					} 
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

					//writeInFile(iter,(int) threshold, targetsize/(iter+1), defexp, contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));

					//writeInFile(iter,(int) threshold, targetsize/(iter+1), sumsol/(iter+1), contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));
					//writeInFile(iter,(int) threshold, targetsize/(ITER), sumsol/(ITER), contractiontime/(ITER), solvingtime/(ITER), revmaptime/(ITER));


				}
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				//String x = df.format(revtime);
				writeInFile( targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}
		}


		//System.out.println("Final Targets size "+ targets.size());


	}





	public static void instantContractionTest(double[][] density, int ITER, int nrow, int ncol, 
			int[] percentages, int[] thresholds, double dmax, int nRes)
	{



		//int nRes = 2;
		//double perc = .20;
		int nTargets = nrow*ncol;
		//double thresholds[]={2,3,4,5};
		//int ITER=4;

		//double percentages[] = {10};
		double[] result = new double[percentages.length];
		int rindex=0;
		double defexp=0;

		//int lstart=1,  lend=4, hstart=8, hend=10;
		//double[][] density=generateRandomDensity( percentages[0], ITER, lstart, lend,  hstart, hend, nTargets);


		for(double perc: percentages)
		{

			for(double threshold: thresholds)
			{
				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				int targetsize=0;


				for(int iter=0; iter<ITER; iter++)
				{

					targets.clear();

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

					//makeZeroSum(gamedata,nTargets);
					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);

					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
					//printEdges(targets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					selectDominatedTargets(targets, domindatednodes, threshold);
					//printNodesWithNeighborsAndPath(domindatednodes, targets);

					Date start = new Date();
					long l1 = start.getTime();


					Date tstart = new Date();
					long tl1 = tstart.getTime();


					instantContraction(domindatednodes, targets, dmax);
					//instantContractionWithAPSP(domindatednodes, targets, dmax);


					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;

					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					targetsize+= targets.size();
					int[][] p = new int[targets.size()][];
					try 
					{
						start = new Date();
						l1 = start.getTime();
						ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
						ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
						/**
						 * map has present id
						 * mapback gives the original ids
						 */
						HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
						HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
						makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
						System.out.println("Total path with duplicates "+pathseq.size());
						//printPaths(pathseq);
						pathseq = removeDuplicatePathSimple(pathseq);
						//printPaths(pathseq);
						System.out.println("Total path without duplicates "+pathseq.size());

						//int k = 2;
						Integer[] input = new Integer[pathseq.size()];
						int[] branch = new int[nRes];//{0,0};//new char[k];


						for(int i=0; i<input.length; i++)
						{
							input[i] = i;
						}
						HashSet jSet=new HashSet();
						if(pathseq.size()==0)
						{
							//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
							//choose the worst payoff for defender

							Double mAxpayoff = Double.MIN_VALUE;
							Double defpayoff = 0.0;
							for(int i=0; i<domindatednodes.size(); i++)
							{
								targets.add(domindatednodes.get(i));
							}
							for(TargetNode x: targets)
							{
								if(x.attackerreward>mAxpayoff)
								{
									mAxpayoff= x.attackerreward;
									defpayoff = x.defenderpenalty;
								}
							}

							sumsol += defpayoff;
							System.out.println("Defender expected payoff "+ defpayoff);
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;


						}
						else
						{
							//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
							if(pathseq.size()<nRes)
							{

								branch = new int[pathseq.size()];
								jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
							}
							else
							{
								jSet=combine(input, nRes, 0, branch, 0, jSet);
							}

							List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
							/**
							 * columns will be combination of paths for each resources. 
							 */
							/**
							 * pmat, where columns will be combination of paths. 
							 * rows are targets. 
							 * each entry will say whether the target is in the joint schedule
							 */
							//jSet.

							//printJointSchedule(jset);

							//printNodesAsNeighbors(dominatednodes);

							p = makePmat(pathseq, jset, mapback, targets);
							//printPathMat(p);

							/**
							 * remove duplicates from p
							 */
							//removeDuplicatesFromP(p);
							//System.out.println();
							//printPathMat(p);

							//System.out.println("Number of targets after contraction "+ targets.size());
							//System.out.println("mip in... ");
							System.out.println("Iter "+iter+", mip in... ");
							double[] probdistribution = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;
							if(probdistribution.equals(null))
							{
								throw new Exception("Prob null...");
							}
							System.out.println("Iter "+iter+", mip out... ");
							//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


							start = new Date();
							l1 = start.getTime();
							int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
							//removeDuplicatesFromP(origpmat);
							//printPathMat(origpmat);
							//System.out.println("\n after mapping back");
							//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);


							Date tstop = new Date();
							long tl2 = tstop.getTime();
							long tdiff = tl2 - tl1;

							sumtime += tdiff;




							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							if(threshold>0)
							{
								revmaptime += diff;
							}



							int maxtargetforattacker = findAttackTarget(origpmat, probdistribution, gamedata);

							double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);
							defexp = defexpectedpayoff;
							System.out.println("Attacked target is "+ maxtargetforattacker);
							System.out.println("Defender expected payoff "+ defexpectedpayoff);


							Logger.logit("instant contraction :  \n");
							for(TargetNode n: domindatednodes)
							{
								Logger.logit(" dominated node "+ n.getTargetid()+"\n");
							}
							Logger.logit("expected payoff "+ defexpectedpayoff+ "\n" + "threshold "+ threshold+"\n");



							sumsol += defexp;
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defexpectedpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/

						}

					} 
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

					//writeInFile(iter,(int) threshold, targetsize/(iter+1), defexp, contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));

					//writeInFile(iter,(int) threshold, targetsize/(iter+1), sumsol/(iter+1), contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));
					//writeInFile(iter,(int) threshold, targetsize/(ITER), sumsol/(ITER), contractiontime/(ITER), solvingtime/(ITER), revmaptime/(ITER));


				}
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				//String x = df.format(revtime);
				writeInFile( targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}
		}


		//System.out.println("Final Targets size "+ targets.size());


	}
	
	

	public static double[][] generateRandomDensity(double d, int iTER, int lowstart,
			int lowlimit, int highstart, int highlimit, int nTargets, boolean randomunif)
	{

		//return null;

		double[][]density = new double[iTER][nTargets];

		if(randomunif)
		{
			for(int i=0; i< iTER; i++)
			{

				for(int j=0; j<nTargets; j++)
				{
					density[i][j]=randInt(lowstart, highlimit);;
				}
			}
			return density;
		}


		int n = (int)Math.ceil(nTargets*(d/100.00));
		int l=0;
		//ArrayList<Integer> added = new ArrayList<Integer>();
		int dominatedcount = 0;
		Random rand = new Random(targets.size());
		int settozero = nTargets - n;


		for(int i=0; i< iTER; i++)
		{

			for(int j=0; j<nTargets; j++)
			{
				density[i][j]=-1;
			}



			ArrayList<Integer> added = new ArrayList<Integer>();

			dominatedcount = 0;





			while(true)
			{
				//int randomNum = rand.nextInt((targets.size() - 1-1) + 1) + 1;
				int randomNum = rand1.nextInt((nTargets - 1-1) + 1) + 1;

				int lowutility =  randInt(lowstart, lowlimit);

				if(!added.contains(randomNum) && randomNum != 0)
				{

					density[i][randomNum]=lowutility;
					added.add(randomNum);
					dominatedcount++;

					if(dominatedcount==settozero)
						break;
				}

			}
			int k=0;

			for(int j=0; j<nTargets; j++)
			{

				if(j==0)
				{
					density[i][j]=highlimit;
				}
				else if(density[i][j]==-1)
				{
					int highutility =  randInt(highstart, highlimit);
					density[i][j]=highutility;

				}
			}



		}

		return density;


	}


	
	



	public static double[][] generateRandomDensityV2(int ncat, int iTER, int[][] ranges, int nTargets, int[] targetsincat)
	{

		//return null;


		double[][]density = new double[iTER][nTargets];

		//int n = (int)Math.ceil(nTargets*(d/100.00));
		int l=0;
		//ArrayList<Integer> added = new ArrayList<Integer>();
		//int dominatedcount = 0;
		//Random rand = new Random(targets.size());
		//int settozero = nTargets - n;


		for(int i=0; i< iTER; i++)
		{
			ArrayList<Integer> notadded = new ArrayList<Integer>();

			for(int j=0; j<nTargets; j++)
			{
				notadded.add(j);
			}

			for(int c=0; c<ncat; c++)
			{
				int n = targetsincat[c];
				int dominatedcount = 0;
				while(true)
				{
					//int randomNum = rand.nextInt((targets.size() - 1-1) + 1) + 1;

					if(notadded.size()>1)
					{

						int randomNum = rand1.nextInt(notadded.size());

						int lowutility =  randInt(ranges[c][0], ranges[c][1]);

						if( randomNum != 0)
						{

							density[i][notadded.get(randomNum)]=lowutility;
							notadded.remove(randomNum);
							dominatedcount++;

							if(dominatedcount==n)
								break;
						}
					}
					if(notadded.size()==1)
						break;

				}



			}
			density[i][0] = 10;
		}







		/*int k=0;

			for(int j=0; j<nTargets; j++)
			{

				if(j==0)
				{
					density[i][j]=highlimit;
				}
				else if(density[i][j]==-1)
				{
					int highutility =  randInt(highstart, highlimit);
					density[i][j]=highutility;

				}
			}

		 */



		return density;


	}



	private static void instantContractionWithExtreamPrunning(
			ArrayList<TargetNode> domindatednodes,
			ArrayList<TargetNode> targets, double dmax) 
	{


		/*
		 * find neighbors of dominated nodes which are not dominated
		 */
		ArrayList<TargetNode> neighbornodes = findNeighborsOfDominated(domindatednodes, targets);

		TargetNode basenode = getTargetNode(0, targets);

		/**
		 * now try to discover shortest path between pair of neighbor nodes through diminated nodes. 
		 */

		for(int i=0; i<neighbornodes.size()-1; i++)
		{
			for(int j=i+1; j<neighbornodes.size(); j++)
			{
				TargetNode nodei = neighbornodes.get(i);
				TargetNode nodej = neighbornodes.get(j);


				if(!nodei.getNeighbors().contains(nodej)) // i j not adjacent
				{
					ArrayList<Integer> pathnodes = new ArrayList<Integer>();
					double shortestdist = findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);

					/**
					 * find the shortest dis from base to nodei and nodej
					 * 
					 *  then if (di+dj+shortestdist<= dmax) then add edge 
					 */

					//ArrayList<Integer> pathnodesdi = new ArrayList<Integer>();
					double di = nodei.getDistfrombase(); //findShortestPath(basenode, nodei, targets, domindatednodes, pathnodesdi);


					//ArrayList<Integer> pathnodesdj = new ArrayList<Integer>();
					//pathnodesdi.clear();
					double dj = nodej.getDistfrombase();//findShortestPath(basenode, nodej, targets, domindatednodes, pathnodesdi);





					if(pathnodes.size()>0 && ((shortestdist+di+dj) <= dmax))
					{
						System.out.println("Pathnodes from node "+ nodei.getTargetid() + " to node "+ nodej.getTargetid()+ ", dist "+ shortestdist);
						System.out.print(neighbornodes.get(i).getTargetid()+"->");
						ArrayList<TargetNode> path = new ArrayList<TargetNode>();
						/*for(Integer x: pathnodes)
						{
							System.out.print(x+"->");
							path.add(getTargetNode(x, targets)); // path nodes from i to j
						}*/

						for(int k=0; k<pathnodes.size(); k++)
						{
							path.add(getTargetNode(pathnodes.get(pathnodes.size()-k-1),targets));
							System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
						}
						System.out.println(neighbornodes.get(j).getTargetid());
						/**
						 * add nodes as neighbor
						 * 
						 * path from node i to j is: path
						 * 
						 * add path from i to j 
						 * 
						 * also reverse the path and add the path from j to i
						 * 
						 */

						nodei.addNeighbor(nodej);
						nodei.addDistance(nodej, shortestdist);

						nodei.setPath(nodej, path);
						/**
						 * compute path utility
						 */
						double pathutility = computePathUtility(nodei, path, nodej);
						nodei.setPathUtility(nodej, pathutility);

						//printNeighbors(dest);
						//printDistances(dest);
						nodej.addNeighbor(nodei);
						nodej.addDistance(nodei, shortestdist);
						nodej.setPathUtility(nodei, pathutility);
						//printNeighbors(dest);
						//printDistances(dest);

						/**
						 * rev the path to dest, add to dest
						 */
						ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
						for(int k=0; k<path.size(); k++)
						{
							pathtosrc.add(path.get(path.size()-k-1));
						}
						//System.out.print(nodej.getTargetid()+"->");
						//ArrayList<TargetNode> pathto = new ArrayList<TargetNode>();
						/*for(TargetNode x: pathtosrc)
						{
							System.out.print(x.getTargetid()+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}
						System.out.println(nodei.getTargetid());*/
						nodej.setPath(nodei, pathtosrc);

						//System.out.println("dmax ");

					}
				}

			}

		}

	}
	
	
	
	private static void instantContraction(
			ArrayList<TargetNode> domindatednodes,
			ArrayList<TargetNode> targets, double dmax) 
	{


		/*
		 * find neighbors of dominated nodes which are not dominated
		 */
		ArrayList<TargetNode> neighbornodes = findNeighborsOfDominated(domindatednodes, targets);

		/**
		 * now try to discover shortest path between pair of neighbor nodes through diminated nodes. 
		 */

		for(int i=0; i<neighbornodes.size()-1; i++)
		{
			for(int j=i+1; j<neighbornodes.size(); j++)
			{
				TargetNode nodei = neighbornodes.get(i);
				TargetNode nodej = neighbornodes.get(j);


				if(!nodei.getNeighbors().contains(nodej)) // i j not adjacent
				{

					ArrayList<Integer> pathnodes = new ArrayList<Integer>();
					double shortestdist = findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);

					//if(pathnodes.size()>0)
					if(pathnodes.size()>0 && shortestdist<dmax)
					{
						//System.out.println("Pathnodes from node "+ nodei.getTargetid() + " to node "+ nodej.getTargetid()+ ", dist "+ shortestdist);
						//System.out.print(neighbornodes.get(i).getTargetid()+"->");
						ArrayList<TargetNode> path = new ArrayList<TargetNode>();
						/*for(Integer x: pathnodes)
						{
							System.out.print(x+"->");
							path.add(getTargetNode(x, targets)); // path nodes from i to j
						}*/

						for(int k=0; k<pathnodes.size(); k++)
						{
							path.add(getTargetNode(pathnodes.get(pathnodes.size()-k-1),targets));
							//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
						}
						//System.out.println(neighbornodes.get(j).getTargetid());
						/**
						 * add nodes as neighbor
						 * 
						 * path from node i to j is: path
						 * 
						 * add path from i to j 
						 * 
						 * also reverse the path and add the path from j to i
						 * 
						 */

						nodei.addNeighbor(nodej);
						nodei.addDistance(nodej, shortestdist);

						nodei.setPath(nodej, path);
						/**
						 * compute path utility
						 */
						double pathutility = computePathUtility(nodei, path, nodej);
						nodei.setPathUtility(nodej, pathutility);

						//printNeighbors(dest);
						//printDistances(dest);
						nodej.addNeighbor(nodei);
						nodej.addDistance(nodei, shortestdist);
						nodej.setPathUtility(nodei, pathutility);
						//printNeighbors(dest);
						//printDistances(dest);

						/**
						 * rev the path to dest, add to dest
						 */
						ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
						for(int k=0; k<path.size(); k++)
						{
							pathtosrc.add(path.get(path.size()-k-1));
						}
						//System.out.print(nodej.getTargetid()+"->");
						//ArrayList<TargetNode> pathto = new ArrayList<TargetNode>();
						/*for(TargetNode x: pathtosrc)
						{
							System.out.print(x.getTargetid()+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}
						System.out.println(nodei.getTargetid());*/
						nodej.setPath(nodei, pathtosrc);

						//System.out.println("dmax ");

					}
				}
				else if(nodei.getDistance(nodej)>1) // if they are neighbors
				{
					//System.out.println(" Considering neighbors:  node "+ nodei.getTargetid() + " and node  "+ nodej.getTargetid());



					ArrayList<Integer> pathnodes = new ArrayList<Integer>();
					double shortestdist = findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);

					//if(pathnodes.size()>0)
					if(pathnodes.size()>0 && shortestdist<dmax && nodei.getDistance(nodej)>shortestdist)  // sgortest than the current path
					{
						//System.out.println("Pathnodes from node "+ nodei.getTargetid() + " to node "+ nodej.getTargetid()+ ", dist "+ shortestdist);
						//System.out.print(neighbornodes.get(i).getTargetid()+"->");
						ArrayList<TargetNode> path = new ArrayList<TargetNode>();
						/*for(Integer x: pathnodes)
						{
							System.out.print(x+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}*/

						for(int k=0; k<pathnodes.size(); k++)
						{
							path.add(getTargetNode(pathnodes.get(pathnodes.size()-k-1),targets));
							//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
						}
						//System.out.println(neighbornodes.get(j).getTargetid());
						/**
						 * update the current path because they are already neighbor
						 * 
						 * path from node i to j is: path
						 * 
						 * add path from i to j 
						 * 
						 * also reverse the path and add the path from j to i
						 * 
						 */

						//nodei.addNeighbor(nodej);
						nodei.removeDistance(nodej);
						nodei.addDistance(nodej, shortestdist);

						nodei.removePath(nodej);
						nodei.setPath(nodej, path);
						/**
						 * compute path utility
						 */
						double pathutility = computePathUtility(nodei, path, nodej);
						nodei.removePathUtility(nodej);
						nodei.setPathUtility(nodej, pathutility);

						//printNeighbors(dest);
						//printDistances(dest);
						//nodej.addNeighbor(nodei);


						nodej.removeDistance(nodei);
						nodej.addDistance(nodei, shortestdist);

						nodej.removePathUtility(nodei);
						nodej.setPathUtility(nodei, pathutility);
						//printNeighbors(dest);
						//printDistances(dest);

						/**
						 * rev the path to dest, add to dest
						 */
						ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
						//System.out.print(nodej.getTargetid()+"->");
						for(int k=0; k<path.size(); k++)
						{

							pathtosrc.add(path.get(path.size()-k-1));
							//System.out.print(path.get(path.size()-k-1).getTargetid()+"->");
						}
						//System.out.print(nodei.getTargetid()+"->");
						//ArrayList<TargetNode> pathto = new ArrayList<TargetNode>();
						/*for(TargetNode x: pathtosrc)
						{
							System.out.print(x.getTargetid()+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}
						System.out.println(nodei.getTargetid());*/
						nodej.removePath(nodei);
						nodej.setPath(nodei, pathtosrc);

						//System.out.println("dmax ");

					}

				}


			}

		}

	}

	
	

	public static void instantContractionWithAPSP(
			ArrayList<TargetNode> domindatednodes,
			ArrayList<TargetNode> targets, double dmax) 
	{


		/*
		 * find neighbors of dominated nodes which are not dominated
		 */
		ArrayList<TargetNode> neighbornodes = findNeighborsOfDominated(domindatednodes, targets);
		
		
		/**
		 * 1. Create map
		 */
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		
		for(int i=1; i<=targets.size(); i++)
		{
			map.put(targets.get(i-1).getTargetid(), i);
			mapback.put(i, targets.get(i-1).getTargetid());
		}
		
		
		ArrayList<Integer> domint = new ArrayList<Integer>();
		
		
		//System.out.println("Dominated nodes : ");
		
		for(TargetNode dn: domindatednodes)
		{
			domint.add(map.get(dn.getTargetid()));
			//System.out.println(dn.getTargetid()+"->"+map.get(dn.getTargetid()));
		}
		
		
		
		int[][] adjacencymatrix = new int[targets.size()+1][targets.size()+1];
		makeAdjacencyMatrix(adjacencymatrix , targets, targets.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(targets.size());
		int[][] apsp = allPairShortestPath.allPairShortestPathWithPathConstruction(adjacencymatrix, domint);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPMatrixZero(apsp, targets, targets.size(), map, mapback);
		
		
		
		
		
		

		/**
		 * now try to discover shortest path between pair of neighbor nodes through diminated nodes. 
		 */

		for(int i=0; i<neighbornodes.size()-1; i++)
		{
			for(int j=i+1; j<neighbornodes.size(); j++)
			{
				TargetNode nodei = neighbornodes.get(i);
				TargetNode nodej = neighbornodes.get(j);


				if(!nodei.getNeighbors().contains(nodej)) // i j not adjacent
				{

					//ArrayList<Integer> pathnodes = new ArrayList<Integer>();
					//double shortestdist = findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);
					
					int src = nodei.getTargetid();
					int dest = nodej.getTargetid();
					
					ArrayList<Integer> pathnodes = allPairShortestPath.getPath(src, dest, map, mapback);
					double shortestdist = apsp[map.get(src)][map.get(dest)];//findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);
					
					

					//if(pathnodes.size()>0)
					if(pathnodes.size()>0 && shortestdist<dmax)
					{
						//System.out.println("Pathnodes from node "+ nodei.getTargetid() + " to node "+ nodej.getTargetid()+ ", dist "+ shortestdist);
						//System.out.print(neighbornodes.get(i).getTargetid()+"->");
						ArrayList<TargetNode> path = new ArrayList<TargetNode>();
						/*for(Integer x: pathnodes)
						{
							System.out.print(x+"->");
							path.add(getTargetNode(x, targets)); // path nodes from i to j
						}*/

						for(int k=0; k<pathnodes.size(); k++)
						{
							path.add(getTargetNode(pathnodes.get(pathnodes.size()-k-1),targets));
							//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
						}
						//System.out.println(neighbornodes.get(j).getTargetid());
						/**
						 * add nodes as neighbor
						 * 
						 * path from node i to j is: path
						 * 
						 * add path from i to j 
						 * 
						 * also reverse the path and add the path from j to i
						 * 
						 */

						nodei.addNeighbor(nodej);
						nodei.addDistance(nodej, shortestdist);

						nodei.setPath(nodej, path);
						/**
						 * compute path utility
						 */
						double pathutility = computePathUtility(nodei, path, nodej);
						nodei.setPathUtility(nodej, pathutility);

						//printNeighbors(dest);
						//printDistances(dest);
						nodej.addNeighbor(nodei);
						nodej.addDistance(nodei, shortestdist);
						nodej.setPathUtility(nodei, pathutility);
						//printNeighbors(dest);
						//printDistances(dest);

						/**
						 * rev the path to dest, add to dest
						 */
						ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
						for(int k=0; k<path.size(); k++)
						{
							pathtosrc.add(path.get(path.size()-k-1));
						}
						//System.out.print(nodej.getTargetid()+"->");
						//ArrayList<TargetNode> pathto = new ArrayList<TargetNode>();
						/*for(TargetNode x: pathtosrc)
						{
							System.out.print(x.getTargetid()+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}
						System.out.println(nodei.getTargetid());*/
						nodej.setPath(nodei, pathtosrc);

						//System.out.println("dmax ");

					}
				}
				else if(/*(nodei.getDistance(nodej)>1) */ (nodei.getPath(nodej).size()>0)) // if they are neighbors
				{
					//System.out.println(" Considering neighbors:  node "+ nodei.getTargetid() + " and node  "+ nodej.getTargetid());



					//ArrayList<Integer> pathnodes = new ArrayList<Integer>();
					//double shortestdist = findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);
					
					int src = nodei.getTargetid();
					int dest = nodej.getTargetid();
					
					ArrayList<Integer> pathnodes = allPairShortestPath.getPath(src, dest, map, mapback);
					double shortestdist = apsp[map.get(src)][map.get(dest)];//findShortestPath(nodei, nodej, targets, domindatednodes, pathnodes);
					
					
					
					

					//if(pathnodes.size()>0)
					if(pathnodes.size()>0 && shortestdist<dmax && nodei.getDistance(nodej)>shortestdist)  // sgortest than the current path
					{
						//System.out.println("Pathnodes from node "+ nodei.getTargetid() + " to node "+ nodej.getTargetid()+ ", dist "+ shortestdist);
						//System.out.print(neighbornodes.get(i).getTargetid()+"->");
						ArrayList<TargetNode> path = new ArrayList<TargetNode>();
						/*for(Integer x: pathnodes)
						{
							System.out.print(x+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}*/

						for(int k=0; k<pathnodes.size(); k++)
						{
							path.add(getTargetNode(pathnodes.get(pathnodes.size()-k-1),targets));
							//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
						}
						//System.out.println(neighbornodes.get(j).getTargetid());
						/**
						 * update the current path because they are already neighbor
						 * 
						 * path from node i to j is: path
						 * 
						 * add path from i to j 
						 * 
						 * also reverse the path and add the path from j to i
						 * 
						 */

						//nodei.addNeighbor(nodej);
						nodei.removeDistance(nodej);
						nodei.addDistance(nodej, shortestdist);

						nodei.removePath(nodej);
						nodei.setPath(nodej, path);
						/**
						 * compute path utility
						 */
						double pathutility = computePathUtility(nodei, path, nodej);
						nodei.removePathUtility(nodej);
						nodei.setPathUtility(nodej, pathutility);

						//printNeighbors(dest);
						//printDistances(dest);
						//nodej.addNeighbor(nodei);


						nodej.removeDistance(nodei);
						nodej.addDistance(nodei, shortestdist);

						nodej.removePathUtility(nodei);
						nodej.setPathUtility(nodei, pathutility);
						//printNeighbors(dest);
						//printDistances(dest);

						/**
						 * rev the path to dest, add to dest
						 */
						ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
						//System.out.print(nodej.getTargetid()+"->");
						for(int k=0; k<path.size(); k++)
						{

							pathtosrc.add(path.get(path.size()-k-1));
							//System.out.print(path.get(path.size()-k-1).getTargetid()+"->");
						}
						//System.out.print(nodei.getTargetid()+"->");
						//ArrayList<TargetNode> pathto = new ArrayList<TargetNode>();
						/*for(TargetNode x: pathtosrc)
						{
							System.out.print(x.getTargetid()+"->");
							//path.add(getTargetNode(x, targets)); // path nodes from i to j
						}
						System.out.println(nodei.getTargetid());*/
						nodej.removePath(nodei);
						nodej.setPath(nodei, pathtosrc);

						//System.out.println("dmax ");

					}

				}


			}

		}

	}


	private static double findShortestPathThrougGraph(TargetNode src,
			TargetNode dest, ArrayList<TargetNode> targets,
			ArrayList<Integer> pathnodes) {


		//ArrayList<Integer> pathnodes= new ArrayList<Integer>();
		/**
		 * we have to build a new tree
		 */

		TargetNode start = new TargetNode(src);
		double shortestdistancecoveredyet = 0.0;

		/**
		 * only tree nodes will be in the lists
		 */
		ArrayList<TargetNode> potentialnodes = new ArrayList<TargetNode>(); 
		ArrayList<TargetNode> closednodes = new ArrayList<TargetNode>();
		potentialnodes.add(start);
		//System.out.println("Start node "+ src.getTargetid());
		TargetNode goal = new TargetNode();

		while( !potentialnodes.isEmpty() )
		{
			//get the node with shortest neighbor distance
			TargetNode tmpnodetoexpand = getNearestNode(potentialnodes);
			//System.out.println("Popping node "+ tmpnodetoexpand.getTargetid()+", distfromstart "+ tmpnodetoexpand.getDistfromstart());
			/**
			 * get the node from the original grid. becasue we need to know the neighbors
			 */
			TargetNode orignode = getTargetNode(tmpnodetoexpand.getTargetid(), targets);
			shortestdistancecoveredyet = tmpnodetoexpand.getDistfromstart();
			//System.out.println("Shortest distance covered yet "+ shortestdistancecoveredyet);
			//System.out.println("Remving tmpnode "+ tmpnodetoexpand.getTargetid());
			//System.out.println("queue size before removing "+ potentialnodes.size());
			potentialnodes.remove(tmpnodetoexpand);
			//System.out.println("queue size after removing "+ potentialnodes.size());
			closednodes.add(tmpnodetoexpand);
			if(tmpnodetoexpand.getTargetid()==dest.getTargetid())
			{
				//System.out.println("Dest node "+ dest.getTargetid()+" found, distance from node "+ src.getTargetid()+ " is "+ tmpnodetoexpand.getDistfromstart());
				goal=tmpnodetoexpand;
				instantcontractionshortestdist= tmpnodetoexpand.getDistfromstart();
				shortestdistancecoveredyet= tmpnodetoexpand.getDistfromstart();
				break;

			}

			//System.out.println("Expanding node "+ orignode.getTargetid());
			for(TargetNode successor: orignode.getNeighbors())
			{
				//if( successor.getTargetid()==dest.getTargetid())
				{
					TargetNode suc = new TargetNode(successor);
					if( (!ifVisited(closednodes,suc)))
					{
						//System.out.println("adding sucessor node "+ suc.getTargetid() + " to "+ orignode.getTargetid());

						suc.setDistfromstart(tmpnodetoexpand.getDistfromstart() + orignode.getDistance(successor));
						//System.out.println("setting disst from start " + suc.getDistfromstart());
						suc.setParent(tmpnodetoexpand);
						//System.out.println("setting parent as " + suc.getParent().getTargetid() + " \n\n");
						potentialnodes.add(suc);
					}
				}
			}
		}
		if(goal.parent==null)
			return -1;


		//ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
		TargetNode tmpgoal = goal;
		TargetNode tmpstart = goal.parent;

		//pathnodes.add(tmpgoal.getTargetid());
		while(tmpstart.parent!=null)
		{

			pathnodes.add(tmpstart.getTargetid());
			tmpstart = tmpstart.parent;
			if(tmpstart.getTargetid()==src.getTargetid())
			{
				//tmppathseq.add(tmpgoal.getTargetid());
				break;
			}

		}


		return shortestdistancecoveredyet;

	}




	private static double findShortestPathThrougGraphWDlimit(TargetNode src,
			TargetNode dest, ArrayList<TargetNode> targets,
			ArrayList<Integer> pathnodes, double dmax) {


		//ArrayList<Integer> pathnodes= new ArrayList<Integer>();
		/**
		 * we have to build a new tree
		 */

		TargetNode start = new TargetNode(src);
		double shortestdistancecoveredyet = 0.0;

		/**
		 * only tree nodes will be in the lists
		 */
		ArrayList<TargetNode> potentialnodes = new ArrayList<TargetNode>(); 
		ArrayList<TargetNode> closednodes = new ArrayList<TargetNode>();
		potentialnodes.add(start);
		//System.out.println("Start node "+ src.getTargetid());
		TargetNode goal = new TargetNode();

		while( !potentialnodes.isEmpty() )
		{
			//get the node with shortest neighbor distance
			TargetNode tmpnodetoexpand = getNearestNode(potentialnodes);
			//System.out.println("Popping node "+ tmpnodetoexpand.getTargetid()+", distfromstart "+ tmpnodetoexpand.getDistfromstart());
			/**
			 * get the node from the original grid. becasue we need to know the neighbors
			 */
			TargetNode orignode = getTargetNode(tmpnodetoexpand.getTargetid(), targets);
			shortestdistancecoveredyet = tmpnodetoexpand.getDistfromstart();
			//System.out.println("Shortest distance covered yet "+ shortestdistancecoveredyet);
			//System.out.println("Remving tmpnode "+ tmpnodetoexpand.getTargetid());
			//System.out.println("queue size before removing "+ potentialnodes.size());
			potentialnodes.remove(tmpnodetoexpand);
			//System.out.println("queue size after removing "+ potentialnodes.size());
			closednodes.add(tmpnodetoexpand);
			if(tmpnodetoexpand.getTargetid()==dest.getTargetid())
			{
				//System.out.println("Dest node "+ dest.getTargetid()+" found, distance from node "+ src.getTargetid()+ " is "+ tmpnodetoexpand.getDistfromstart());
				goal=tmpnodetoexpand;
				instantcontractionshortestdist= tmpnodetoexpand.getDistfromstart();
				shortestdistancecoveredyet= tmpnodetoexpand.getDistfromstart();
				break;

			}

			//System.out.println("Expanding node "+ orignode.getTargetid());
			for(TargetNode successor: orignode.getNeighbors())
			{
				//if( successor.getTargetid()==dest.getTargetid())
				{
					TargetNode suc = new TargetNode(successor);
					if( (!ifVisited(closednodes,suc)))
					{
						//System.out.println("adding sucessor node "+ suc.getTargetid() + " to "+ orignode.getTargetid());

						suc.setDistfromstart(tmpnodetoexpand.getDistfromstart() + orignode.getDistance(successor));
						//System.out.println("setting disst from start " + suc.getDistfromstart());
						suc.setParent(tmpnodetoexpand);
						//System.out.println("setting parent as " + suc.getParent().getTargetid() + " \n\n");
						potentialnodes.add(suc);
					}
				}
			}
		}
		if(goal.parent==null)
			return -1;


		//ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
		TargetNode tmpgoal = goal;
		TargetNode tmpstart = goal.parent;

		//pathnodes.add(tmpgoal.getTargetid());
		while(tmpstart.parent!=null)
		{

			pathnodes.add(tmpstart.getTargetid());
			tmpstart = tmpstart.parent;
			if(tmpstart.getTargetid()==src.getTargetid())
			{
				//tmppathseq.add(tmpgoal.getTargetid());
				break;
			}

		}


		return shortestdistancecoveredyet;

	}


	private static double findShortestPath(TargetNode src,
			TargetNode dest, ArrayList<TargetNode> targets,
			ArrayList<TargetNode> domindatednodes, ArrayList<Integer> pathnodes) {


		//ArrayList<Integer> pathnodes= new ArrayList<Integer>();
		/**
		 * we have to build a new tree
		 */

		TargetNode start = new TargetNode(src);
		double shortestdistancecoveredyet = 0.0;

		/**
		 * only tree nodes will be in the lists
		 */
		ArrayList<TargetNode> potentialnodes = new ArrayList<TargetNode>(); 
		ArrayList<TargetNode> closednodes = new ArrayList<TargetNode>();
		potentialnodes.add(start);
		//System.out.println("Start node "+ src.getTargetid());
		TargetNode goal = new TargetNode();

		while( !potentialnodes.isEmpty() )
		{
			//get the node with shortest neighbor distance
			TargetNode tmpnodetoexpand = getNearestNode(potentialnodes);
			//System.out.println("Popping node "+ tmpnodetoexpand.getTargetid()+", distfromstart "+ tmpnodetoexpand.getDistfromstart());
			/**
			 * get the node from the original grid. becasue we need to know the neighbors
			 */
			TargetNode orignode = getTargetNode(tmpnodetoexpand.getTargetid(), targets);
			shortestdistancecoveredyet = tmpnodetoexpand.getDistfromstart();
			//System.out.println("Shortest distance covered yet "+ shortestdistancecoveredyet);
			//System.out.println("Remving tmpnode "+ tmpnodetoexpand.getTargetid());
			//System.out.println("queue size before removing "+ potentialnodes.size());
			potentialnodes.remove(tmpnodetoexpand);
			//System.out.println("queue size after removing "+ potentialnodes.size());
			closednodes.add(tmpnodetoexpand);
			if(tmpnodetoexpand.getTargetid()==dest.getTargetid())
			{
				//System.out.println("Dest node "+ dest.getTargetid()+" found, distance from node "+ src.getTargetid()+ " is "+ tmpnodetoexpand.getDistfromstart());
				goal=tmpnodetoexpand;
				instantcontractionshortestdist= tmpnodetoexpand.getDistfromstart();
				shortestdistancecoveredyet= tmpnodetoexpand.getDistfromstart();
				break;

			}

			//System.out.println("Expanding node "+ orignode.getTargetid());
			for(TargetNode successor: orignode.getNeighbors())
			{
				if(domindatednodes.contains(successor) || successor.getTargetid()==dest.getTargetid())
				{
					TargetNode suc = new TargetNode(successor);
					if( (!ifVisited(closednodes,suc)))
					{
						//System.out.println("adding sucessor node "+ suc.getTargetid() + " to "+ orignode.getTargetid());

						suc.setDistfromstart(tmpnodetoexpand.getDistfromstart() + orignode.getDistance(successor));
						//System.out.println("setting disst from start " + suc.getDistfromstart());
						suc.setParent(tmpnodetoexpand);
						//System.out.println("setting parent as " + suc.getParent().getTargetid() + " \n\n");
						potentialnodes.add(suc);
					}
				}
			}
		}
		if(goal.parent==null)
			return -1;


		//ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
		TargetNode tmpgoal = goal;
		TargetNode tmpstart = goal.parent;

		//pathnodes.add(tmpgoal.getTargetid());
		while(tmpstart.parent!=null)
		{

			pathnodes.add(tmpstart.getTargetid());
			tmpstart = tmpstart.parent;
			if(tmpstart.getTargetid()==src.getTargetid())
			{
				//tmppathseq.add(tmpgoal.getTargetid());
				break;
			}

		}


		return shortestdistancecoveredyet;

	}

	private static ArrayList<TargetNode> findNeighborsOfDominated(
			ArrayList<TargetNode> domindatednodes,
			ArrayList<TargetNode> targets) {

		ArrayList<TargetNode> nn= new ArrayList<TargetNode>();
		for(TargetNode t: domindatednodes)
		{
			for(TargetNode n: t.getNeighbors())
			{
				if(!domindatednodes.contains(n) && !nn.contains(n))
				{
					nn.add(n);
					//System.out.println("Node "+n.getTargetid()+" is added to neighbor list");
				}
			}
		}


		return nn;
	}


	public static double[] basicInstantAbstraction(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;

		int finalsize=-1;



		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus instant :  \n");
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();

		long contracttime = 0;
		long solvtime = 0;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();



			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);




			/*System.out.print("\n original targets \n ");

			for(TargetNode s: targets)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.print("\n dupli targets \n ");

			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			/*System.out.print("\n dom targets \n ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date start = new Date();
			long l1 = start.getTime();

			instantContraction(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contracttime += diff;






			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			/*System.out.print("\n dupli targets after contraction \n ");
			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);

			//printPaths(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();





				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes );


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);
				double defp = expectedDefenderPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);



				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution,originalmap);



					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus instant expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					/*for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");*/
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);

				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}

				//sgc.contractGraph(domindatednodes, targets,dmax);
				instantContraction(domindatednodes, targets,dmax);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contracttime, solvtime, finalsize};





	}
	
	


	public static double[] basicGCSingleInstantAllPathsLP(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;

		int finalsize=-1;



		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus instant :  \n");
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();

		long contracttime = 0;
		long solvtime = 0;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();



			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);




			/*System.out.print("\n original targets \n ");

			for(TargetNode s: targets)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.print("\n dupli targets \n ");

			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			/*System.out.print("\n dom targets \n ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date start = new Date();
			long l1 = start.getTime();

			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contracttime += diff;






			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			/*System.out.print("\n dupli targets after contraction \n ");
			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			
			/*ArrayList<Integer> cur = new ArrayList<Integer>();
			for(TargetNode x: tmpgraph)
			{
				cur.add(x.getTargetid());
			}
			
			ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, cur, nRes);
			
			*/
			
			
			
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			
			
			/*pathseq = buildGreedyPathMultRes(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			
			
			int icount =0;
			for(int i=0; i<targets.size(); i++)
			{

				map.put(targets.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;

			}*/
			
			
			
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);

			//printPaths(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();





				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,astrategy );


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;


				//findAttackTarget(p, probdistribution, gamedata);

				int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
				attackedtarget = mapback.get(attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				//double defp = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);



				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution,originalmap);



					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus instant expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					/*for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");*/
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);

				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}

				//sgc.contractGraph(domindatednodes, targets,dmax);
				instantContractionWithAPSP(domindatednodes, targets,dmax);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contracttime, solvtime, finalsize};





	}
	
	

	public static double[] basicGCMultiInstantAllPathsLP(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;

		int finalsize=-1;



		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus instant :  \n");
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();

		long contracttime = 0;
		long solvtime = 0;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();



			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);




			/*System.out.print("\n original targets \n ");

			for(TargetNode s: targets)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.print("\n dupli targets \n ");

			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			/*System.out.print("\n dom targets \n ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date start = new Date();
			long l1 = start.getTime();

			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contracttime += diff;






			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			/*System.out.print("\n dupli targets after contraction \n ");
			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			
			/*ArrayList<Integer> cur = new ArrayList<Integer>();
			for(TargetNode x: tmpgraph)
			{
				cur.add(x.getTargetid());
			}
			
			ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, cur, nRes);
			
			*/
			
			
			
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			
			
			/*pathseq = buildGreedyPathMultRes(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			
			
			int icount =0;
			for(int i=0; i<targets.size(); i++)
			{

				map.put(targets.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;

			}*/
			
			
			
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);

			//printPaths(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();





				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,astrategy );


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;


				//findAttackTarget(p, probdistribution, gamedata);

				int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
				attackedtarget = mapback.get(attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				//double defp = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);



				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution,originalmap);



					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus instant expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					/*for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");*/
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);

				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}

				//sgc.contractGraph(domindatednodes, targets,dmax);
				instantContractionWithAPSP(domindatednodes, targets,dmax);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contracttime, solvtime, finalsize};





	}

	
	


	public static double[] basicAPSPGCMultiInstantAllPathsLP(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;

		int finalsize=-1;



		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus instant :  \n");
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();

		long contracttime = 0;
		long solvtime = 0;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();



			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);




			/*System.out.print("\n original targets \n ");

			for(TargetNode s: targets)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.print("\n dupli targets \n ");

			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			/*System.out.print("\n dom targets \n ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date start = new Date();
			long l1 = start.getTime();

			//instantContraction(domindatednodes, tmpgraph, dmax);
			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contracttime += diff;






			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			/*System.out.print("\n dupli targets after contraction \n ");
			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			
			/*ArrayList<Integer> cur = new ArrayList<Integer>();
			for(TargetNode x: tmpgraph)
			{
				cur.add(x.getTargetid());
			}
			
			ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, cur, nRes);
			
			*/
			
			
			
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			
			
			/*pathseq = buildGreedyPathMultRes(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			
			
			int icount =0;
			for(int i=0; i<targets.size(); i++)
			{

				map.put(targets.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;

			}*/
			
			
			
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);

			//printPaths(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();





				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,astrategy );


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;


				//findAttackTarget(p, probdistribution, gamedata);

				int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
				attackedtarget = mapback.get(attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				//double defp = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);



				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution,originalmap);



					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus instant expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					/*for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");*/
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);

				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}

				//sgc.contractGraph(domindatednodes, targets,dmax);
				instantContraction(domindatednodes, targets,dmax);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contracttime, solvtime, finalsize};





	}




	public static double[] basicGCMultiInstantGP3LP(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;

		int finalsize=-1;



		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus instant :  \n");
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();

		long contracttime = 0;
		long solvtime = 0;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();



			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);




			/*System.out.print("\n original targets \n ");

			for(TargetNode s: targets)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.print("\n dupli targets \n ");

			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			/*System.out.print("\n dom targets \n ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date start = new Date();
			long l1 = start.getTime();

			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contracttime += diff;






			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			/*System.out.print("\n dupli targets after contraction \n ");
			for(TargetNode s: tmpgraph)
			{
				System.out.print(s.getTargetid()+" ");
			}*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			//ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			
			ArrayList<Integer> cur = new ArrayList<Integer>();
			for(TargetNode x: tmpgraph)
			{
				cur.add(x.getTargetid());
			}
			
			
			
			
			
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			
			pathseq = generatePathsGreedy3(dmax, gamedata, tmpgraph, cur, nRes);
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			
			
			//pathseq = buildGreedyPathMultRes(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			
			
			int icount =0;
			for(int i=0; i<targets.size(); i++)
			{

				map.put(targets.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;

			}
			
			
			
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);

			//printPaths(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();





				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,astrategy );


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;


				//findAttackTarget(p, probdistribution, gamedata);

				int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
				attackedtarget = mapback.get(attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				//double defp = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);



				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution,originalmap);



					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus instant expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					/*for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");*/
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);

				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}

				//sgc.contractGraph(domindatednodes, targets,dmax);
				instantContractionWithAPSP(domindatednodes, targets,dmax);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contracttime, solvtime, finalsize};





	}


	public static double[] basicInstantAbstractionWithExtreamPruning(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;





		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus instant :  \n");
		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);

			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			//sgc.contractGraph(domindatednodes, tmpgraph);
			instantContractionWithExtreamPrunning(domindatednodes, tmpgraph, dmax);



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			//System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");
				double[] probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution,originalmap);

					Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");

					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					System.out.println("Final #of targets "+ targets.size());
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);

				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
				}

				//sgc.contractGraph(domindatednodes, targets);
				instantContractionWithExtreamPrunning(domindatednodes, targets,dmax);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr};





	}




	public static double[] basicSeqAbstraction(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;





		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();
		int itr=1;
		Logger.logit("basic plus seq :  \n");

		long contractime =0;
		long solvtime = 0;
		int finalsize=-1;;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);


			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();
			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove  );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);



			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			System.out.print("Iter "+itr+ " dominated nodes \n");

			for(TargetNode n: domindatednodes)
			{
				System.out.print(n.getTargetid()+ " ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			Date start = new Date();
			long l1 = start.getTime();

			sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractime += diff;
			//instantContraction(domindatednodes, tmpgraph);



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			System.out.println("Considering targets in temporary graph "+ tmpgraph.size());

			printtargets(tmpgraph);


			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();





				double[] probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution, originalmap);

					//Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus seq expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");



					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);
				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
				}

				for(TargetNode n: domindatednodes)
				{
					//Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}



				sgc.contractGraph(domindatednodes, targets, dmax);
				//instantContraction(domindatednodes, targets);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contractime, solvtime, finalsize};





	}

	
	public static double[] noContractionNoColumnGeneration(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;





		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();
		int itr=1;
		Logger.logit("basic plus seq :  \n");

		long contractime =0;
		long solvtime = 0;
		int finalsize=-1;;

		//while(true)
		//	{

		//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);


		//ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove  );
		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);



		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		selectDominatedTargets(targets, domindatednodes, threshold);

		System.out.print("Iter "+itr+ " dominated nodes \n");

		for(TargetNode n: domindatednodes)
		{
			System.out.print(n.getTargetid()+ " ");
		}
		System.out.println();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		Date start = new Date();
		long l1 = start.getTime();

		//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractime += diff;
		//instantContraction(domindatednodes, tmpgraph);



		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		System.out.println("Considering targets in graph "+ targets.size());

		printtargets(targets);


		int[][] p = new int[targets.size()][]; // p matrix
		ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		/**
		 * map has present id
		 * mapback gives the original ids
		 */
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
		makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets);
		//makeOriginalMapping(originalmap, originalmapback, targets);
		//printPaths(pathseq);
		System.out.println("Total path with duplicates "+pathseq.size());
		//pathseq = removeDuplicatePathSimple(pathseq);
		System.out.println("Total path without duplicates "+pathseq.size()+"\n");
		//printPaths(pathseq);

		Integer[] input = new Integer[pathseq.size()];
		int[] branch = new int[nRes];//{0,0};//new char[k];


		for(int i=0; i<input.length; i++)
		{
			input[i] = i;
		}
		HashSet jSet=new HashSet();


		if(pathseq.size()==0)
		{
			//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
			//choose the worst payoff for defender

			Double mAxpayoff = Double.MIN_VALUE;
			Double defpayoff = 0.0;
			for(int i=0; i<domindatednodes.size(); i++)
			{
				targets.add(domindatednodes.get(i));
			}
			for(TargetNode x: targets)
			{
				if(x.attackerreward>mAxpayoff)
				{
					mAxpayoff= x.attackerreward;
					defpayoff = x.defenderpenalty;
				}
			}


			//System.out.println("Defender expected payoff "+ defpayoff);
			/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



		}
		else
		{
			//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
			if(pathseq.size()<nRes)
			{

				branch = new int[pathseq.size()];
				jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
			}
			else
			{
				jSet=combine(input, nRes, 0, branch, 0, jSet);
			}

			List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
			/**
			 * columns will be combination of paths for each resources. 
			 */
			/**
			 * pmat, where columns will be combination of paths. 
			 * rows are targets. 
			 * each entry will say whether the target is in the joint schedule
			 */
			//jSet.

			//printJointSchedule(jset);

			//printNodesAsNeighbors(dominatednodes);

			p = makePmat(pathseq, jset, mapback, targets);


			start = new Date();
			l1 = start.getTime();

			HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
			double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, targets, nRes, astrategy );


			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			solvtime += diff;


			int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);

			if(probdistribution.equals(null))
			{
				throw new Exception("Prob null...");
			}



			defexp = expectedPayoffDefWMapping(attackedtarget, p, gamedata, probdistribution, map);

			//Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			System.out.println("Final #of targets "+ targets.size());
			finalsize = targets.size();
			Logger.logit(" basic plus seq expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			/*for (TargetNode n: targets)
				{
					Logger.logit(n.getTargetid()+" ");
				}
				Logger.logit("\n ");
			 */
		}


		return new double[]{defexp, threshold, itr, contractime, solvtime, finalsize};





	}
	
	
	private static void printStrategyPaths(
			double[] probdistribution, List<ArrayList<Integer>> superjset, 

			ArrayList<ArrayList<Integer>> finalpaths) {

		// run bfs for every path

		//ArrayList<SuperTarget> goals = new ArrayList<SuperTarget>();

		//HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		System.out.print("DEfender patrollling path :\n");
		for(ArrayList<Integer> jschedule : superjset)
		{

			if(probdistribution[superjset.indexOf(jschedule)]>0)
			{
				for(Integer pinddex: jschedule)
				{


							for(Integer n: finalpaths.get(pinddex))
							{
								if(finalpaths.get(pinddex).indexOf(n)>0)
									System.out.print((int)n+"->");
								else
									System.out.print((int)n+"->");

							}
						
					
					System.out.print("\n");
					
				}
			}

		}
//return goals;
	}
	
	
	public static double[] baseline(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;





		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();
		int itr=1;
		Logger.logit("basic plus seq :  \n");

		long contractime =0;
		long solvtime = 0;
		int finalsize=-1;;

		//while(true)
		//	{

		//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);


		//ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove  );
		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);



		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		selectDominatedTargets(targets, domindatednodes, threshold);

		System.out.print("Iter "+itr+ " dominated nodes \n");

		for(TargetNode n: domindatednodes)
		{
			System.out.print(n.getTargetid()+ " ");
		}
		System.out.println();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		Date start = new Date();
		long l1 = start.getTime();

		//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractime += diff;
		//instantContraction(domindatednodes, tmpgraph);



		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		System.out.println("Considering targets in graph "+ targets.size());

		printtargets(targets);


		int[][] p = new int[targets.size()][]; // p matrix
		ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		/**
		 * map has present id
		 * mapback gives the original ids
		 */
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
		makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets);
		//makeOriginalMapping(originalmap, originalmapback, targets);
		//printPaths(pathseq);
		System.out.println("Total path with duplicates "+pathseq.size());
		pathseq = removeDuplicatePathSimple(pathseq);
		System.out.println("Total path without duplicates "+pathseq.size()+"\n");
		printPaths(pathseq);

		Integer[] input = new Integer[pathseq.size()];
		int[] branch = new int[nRes];//{0,0};//new char[k];


		for(int i=0; i<input.length; i++)
		{
			input[i] = i;
		}
		HashSet jSet=new HashSet();


		if(pathseq.size()==0)
		{
			//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
			//choose the worst payoff for defender

			Double mAxpayoff = Double.MIN_VALUE;
			Double defpayoff = 0.0;
			for(int i=0; i<domindatednodes.size(); i++)
			{
				targets.add(domindatednodes.get(i));
			}
			for(TargetNode x: targets)
			{
				if(x.attackerreward>mAxpayoff)
				{
					mAxpayoff= x.attackerreward;
					defpayoff = x.defenderpenalty;
				}
			}


			//System.out.println("Defender expected payoff "+ defpayoff);
			/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



		}
		else
		{
			//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
			if(pathseq.size()<nRes)
			{

				branch = new int[pathseq.size()];
				jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
			}
			else
			{
				jSet=combine(input, nRes, 0, branch, 0, jSet);
			}

			List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
			/**
			 * columns will be combination of paths for each resources. 
			 */
			/**
			 * pmat, where columns will be combination of paths. 
			 * rows are targets. 
			 * each entry will say whether the target is in the joint schedule
			 */
			//jSet.

			//printJointSchedule(jset);

			//printNodesAsNeighbors(dominatednodes);

			p = makePmat(pathseq, jset, mapback, targets);


			start = new Date();
			l1 = start.getTime();

			HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
			double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, targets, nRes, astrategy );


			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			solvtime += diff;


			int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);

			if(probdistribution.equals(null))
			{
				throw new Exception("Prob null...");
			}

			printStrategyPaths(probdistribution, jset, pathseq);

			defexp = expectedPayoffDefWMapping(attackedtarget, p, gamedata, probdistribution, map);

			//Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			System.out.println("Final #of targets "+ targets.size());
			finalsize = targets.size();
			//Logger.logit(" basic plus seq expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			/*for (TargetNode n: targets)
				{
					Logger.logit(n.getTargetid()+" ");
				}
				Logger.logit("\n ");
			 */
			
			System.out.println("attacked target "+attackedtarget + " \n Def exp payoff: "+ defexp);
		}


		return new double[]{defexp, threshold, itr, contractime, solvtime, finalsize};





	}

	




	




	public static double[] basicSeqAbstractionLP(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;





		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();
		int itr=1;
		Logger.logit("basic plus seq :  \n");

		long contractime =0;
		long solvtime = 0;
		int finalsize=-1;;

		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);


			ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();
			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove  );
			//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);



			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			System.out.print("Iter "+itr+ " dominated nodes \n");

			for(TargetNode n: domindatednodes)
			{
				System.out.print(n.getTargetid()+ " ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			Date start = new Date();
			long l1 = start.getTime();

			sgc.contractGraph(domindatednodes, tmpgraph, dmax);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractime += diff;
			//instantContraction(domindatednodes, tmpgraph);



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			System.out.println("Considering targets in temporary graph "+ tmpgraph.size());

			//printtargets(tmpgraph);


			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<Integer> currenttargets = new ArrayList<Integer>();
			
			for(TargetNode x: tmpgraph)
			{
				currenttargets.add(x.getTargetid());
			}
			
			
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets , nRes) ;//generatePaths(dmax, gamedata, tmpgraph);
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			//pathseq = generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets , nRes) ;
			//pathseq = buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			
			
			
			/*int icount =0;
			for(int i=0; i<targets.size(); i++)
			{
				map.put(targets.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;
			}*/
			
			
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");

				start = new Date();
				l1 = start.getTime();




				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,astrategy);


				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvtime += diff;

				int attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
				attackedtarget = mapback.get(attackedtarget);



				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution, originalmap);

					//Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					System.out.println("Final #of targets "+ tmpgraph.size());
					finalsize = tmpgraph.size();
					Logger.logit(" basic plus seq expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					for (TargetNode n: targets)
					{
						Logger.logit(n.getTargetid()+" ");
					}
					Logger.logit("\n ");



					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);
				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
				}

				for(TargetNode n: domindatednodes)
				{
					//Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
					permanentdomindatednodes.add(n);
				}



				sgc.contractGraph(domindatednodes, targets, dmax);
				//instantContraction(domindatednodes, targets);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr, contractime, solvtime, finalsize};





	}
	
	
	
	public static double[] solvingSubgraphBaseline(ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps,
			double dmax, int[][] gamedata, int nRes, int base, int dest, int[][] apspmat, 
			HashMap<Integer,Integer> apspmap, HashMap<String,HashMap<ArrayList<Integer>,Double>> subgamepaths/*, double r*/ ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;
		double[] attackedtarget = {-1, -1};
		
		String key = base+","+dest+","+(int)dmax;
		HashMap<ArrayList<Integer>,Double> spaths = new HashMap<ArrayList<Integer>, Double>();
		ArrayList<Integer> spath = new ArrayList<Integer>();


		if((base != dest)  &&  (apspmat[apspmap.get(base)][apspmap.get(dest)]>dmax) )
		{
			
			//System.out.println("@@@@@@@@@@@Invalid d<dmax \n");
			/*
			// target can't be protected
			// invalid subgame
			// attacker gets most
			double max = Double.MIN_VALUE;
			for(TargetNode n: targetmaps.values())
			{
				if(max<n.attackerreward)
				{
					max = n.attackerreward;
				}
			}*/
			//subgamepaths.put(key, value)
			spaths.put(spath,0.0 );
			//subgamepaths.put(Integer.parseInt(key), spaths);
			
			return new double[]{777, -777, -1};
			
		}
		if(dmax==0 && base==dest)
		{
			
			if(targets.size()==1)
			{
				
				spath.add(targets.get(0).getTargetid());
				spaths.put(spath,1.0 );
				subgamepaths.put((key), spaths);
				return new double[]{0, 0 , targets.get(0).getTargetid()};
			}
			//System.out.println("@@@@@@@@@@@only entry can be covered, src=dest, d=0 \n");
			// target can't be protected
			// invalid subgame
			// attacker gets most
			/*double max = 0;
			for(TargetNode n: targetmaps.values())
			{
				if(max<n.attackerreward && n.getTargetid() != base)
				{
					max = n.attackerreward;
				}
			}*/
			//ArrayList<Integer> spath = new ArrayList<Integer>();
			spaths.put(spath,0.0 );
			//subgamepaths.put(Integer.parseInt(key), spaths);
			
			return new double[]{777, -777,-1};
		}



		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		//ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();
		//int itr=1;
		//Logger.logit("basic plus seq :  \n");

		long contractime =0;
		long solvtime = 0;
		int finalsize=-1;;

		//while(true)
		//	{

		//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);


		//ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove  );
		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);



		/*ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		selectDominatedTargets(targets, domindatednodes, threshold);*/

		/*System.out.print("Iter "+itr+ " dominated nodes \n");

		for(TargetNode n: domindatednodes)
		{
			System.out.print(n.getTargetid()+ " ");
		}*/
		//System.out.println();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		Date start = new Date();
		long l1 = start.getTime();

		//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractime += diff;
		//instantContraction(domindatednodes, tmpgraph);



		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		//System.out.println("Considering targets in graph "+ targets.size());

		//printtargets(targets);


		int[][] p = new int[targets.size()][]; // p matrix
		ArrayList<TargetNode> goals = generatePathsWithSrcDest (dmax, gamedata, targets, targetmaps, base, dest);
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		/**
		 * map has present id
		 * mapback gives the original ids
		 */
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
		makePathSeqSrcDest(pathseq, goals, goals.size(), targets.size(), map, mapback, targets);
		//makeOriginalMapping(originalmap, originalmapback, targets);
		printPaths(pathseq);
		//System.out.println("Total path with duplicates "+pathseq.size());
		pathseq = removeDuplicatePathSimple(pathseq);
		//System.out.println("Total path without duplicates "+pathseq.size()+"\n");
		printPaths(pathseq);

		Integer[] input = new Integer[pathseq.size()];
		int[] branch = new int[nRes];//{0,0};//new char[k];


		for(int i=0; i<input.length; i++)
		{
			input[i] = i;
		}
		HashSet jSet=new HashSet();


		if(pathseq.size()==0)
		{
			
			
			//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
			//choose the worst payoff for defender

			Double mAxpayoff = Double.MIN_VALUE;
			Double defpayoff = 0.0;
			/*for(int i=0; i<domindatednodes.size(); i++)
			{
				targets.add(domindatednodes.get(i));
			}*/
			for(TargetNode x: targets)
			{
				if(x.attackerreward>mAxpayoff && x.getTargetid()!=base)
				{
					mAxpayoff= x.attackerreward;
					defpayoff = x.defenderpenalty;
					attackedtarget[0] = x.getTargetid();
					attackedtarget[1] = mAxpayoff;
				}
			}
			//throw new Exception("Nothing should be here");
			return new double[]{777, -777, -1 };


			//System.out.println("Defender expected payoff "+ defpayoff);
			/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



		}
		else
		{
			//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
			if(pathseq.size()<nRes)
			{

				branch = new int[pathseq.size()];
				jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
			}
			else
			{
				jSet=combine(input, nRes, 0, branch, 0, jSet);
			}

			List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
			/**
			 * columns will be combination of paths for each resources. 
			 */
			/**
			 * pmat, where columns will be combination of paths. 
			 * rows are targets. 
			 * each entry will say whether the target is in the joint schedule
			 */
			//jSet.

			//printJointSchedule(jset);

			//printNodesAsNeighbors(dominatednodes);

			p = makePmat(pathseq, jset, mapback, targets);


			start = new Date();
			l1 = start.getTime();

			HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
			double[] probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, targets, nRes, astrategy );


			stop = new Date();
			l2 = stop.getTime();
			diff = l2 - l1;

			solvtime += diff;


			attackedtarget = findAttackTargetWHashMap(p, probdistribution, gamedata, map, mapback);

			if(probdistribution.equals(null))
			{
				throw new Exception("Prob null...");
			}



			defexp = expectedPayoffDefWMapping((int)attackedtarget[0], p, gamedata, probdistribution, map);
			
			addSubGamePaths(key, subgamepaths, spaths, probdistribution, jset, pathseq);

			//Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			//System.out.println("Final #of targets "+ targets.size());
			finalsize = targets.size();
			//Logger.logit(" basic plus seq expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			/*for (TargetNode n: targets)
				{
					Logger.logit(n.getTargetid()+" ");
				}
				Logger.logit("\n ");
			 */
		}

		
		return new double[]{defexp, attackedtarget[1], attackedtarget[0]};





	}


	
	
	public static double[] solvingSubgraphBaseline2(ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps,
			double dmax, int[][] gamedata, int nRes, int base, int dest, int[][] apspmat, 
			HashMap<Integer,Integer> apspmap, HashMap<String,HashMap<ArrayList<Integer>,Double>> subgamepaths/*, double r*/, int[] disallocation ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;
		double[] attackedtarget = {-1, -1};
		
		
		HashMap<ArrayList<Integer>,Double> spaths = new HashMap<ArrayList<Integer>, Double>();
		ArrayList<Integer> spath = new ArrayList<Integer>();


		
		if(base==dest && targets.size()>1)
		{
			return new double[]{0, 0 , -1}; 
		}


		if(targets.size()==1)
		{
			String key = base+","+dest+","+0;
			spath.add(targets.get(0).getTargetid());
			spaths.put(spath,1.0 );
			subgamepaths.put((key), spaths);
			disallocation[0] = 0;
			return new double[]{0, 0 , targets.get(0).getTargetid()};
		}
			//System.out.println("@@@@@@@@@@@only entry can be covered, src=dest, d=0 \n");
			// target can't be protected
			// invalid subgame
			// attacker gets most
			/*double max = 0;
			for(TargetNode n: targetmaps.values())
			{
				if(max<n.attackerreward && n.getTargetid() != base)
				{
					max = n.attackerreward;
				}
			}*/
			//ArrayList<Integer> spath = new ArrayList<Integer>();
			
			//subgamepaths.put(Integer.parseInt(key), spaths);
			
			

		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		//ArrayList<TargetNode> permanentdomindatednodes = new ArrayList<TargetNode>();
		//int itr=1;
		//Logger.logit("basic plus seq :  \n");

		long contractime =0;
		long solvtime = 0;
		int finalsize=-1;;

		//while(true)
		//	{

		//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);


		//ArrayList<TargetNode> domindatednodestoremove = new ArrayList<TargetNode>();
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets, permanentdomindatednodes,domindatednodestoremove  );
		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodestoremove, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodestoremove, tmpgraph);



		/*ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		selectDominatedTargets(targets, domindatednodes, threshold);*/

		/*System.out.print("Iter "+itr+ " dominated nodes \n");

		for(TargetNode n: domindatednodes)
		{
			System.out.print(n.getTargetid()+ " ");
		}*/
		//System.out.println();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		Date start = new Date();
		long l1 = start.getTime();

		//sgc.contractGraph(domindatednodes, tmpgraph, dmax);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractime += diff;
		//instantContraction(domindatednodes, tmpgraph);



		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		//System.out.println("Considering targets in graph "+ targets.size());

		//printtargets(targets);
		
		
		
		ArrayList<Integer> onepath = buildOneGreedyPathWithSrcDest(targets, dmax, targets.size(), base, dest, nRes, disallocation);


		//int[][] p = new int[targets.size()][]; // p matrix
		//ArrayList<TargetNode> goals = generateOnePathWithSrcDest (dmax, gamedata, targets, targetmaps, base, dest, disallocation);
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		
		/*if(goals.size()==0)
		{
			return new double[]{0, 0 , -1}; 
		}
		*/
		/**
		 * map has present id
		 * mapback gives the original ids
		 */
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
		//HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
		//makePathSeqSrcDest(pathseq, goals, goals.size(), targets.size(), map, mapback, targets);
		//makeOriginalMapping(originalmap, originalmapback, targets);
		pathseq.add(onepath);
		printPaths(pathseq);
		//System.out.println("Total path with duplicates "+pathseq.size());
		pathseq = removeDuplicatePathSimple(pathseq);
		//System.out.println("Total path without duplicates "+pathseq.size()+"\n");
		printPaths(pathseq);

		String key = base+","+dest+","+ disallocation[0];
		addSubGamePaths(key, subgamepaths, spaths, pathseq);

			//Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			//System.out.println("Final #of targets "+ targets.size());
			finalsize = targets.size();
			//Logger.logit(" basic plus seq expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
			/*for (TargetNode n: targets)
				{
					Logger.logit(n.getTargetid()+" ");
				}
				Logger.logit("\n ");
			 */
			
			// TODO 
			// find defender expected payoff
			
		if(pathseq.size()>1)
		{
			throw new Exception("More than one paths");
		}
		
		
		int maxtarget = -1;
		double maxattackerpayoff = Double.NEGATIVE_INFINITY;
		for(TargetNode t: targetmaps.values())
		{
			// if  covered then 
			double tmppayoff = -1;
			if(pathseq.get(0).contains(t.getTargetid()))
			{
				tmppayoff = t.attackerpenalty;
			}
			else
			{
				tmppayoff = t.attackerreward;
			}
			
			if(maxattackerpayoff < tmppayoff)
			{
				maxattackerpayoff = tmppayoff;
				maxtarget = t.getTargetid();
			}
		}
		

		
		return new double[]{-maxattackerpayoff, maxattackerpayoff, maxtarget};





	}





	private static void addSubGamePaths(String key, HashMap<String, HashMap<ArrayList<Integer>,Double>> subgamepaths,
			HashMap<ArrayList<Integer>,Double> spaths, double[] probdistribution,
			List<ArrayList<Integer>> jset, ArrayList<ArrayList<Integer>> pathseq) {
		
		
		
		for(int probindex=0; probindex<probdistribution.length; probindex++)
		{
			//if(probdistribution[probindex]>0)
			{
				// get the path index from jset
				for(int pathindex: jset.get(probindex))
				{
					// get the path 
					ArrayList<Integer> spath = new ArrayList<Integer>();
					for(int pathnode: pathseq.get(pathindex))
					{
						// add nodes of path
						spath.add(pathnode);
						
					}
					spaths.put(spath,probdistribution[probindex] );
				}
			}
		}
		subgamepaths.put((key), spaths);
		
		
	}
	
	

	private static void addSubGamePaths(String key, HashMap<String, HashMap<ArrayList<Integer>,Double>> subgamepaths,
			HashMap<ArrayList<Integer>,Double> spaths, ArrayList<ArrayList<Integer>> pathseq) {
		
		
		
		
		spaths.put(pathseq.get(0), 1.0);
		subgamepaths.put((key), spaths);
		
		
	}

	public static double[] basicSeqAbstractionExtreamPruning(ArrayList<TargetNode> targets, double threshold,
			double dmax, int[][] gamedata, SecurityGameContraction sgc, int nRes ) throws Exception
	{
		double attackeru=0;
		double attackerv=0;
		double defexp=0;





		// repeat untill u==v
		//ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

		int itr=1;
		Logger.logit("basic plus seq :  \n");
		while(true)
		{
			if(itr==100)
				break;
			//System.out.println("\n######Itr: "+itr+", threshold = "+threshold);
			ArrayList<TargetNode> tmpgraph = getDuplicateGraph(targets);
			ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

			selectDominatedTargets(tmpgraph, domindatednodes, threshold);

			SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			sgc.contractGraphWithExtreamPruning(domindatednodes, tmpgraph, dmax);
			//instantContraction(domindatednodes, tmpgraph);



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

			int[][] p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			makeOriginalMapping(originalmap, originalmapback, targets);
			//printPaths(pathseq);
			//System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];


			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);

				/**
				 * remove duplicates from p
				 */
				//removeDuplicatesFromP(p);
				//System.out.println();
				//printPathMat(p);

				//System.out.println("Number of targets after contraction "+ targets.size());
				//System.out.println("mip in... ");
				double[] probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}
				//System.out.println("mip out... ");
				//printPathWithPositiveCoverage(p, probdistribution, jset, pathseq, map);



				int[][] origpmat = makeOrigPMat(p, pathseq, jset, targets.size(), domindatednodes, originalmap, originalmapback, tmpgraph);
				//removeDuplicatesFromP(origpmat);
				//printPathMat(origpmat);
				//System.out.println("\n after mapping back");
				//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



				int xx= findAttackTargetWMapping(origpmat, probdistribution, gamedata, originalmap, originalmapback);

				int maxtargetforattacker = originalmapback.get(xx); // add map back


				//double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);

				//System.out.println("Attacked target is "+ maxtargetforattacker);
				//int v= getTargetNode(maxtargetforattacker, targets).getTargetid();
				attackerv = expectedAttackerPayoff(maxtargetforattacker, origpmat, probdistribution, gamedata, originalmap);
				//System.out.println("v= "+attackerv);
				//System.out.println("Attacked expected payoff "+ getTargetNode(maxtargetforattacker, targets).attackerreward);
				//System.out.println("Defender expected payoff "+ defexpectedpayoff);

				System.out.println("iter "+itr+", u= "+attackeru +", v= "+attackerv );
				DecimalFormat df = new DecimalFormat("#.####");
				df.setRoundingMode(RoundingMode.CEILING);
				double du=Double.parseDouble(df.format(attackeru));
				double dv=Double.parseDouble(df.format(attackerv));

				if(du==dv)
				{
					/*domindatednodes = new ArrayList<TargetNode>();

					selectDominatedTargets(targets, domindatednodes, attackeru);



					sgc.contractGraphV2(domindatednodes, targets);
					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);*/
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					defexp = expectedPayoffDefWMapping(maxtargetforattacker, origpmat, gamedata, probdistribution, originalmap);

					Logger.logit("expected payoff "+ defexp+ "\n" + "threshold "+ threshold+"\n");
					System.out.println("Final #of targets "+ targets.size());
					break;
				}


				domindatednodes = new ArrayList<TargetNode>();

				selectDominatedTargets(targets, domindatednodes, attackeru);
				for(TargetNode n: domindatednodes)
				{
					Logger.logit("iter "+itr+", dominated node "+ n.getTargetid()+"\n");
				}


				sgc.contractGraphWithExtreamPruning(domindatednodes, targets, dmax);
				//instantContraction(domindatednodes, targets);


				SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
				SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);




				threshold = attackerv;

			}
			//break;

			itr++;

		} // outer while loop

		return new double[]{defexp, threshold, itr};





	}



	private static void setDensity(int[][] gamedata, double perc, int nTargets) {





	}

	private static void makeZeroSum(int[][] gamedata, int n) {

		for(int i=0;i<n; i++)
		{
			gamedata[i][1]=0;
			gamedata[i][3]=0;
			gamedata[i][2]=gamedata[i][0];


		}


	}

	private static void makeOriginalMapping(
			HashMap<Integer, Integer> originalmap,
			HashMap<Integer, Integer> originalmapback,
			ArrayList<TargetNode> targets) {


		int icount =0;
		for(int i=0; i<targets.size(); i++)
		{

			originalmap.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			originalmapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}

	}

	private static void makeStarGraph(int[][] gamedata, int nTargets) {
		targets.clear();







		TargetNode x[] = new TargetNode[nTargets];




		for(int i=0; i<nTargets; i++)
		{
			/*gamedata[i][1]=0;
			gamedata[i][0]=0;
			gamedata[i][3]=0;
			gamedata[i][2]=0;
			 */

			x[i] = new TargetNode(i, gamedata[i][2]);
			x[i].attackerpenalty =gamedata[i][3];
			x[i].attackerreward=gamedata[i][2];
			x[i].defenderpenalty=gamedata[i][1];
			x[i].defenderreward=gamedata[i][0];

			targets.add(x[i]);
			if(i>0)
			{
				x[i].addNeighbor(x[0]);
				x[0].addNeighbor(x[i]);
				ArrayList<TargetNode> p = new ArrayList<TargetNode>();
				x[i].setPath(x[0], p);
				x[0].setPath(x[i], p);
				x[i].setPathUtility(x[0], 0.0);
				x[0].setPathUtility(x[i], 0.0);
				x[i].addDistance(x[0], 1.0);
				x[0].addDistance(x[i], 1.0);
			}

		}




		/*gamedata[0][1]=0;
		gamedata[0][0]=10;
		gamedata[0][3]=0;
		gamedata[0][2]=10;

		gamedata[1][1]=0;
		gamedata[1][0]=10;
		gamedata[1][3]=0;
		gamedata[1][2]=10;

		gamedata[2][1]=0;
		gamedata[2][0]=10;
		gamedata[2][3]=0;
		gamedata[2][2]=10;

		gamedata[3][1]=0;
		gamedata[3][0]=9;
		gamedata[3][3]=0;
		gamedata[3][2]=9;


		gamedata[4][1]=0;
		gamedata[4][0]=0;
		gamedata[4][3]=0;
		gamedata[4][2]=0;

		gamedata[5][1]=0;
		gamedata[5][0]=0;
		gamedata[5][3]=0;
		gamedata[5][2]=0;

		gamedata[6][1]=0;
		gamedata[6][0]=0;
		gamedata[6][3]=0;
		gamedata[6][2]=0;

		gamedata[7][1]=0;
		gamedata[7][0]=0;
		gamedata[7][3]=0;
		gamedata[7][2]=0;




		//for(int i=0; i<n; i++)
		{
			TargetNode t = new TargetNode(0, 10);
			t.attackerpenalty =0;
			t.attackerreward=10;
			t.defenderpenalty=0;
			t.defenderreward=10;
			TargetNode t1 = new TargetNode(1, 10);
			t1.attackerpenalty =0;
			t1.attackerreward=10;
			t1.defenderpenalty=0;
			t1.defenderreward=10;
			TargetNode t2 = new TargetNode(2, 10);
			t2.attackerpenalty =0;
			t2.attackerreward=10;
			t2.defenderpenalty=0;
			t2.defenderreward=10;
			TargetNode t3 = new TargetNode(3, 9);
			t3.attackerpenalty =0;
			t3.attackerreward=9;
			t3.defenderpenalty=0;
			t3.defenderreward=9;
			TargetNode t4 = new TargetNode(4, 0);
			t4.attackerpenalty =0;
			t4.attackerreward=0;
			t4.defenderpenalty=0;
			t4.defenderreward=0;
			TargetNode t5 = new TargetNode(5, 0);
			t5.attackerpenalty =0;
			t5.attackerreward=0;
			t5.defenderpenalty=0;
			t5.defenderreward=0;
			TargetNode t6 = new TargetNode(6, 0);
			t6.attackerpenalty =0;
			t6.attackerreward=0;
			t6.defenderpenalty=0;
			t6.defenderreward=0;

			t1.addNeighbor(t);
			t.addNeighbor(t1);
			ArrayList<TargetNode> p = new ArrayList<TargetNode>();
			t1.setPath(t, p);
			t.setPath(t1, p);
			t1.setPathUtility(t, 0.0);
			t.setPathUtility(t1, 0.0);
			t1.addDistance(t, 1.0);
			t.addDistance(t1, 1.0);



			t2.addNeighbor(t);
			t.addNeighbor(t2);
			p = new ArrayList<TargetNode>();
			t2.setPath(t, p);
			t.setPath(t2, p);
			t2.setPathUtility(t, 0.0);
			t.setPathUtility(t2, 0.0);
			t2.addDistance(t, 1.0);
			t.addDistance(t2, 1.0);




			t3.addNeighbor(t);
			t.addNeighbor(t3);
			p = new ArrayList<TargetNode>();
			t3.setPath(t, p);
			t.setPath(t3, p);
			t3.setPathUtility(t, 0.0);
			t.setPathUtility(t3, 0.0);
			t3.addDistance(t, 1.0);
			t.addDistance(t3, 1.0);



			t4.addNeighbor(t);
			t.addNeighbor(t4);
			p = new ArrayList<TargetNode>();
			t4.setPath(t, p);
			t.setPath(t4, p);
			t4.setPathUtility(t, 0.0);
			t.setPathUtility(t4, 0.0);
			t4.addDistance(t, 1.0);
			t.addDistance(t4, 1.0);



			t5.addNeighbor(t);
			t.addNeighbor(t5);
			p = new ArrayList<TargetNode>();
			t5.setPath(t, p);
			t.setPath(t5, p);
			t5.setPathUtility(t, 0.0);
			t.setPathUtility(t5, 0.0);
			t5.addDistance(t, 1.0);
			t.addDistance(t5, 1.0);





			t6.addNeighbor(t);
			t.addNeighbor(t6);
			t6.setPath(t, p);
			t.setPath(t6, p);
			t6.setPathUtility(t, 0.0);
			t.setPathUtility(t6, 0.0);
			t6.addDistance(t, 1.0);
			t.addDistance(t6, 1.0);



		 */





		/*targets.add(t);
			targets.add(t1);
			targets.add(t2);
			targets.add(t3);
			targets.add(t4);
			targets.add(t5);
			targets.add(t6);*/

		//}



	}

	private static void selectDominatedTargets(ArrayList<TargetNode> tmpgraph,
			ArrayList<TargetNode> domindatednodes, double threshold) {

		for(TargetNode t: tmpgraph)
		{
			System.out.println("Threshld : "+ threshold + ", target "+ t.getTargetid() +", areward: " +t.attackerreward);
			if(t.getAnimaldensity()<threshold && t.getTargetid()!=0)
			{
				System.out.println("Target "+t.getTargetid()+ " is added to the dominated list, u: "+ t.getAnimaldensity());
				domindatednodes.add(t);
			}
		}


	}



	private static ArrayList<TargetNode> getDuplicateGraph(
			ArrayList<TargetNode> src, ArrayList<TargetNode> permanentdomindatednodes, 
			ArrayList<TargetNode> domindatednodestoremove) {


		///
		ArrayList<TargetNode> duplicategraph = new ArrayList<TargetNode>();

		for(TargetNode s: src)
		{
			TargetNode tmp = new TargetNode(s.getTargetid(), s.getAnimaldensity());
			tmp.attackerpenalty = s.attackerpenalty;
			tmp.attackerreward = s.attackerreward;
			tmp.defenderpenalty= s.defenderpenalty;
			tmp.defenderreward = s.defenderreward;
			tmp.setDistfrombase(s.getDistfrombase());
			duplicategraph.add(tmp);
		}


		for(TargetNode s: permanentdomindatednodes)
		{
			TargetNode tmp = new TargetNode(s.getTargetid(), s.getAnimaldensity());
			tmp.attackerpenalty = s.attackerpenalty;
			tmp.attackerreward = s.attackerreward;
			tmp.defenderpenalty= s.defenderpenalty;
			tmp.defenderreward = s.defenderreward;
			//tmp.setDistfrombase(s.getDistfrombase());
			duplicategraph.add(tmp);
			domindatednodestoremove.add(tmp);
		}


		/*System.out.print("\n original targets \n ");

		for(TargetNode s: src)
		{
			System.out.print(s.getTargetid()+" ");
		}*/

		/*System.out.print("\n dominated targets \n ");

		for(TargetNode s: permanentdomindatednodes)
		{
			System.out.print(s.getTargetid()+" ");
		}*/


		/*System.out.print("\n Duplicate targets \n ");

		for(TargetNode s: duplicategraph)
		{
			System.out.print(s.getTargetid()+" ");
		}

		System.out.println();*/







		int tindex = 0;
		for(TargetNode s: src)
		{
			TargetNode t = duplicategraph.get(tindex);
			// add all neighbors

			for(TargetNode nei: s.getNeighbors())
			{
				for(TargetNode t2: duplicategraph)
				{
					if(t.getTargetid()!= t2.getTargetid() && nei.getTargetid()==t2.getTargetid())
					{
						t.addNeighbor(t2);
						//System.out.println("Adding target "+t2.getTargetid()+" as node "+ t.getTargetid()+"'s neighbor");
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//System.out.println("Pathnodes to target "+ t2.getTargetid());



						/*for(TargetNode pnode: s.getPath(nei))
						{
							System.out.print(pnode.getTargetid()+"->");
						}
						System.out.println();
						 */
						int f=1;


						for(TargetNode pnode: s.getPath(nei))
						{
							for(TargetNode t3: duplicategraph)
							{
								if(t3.getTargetid()==pnode.getTargetid())
								{
									/**
									 * modification
									 */
									if(pathnodes.size()>0)
									{
										/**
										 * get the last index node add them as neighbors with new one
										 */
										TargetNode lastnodeyet = pathnodes.get(pathnodes.size()-1);
										lastnodeyet.addNeighbor(t3);

									}
									pathnodes.add(t3);
									break;
								}
							}
						}

						/*System.out.println("recreated Pathnodes to target "+ t2.getTargetid());
						for(TargetNode pnode: pathnodes)
						{
							System.out.print(pnode.getTargetid()+"->");
						}
						System.out.println();
						 */

						t.setPath(t2, pathnodes);
						t.setPathUtility(t2, 0.0);

						t.addDistance(t2, s.getDistance(nei));
						break;
					}
				}
			}

			tindex++;
		}








		return duplicategraph;
	}


	
	
	


	public static ArrayList<TargetNode> getDuplicateGraph(
			ArrayList<TargetNode> src) {


		///
		ArrayList<TargetNode> duplicategraph = new ArrayList<TargetNode>();

		for(TargetNode s: src)
		{
			TargetNode tmp = new TargetNode(s.getTargetid(), s.getAnimaldensity());
			tmp.attackerpenalty = s.attackerpenalty;
			tmp.attackerreward = s.attackerreward;
			tmp.defenderpenalty= s.defenderpenalty;
			tmp.defenderreward = s.defenderreward;
			tmp.setDistfrombase(s.getDistfrombase());
			duplicategraph.add(tmp);
		}




		/*System.out.print("\n original targets \n ");

		for(TargetNode s: src)
		{
			System.out.print(s.getTargetid()+" ");
		}*/


		int tindex = 0;
		for(TargetNode s: src)
		{
			TargetNode t = duplicategraph.get(tindex);
			// add all neighbors

			for(TargetNode nei: s.getNeighbors())
			{
				for(TargetNode t2: duplicategraph)
				{
					if(t.getTargetid()!= t2.getTargetid() && nei.getTargetid()==t2.getTargetid())
					{
						t.addNeighbor(t2);
						//System.out.println("Adding target "+t2.getTargetid()+" as node "+ t.getTargetid()+"'s neighbor");
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//System.out.println("Pathnodes to target "+ t2.getTargetid());



						/*for(TargetNode pnode: s.getPath(nei))
						{
							System.out.print(pnode.getTargetid()+"->");
						}
						System.out.println();*/

						int f=1;


						for(TargetNode pnode: s.getPath(nei))
						{
							for(TargetNode t3: duplicategraph)
							{
								if(t3.getTargetid()==pnode.getTargetid())
								{
									/**
									 * modification
									 */
									if(pathnodes.size()>0)
									{
										/**
										 * get the last index node add them as neighbors with new one
										 */
										TargetNode lastnodeyet = pathnodes.get(pathnodes.size()-1);
										lastnodeyet.addNeighbor(t3);

									}
									pathnodes.add(t3);
									break;
								}
							}
						}

						//System.out.println("recreated Pathnodes to target "+ t2.getTargetid());
						/*for(TargetNode pnode: pathnodes)
						{
							System.out.print(pnode.getTargetid()+"->");
						}
						System.out.println();*/


						t.setPath(t2, pathnodes);
						t.setPathUtility(t2, 0.0);

						t.addDistance(t2, s.getDistance(nei));
						break;
					}
				}
			}

			tindex++;
		}








		return duplicategraph;
	}

	private static double initializeThreshold(ArrayList<Integer> greedypath) {

		double maxutility =Double.NEGATIVE_INFINITY;
		int maxtarget=-1;
		for(TargetNode t: targets )
		{
			//System.out.println("Target "+ t.getTargetid()+ ", utility "+ t.attackerreward);
			if(!greedypath.contains(t.getTargetid()) && t.attackerreward>maxutility)
			{
				maxutility=t.attackerreward;
				maxtarget=t.getTargetid();

			}
			//System.out.println("max utility "+ maxutility+" max Target "+ maxtarget);

		}

		return maxutility;
	}

	private static ArrayList<Integer> greedyFirstRoute(double dmax, int[][] gamedata, ArrayList<TargetNode> targets) throws Exception {


		ArrayList<TargetNode> goals = generatePathsGreedy(dmax, gamedata, targets);
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		/**
		 * map has present id
		 * mapback gives the original ids
		 */
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
		//printPaths(pathseq);
		System.out.println("Total path with duplicates "+pathseq.size());
		pathseq = removeDuplicatePathSimple(pathseq);
		System.out.println("Total path without duplicates "+pathseq.size());

		//printPaths(pathseq);

		double maxutility = Double.NEGATIVE_INFINITY;
		int maxpathindex=0;
		for(ArrayList<Integer> path: pathseq)
		{
			double tmputility = calculatePathUtility(path, gamedata);
			System.out.println("PAth "+ pathseq.indexOf(path) + " utility : "+ tmputility);
			if(tmputility>maxutility)
			{
				maxutility= tmputility;
				maxpathindex = pathseq.indexOf(path);

			}
			System.out.println("Max path index "+ maxpathindex+ ", utility : "+ maxutility);
		}
		return pathseq.get(maxpathindex);


	}

	private static double calculatePathUtility(ArrayList<Integer> path, int [][] gamedata) {
		double utility=0;
		ArrayList<Integer> considered= new ArrayList<Integer>();
		for(Integer n: path)
		{
			if(!considered.contains(n))
			{
				utility += targets.get(n).attackerreward; 
				considered.add(n);
			}
		}
		return utility;
	}


	public static void testPathGeneration()
	{
		int[] contractionsizes = {0,2,5,8,10};
		int ITER = 100;
		double[] result = new double[contractionsizes.length];
		int rindex=0;
		double perc = 80;
		for(int contractionsize: contractionsizes)
		{
			double sumsol = 0;
			//long sumtime = 0;
			long contractiontime=0;
			long solvingtime=0;
			long revmaptime=0;


			for(int iter=0; iter<ITER; iter++)
			{
				targets.clear();

				int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
				int nrow= 5;
				int ncol = 5;
				int dmax = 6;
				int nUnaccesstargets = contractionsize;
				int nRes = 3;
				int nTargets = nrow*ncol;
				System.out.println("\n Iter "+ iter);
				System.out.println("Number of targets "+ nrow*ncol);
				System.out.println("dmax "+ dmax);
				System.out.println("Unnecessary targets "+ nUnaccesstargets);
				System.out.println("nRes "+ nRes);
				SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

				//chooseDummyNodes(nUnaccesstargets);
				assignRandomDensity(perc/100, gamedata, targets);
				Date start = new Date();
				long l1 = start.getTime();
				ArrayList<TargetNode> dominatednodes = chooseContractedNodes(nUnaccesstargets, gamedata);




				System.out.println("Contracted targets "+ dominatednodes.size());
				//dominatednodes.add(targets.get(1));
				//buildOrigGraph(nrow, ncol, gamedata);
				//ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes);
				ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes, targets, dmax);
				Date stop = new Date();
				long l2 = stop.getTime();
				long diff = l2 - l1;
				contractiontime += diff;

				SecurityGameContraction.removePathsToDominatedNodes(contractednodes, targets);
				SecurityGameContraction.removeDominatedTargets(contractednodes, targets);
				//SecurityGameContraction.printNodesWithNeighborsAndPath(dominatednodes);
				int[][] p = new int[targets.size()][];
				try 
				{
					start = new Date();
					l1 = start.getTime();
					ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
					ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
					/**
					 * map has present id
					 * mapback gives the original ids
					 */
					HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
					HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
					makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
					System.out.println("Total path with duplicates "+pathseq.size());
					pathseq = removeDuplicatePathSimple(pathseq);
					//printPaths(pathseq);
					System.out.println("Total path without duplicates "+pathseq.size());

					//int k = 2;
					Integer[] input = new Integer[pathseq.size()];
					int[] branch = new int[nRes];//{0,0};//new char[k];


					for(int i=0; i<input.length; i++)
					{
						input[i] = i;
					}
					HashSet jSet=new HashSet();
					if(pathseq.size()==0)
					{
						//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
						//choose the worst payoff for defender

						Double mAxpayoff = Double.MIN_VALUE;
						Double defpayoff = 0.0;
						for(int i=0; i<contractednodes.size(); i++)
						{
							targets.add(contractednodes.get(i));
						}
						for(TargetNode x: targets)
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}

						sumsol += defpayoff;
						System.out.println("Defender expected payoff "+ defpayoff);
						/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/
						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						solvingtime += diff;


					}
					else
					{
						//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
						if(pathseq.size()<nRes)
						{

							branch = new int[pathseq.size()];
							jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
						}
						else
						{
							jSet=combine(input, nRes, 0, branch, 0, jSet);
						}

						List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
						/**
						 * columns will be combination of paths for each resources. 
						 */
						/**
						 * pmat, where columns will be combination of paths. 
						 * rows are targets. 
						 * each entry will say whether the target is in the joint schedule
						 */
						//jSet.

						//printJointSchedule(jset);

						//printNodesAsNeighbors(dominatednodes);

						p = makePmat(pathseq, jset, mapback, targets);
						//printPathMat(p);

						/**
						 * remove duplicates from p
						 */
						//removeDuplicatesFromP(p);
						//System.out.println();
						//printPathMat(p);

						System.out.println("Number of targets after contraction "+ targets.size());
						System.out.println("mip in... ");
						double[] coverage = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						solvingtime += diff;
						if(coverage.equals(null))
						{
							throw new Exception("Prob null...");
						}
						System.out.println("mip out... ");
						//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


						start = new Date();
						l1 = start.getTime();
						int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, dominatednodes, map, mapback, targets);
						//removeDuplicatesFromP(origpmat);
						//printPathMat(origpmat);
						//System.out.println("\n after mapping back");
						//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);

						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						if(contractionsize>0)
						{
							revmaptime += diff;
						}


						/*for(int i=0; i<coverage.length; i++)
			{
				for(int j=0; j<coverage[i].length; j++)
				{
					if(coverage[i][j]>0)
					{
						System.out.println("selected path : " + j);

						for(int k=0; k<origpmat.length; k++)
						{
							System.out.print(origpmat[k][j] + " ");
						}
						System.out.println();

					}
				}
			}*/			

						int maxtargetforattacker = findAttackTarget(origpmat, coverage, gamedata);

						double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, coverage);

						System.out.println("Attacked target is "+ maxtargetforattacker);
						System.out.println("Defender expected payoff "+ defexpectedpayoff);



						sumsol += defexpectedpayoff;
						/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defexpectedpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/

					}

				} 
				catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				targets.clear();
			}
			//double avgtime = sumtime/(ITER);
			double avgsol = sumsol/ITER;
			result[rindex++] = avgsol;
			DecimalFormat df = new DecimalFormat("#.#######");
			double revtime = revmaptime/ITER;
			String x = df.format(revtime);
			writeInFile(contractionsize, avgsol, contractiontime/ITER, solvingtime/ITER, Long.parseLong(x));

		}
		writeInFile(result,contractionsizes);

	}


	public static void testPathGenerationV3(double[][] density, int ITER, int nrow, int ncol,
			int[] percentages, int[] thresholds, double dmax, int nRes)
	{
		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;

		double defexp=0;
		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();

					int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);


					//ArrayList<TargetNode> dominatednodes = chooseContractedNodes(nUnaccesstargets, gamedata);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					selectDominatedTargets(targets, domindatednodes, threshold);

					System.out.println("Contracted targets "+ domindatednodes.size());
					//dominatednodes.add(targets.get(1));
					//buildOrigGraph(nrow, ncol, gamedata);
					//ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes);
					Date start = new Date();
					long l1 = start.getTime();


					Date tstart = new Date();
					long tl1 = tstart.getTime();

					sgc.contractGraph(domindatednodes, targets, dmax);
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;

					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);
					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					targetsize+=targets.size();
					int[][] p = new int[targets.size()][];
					try 
					{
						start = new Date();
						l1 = start.getTime();
						ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
						ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
						/**
						 * map has present id
						 * mapback gives the original ids
						 */
						HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
						HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
						makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
						System.out.println("Total path with duplicates "+pathseq.size());
						pathseq = removeDuplicatePathSimple(pathseq);
						//printPaths(pathseq);
						System.out.println("Total path without duplicates "+pathseq.size());

						//int k = 2;
						Integer[] input = new Integer[pathseq.size()];
						int[] branch = new int[nRes];//{0,0};//new char[k];


						for(int i=0; i<input.length; i++)
						{
							input[i] = i;
						}
						HashSet jSet=new HashSet();
						if(pathseq.size()==0)
						{
							//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
							//choose the worst payoff for defender

							Double mAxpayoff = Double.MIN_VALUE;
							Double defpayoff = 0.0;
							for(int i=0; i<domindatednodes.size(); i++)
							{
								targets.add(domindatednodes.get(i));
							}
							for(TargetNode x: targets)
							{
								if(x.attackerreward>mAxpayoff)
								{
									mAxpayoff= x.attackerreward;
									defpayoff = x.defenderpenalty;
								}
							}

							sumsol += defpayoff;
							System.out.println("Defender expected payoff "+ defpayoff);
							/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;


						}
						else
						{
							//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
							if(pathseq.size()<nRes)
							{

								branch = new int[pathseq.size()];
								jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
							}
							else
							{
								jSet=combine(input, nRes, 0, branch, 0, jSet);
							}

							List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
							/**
							 * columns will be combination of paths for each resources. 
							 */
							/**
							 * pmat, where columns will be combination of paths. 
							 * rows are targets. 
							 * each entry will say whether the target is in the joint schedule
							 */
							//jSet.

							//printJointSchedule(jset);

							//printNodesAsNeighbors(dominatednodes);

							p = makePmat(pathseq, jset, mapback, targets);
							//printPathMat(p);

							/**
							 * remove duplicates from p
							 */
							//removeDuplicatesFromP(p);
							//System.out.println();
							//printPathMat(p);

							System.out.println("Number of targets after contraction "+ targets.size());
							//System.out.println("mip in... ");
							System.out.println("Iter "+iter+", mip in... ");
							double[] coverage = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;
							if(coverage.equals(null))
							{
								throw new Exception("Prob null...");
							}
							System.out.println("mip out... ");
							//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


							start = new Date();
							l1 = start.getTime();
							int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
							//removeDuplicatesFromP(origpmat);
							//printPathMat(origpmat);
							//System.out.println("\n after mapping back");
							//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);



							Date tstop = new Date();
							long tl2 = tstop.getTime();
							long tdiff = tl2 - tl1;

							sumtime += tdiff;



							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							if(threshold>0)
							{
								revmaptime += diff;
							}


							/*for(int i=0; i<coverage.length; i++)
			{
				for(int j=0; j<coverage[i].length; j++)
				{
					if(coverage[i][j]>0)
					{
						System.out.println("selected path : " + j);

						for(int k=0; k<origpmat.length; k++)
						{
							System.out.print(origpmat[k][j] + " ");
						}
						System.out.println();

					}
				}
			}*/			

							int maxtargetforattacker = findAttackTarget(origpmat, coverage, gamedata);

							double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, coverage);
							defexp=defexpectedpayoff;
							System.out.println("Attacked target is "+ maxtargetforattacker);
							System.out.println("Defender expected payoff "+ defexpectedpayoff);



							sumsol += defexp;
							/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defexpectedpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/

						}

					} 
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					//targets.clear();

				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile( (int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);

	}





	public static void testPathGenerationWithExtreamPruningV3(double[][] density, int ITER, int nrow, int ncol,
			int[] percentages, int[] thresholds, double dmax)
	{
		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;

		double defexp=0;
		//int nUnaccesstargets = percentage;
		int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				//long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();

					int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);


					//ArrayList<TargetNode> dominatednodes = chooseContractedNodes(nUnaccesstargets, gamedata);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					selectDominatedTargets(targets, domindatednodes, threshold);

					System.out.println("Contracted targets "+ domindatednodes.size());
					//dominatednodes.add(targets.get(1));
					//buildOrigGraph(nrow, ncol, gamedata);
					//ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes);

					preComputeShortestPaths(0, targets, domindatednodes);

					Date start = new Date();
					long l1 = start.getTime();
					sgc.contractGraphWithExtreamPruning(domindatednodes, targets, dmax);
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;

					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);
					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					targetsize+=targets.size();
					int[][] p = new int[targets.size()][];
					try 
					{
						start = new Date();
						l1 = start.getTime();
						ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
						ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
						/**
						 * map has present id
						 * mapback gives the original ids
						 */
						HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
						HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
						makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
						System.out.println("Total path with duplicates "+pathseq.size());
						pathseq = removeDuplicatePathSimple(pathseq);
						//printPaths(pathseq);
						System.out.println("Total path without duplicates "+pathseq.size());

						//int k = 2;
						Integer[] input = new Integer[pathseq.size()];
						int[] branch = new int[nRes];//{0,0};//new char[k];


						for(int i=0; i<input.length; i++)
						{
							input[i] = i;
						}
						HashSet jSet=new HashSet();
						if(pathseq.size()==0)
						{
							//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
							//choose the worst payoff for defender

							Double mAxpayoff = Double.MIN_VALUE;
							Double defpayoff = 0.0;
							for(int i=0; i<domindatednodes.size(); i++)
							{
								targets.add(domindatednodes.get(i));
							}
							for(TargetNode x: targets)
							{
								if(x.attackerreward>mAxpayoff)
								{
									mAxpayoff= x.attackerreward;
									defpayoff = x.defenderpenalty;
								}
							}

							sumsol += defpayoff;
							System.out.println("Defender expected payoff "+ defpayoff);
							/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;


						}
						else
						{
							//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
							if(pathseq.size()<nRes)
							{

								branch = new int[pathseq.size()];
								jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
							}
							else
							{
								jSet=combine(input, nRes, 0, branch, 0, jSet);
							}

							List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
							/**
							 * columns will be combination of paths for each resources. 
							 */
							/**
							 * pmat, where columns will be combination of paths. 
							 * rows are targets. 
							 * each entry will say whether the target is in the joint schedule
							 */
							//jSet.

							//printJointSchedule(jset);

							//printNodesAsNeighbors(dominatednodes);

							p = makePmat(pathseq, jset, mapback, targets);
							//printPathMat(p);

							/**
							 * remove duplicates from p
							 */
							//removeDuplicatesFromP(p);
							//System.out.println();
							//printPathMat(p);

							System.out.println("Number of targets after contraction "+ targets.size());
							//System.out.println("mip in... ");
							System.out.println("Iter "+iter+", mip in... ");
							double[] coverage = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;
							if(coverage.equals(null))
							{
								throw new Exception("Prob null...");
							}
							System.out.println("mip out... ");
							//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


							start = new Date();
							l1 = start.getTime();
							int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);
							//removeDuplicatesFromP(origpmat);
							//printPathMat(origpmat);
							//System.out.println("\n after mapping back");
							//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);

							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							if(threshold>0)
							{
								revmaptime += diff;
							}


							/*for(int i=0; i<coverage.length; i++)
			{
				for(int j=0; j<coverage[i].length; j++)
				{
					if(coverage[i][j]>0)
					{
						System.out.println("selected path : " + j);

						for(int k=0; k<origpmat.length; k++)
						{
							System.out.print(origpmat[k][j] + " ");
						}
						System.out.println();

					}
				}
			}*/			

							int maxtargetforattacker = findAttackTarget(origpmat, coverage, gamedata);

							double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, coverage);
							defexp=defexpectedpayoff;
							System.out.println("Attacked target is "+ maxtargetforattacker);
							System.out.println("Defender expected payoff "+ defexpectedpayoff);



							sumsol += defexp;
							/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defexpectedpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/

						}

					} 
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					//targets.clear();

				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile( (int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, Long.parseLong(x));
			}

		}
		//writeInFile(result,contractionsizes);

	}


	private static void preComputeShortestPaths(int basetargetid,
			ArrayList<TargetNode> targets, ArrayList<TargetNode> domindatednodes) {


		ArrayList<Integer> pathnodes = new ArrayList<Integer>();
		TargetNode base = getTargetNode(basetargetid, targets);
		for(TargetNode n: targets)
		{
			double dist = findShortestPathThrougGraph(base, n, targets, pathnodes);
			System.out.println("Distance from base to target "+n.getTargetid() + " : "+dist);
			n.setDistfrombase(dist);
		}



	}

	public static void testPathGenerationv2()
	{
		int[] contractionsizes = {8};
		int ITER = 1;
		double[] result = new double[contractionsizes.length];
		int rindex=0;
		double perc = 80;
		for(int contractionsize: contractionsizes)
		{
			double sumsol = 0;
			//long sumtime = 0;
			long contractiontime=0;
			long solvingtime=0;
			long revmaptime=0;


			for(int iter=0; iter<ITER; iter++)
			{
				targets.clear();

				int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
				int nrow= 4;
				int ncol = 4;
				int dmax = 6;
				int nUnaccesstargets = contractionsize;
				int nRes = 2;
				int nTargets = nrow*ncol;
				System.out.println("\n Iter "+ iter);
				System.out.println("Number of targets "+ nrow*ncol);
				System.out.println("dmax "+ dmax);
				System.out.println("Unnecessary targets "+ nUnaccesstargets);
				System.out.println("nRes "+ nRes);
				SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

				//chooseDummyNodes(nUnaccesstargets);
				//assignSimilarUtility(perc/100, gamedata);
				Date start = new Date();
				long l1 = start.getTime();
				ArrayList<TargetNode> dominatednodes = chooseContractedNodes(nUnaccesstargets, gamedata);

				System.out.println("Contracted targets "+ dominatednodes.size());
				//dominatednodes.add(targets.get(1));
				//buildOrigGraph(nrow, ncol, gamedata);
				//ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes);
				ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes, targets, dmax);
				Date stop = new Date();
				long l2 = stop.getTime();
				long diff = l2 - l1;
				contractiontime += diff;

				SecurityGameContraction.removePathsToDominatedNodes(contractednodes, targets);
				SecurityGameContraction.removeDominatedTargets(contractednodes, targets);
				//SecurityGameContraction.printNodesWithNeighborsAndPath(dominatednodes);
				int[][] p = new int[targets.size()][];
				try 
				{
					start = new Date();
					l1 = start.getTime();
					ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);
					ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
					/**
					 * map has present id
					 * mapback gives the original ids
					 */
					HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
					HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
					makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets );
					System.out.println("Total path with duplicates "+pathseq.size());
					pathseq = removeDuplicatePathSimple(pathseq);
					printPaths(pathseq);
					System.out.println("Total path without duplicates "+pathseq.size());

					//int k = 2;
					Integer[] input = new Integer[pathseq.size()];
					int[] branch = new int[nRes];//{0,0};//new char[k];


					for(int i=0; i<input.length; i++)
					{
						input[i] = i;
					}
					HashSet jSet=new HashSet();
					if(pathseq.size()==0)
					{
						//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
						//choose the worst payoff for defender

						Double mAxpayoff = Double.MIN_VALUE;
						Double defpayoff = 0.0;
						for(int i=0; i<contractednodes.size(); i++)
						{
							targets.add(contractednodes.get(i));
						}
						for(TargetNode x: targets)
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
						}

						sumsol += defpayoff;
						System.out.println("Defender expected payoff "+ defpayoff);
						/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/
						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						solvingtime += diff;


					}
					else
					{
						//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
						if(pathseq.size()<nRes)
						{

							branch = new int[pathseq.size()];
							jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
						}
						else
						{
							jSet=combine(input, nRes, 0, branch, 0, jSet);
						}

						List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
						/**
						 * columns will be combination of paths for each resources. 
						 */
						/**
						 * pmat, where columns will be combination of paths. 
						 * rows are targets. 
						 * each entry will say whether the target is in the joint schedule
						 */
						//jSet.

						printJointSchedule(jset);

						//printNodesAsNeighbors(dominatednodes);

						p = makePmat(pathseq, jset, mapback, targets);
						//printPathMat(p);

						/**
						 * remove duplicates from p
						 */
						//removeDuplicatesFromP(p);
						//System.out.println();
						//printPathMat(p);

						System.out.println("Number of targets after contraction "+ targets.size());
						System.out.println("mip in... ");
						double[] coverage = MIPSolver4.solve(p, gamedata, SecurityGameContraction.targets, nRes);
						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						solvingtime += diff;
						if(coverage.equals(null))
						{
							throw new Exception("Prob null...");
						}
						System.out.println("mip out... ");
						//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


						start = new Date();
						l1 = start.getTime();
						int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, dominatednodes, map, mapback, targets);
						//removeDuplicatesFromP(origpmat);
						//printPathMat(origpmat);
						//System.out.println("\n after mapping back");
						//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);

						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						if(contractionsize>0)
						{
							revmaptime += diff;
						}


						/*for(int i=0; i<coverage.length; i++)
			{
				for(int j=0; j<coverage[i].length; j++)
				{
					if(coverage[i][j]>0)
					{
						System.out.println("selected path : " + j);

						for(int k=0; k<origpmat.length; k++)
						{
							System.out.print(origpmat[k][j] + " ");
						}
						System.out.println();

					}
				}
			}*/			

						int maxtargetforattacker = findAttackTarget(origpmat, coverage, gamedata);

						double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, coverage);

						System.out.println("Attacked target is "+ maxtargetforattacker);
						System.out.println("Defender expected payoff "+ defexpectedpayoff);



						sumsol += defexpectedpayoff;
						/*try
						{
							PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
							pw.append(iter+ "," + defexpectedpayoff+"\n");
							pw.close();

						}
						catch(Exception e)
						{

						}*/

					}

				} 
				catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				targets.clear();
			}
			//double avgtime = sumtime/(ITER);
			double avgsol = sumsol/ITER;
			result[rindex++] = avgsol;
			DecimalFormat df = new DecimalFormat("#.#######");
			double revtime = revmaptime/ITER;
			String x = df.format(revtime);
			writeInFile(contractionsize, avgsol, contractiontime/ITER, solvingtime/ITER, Long.parseLong(x));

		}
		writeInFile(result,contractionsizes);

	}




	public static int randInt(int min, int max) {

		// Usually this should be a field rather than a method variable so
		// that it is not re-seeded every call.


		// nextInt is normally exclusive of the top value,
		// so add 1 to make it inclusive
		int randomNum = rand1.nextInt((max - min) + 1) + min;

		return randomNum;
	}


	public static void assignRandomDensityZeroSum(double[][] density, int[][] gamedata, ArrayList<TargetNode> targets, int iter) {

		for(int i=0; i<targets.size(); i++)
		{
			targets.get(i).defenderreward=0;//density[iter][i];
			targets.get(i).setAnimaldensity(density[iter][i]);
			targets.get(i).defenderpenalty=-density[iter][i];
			targets.get(i).attackerpenalty=0;
			targets.get(i).attackerreward=density[iter][i];



			gamedata[i][0] = 0;//(int)density[iter][i] ;
			gamedata[i][1] = -(int)density[iter][i];
			gamedata[i][2] = (int)density[iter][i];  // uncovered
			gamedata[i][3] = 0; //covered

		}







	}

	
	
	private static void assignRandomDensityZeroSum(HashMap<Integer,TargetNode> targetmaps, int[][] gamedata, ArrayList<TargetNode> targets, int iter) {

		for(int i=0; i<targets.size(); i++)
		{
			/*targets.get(i).defenderreward=0;//density[iter][i];
			targets.get(i).setAnimaldensity(density[iter][i]);
			targets.get(i).defenderpenalty=-density[iter][i];
			targets.get(i).attackerpenalty=0;
			targets.get(i).attackerreward=density[iter][i];*/
			
			
			//targets.add(targetmaps.get(i));



			gamedata[i][0] = 0;//(int)density[iter][i] ;
			gamedata[i][1] = (int)targetmaps.get(i).defenderpenalty;
			gamedata[i][2] = (int)targetmaps.get(i).attackerreward;  // uncovered
			gamedata[i][3] = 0; //covered

		}







	}


	private static void assignRandomDensity(double d, int[][] gamedata, ArrayList<TargetNode> targets) {

		int lowstart=1; // dont put 0
		int lowlimit=1;
		int highstart=8;
		int highlimit=10;


		int n = (int)Math.ceil(targets.size()*d);
		int l=0;
		ArrayList<Integer> added = new ArrayList<Integer>();
		int dominatedcount = 0;
		Random rand = new Random();
		int settozero = targets.size() - n;
		int g=4;
		while(true)
		{
			//int randomNum = rand.nextInt((targets.size() - 1-1) + 1) + 1;
			int randomNum = rand.nextInt((targets.size() - 1-1) + 1) + 1;

			int lowutility =  randInt(lowstart, lowlimit);

			if(!added.contains(randomNum))
			{
				added.add(randomNum);
				targets.get(randomNum).defenderreward=lowutility;
				targets.get(randomNum).setAnimaldensity(lowutility);
				targets.get(randomNum).defenderpenalty=-lowutility;
				targets.get(randomNum).attackerpenalty=-lowutility;
				targets.get(randomNum).attackerreward=lowutility;


				gamedata[randomNum][0] = lowutility;
				gamedata[randomNum][1] = -lowutility;
				gamedata[randomNum][2] = lowutility;
				gamedata[randomNum][3] = -lowutility;
				dominatedcount++;
				g++;
				System.out.println("target "+ g + " has 0 payoff");
				//}
				if(dominatedcount==settozero)
					break;
			}

		}
		int k=0;
		for(TargetNode x: targets)
		{

			int highutility =  randInt(highstart, highlimit);

			if(x.getAnimaldensity()==0)
			{
				/*targets.get(x.getTargetid()).defenderreward=highutility;
				targets.get(x.getTargetid()).setAnimaldensity(highutility);
				targets.get(x.getTargetid()).defenderpenalty=-highutility;
				targets.get(x.getTargetid()).attackerpenalty=-highutility;
				targets.get(x.getTargetid()).attackerreward=highutility;


				gamedata[x.getTargetid()][0] =highutility ;
				gamedata[x.getTargetid()][1] = -highutility;
				gamedata[x.getTargetid()][2] = highutility;
				gamedata[x.getTargetid()][3] = -highutility;*/
				targets.get(k).defenderreward=highutility;
				targets.get(k).setAnimaldensity(highutility);
				targets.get(k).defenderpenalty=-highutility;
				targets.get(x.getTargetid()).attackerpenalty=-highutility;
				targets.get(x.getTargetid()).attackerreward=highutility;


				gamedata[k][0] =highutility ;
				gamedata[k][1] = -highutility;
				gamedata[k][2] = highutility;
				gamedata[k][3] = -highutility;
				//System.out.println("target "+ x.getTargetid() + " has positive payoff");
				k=k+1;
			}
		}


	}

	private static void writeInFile(double[] result, int[] contractionsizes) {

		double[] epsilon = new double[result.length];
		int i=0;
		for(int csize: contractionsizes)
		{
			epsilon[i] = (result[i]- result[0])/(result[0]);
			epsilon[i] = epsilon[i]*100;
			try
			{
				PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
				pw.append(csize+","+epsilon[i]+"\n");
				pw.close();

			}
			catch(Exception e)
			{

			}
			i++;

		}


	}

	private static void chooseDummyNodes(int nUnaccesstargets) {




	}



	private static void writeInFile(String filename, int finalsize, double avgsol, double contracttime,
			double solvingtime, double revmaptime) 
	{

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File(filename),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(finalsize+","+avgsol+ ","+contracttime+"," + solvingtime+ ","+revmaptime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}




	private static void writeInFile(double contractionsize, double avgsol, double contracttime,
			long solvingtime, long revmaptime) 
	{

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(contractionsize+ "," + avgsol+ ","+contracttime+"," + solvingtime+ ","+revmaptime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}

	private static void printPmatWithPositiveCoverage(int[][] origpmat,
			double[] coverage, List<ArrayList<Integer>> jset,
			ArrayList<ArrayList<Integer>> pathseq, HashMap<Integer, Integer> map) {


		for(int i=0; i<coverage.length; i++)
		{
			if(coverage[i]>0)
			{
				System.out.println("Joint schedule "+ i+ " is selected \n Paths in this schedule are : \n");
				for(Integer x: jset.get(i))
				{
					System.out.println("path "+x+":");
					System.out.println("\ncovered targets: ");
					//ArrayList<Integer> nodes = pathseq.get(x);
					for(int j=0; j<origpmat.length; j++)
					{
						//if(origpmat[j][i]>0)
						//int id = map.get(nodes.get(j));
						System.out.println(j+" is covered "+ origpmat[j][i]);
					}


				}
			}
			//System.out.println();
		}

	}

	public static void printPaths(ArrayList<ArrayList<Integer>> pathseq) {


		System.out.println("Paths ");
		Logger.logit("paths : \n");

		int i=0;
		for(ArrayList<Integer> path: pathseq)
		{
			System.out.print(i+" : ");

			for(Integer x: path)
			{
				System.out.print(x+"->");
				//Logger.logit(x+"->");
			}
			System.out.println();
			//Logger.logit("\n");
			i++;

		}



	}

	private static void printPathWithPositiveCoverage(int[][] origpmat,
			double[] coverage, List<ArrayList<Integer>> jset, ArrayList<ArrayList<Integer>> pathseq,
			HashMap<Integer,Integer> map) {


		for(int i=0; i<coverage.length; i++)
		{
			if(coverage[i]>0)
			{
				System.out.println("Joint schedule "+ i+ " is selected \n Paths in this schedule are : \n");
				for(Integer x: jset.get(i))
				{
					System.out.println("path "+x+":");
					System.out.println("\ncovered targets: ");
					ArrayList<Integer> nodes = pathseq.get(x);
					for(int j=0; j<nodes.size(); j++)
					{
						//if(origpmat[j][i]>0)
						int id = map.get(nodes.get(j));
						System.out.println(nodes.get(j)+" is covered ");
					}


				}
			}
			//System.out.println();
		}


	}

	public static void printJointSchedule(List<ArrayList<Integer>> jset) {

		System.out.println();
		for(int i=0; i<jset.size(); i++)
		{
			System.out.print("[");
			for(int j=0; j<jset.get(i).size(); j++)
			{
				System.out.print(jset.get(i).get(j)+",");
			}
			System.out.println("]");
		}
		System.out.println();


	}

	private static void removeDuplicatesFromP(int[][] p) {


		for(int i=0; i<p.length; i++)
		{
			for(int j=0; j<p[i].length; j++)
			{
				/**
				 * see if p[i][j] is 1
				 * id so, make any 1 after that 0
				 */
				if(p[i][j]==1)
				{
					for(int k=j+1; k<p[i].length; k++)
					{
						if(p[i][k]==1)
						{
							p[i][k]=0;
						}
					}
				}
			}
		}



	}

	public static int[][] makePmat(ArrayList<ArrayList<Integer>> pathseq,
			List<ArrayList<Integer>> jset, HashMap<Integer,Integer> mapback, ArrayList<TargetNode> targets) {





		int[][] pmat = new int[targets.size()][jset.size()];



		/**
		 * mapping happened
		 */
		for(int t=0; t<targets.size(); t++)
		{
			for(int j=0; j<jset.size(); j++)
			{
				/**
				 * check if target t is in j schedule
				 */
				int targetid = mapback.get(t);//targets.get(t).getTargetid();
				int isinj = isInJointSchedule(targetid, jset.get(j), pathseq);
				pmat[t][j] = isinj;

			}
		}


		return pmat;
	}
	
	
	
	public static int[][] makeSuperPmat(ArrayList<ArrayList<Integer>> pathseq,
			List<ArrayList<Integer>> jset, HashMap<Integer,Integer> mapback, HashMap<Integer, 
			SuperTarget> supertargets, HashMap<Integer,Integer> map) {





		int[][] pmat = new int[supertargets.size()][jset.size()];



		/**
		 * mapping happened
		 */
		for(Integer stid: supertargets.keySet())
		{
			for(int j=0; j<jset.size(); j++)
			{
				/**
				 * check if target t is in j schedule
				 */
				int targetid = stid;//targets.get(t).getTargetid();
				int isinj = isInJointSchedule(targetid, jset.get(j), pathseq);
				pmat[map.get(stid)][j] = isinj;

			}
		}


		return pmat;
	}
	
	

	private static int isInJointSchedule(int target, ArrayList<Integer> jointschedule,
			ArrayList<ArrayList<Integer>> pathseq) {


		for(Integer j: jointschedule)
		{
			ArrayList<Integer> path = pathseq.get(j);
			for(Integer tt: path)
			{
				if(target==tt)
					return 1;
			}
		}


		return 0;
	}

	static HashSet combine(Integer[] arr, int k, int startId, int[] branch, int numElem,HashSet arrSet)
	{
		if (numElem == k)
		{
			//System.out.println("k: "+k+(Arrays.toString(branch)));
			ArrayList<Integer> mySet = new ArrayList<Integer>();
			for(int i=0;i<branch.length;i++)
			{
				mySet.add(branch[i]);
			}
			arrSet.add(mySet);
			return arrSet;
		}

		for (int i = startId; i < arr.length; ++i)
		{
			branch[numElem++]=arr[i];
			combine(arr, k, ++startId, branch, numElem, arrSet);
			--numElem;
		}
		return arrSet;
	}



	public static ArrayList<ArrayList<Integer>> removeDuplicatePathSimple(
			ArrayList<ArrayList<Integer>> pathseq) {

		ArrayList<Integer> duplicatepaths = new ArrayList<Integer>();


		for(int i=0; i<pathseq.size()-1; i++)
		{
			if(!duplicatepaths.contains(i))
			{

				for(int j=i+1; j<pathseq.size(); j++)
				{
					//if(pathseq.get(i).size()==pathseq.get(j).size())
					{
						boolean same = true;

						for(int k=0; k<pathseq.get(j).size(); k++)
						{
							boolean isinj = isInSuperJ(pathseq.get(j).get(k), pathseq.get(i));
							if(!isinj)
							{
								same=false;
								break;
							}

						}
						if(same)
						{
							if(!duplicatepaths.contains(j))
							{
								//System.out.println("Path "+ j +" is same as path "+ i);
								duplicatepaths.add(j);
							}
						}
					}


				}
			}



		}

		ArrayList<ArrayList<Integer>> pathseqcopy = new ArrayList<ArrayList<Integer>>();

		int index = 0;
		for(ArrayList<Integer> p: pathseq)
		{
			if(!duplicatepaths.contains(index))
			{
				pathseqcopy.add(p);
			}

			index++;
		}

		return pathseqcopy;

	}


	private static ArrayList<ArrayList<Integer>> removeDuplicatePath(
			ArrayList<ArrayList<Integer>> pathseq) {

		ArrayList<Integer> duplicatepaths = new ArrayList<Integer>();
		//printPaths(pathseq);

		for(int i=0; i<pathseq.size()-1; i++)
		{
			if(!duplicatepaths.contains(i))
			{

				for(int j=i+1; j<pathseq.size(); j++)
				{
					boolean same = true;


					//if(pathseq.get(j).size()>=pathseq.get(i).size())
					{
						int count=0;

						for(int k=0; k<pathseq.get(i).size(); k++)
						{
							boolean isinj = isInJ(pathseq.get(i).get(k), pathseq.get(j));
							if(!isinj)
							{
								same=false;
								break;
							}
							count++;

						}
						if(same)
						{

							Set<Integer> seti = new HashSet<Integer>(pathseq.get(i));
							Set<Integer> setj = new HashSet<Integer>(pathseq.get(j));



							if(pathseq.get(i).size()>=pathseq.get(j).size() && (setj.size()<=seti.size()))
							{
								if(!duplicatepaths.contains(j))
								{
									//System.out.println("\nPath "+ j +" is same as path "+ i+ ", path "+ j + " added to duplicate list");

									//System.out.print("Path "+ j + " : ");

									/*for(Integer p: pathseq.get(j))
									{
										System.out.print(p+"-> ");
									}*/


									/*System.out.print("\nPath "+ i + " : ");

									for(Integer p: pathseq.get(i))
									{
										System.out.print(p+"-> ");
									}

									System.out.println();*/


									duplicatepaths.add(j);
								}
							}
							/*else
							{
								if(!duplicatepaths.contains(i))
								{
									System.out.println("\nPath "+ j +" is same as path "+ i+ ", path "+ i + " added to duplicate list");

									System.out.print("Path "+ j + " : ");

									for(Integer p: pathseq.get(j))
									{
										System.out.print(p+"-> ");
									}


									System.out.print("\nPath "+ i + " : ");

									for(Integer p: pathseq.get(i))
									{
										System.out.print(p+"-> ");
									}

									System.out.println();


									duplicatepaths.add(i);
								}

							}*/
						}
					}
					/*else
					{
						int count=0;
						for(int k=0; k<pathseq.get(j).size(); k++)
						{
							boolean isinj = isInJ(pathseq.get(j).get(k), pathseq.get(i));
							if(!isinj)
							{
								same=false;
								break;
							}
							count++;

						}
						if(same)
						{
							if(!duplicatepaths.contains(j))
							{
								//System.out.println("Path "+ j +" is same as path "+ i);
								duplicatepaths.add(j);
							}
						}
					}*/


				}
			}



		}

		ArrayList<ArrayList<Integer>> pathseqcopy = new ArrayList<ArrayList<Integer>>();

		int index = 0;
		for(ArrayList<Integer> p: pathseq)
		{
			if(!duplicatepaths.contains(index))
			{
				//System.out.println("Adding path "+ index);
				pathseqcopy.add(p);
			}

			index++;
		}

		//printPaths(pathseqcopy);

		return pathseqcopy;

	}

	private static boolean isInJ(Integer integer, ArrayList<Integer> arrayList) {


		for(Integer x: arrayList)
		{
			if(integer.equals(x))
				return true;
		}


		return false;
	}
	
	
	public static boolean isInSuperJ(Integer t, ArrayList<Integer> path) {


		for(Integer p: path)
		{
			if(p.equals(t))
			{
				//System.out.println( p+ " is equal to "+ t);
				return true;
			}
			else
			{
				//System.out.println( p+ " is not equal to "+ t);
			}
		}


		return false;
	}

	public static double expectedPayoffDef(int target,
			int[][] origpmat, int[][] gamedata, double[] coverage) {

		double prob =0;

		for(int path=0; path<coverage.length; path++)
		{
			prob += origpmat[target][path]*coverage[path];
		}

		double payoff = prob*(gamedata[target][0]) + (1-prob)*gamedata[target][1];
		payoff = Math.round(payoff*100)/100.0;




		return payoff;
	}
	
	
	public static double expectedPayoffDef(int target,
			int[][] origpmat, HashMap<Integer, TargetNode> targetmaps, double[] coverage) {

		double prob =0;

		for(int path=0; path<coverage.length; path++)
		{
			prob += origpmat[target][path]*coverage[path];
		}

		double payoff = prob*(targetmaps.get(target).defenderreward) + (1-prob)*targetmaps.get(target).defenderpenalty;
		payoff = Math.round(payoff*100)/100.0;




		return payoff;
	}




	public static double expectedPayoffAtt(int target,
			int[][] origpmat, int[][] gamedata, double[] coverage) {

		double prob =0;

		for(int path=0; path<coverage.length; path++)
		{
			prob += origpmat[target][path]*coverage[path];
		}

		double payoff = prob*(gamedata[target][3]) + (1-prob)*gamedata[target][2];
		payoff = Math.round(payoff*100)/100.0;




		return payoff;
	}
	
	
	public static double expectedPayoffAtt(int target,
			int[][] origpmat, HashMap<Integer, TargetNode> targetmaps, double[] coverage) {

		double prob =0;

		for(int path=0; path<coverage.length; path++)
		{
			prob += origpmat[target][path]*coverage[path];
		}

		double payoff = prob*(targetmaps.get(target).attackerpenalty) + (1-prob)*targetmaps.get(target).attackerreward;
		payoff = Math.round(payoff*100)/100.0;




		return payoff;
	}



	private static double expectedPayoffDefWMapping(int target,
			int[][] origpmat, int[][] gamedata, double[] coverage, HashMap<Integer,Integer> originalmap) {

		double prob =0;

		for(int path=0; path<coverage.length; path++)
		{
			prob += origpmat[originalmap.get(target)][path]*coverage[path];
		}

		double payoff = prob*(gamedata[target][0]) + (1-prob)*gamedata[target][1];
		payoff = Math.round(payoff*100)/100.0;




		return payoff;
	}


	private static void buildOrigGraph(int numRow, int numCol, int[][] gamedata) {



		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, gamedata[target][0]);
			node.defenderreward = gamedata[target][0];
			node.defenderpenalty = gamedata[target][1];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = gamedata[target][3];

			origtargets.add(node);
			if(target==0)
			{
				//graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;

		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				origtargets.get(targetid).setRowCol(row, col);
				origtargets.get(targetid).setCoinvalue(gamedata[targetid][0]);
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						origtargets.get(targetid).addNeighbor(origtargets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						origtargets.get(targetid).setPath(origtargets.get(neighborid), pathnodes);
						origtargets.get(targetid).setPathUtility(origtargets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						Double distance = 1.0;//rand.nextDouble()*10+5;
						origtargets.get(targetid).addDistance(origtargets.get(neighborid), Math.floor(distance));
						origtargets.get(neighborid).addDistance(origtargets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}


	}
	
	public static int findAttackSuperTargetWMapping(int[][] pmat, double[] coverage, HashMap<Integer, SuperTarget> currentst, 
			HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		double[] expectedpayoffs = expectedAttackerSTPayoffsWithMapping(pmat, coverage,  currentst, map, mapback);

		/*for(int i=0; i<expectedpayoffs.length; i++)
		{
			System.out.println("target "+ mapback.get(i) +", attkr payoff "+ expectedpayoffs[i]);
		}*/
		double max = Double.NEGATIVE_INFINITY;
		int maxtarget = -1;

		for(int i=0; i<expectedpayoffs.length; i++)
		{
			if(max<expectedpayoffs[i])
			{
				max = expectedpayoffs[i];
				maxtarget = i;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(int i=0; i<pmat.length; i++)
		{
			if((i!=maxtarget) && (expectedpayoffs[i] == max))
			{
				//System.out.println("tied target "+mapback.get(i) +", attkr payoff "+ expectedpayoffs[i]);

				tiedtargets.add(i);
			}
		}

		int selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{

			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			double defenderexpectedpayoffs[] =  defenderExpectedSTPayoffsWithMapping(tiedtargets, coverage, currentst, pmat, map, mapback);
			/*for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				System.out.println("target "+mapback.get(tiedtargets.get(i)) +", defender payoff "+ defenderexpectedpayoffs[i]);
			}*/
			double defmax = Double.NEGATIVE_INFINITY;
			int defmaxtarget = -1;

			for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				if(defmax<defenderexpectedpayoffs[i])
				{
					defmax = defenderexpectedpayoffs[i];
					defmaxtarget = i;
				}
			}
			selectedtarget = defmaxtarget;
			return tiedtargets.get(selectedtarget);
		}



		return selectedtarget;
	}
	
	
	
	public static int findAttackTargetWMapping(int[][] pmat, double[] coverage, int[][] gamedata, 
			HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		double[] expectedpayoffs = expectedAttackerPayoffsWithMapping(pmat, coverage, gamedata, map, mapback);

		/*for(int i=0; i<expectedpayoffs.length; i++)
		{
			System.out.println("target "+ mapback.get(i) +", attkr payoff "+ expectedpayoffs[i]);
		}*/
		double max = Double.NEGATIVE_INFINITY;
		int maxtarget = -1;

		for(int i=0; i<expectedpayoffs.length; i++)
		{
			if(max<expectedpayoffs[i])
			{
				max = expectedpayoffs[i];
				maxtarget = i;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(int i=0; i<pmat.length; i++)
		{
			if((i!=maxtarget) && (expectedpayoffs[i] == max))
			{
				//System.out.println("tied target "+mapback.get(i) +", attkr payoff "+ expectedpayoffs[i]);

				tiedtargets.add(i);
			}
		}

		int selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{

			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			double defenderexpectedpayoffs[] =  defenderExpectedPayoffsWithMapping(tiedtargets, coverage, gamedata, pmat, map, mapback);
			/*for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				System.out.println("target "+mapback.get(tiedtargets.get(i)) +", defender payoff "+ defenderexpectedpayoffs[i]);
			}*/
			double defmax = Double.NEGATIVE_INFINITY;
			int defmaxtarget = -1;

			for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				if(defmax<defenderexpectedpayoffs[i])
				{
					defmax = defenderexpectedpayoffs[i];
					defmaxtarget = i;
				}
			}
			selectedtarget = defmaxtarget;
			return tiedtargets.get(selectedtarget);
		}



		return selectedtarget;
	}
	
	
	

	private static double[] findAttackTargetWHashMap(int[][] pmat, double[] coverage, int[][] gamedata, 
			HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		HashMap<Integer,Double> expectedpayoffs = expectedAttackerPayoffsWithHashMap(pmat, coverage, gamedata, map, mapback);

		/*for(int i=0; i<expectedpayoffs.length; i++)
		{
			System.out.println("target "+ mapback.get(i) +", attkr payoff "+ expectedpayoffs[i]);
		}*/
		double max = Double.NEGATIVE_INFINITY;
		int maxtarget = -1;

		for(int i: expectedpayoffs.keySet())
		{
			if(max<expectedpayoffs.get(i))
			{
				max = expectedpayoffs.get(i);
				maxtarget = i;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(int i: expectedpayoffs.keySet())
		{
			if((i!=maxtarget) && (expectedpayoffs.get(i) == max))
			{
				//System.out.println("tied target "+mapback.get(i) +", attkr payoff "+ expectedpayoffs[i]);

				tiedtargets.add(i);
			}
		}

		int selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{

			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			HashMap<Integer,Double> defenderexpectedpayoffs =  defenderExpectedPayoffsWithHashMap(tiedtargets, coverage, gamedata, pmat, map, mapback);
			/*for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				System.out.println("target "+mapback.get(tiedtargets.get(i)) +", defender payoff "+ defenderexpectedpayoffs[i]);
			}*/
			double defmax = Double.NEGATIVE_INFINITY;
			int defmaxtarget = -1;

			for(int i : defenderexpectedpayoffs.keySet())
			{
				if(defmax<defenderexpectedpayoffs.get(i))
				{
					defmax = defenderexpectedpayoffs.get(i);
					defmaxtarget = i;
				}
			}
			selectedtarget = defmaxtarget;
			return new double[]{selectedtarget, max};
		}



		return new double[]{selectedtarget, max};
	}

	public static int findAttackTarget(int[][] pmat, double[] coverage, int[][] gamedata) {


		double[] expectedpayoffs = expectedAttackerPayoffs(pmat, coverage, gamedata);

		/*for(int i=0; i<expectedpayoffs.length; i++)
		{
			System.out.println("target "+i+", attkr payoff "+ expectedpayoffs[i]);
		}*/
		double max = Double.NEGATIVE_INFINITY;
		int maxtarget = -1;

		for(int i=0; i<expectedpayoffs.length; i++)
		{
			if(max<expectedpayoffs[i])
			{
				max = expectedpayoffs[i];
				maxtarget = i;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(int i=0; i<pmat.length; i++)
		{
			if((i!=maxtarget) && (expectedpayoffs[i] == max))
			{
				//System.out.println("tied target "+i+", attkr payoff "+ expectedpayoffs[i]);

				tiedtargets.add(i);
			}
		}

		int selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{
			//System.out.println("target "+i+", attkr payoff "+ expectedpayoffs[i]);
			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			double defenderexpectedpayoffs[] =  defenderExpectedPayoffs(tiedtargets, coverage, gamedata, pmat);
			/*for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				System.out.println("target "+tiedtargets.get(i)+", defender payoff "+ defenderexpectedpayoffs[i]);
			}*/
			double defmax = Double.NEGATIVE_INFINITY;
			int defmaxtarget = -1;

			for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				if(defmax<defenderexpectedpayoffs[i])
				{
					defmax = defenderexpectedpayoffs[i];
					defmaxtarget = i;
				}
			}
			selectedtarget = defmaxtarget;
			return tiedtargets.get(selectedtarget);
		}



		return selectedtarget;
	}
	
	
	public static int findAttackTarget(int[][] pmat, double[] coverage, HashMap<Integer, TargetNode> targetmaps) {


		double[] expectedpayoffs = expectedAttackerPayoffs(pmat, coverage, targetmaps);

		/*for(int i=0; i<expectedpayoffs.length; i++)
		{
			System.out.println("target "+i+", attkr payoff "+ expectedpayoffs[i]);
		}*/
		double max = Double.NEGATIVE_INFINITY;
		int maxtarget = -1;

		for(int i=0; i<expectedpayoffs.length; i++)
		{
			if(max<expectedpayoffs[i])
			{
				max = expectedpayoffs[i];
				maxtarget = i;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(int i=0; i<pmat.length; i++)
		{
			if((i!=maxtarget) && (expectedpayoffs[i] == max))
			{
				//System.out.println("tied target "+i+", attkr payoff "+ expectedpayoffs[i]);

				tiedtargets.add(i);
			}
		}

		int selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{
			//System.out.println("target "+i+", attkr payoff "+ expectedpayoffs[i]);
			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			double defenderexpectedpayoffs[] =  defenderExpectedPayoffs(tiedtargets, coverage, targetmaps, pmat);
			/*for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				System.out.println("target "+tiedtargets.get(i)+", defender payoff "+ defenderexpectedpayoffs[i]);
			}*/
			double defmax = Double.NEGATIVE_INFINITY;
			int defmaxtarget = -1;

			for(int i=0; i<defenderexpectedpayoffs.length; i++)
			{
				if(defmax<defenderexpectedpayoffs[i])
				{
					defmax = defenderexpectedpayoffs[i];
					defmaxtarget = i;
				}
			}
			selectedtarget = defmaxtarget;
			return tiedtargets.get(selectedtarget);
		}



		return selectedtarget;
	}
	
	
	private static HashMap<Integer,Double> defenderExpectedPayoffsWithHashMap(
			ArrayList<Integer> targets, double[] coverage,
			int[][] gamedata, int [][] origpmat, 
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback) {

		HashMap<Integer,Double> expectedpayoff = new HashMap<Integer,Double>();

		double c = 0.0;



		for(int target : targets)
		{
			//int target = targets.get(t);

			c=0.0;
			for(int path=0; path<coverage.length; path++)
			{
				c += origpmat[originalmap.get(target)][path]*coverage[path];


				/*double defreward = origpmat[target][path]*coverage[path]*gamedata[target][0];
					double defpenalty = (1-origpmat[target][path])*(1-coverage[path])*gamedata[target][1];
					sum = sum + defpenalty + defreward;*/
			}


			double tmpexp = c*gamedata[target][0] + (1-c)*gamedata[target][1];
			expectedpayoff.put(target, Math.round(tmpexp*100)/100.00);

		}


		return expectedpayoff;
	}

	
	

	private static double[] defenderExpectedPayoffsWithMapping(
			ArrayList<Integer> tiedtargets, double[] coverage,
			int[][] gamedata, int [][] origpmat, 
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback) {

		double[] expectedpayoff = new double[tiedtargets.size()];

		double[] c = new double[tiedtargets.size()];



		for(int t=0; t<expectedpayoff.length; t++)
		{
			int target = tiedtargets.get(t);

			for(int path=0; path<coverage.length; path++)
			{
				c[t] += origpmat[target][path]*coverage[path];


				/*double defreward = origpmat[target][path]*coverage[path]*gamedata[target][0];
					double defpenalty = (1-origpmat[target][path])*(1-coverage[path])*gamedata[target][1];
					sum = sum + defpenalty + defreward;*/
			}


			expectedpayoff[t] = c[t]*gamedata[originalmapback.get(target)][0] + (1-c[t])*gamedata[originalmapback.get(target)][1];
			expectedpayoff[t] = Math.round(expectedpayoff[t]*100)/100.00;

		}


		return expectedpayoff;
	}
	
	private static double[] defenderExpectedSTPayoffsWithMapping(
			ArrayList<Integer> tiedtargets, double[] coverage,
			HashMap<Integer, SuperTarget> currentst, int [][] origpmat, 
			HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {

		double[] expectedpayoff = new double[tiedtargets.size()];

		double[] c = new double[tiedtargets.size()];



		for(int t=0; t<expectedpayoff.length; t++)
		{
			int target = tiedtargets.get(t);

			for(int path=0; path<coverage.length; path++)
			{
				c[t] += origpmat[target][path]*coverage[path];


				/*double defreward = origpmat[target][path]*coverage[path]*gamedata[target][0];
					double defpenalty = (1-origpmat[target][path])*(1-coverage[path])*gamedata[target][1];
					sum = sum + defpenalty + defreward;*/
			}


			expectedpayoff[t] = c[t]*currentst.get(mapback.get(target)).defenderreward + (1-c[t])*currentst.get(mapback.get(target)).defenderpenalty;
			expectedpayoff[t] = Math.round(expectedpayoff[t]*100)/100.00;

		}


		return expectedpayoff;
	}


	private static double[] defenderExpectedPayoffs(
			ArrayList<Integer> tiedtargets, double[] coverage,
			int[][] gamedata, int [][] origpmat) {

		double[] expectedpayoff = new double[tiedtargets.size()];

		double[] c = new double[tiedtargets.size()];



		for(int t=0; t<expectedpayoff.length; t++)
		{
			int target = tiedtargets.get(t);

			for(int path=0; path<coverage.length; path++)
			{
				c[t] += origpmat[target][path]*coverage[path];


				/*double defreward = origpmat[target][path]*coverage[path]*gamedata[target][0];
					double defpenalty = (1-origpmat[target][path])*(1-coverage[path])*gamedata[target][1];
					sum = sum + defpenalty + defreward;*/
			}


			expectedpayoff[t] = c[t]*gamedata[target][0] + (1-c[t])*gamedata[target][1];
			expectedpayoff[t]= Math.round(expectedpayoff[t]*100)/100.0;

		}


		return expectedpayoff;
	}
	
	
	private static double[] defenderExpectedPayoffs(
			ArrayList<Integer> tiedtargets, double[] coverage,
			HashMap<Integer, TargetNode> targetmaps, int [][] origpmat) {

		double[] expectedpayoff = new double[tiedtargets.size()];

		double[] c = new double[tiedtargets.size()];



		for(int t=0; t<expectedpayoff.length; t++)
		{
			int target = tiedtargets.get(t);

			for(int path=0; path<coverage.length; path++)
			{
				c[t] += origpmat[target][path]*coverage[path];


				/*double defreward = origpmat[target][path]*coverage[path]*gamedata[target][0];
					double defpenalty = (1-origpmat[target][path])*(1-coverage[path])*gamedata[target][1];
					sum = sum + defpenalty + defreward;*/
			}


			expectedpayoff[t] = c[t]*targetmaps.get(target).defenderreward + (1-c[t])*targetmaps.get(target).defenderpenalty;
			expectedpayoff[t]= Math.round(expectedpayoff[t]*100)/100.0;

		}


		return expectedpayoff;
	}


	public static double expectedAttackerPayoff(int target, int[][] p,
			double[] probdistribution, int[][] gamedata, HashMap<Integer, Integer> map) {


		double expectedpayoffs = 0.0;

		double c = 0;
		//for(int target=0; target<expectedpayoffs.length; target++)
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{
				c += p[map.get(target)][jointschedule]* probdistribution[jointschedule];
				//System.out.println(" attk: target "+ target+ ", probj["+jointschedule+"] "+ probdistribution[jointschedule]+ ", coverage c:"+c);

			}

			double x = c*gamedata[target][3];
			double y = (1-c)*gamedata[target][2];
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs = Math.round((x+y)*100)/100.00; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}
	
	
	public static double expectedAttackerSTPayoff(int target, int[][] p,
			double[] probdistribution, HashMap<Integer,SuperTarget> currentst,
			 HashMap<Integer, Integer> map) {


		double expectedpayoffs = 0.0;

		double c = 0;
		//for(int target=0; target<expectedpayoffs.length; target++)
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{
				c += p[map.get(target)][jointschedule]* probdistribution[jointschedule];
				//System.out.println(" attk: target "+ target+ ", probj["+jointschedule+"] "+ probdistribution[jointschedule]+ ", coverage c:"+c);

			}

			double x = c*currentst.get(target).attackerpenalty;
			double y = (1-c)*currentst.get(target).attackerreward;
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs = Math.round((x+y)*100)/100.00; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}



	private static double expectedDefenderPayoff(int target, int[][] p,
			double[] probdistribution, int[][] gamedata, HashMap<Integer, Integer> map) {


		double expectedpayoffs = 0.0;

		double c = 0;
		//for(int target=0; target<expectedpayoffs.length; target++)
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{
				c += p[map.get(target)][jointschedule]* probdistribution[jointschedule];
				//System.out.println(" attk: target "+ target+ ", probj["+jointschedule+"] "+ probdistribution[jointschedule]+ ", coverage c:"+c);

			}

			double x = c*gamedata[target][0];
			double y = (1-c)*gamedata[target][1];
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs = Math.round((x+y)*100)/100.00; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}



	private static double[] expectedAttackerPayoffsWithMapping(int[][] origpmat,
			double[] coverage, int[][] gamedata, 
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback) {


		double[] expectedpayoffs = new double[origpmat.length];

		double[] c = new double[origpmat.length];
		for(int target=0; target<expectedpayoffs.length; target++)
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<coverage.length; jointschedule++)
			{
				c[target] += origpmat[target][jointschedule]* coverage[jointschedule];

			}

			double x = c[target]*gamedata[originalmapback.get(target)][3];
			double y = (1-c[target])*gamedata[originalmapback.get(target)][2];
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs[target] = Math.round((x+y)*100)/100.00; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}
	
	
	private static double[] expectedAttackerSTPayoffsWithMapping(int[][] pmat,
			double[] coverage, HashMap<Integer, SuperTarget> currentst, 
			HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		double[] expectedpayoffs = new double[pmat.length];

		double[] c = new double[pmat.length];
		for(int stid: currentst.keySet())
		{
			//System.out.println("\n\nTarget "+ target);


			int tid = map.get(stid);

			for(int jointschedule=0; jointschedule<coverage.length; jointschedule++)
			{
				c[tid] += pmat[tid][jointschedule]* coverage[jointschedule];

			}

			//double x = c[tid]*gamedata[originalmapback.get(target)][3];
			double x = c[tid]*currentst.get(stid).attackerpenalty;
			//double y = (1-c[target])*gamedata[originalmapback.get(target)][2];
			double y = (1-c[tid])*currentst.get(stid).attackerreward;
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs[tid] = Math.round((x+y)*100)/100.00; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}
	
	
	
	private static HashMap<Integer,Double> expectedAttackerPayoffsWithHashMap(int[][] origpmat,
			double[] coverage, int[][] gamedata, 
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback) {

		HashMap<Integer,Double> exppayoff = new HashMap<Integer,Double>();
		//double[] expectedpayoffs = new double[origpmat.length];

		double[] c = new double[origpmat.length];
		for(int target: originalmap.values())
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<coverage.length; jointschedule++)
			{
				c[target] += origpmat[target][jointschedule]* coverage[jointschedule];

			}

			double x = c[target]*gamedata[originalmapback.get(target)][3];
			double y = (1-c[target])*gamedata[originalmapback.get(target)][2];
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			//expectedpayoffs[target] = Math.round((x+y)*100)/100.00; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);
			exppayoff.put(originalmapback.get(target),  Math.round((x+y)*100)/100.00);

		}


		return exppayoff;
	}


	private static double[] expectedAttackerPayoffs(int[][] origpmat,
			double[] coverage, HashMap<Integer, TargetNode> targetmaps) {


		double[] expectedpayoffs = new double[origpmat.length];

		double[] c = new double[origpmat.length];
		for(int target=0; target<expectedpayoffs.length; target++)
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<coverage.length; jointschedule++)
			{
				c[target] += origpmat[target][jointschedule]* coverage[jointschedule];

			}

			double x = c[target]*targetmaps.get(target).attackerpenalty;
			double y = (1-c[target])*targetmaps.get(target).attackerreward;
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs[target] =  Math.round(((x+y)*100))/100.0; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}
	
	private static double[] expectedAttackerPayoffs(int[][] origpmat,
			double[] coverage, int[][] gamedata) {


		double[] expectedpayoffs = new double[origpmat.length];

		double[] c = new double[origpmat.length];
		for(int target=0; target<expectedpayoffs.length; target++)
		{
			//System.out.println("\n\nTarget "+ target);




			for(int jointschedule=0; jointschedule<coverage.length; jointschedule++)
			{
				c[target] += origpmat[target][jointschedule]* coverage[jointschedule];

			}

			double x = c[target]*gamedata[target][3];
			double y = (1-c[target])*gamedata[target][2];
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs[target] =  Math.round(((x+y)*100))/100.0; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}
	
	
	
	public static int[][] makeOrigPMatWOMap(int[][] p,
			ArrayList<ArrayList<Integer>> pathseq, 
			List<ArrayList<Integer>> jset, int nTargets, ArrayList<TargetNode> dominatednodes,
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback, 
			ArrayList<TargetNode> targets) throws Exception {



		int[][] origpmat = new int[nTargets][jset.size()];

		int jindex = 0;
		for(ArrayList<Integer> j: jset)
		{

			ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
			for(Integer jpath: j)
			{
				paths.add(pathseq.get(jpath));
			}


			// 0 11 18 2 5 9 13 19 6 7 12 14 15 16 20 21 22 24 

			for(ArrayList<Integer> path: paths)
			{
				//printNodesWithNeighborsAndPath(dominatednodes, targets);
				/*System.out.println("Considering path :");
				printGreedyPath(path);
				System.out.println();*/

				for(int targetindex=0; targetindex<path.size()-1; targetindex++)
				{

					int target = path.get(targetindex);
					int nexttarget = path.get(targetindex+1);
					
					if(target==0 && nexttarget==4)
					{
						System.out.println();
					}

					TargetNode targetnode = getTargetNode(target, targets);
					TargetNode nexttargetnode = getTargetNode(nexttarget, targets);

					origpmat[target][jindex] = 1; 
					origpmat[nexttarget][jindex] = 1; 


					ArrayList<TargetNode> inbetweennodes = targetnode.getPath(nexttargetnode);
					try{
						for(TargetNode t: inbetweennodes)
						{
							//System.out.println("Nnode "+t.getTargetid() + " is used as inbetween node in joint schedule "+ jindex);
							//origpmat[originalmap.get(t.getTargetid())][jindex] = 1; // added mapping : modified for contraction 
							origpmat[t.getTargetid()][jindex] = 1; // added mapping : modified for instant contraction 
						}
					}
					catch(Exception ex)
					{
						//printNodesWithNeighborsAndPath(dominatednodes, targets);
						throw new Exception("Mapping problem in orig map");
						
					}


				}

			}
			jindex++;
		}


		return origpmat;
	}



	public static int[][] makeSuperOrigPMatWOMap(int[][] p,
			ArrayList<ArrayList<Integer>> pathseq, 
			List<ArrayList<Integer>> jset, int nTargets,
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback, 
			HashMap<Integer,TargetNode> targetmaps, HashMap<Integer,SuperTarget> currentst, HashMap<Integer,ArrayList<Integer>> stpaths) throws Exception {



		int[][] origpmat = new int[nTargets][jset.size()];

		int jindex = 0;
		for(ArrayList<Integer> j: jset)
		{

			ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
			for(Integer jpath: j)
			{
				paths.add(pathseq.get(jpath));
			}


			// 0 11 18 2 5 9 13 19 6 7 12 14 15 16 20 21 22 24 

			for(ArrayList<Integer> path: paths)
			{
				//printNodesWithNeighborsAndPath(dominatednodes, targets);
				/*System.out.println("Considering path :");
				printGreedyPath(path);
				System.out.println();*/

				for(int targetindex=0; targetindex<path.size()-1; targetindex++)
				{

					// get the supertarget
					int st = path.get(targetindex);
					
					// for all the nodes in the supertarget assign 1
					
					for(TargetNode t: currentst.get(st).nodes.values())
					{
						
						origpmat[t.getTargetid()][jindex] = 1;
					}
					
					
					/*if(currentst.get(st).nodes.size()==1)
					{
						for(TargetNode t: currentst.get(st).nodes.values())
						{
							
							origpmat[t.getTargetid()][jindex] = 1;
						}
						
					}
					else if(currentst.get(st).nodes.size() > 1)
					{
						ArrayList<Integer> pathnodes = stpaths.get(st);
						for(TargetNode t: currentst.get(st).nodes.values())
						{
							if(pathnodes.contains(t.getTargetid()))
							{
							
								origpmat[t.getTargetid()][jindex] = 1;
							}
						}
						
					}*/
					
					
				}

			}
			jindex++;
		}


		return origpmat;
	}



	private static int[][] makeOrigPMat(int[][] p,
			ArrayList<ArrayList<Integer>> pathseq, 
			List<ArrayList<Integer>> jset, int nTargets, ArrayList<TargetNode> dominatednodes,
			HashMap<Integer,Integer> originalmap, HashMap<Integer,Integer> originalmapback, ArrayList<TargetNode> targets) {



		int[][] origpmat = new int[nTargets][jset.size()];

		int jindex = 0;
		for(ArrayList<Integer> j: jset)
		{

			ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
			for(Integer jpath: j)
			{
				paths.add(pathseq.get(jpath));
			}


			for(ArrayList<Integer> path: paths)
			{
				for(int targetindex=0; targetindex<path.size()-1; targetindex++)
				{

					int target = path.get(targetindex);
					int nexttarget = path.get(targetindex+1);

					TargetNode targetnode = getTargetNode(target, targets);
					TargetNode nexttargetnode = getTargetNode(nexttarget, targets);

					origpmat[originalmap.get(target)][jindex] = 1; //mapping added 
					origpmat[originalmap.get(nexttarget)][jindex] = 1; //mapping added


					ArrayList<TargetNode> inbetweennodes = targetnode.getPath(nexttargetnode);
					for(TargetNode t: inbetweennodes)
					{

						if(isTInDomNodes(t, dominatednodes))
						{

							//System.out.println("Nnode "+t.getTargetid() + " is used as inbetween node in joint schedule "+ jindex);

							int tid = t.getTargetid();

							int mappedid = originalmap.get(tid);


							origpmat[mappedid][jindex] = 1; // added mapping : modified for contraction 
							//origpmat[t.getTargetid()][jindex] = 1; // added mapping : modified for instant contraction 
						}
					}


				}

			}
			jindex++;
		}


		return origpmat;
	}

	private static boolean isTInDomNodes(TargetNode t,
			ArrayList<TargetNode> dominatednodes) {


		for(TargetNode n: dominatednodes)
		{
			if(t.getTargetid()==n.getTargetid())
				return true;
		}


		return false;
	}

	public static TargetNode getTargetNode(int targetid, ArrayList<TargetNode> targets) {

		for(TargetNode t: targets)
		{
			if(t.getTargetid()==targetid)
			{
				return t;
			}
		}

		return null;


	}



	private static void makeOnlyPathSeq(ArrayList<ArrayList<Integer>> pathseq
			,ArrayList<TargetNode> goals, int pathcounter, int nTargets, 
			ArrayList<TargetNode> targets) {





		//System.out.println("Total path "+ pathcounter);

		//int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			//Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			tmppathseq.add(tmpgoal.getTargetid());
			while(tmpstart.parent!=null)
			{
				//Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmppathseq.add(tmpstart.getTargetid());
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					//tmppathseq.add(tmpgoal.getTargetid());
					//break;
				}

			}
			pathindex++;
			pathseq.add(tmppathseq);
		}

		//return Pmat;





	}



	private static void makePathSeq(ArrayList<ArrayList<Integer>> pathseq
			,ArrayList<TargetNode> goals, int pathcounter, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback,
			ArrayList<TargetNode> targets) {


		int icount =0;
		//map = new HashMap<Integer, Integer>();
		//mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		//System.out.println("Total path "+ pathcounter);

		//int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			//Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			tmppathseq.add(tmpgoal.getTargetid());
			while(tmpstart.parent!=null)
			{
				//Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmppathseq.add(tmpstart.getTargetid());
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					//tmppathseq.add(tmpgoal.getTargetid());
					//break;
				}

			}
			pathindex++;
			pathseq.add(tmppathseq);
		}

		//return Pmat;





	}
	
	
	
	public static void makePathSeqSrcDest(ArrayList<ArrayList<Integer>> pathseq
			,ArrayList<TargetNode> goals, int pathcounter, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback,
			ArrayList<TargetNode> targets) {


		int icount =0;
		//map = new HashMap<Integer, Integer>();
		//mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		//System.out.println("Total path "+ pathcounter);

		//int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			
			ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
			buildpathSeqq(n, tmppathseq);
			
			pathindex++;
			pathseq.add(tmppathseq);
		}

		//return Pmat;





	}
	
	
	

	public static void makeClusterPathSeq(ArrayList<TargetNode> goals, ArrayList<Integer> spath) {

		int pathindex = 0;
		for(TargetNode n: goals)
		{
			
			//ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
			buildpathSeqq(n, spath);
			
			
		}

	}




	private static ArrayList<Integer> buildpathSeqq(TargetNode n, ArrayList<Integer> tmppathseq) {
		
		if(n.parent==null)
			return null;
		//System.out.println("parent  "+ n.parent.getTargetid());
		buildpathSeqq(n.parent, tmppathseq);
		tmppathseq.add(n.getTargetid());
		//System.out.println("Adding "+ n.getTargetid());
		
		
		
		return tmppathseq;
	}

	private static void makeSlavePathSeq(ArrayList<ArrayList<Integer>> pathseq
			,ArrayList<TargetNode> goals
			) {


		int icount =0;
		//map = new HashMap<Integer, Integer>();
		//mapback = new HashMap<Integer, Integer>();




		//System.out.println("Total path "+ pathcounter);

		//int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			ArrayList<Integer> tmppathseq = new ArrayList<Integer>();
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			//Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			tmppathseq.add(tmpgoal.getTargetid());
			while(tmpstart.parent!=null)
			{
				//Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmppathseq.add(tmpstart.getTargetid());
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					//tmppathseq.add(tmpgoal.getTargetid());
					//break;
				}

			}
			pathindex++;
			pathseq.add(tmppathseq);
		}

		//return Pmat;





	}

	public static ArrayList<TargetNode> generatePaths(double dmax, int[][] gamedata, ArrayList<TargetNode> targets) throws Exception
	{
		ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getPathsByBFS(dmax,paths, gamedata, targets);
		//printPathMat(pathsmat);

		return path;

	}
	
	
	public static ArrayList<TargetNode> generatePathsWithSrcDest(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			HashMap<Integer, TargetNode> targetmaps, int base, int dest) throws Exception
	{
		//ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getPathsByBFSWithSrcDest(dmax, gamedata, targets, targetmaps, base, dest);
		//printPathMat(pathsmat);

		return path;

	}
	

	public static ArrayList<TargetNode> generateOnePathWithSrcDest(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			HashMap<Integer, TargetNode> targetmaps, int base, int dest, int [] distallocation) throws Exception
	{
		//ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getOnePathByBFSWithSrcDest(dmax, gamedata, targets, targetmaps, base, dest, distallocation);
		//printPathMat(pathsmat);

		return path;

	}
	
	
	
	


	public static ArrayList<TargetNode> generatePathsGreedy2SrcDest(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes, int srcid, int destid) throws Exception
			{
		ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getPathsByBFSGreedy2SrcDest(dmax,paths, gamedata, targets, currenttargets, nRes, srcid,destid);
		//printPathMat(pathsmat);

		return path;

			}
	


	public static ArrayList<TargetNode> generatePathsGreedy2(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes) throws Exception
			{
		ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getPathsByBFSGreedy2(dmax,paths, gamedata, targets, currenttargets, nRes);
		//printPathMat(pathsmat);

		return path;

			}
	
	

	public static ArrayList<ArrayList<Integer>> generatePathsGreedy3WithAPSP(int baseid, double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes) throws Exception
			{


		TargetNode base = getTargetNode(baseid, targets);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		for(int i=1; i<=targets.size(); i++)
		{
			map.put(targets.get(i-1).getTargetid(), i);
			mapback.put(i, targets.get(i-1).getTargetid());
			graphint.add(map.get(targets.get(i-1).getTargetid()));
		}
		
		
		
		int[][] adjacencymatrix = new int[targets.size()+1][targets.size()+1];
		makeAdjacencyMatrix(adjacencymatrix , targets, targets.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(targets.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPMatrixZero(apsp, targets, targets.size(), map, mapback);
		
       // System.out.println("FInding  "+ apsp[map.get(11)][map.get(82)]);
     //   System.out.println("FInding  "+ apsp[map.get(82)][map.get(11)]);
		
		
		


		for(TargetNode dest: targets)
		{
			if(dest.getTargetid()!=0)
			{

				//System.out.println("FInding shortest dist for target "+ dest.getTargetid());
				ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
				ArrayList<TargetNode>  pnodes = new ArrayList<TargetNode>();
				if(base.getNeighbors().contains(dest))
				{
					pnodes = base.getPath(dest);
					for(int k=0; k<pnodes.size(); k++)
					{

						//System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
						//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

					}

				}
				else
				{
					//double distcovered1 = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
					//System.out.print("dist covered "+ distcovered1+"\n");


					int src = base.getTargetid();
					int des = dest.getTargetid();



					double distcovered = apsp[map.get(src)][map.get(des)];
					System.out.print("dist covered "+ distcovered+"\n");

					if(distcovered<=dmax/2)
					{
						ArrayList<Integer>	tmppathnodes = allPairShortestPath.getPath(src, des, map, mapback);

						for(int k=0; k<tmppathnodes.size(); k++)
						{
							pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
						}
					}
					else
						continue;


				}
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				tmppath.add(base.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
			//	System.out.print(dest.getTargetid()+"\n");
				tmppath.add(dest.getTargetid());
				System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				System.out.print("complete path : \n");
				for(int k=0; k<tmppath.size(); k++)
				{

					System.out.print(tmppath.get(k)+"->");
				}
				paths.add(tmppath);
				System.out.print("\n");
			}



		}
		return paths;
			}

	
	
	

	public static ArrayList<ArrayList<Integer>> generatePathsGreedy3WithAPSP(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes) throws Exception
			{


		TargetNode base = getTargetNode(0, targets);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		for(int i=1; i<=targets.size(); i++)
		{
			map.put(targets.get(i-1).getTargetid(), i);
			mapback.put(i, targets.get(i-1).getTargetid());
			graphint.add(map.get(targets.get(i-1).getTargetid()));
		}
		
		
		
		int[][] adjacencymatrix = new int[targets.size()+1][targets.size()+1];
		makeAdjacencyMatrix(adjacencymatrix , targets, targets.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(targets.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPMatrixZero(apsp, targets, targets.size(), map, mapback);
		
       // System.out.println("FInding  "+ apsp[map.get(11)][map.get(82)]);
     //   System.out.println("FInding  "+ apsp[map.get(82)][map.get(11)]);
		
		
		


		for(TargetNode dest: targets)
		{
			if(dest.getTargetid()!=0)
			{

				//System.out.println("FInding shortest dist for target "+ dest.getTargetid());
				ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
				ArrayList<TargetNode>  pnodes = new ArrayList<TargetNode>();
				
				
				int src = base.getTargetid();
				int des = dest.getTargetid();



				double distcovered = apsp[map.get(src)][map.get(des)];
				if(base.getNeighbors().contains(dest) && distcovered<=dmax/2)
				{
					
					
					
					pnodes = base.getPath(dest);
					for(int k=0; k<pnodes.size(); k++)
					{

						//System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
						//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

					}

				}
				else
				{
					//double distcovered1 = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
					//System.out.print("dist covered "+ distcovered1+"\n");


					src = base.getTargetid();
					des = dest.getTargetid();



					distcovered = apsp[map.get(src)][map.get(des)];
					//System.out.print("dist covered "+ distcovered+"\n");

					if(distcovered<=dmax/2)
					{
						ArrayList<Integer>	tmppathnodes = allPairShortestPath.getPath(src, des, map, mapback);

						for(int k=0; k<tmppathnodes.size(); k++)
						{
							pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
						}
					}
					else
						continue;


				}
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				tmppath.add(base.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
			//	System.out.print(dest.getTargetid()+"\n");
				tmppath.add(dest.getTargetid());
				//System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				//System.out.print("complete path : \n");
				/*for(int k=0; k<tmppath.size(); k++)
				{

					//System.out.print(tmppath.get(k)+"->");
				}*/
				paths.add(tmppath);
				//System.out.print("\n");
			}



		}
		return paths;
			}
	
	
	
	public static ArrayList<Integer> pathForAT(double dmax, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes, int attackedtarget) throws Exception
			{


		TargetNode base = getTargetNode(0, targets);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		for(int i=1; i<=targets.size(); i++)
		{
			map.put(targets.get(i-1).getTargetid(), i);
			mapback.put(i, targets.get(i-1).getTargetid());
			graphint.add(map.get(targets.get(i-1).getTargetid()));
		}
		
		
		
		int[][] adjacencymatrix = new int[targets.size()+1][targets.size()+1];
		makeAdjacencyMatrix(adjacencymatrix , targets, targets.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(targets.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPMatrixZero(apsp, targets, targets.size(), map, mapback);
		
       // System.out.println("FInding  "+ apsp[map.get(11)][map.get(82)]);
     //   System.out.println("FInding  "+ apsp[map.get(82)][map.get(11)]);
		
		
		TargetNode dest = getTargetNode(attackedtarget, targets);


		//for(TargetNode dest: targets)
		{
			if(dest.getTargetid()!=0)
			{

				//System.out.println("FInding shortest dist for target "+ dest.getTargetid());
				ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
				ArrayList<TargetNode>  pnodes = new ArrayList<TargetNode>();
				
				
				int src = base.getTargetid();
				int des = dest.getTargetid();



				double distcovered = apsp[map.get(src)][map.get(des)];
				if(base.getNeighbors().contains(dest) && distcovered<=dmax/2)
				{
					
					
					
					pnodes = base.getPath(dest);
					for(int k=0; k<pnodes.size(); k++)
					{

						//System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
						//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

					}

				}
				else
				{
					//double distcovered1 = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
					//System.out.print("dist covered "+ distcovered1+"\n");


					src = base.getTargetid();
					des = dest.getTargetid();



					distcovered = apsp[map.get(src)][map.get(des)];
					//System.out.print("dist covered "+ distcovered+"\n");

					if(distcovered<=dmax/2)
					{
						ArrayList<Integer>	tmppathnodes = allPairShortestPath.getPath(src, des, map, mapback);

						for(int k=0; k<tmppathnodes.size(); k++)
						{
							pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
						}
					}
					


				}
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				tmppath.add(base.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
			//	System.out.print(dest.getTargetid()+"\n");
				tmppath.add(dest.getTargetid());
				//System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				//System.out.print("complete path : \n");
				/*for(int k=0; k<tmppath.size(); k++)
				{

					//System.out.print(tmppath.get(k)+"->");
				}*/
				paths.add(tmppath);
				//System.out.print("\n");
			}



		}
		return paths.get(0);
			}
	
	
	
	public static ArrayList<ArrayList<Integer>> generatePathsForSuperTargetsAPSP(double dmax, HashMap<Integer, SuperTarget> sts,
			HashMap<Integer,TargetNode> targetmaps, int nRes, HashMap<Integer,Double> dstravel) throws Exception
			{


		SuperTarget base = sts.get(0);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		
		int i = 1;
		for(Integer stid: sts.keySet())
		{
			map.put(stid, i);
			mapback.put(i, stid);
			//graphint.add(map.get(targets.get(i-1).getTargetid()));
			i++;
		}
		
		
		
		int[][] adjacencymatrix = new int[sts.size()+1][sts.size()+1];
		makeAdjacencyMatrixST(adjacencymatrix , sts, sts.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(sts.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPSTMatrixZero(apsp, sts, sts.size(), map, mapback);
		
       // System.out.println("FInding  "+ apsp[map.get(11)][map.get(82)]);
     //   System.out.println("FInding  "+ apsp[map.get(82)][map.get(11)]);
		
		
		


        for(SuperTarget dest: sts.values())
        {
        	if(dest.stid!=0)
        	{

        		//System.out.println("FInding shortest dist for target "+ dest.stid);
        		ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
        		ArrayList<SuperTarget>  pnodes = new ArrayList<SuperTarget>();
        		if(base.neighbors.containsKey(dest.stid))
				{
					//pnodes = base.getPath(dest);
        			int src = base.stid;
            		int des = dest.stid;


            		//TODO consider distance, intra cluster
            		
            		
            		
            		double distcovered = apsp[map.get(src)][map.get(des)]+dstravel.get(dest.stid);
            		//System.out.print("dist covered "+ distcovered+"\n");
            		
            		if(distcovered > dmax/2)
            			continue;
					/*for(int k=0; k<pnodes.size(); k++)
					{

						//System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
						//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

					}*/

				}
        		else
        		{

        			
        		//double distcovered1 = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
        		//System.out.print("dist covered "+ distcovered1+"\n");


        		int src = base.stid;
        		int des = dest.stid;


        		//TODO consider distance, intra cluster
        		
        		if(!dstravel.containsKey(dest.stid)) // no path exists
        		{
        			continue;
        		}
        		
        		
        		double distcovered = apsp[map.get(src)][map.get(des)]+dstravel.get(dest.stid);
        		//System.out.print("dist covered "+ distcovered+"\n");

        		if(distcovered<=dmax/2)
        		{
        			ArrayList<Integer>	tmppathnodes = allPairShortestPath.getPath(src, des, map, mapback);

        			for(int k=0; k<tmppathnodes.size(); k++)
        			{
        				pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
        			}
        		}
        		else
        			continue;
        		
        		
        		//throw new Exception("Base to not neighbor for initial set of paths **********8");

        		}

        		ArrayList<Integer> tmppath = new ArrayList<Integer>();
        		//System.out.print("\n0->");
        		tmppath.add(base.stid);
        		for(int k=0; k<pathnodes.size(); k++)
        		{
        			tmppath.add(pathnodes.get(pathnodes.size()-k-1));
        			//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
        		}
        		//	System.out.print(dest.getTargetid()+"\n");
        		tmppath.add(dest.stid);
        		//System.out.print("\n");
        		/**
        		 * make rev path
        		 */
        		for(int j=tmppath.size()-2; j>=0; j--)
        		{
        			tmppath.add(tmppath.get(j));
        		}
        		//System.out.print("complete path : \n");
        		/*for(int k=0; k<tmppath.size(); k++)
        		{

        			//System.out.print(tmppath.get(k)+"->");
        		}*/
        		paths.add(tmppath);
        		//System.out.print("\n");
        	}



        }
		return paths;
			}
	
	
	
	
	
	public static void addSuperTargetsAPSP(double dmax, HashMap<Integer, SuperTarget> sts,
			HashMap<Integer,TargetNode> targetmaps, int nRes, HashMap<Integer,Double> dstravel,
			ArrayList<ArrayList<Integer>> newpathseq, ArrayList<Integer> currentattackedsupertargets) throws Exception
			{


		SuperTarget base = sts.get(0);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		
		int i = 1;
		for(Integer stid: sts.keySet())
		{
			map.put(stid, i);
			mapback.put(i, stid);
			//graphint.add(map.get(targets.get(i-1).getTargetid()));
			i++;
		}
		
		
		
		int[][] adjacencymatrix = new int[sts.size()+1][sts.size()+1];
		makeAdjacencyMatrixST(adjacencymatrix , sts, sts.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(sts.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPSTMatrixZero(apsp, sts, sts.size(), map, mapback);
		
       // System.out.println("FInding  "+ apsp[map.get(11)][map.get(82)]);
     //   System.out.println("FInding  "+ apsp[map.get(82)][map.get(11)]);
		
		
		


        for(Integer supertarget: currentattackedsupertargets)
        {
        	
        	SuperTarget dest = sts.get(supertarget);
        	
        	if(dest.stid!=0)
        	{

        		//System.out.println("FInding shortest dist for target "+ dest.stid);
        		ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
        		ArrayList<SuperTarget>  pnodes = new ArrayList<SuperTarget>();
        		if(base.neighbors.containsKey(dest.stid))
				{
					//pnodes = base.getPath(dest);
        			int src = base.stid;
            		int des = dest.stid;


            		//TODO consider distance, intra cluster
            		
            		
            		
            		double distcovered = apsp[map.get(src)][map.get(des)]+dstravel.get(dest.stid);
            		//System.out.print("dist covered "+ distcovered+"\n");
            		
            		if(distcovered > dmax/2)
            			continue;
					/*for(int k=0; k<pnodes.size(); k++)
					{

						//System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
						//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

					}*/

				}
        		else
        		{

        			
        		//double distcovered1 = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
        		//System.out.print("dist covered "+ distcovered1+"\n");


        		int src = base.stid;
        		int des = dest.stid;


        		//TODO consider distance, intra cluster
        		
        		if(!dstravel.containsKey(dest.stid)) // no path exists
        		{
        			continue;
        		}
        		
        		
        		double distcovered = apsp[map.get(src)][map.get(des)]+dstravel.get(dest.stid);
        		//System.out.print("dist covered "+ distcovered+"\n");

        		if(distcovered<=dmax/2)
        		{
        			ArrayList<Integer>	tmppathnodes = allPairShortestPath.getPath(src, des, map, mapback);

        			for(int k=0; k<tmppathnodes.size(); k++)
        			{
        				pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
        			}
        		}
        		else
        			continue;
        		
        		
        		//throw new Exception("Base to not neighbor for initial set of paths **********8");

        		}

        		ArrayList<Integer> tmppath = new ArrayList<Integer>();
        		//System.out.print("\n0->");
        		tmppath.add(base.stid);
        		for(int k=0; k<pathnodes.size(); k++)
        		{
        			tmppath.add(pathnodes.get(pathnodes.size()-k-1));
        			//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
        		}
        		//	System.out.print(dest.getTargetid()+"\n");
        		tmppath.add(dest.stid);
        		//System.out.print("\n");
        		/**
        		 * make rev path
        		 */
        		for(int j=tmppath.size()-2; j>=0; j--)
        		{
        			tmppath.add(tmppath.get(j));
        		}
        		//System.out.print("complete path : \n");
        		/*for(int k=0; k<tmppath.size(); k++)
        		{

        			//System.out.print(tmppath.get(k)+"->");
        		}*/
        		newpathseq.add(tmppath);
        		//System.out.print("\n");
        	}



        }
		//return paths;
			}
	
	
	
	
	
	
	public static boolean getPathNodess(ArrayList<Integer> pathnodes, ArrayList<TargetNode> targets, int srcid, int destid, 
			double dmax, HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, AllPairShortestPath allPairShortestPath, int[][] apsp)
	{
		//System.out.println("FInding shortest dist for target "+ dest.getTargetid());
		//ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
		ArrayList<TargetNode>  pnodes = new ArrayList<TargetNode>();
		
		TargetNode src = getTargetNode(srcid, targets);
		TargetNode dest = getTargetNode(destid, targets);
		
		
		double dist=0;
		if(src.getNeighbors().contains(dest))
		{
			pnodes = src.getPath(dest);
			dist = src.getDistance(dest);
			if(dist<=(dmax/2))
			{
				/*for(int i=0; i<pnodes.size(); i++)
				{
					pathnodes.add(pnodes.get(i).getTargetid());
				}*/
				return true;
			}
			else
				return false;
			/*for(int k=0; k<pnodes.size(); k++)
			{

				//System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
				//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

			}*/

		}
		else
		{
			//double distcovered1 = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
			//System.out.print("dist covered "+ distcovered1+"\n");


			//int srcid = src.getTargetid();
			//int destid = dest.getTargetid();



			double distcovered = apsp[map.get(srcid)][map.get(destid)];
			System.out.print("dist covered "+ distcovered+"\n");

			if(distcovered<=dmax/2)
			{
				ArrayList<Integer>	tmppathnodes = allPairShortestPath.getPath(srcid, destid, map, mapback);

				for(int k=0; k<tmppathnodes.size(); k++)
				{
					pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
				}
				return true;
			}
			else
			{
				return false;
			}
			
			
			


		}
		//return false;
	
	}
	
	


	public static ArrayList<ArrayList<Integer>> generatePathsGreedy3WithSrcDest(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes, int srcid, int destid) throws Exception
			{


		TargetNode base = getTargetNode(srcid, targets);
		TargetNode dest = getTargetNode(destid, targets);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		for(int i=1; i<=targets.size(); i++)
		{
			map.put(targets.get(i-1).getTargetid(), i);
			mapback.put(i, targets.get(i-1).getTargetid());
			graphint.add(map.get(targets.get(i-1).getTargetid()));
		}
		
		
		
		int[][] adjacencymatrix = new int[targets.size()+1][targets.size()+1];
		makeAdjacencyMatrix(adjacencymatrix , targets, targets.size(), map, mapback);
		
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(targets.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);
        
        
        purifyAPSPMatrixZero(apsp, targets, targets.size(), map, mapback);
		
       // System.out.println("FInding  "+ apsp[map.get(11)][map.get(82)]);
     //   System.out.println("FInding  "+ apsp[map.get(82)][map.get(11)]);
		
		
		


		for(TargetNode mid: targets)
		{
			if(mid.getTargetid()!=srcid /*&& mid.getTargetid()!=destid*/)
			{
				ArrayList<Integer> pathnodes = new ArrayList<Integer>();
				boolean ok = getPathNodess(pathnodes, targets, srcid, mid.getTargetid(), dmax, map, mapback, allPairShortestPath, apsp);
				if(!ok)
					continue;
				
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				tmppath.add(base.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
			//	System.out.print(dest.getTargetid()+"\n");
				tmppath.add(mid.getTargetid());
				System.out.print("\n");
				
				/**
				 * make  path from mid to dest 
				 */
				
				
				
				pathnodes = new ArrayList<Integer>();
				ok = getPathNodess(pathnodes, targets,  mid.getTargetid(), destid, dmax, map, mapback, allPairShortestPath, apsp);
				if(!ok)
					continue;
				
				//ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				//tmppath.add(base.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
			//	System.out.print(dest.getTargetid()+"\n");
				tmppath.add(dest.getTargetid());
				
				
				
				
				/*for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}*/
				
				
				System.out.print("complete path : \n");
				for(int k=0; k<tmppath.size(); k++)
				{

					System.out.print(tmppath.get(k)+"->");
				}
				paths.add(tmppath);
				System.out.print("\n");
			}



		}
		return paths;
			}

	

	
	
	


	public static ArrayList<ArrayList<Integer>> generatePathsGreedy3(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			ArrayList<Integer> currenttargets, int nRes) throws Exception
			{


		TargetNode base = getTargetNode(0, targets);
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();


		for(TargetNode dest: targets)
		{
			if(dest.getTargetid()!=0)
			{

				System.out.println("FInding shortest dist for target "+ dest.getTargetid());
				ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
				ArrayList<TargetNode>  pnodes = new ArrayList<TargetNode>();
				double dist = 0;
				if(base.getNeighbors().contains(dest))
				{
					pnodes = base.getPath(dest);
					dist = base.getDistance(dest);
					if(dist>(dmax/2.0))
						continue;
					for(int k=0; k<pnodes.size(); k++)
					{

						System.out.print(pnodes.get(pnodes.size()-k-1).getTargetid()+"->");
						//pathnodes.add(pnodes.get(pnodes.size()-k-1).getTargetid());

					}

				}
				else
				{
					double distcovered = findShortestPathThrougGraphWDlimit(base, dest, targets, pathnodes, dmax/2);
					System.out.print("dist covered "+ distcovered+"\n");
					if(distcovered>(dmax/2))
						continue;
					
				}
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				tmppath.add(base.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
				//System.out.print(dest.getTargetid()+"\n");
				tmppath.add(dest.getTargetid());
				System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				System.out.print("complete path : \n");
				for(int k=0; k<tmppath.size(); k++)
				{

					System.out.print(tmppath.get(k)+"->");
				}
				paths.add(tmppath);
				System.out.print("\n");
			}



		}
		return paths;
			}



	public static ArrayList<TargetNode> generatePathsSlave(double dmax, int[][] gamedata, ArrayList<TargetNode> targets,
			int attackedtarget, int nRes, ArrayList<Integer> currenttargets) throws Exception
			{
		ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getPathsByBFSSlaveGreedy2(dmax,paths, gamedata, targets, attackedtarget, nRes, currenttargets);
		//printPathMat(pathsmat);

		return path;

			}




	public static ArrayList<TargetNode> generatePathsGreedy(double dmax, int[][] gamedata, ArrayList<TargetNode> targets) throws Exception
	{
		ArrayList<ArrayList<TargetNode>> paths = new ArrayList<ArrayList<TargetNode>>();

		ArrayList<TargetNode> path = getPathsByBFSGreedy(dmax,paths, gamedata, targets);
		//printPathMat(pathsmat);

		return path;

	}

	private static void printPathMat(int[][] pathsmat) {

		for(int i=0; i<pathsmat.length; i++)
		{
			for(int j=0; j<pathsmat[i].length; j++)
			{
				System.out.print(pathsmat[i][j]+" ");
			}
			System.out.println();

		}

	}


	private static ArrayList<TargetNode> getPathsByBFSGreedy(double dmax,
			ArrayList<ArrayList<TargetNode>> paths, int[][] gamedata, ArrayList<TargetNode> targets) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		int nTargets = targets.size(); 
		TargetNode start = new TargetNode(targets.get(0));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			//	System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//	System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//	System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==start.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet==dmax)
			{
				//	System.out.println("Adding node "+ node.getTargetid() +" to goals...");
				//printPath(node, node);
				pathcounter++;
				//System.out.println();
				goals.add(node);
			}
			else if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				ArrayList<TargetNode> succs = Expand(node, dmax, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}
	
	
	
	private static ArrayList<TargetNode> getPathsByBFSGreedy2(double dmax,
			ArrayList<ArrayList<TargetNode>> paths, int[][] gamedata,
			ArrayList<TargetNode> targets, ArrayList<Integer> currenttargets, int nRes) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		ArrayList<Integer> donetargets = new ArrayList<Integer>();
		int nTargets = targets.size(); 
		TargetNode start = new TargetNode(targets.get(0));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==start.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet);
				//printPath(node, node);

				//System.out.println();
				//check how many targets are covered

				int coveredcount = coveredCount(donetargets, node, currenttargets);
				/*System.out.print("done targets\n");
				for(Integer n: donetargets)
				{
					System.out.print(n+" ");
				}
				System.out.print("\n");*/

				
				goals.add(node);
				pathcounter++;

				/*if(coveredcount>0)
				{
					goals.add(node);
					pathcounter++;
				}
				else if(donetargets.size()==currenttargets.size())
				{
					goals.add(node);
					pathcounter++;
				}*/

				//if(pathcounter==500)
				//break;

				if(donetargets.size()==currenttargets.size())
					break;
			}
			else if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<TargetNode> succs = Expand(node, dmax, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}
	
	


	private static ArrayList<TargetNode> getPathsByBFSGreedy2SrcDest(double dmax,
			ArrayList<ArrayList<TargetNode>> paths, int[][] gamedata,
			ArrayList<TargetNode> targets, ArrayList<Integer> currenttargets, int nRes, int srcid, int destid) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		ArrayList<Integer> donetargets = new ArrayList<Integer>();
		int nTargets = targets.size(); 
		TargetNode start = new TargetNode(getTargetNode(srcid, targets));
		TargetNode dest = new TargetNode(getTargetNode(destid, targets));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			System.out.println("Polling from queue ");
			System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==dest.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet);
				//printPath(node, node);

				//System.out.println();
				//check how many targets are covered

				int coveredcount = coveredCount(donetargets, node, currenttargets);
				/*System.out.print("done targets\n");
				for(Integer n: donetargets)
				{
					System.out.print(n+" ");
				}
				System.out.print("\n");*/

				
				goals.add(node);
				pathcounter++;

				/*if(coveredcount>0)
				{
					goals.add(node);
					pathcounter++;
				}
				else if(donetargets.size()==currenttargets.size())
				{
					goals.add(node);
					pathcounter++;
				}*/

				//if(pathcounter==500)
				//break;

				if(donetargets.size()==currenttargets.size())
					break;
			}
			else if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<TargetNode> succs = Expand(node, dmax, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}



	private static ArrayList<TargetNode> getPathsByBFSSlaveGreedy2(double dmax,
			ArrayList<ArrayList<TargetNode>> paths, int[][] gamedata,
			ArrayList<TargetNode> targets, int attackedtarget, int nRes, ArrayList<Integer> currenttargets) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		//ArrayList<Integer> donetargets = new ArrayList<Integer>();
		int nTargets = targets.size(); 
		ArrayList<Integer> donetargets = new ArrayList<Integer>();
		TargetNode start = new TargetNode(targets.get(0));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==start.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet);
				//printPath(node, node);

				//System.out.println();
				//check how many targets are covered

				/*boolean covered = isAttackedtargetcovered(node, attackedtarget);
				if(covered)
				{
					System.out.println("covered atatcked target , pcount "+ pathcounter);
					goals.add(node);
					pathcounter++;
					break;
				}*/

				int coveredcount = coveredCount(donetargets, node, currenttargets);
				goals.add(node);
				pathcounter++;

				if(pathcounter==1000)
					break;


				/*if(pathcounter>(10*nRes))
					break;*/
			}
			else if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<TargetNode> succs = Expand(node, dmax, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}



	public static int coveredCount(ArrayList<Integer> donetargets,
			TargetNode node, ArrayList<Integer> currenttargets) 
	{

		int count =0;

		if(currenttargets.contains(node.getTargetid()) && !donetargets.contains(node.getTargetid()))
		{

			donetargets.add(node.getTargetid());
			count++;
		}




		while(node.parent!=null)
		{


			if(currenttargets.contains(node.getTargetid()) && !donetargets.contains(node.getTargetid()))
			{

				donetargets.add(node.getTargetid());
				count++;
			}
			node= node.parent;

		}
		return count;
	}
	
	
	public static int coveredCount(ArrayList<Integer> donetargets,
			TargetNode node, HashMap<Integer, TargetNode> currenttargets) 
	{

		int count =0;

		if(currenttargets.containsKey(node.getTargetid()) && !donetargets.contains(node.getTargetid()))
		{

			donetargets.add(node.getTargetid());
			count++;
		}




		while(node.parent!=null)
		{


			if(currenttargets.containsKey(node.getTargetid()) && !donetargets.contains(node.getTargetid()))
			{

				donetargets.add(node.getTargetid());
				count++;
			}
			node= node.parent;

		}
		return count;
	}
	
	
	

	private static int coveredCount(ArrayList<Integer> donetargets,
			TargetNode node) 
	{

		int count =0;

		if( !donetargets.contains(node.getTargetid()))
		{

			donetargets.add(node.getTargetid());
			count++;
		}




		while(node.parent!=null)
		{


			if(!donetargets.contains(node.getTargetid()))
			{

				donetargets.add(node.getTargetid());
				count++;
			}
			node= node.parent;

		}
		return count;
	}


	private static boolean isAttackedtargetcovered(
			TargetNode node, int attackedtarget) 
	{
		if(node.getTargetid()==attackedtarget)
		{

			return true;

		}
		while(node.parent!=null)
		{
			if(node.getTargetid()==attackedtarget)
			{
				return true;
			}
			node= node.parent;

		}
		return false;
	}


	
	
	private static ArrayList<TargetNode> getOnePathByBFSWithSrcDest(double dmax, int[][] gamedata,
			ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int src, int dest, int[] distallocation) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		int nTargets = targets.size(); 
		TargetNode start = new TargetNode(targetmaps.get(src));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==dest) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//System.out.println();
				//printPath(node);
				ArrayList<Integer> donetargets = new ArrayList<Integer>();
				
				int coveredcount = coveredCount(donetargets, node);
				
				pathcounter++;
				
				int c = (int)(.6*targets.size());
				if(coveredcount >= (c))
				{
				//System.out.println();
					goals.add(node);
					distallocation[0] = (int)node.distancecoveredyet;
					break;
				}
				/*if(pathcounter>5000)
					break;*/
			}
			if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<TargetNode> succs = ExpandTarget(node, dmax, targetmaps, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}
	
	
	private static ArrayList<TargetNode> getPathsByBFSWithSrcDest(double dmax, int[][] gamedata,
			ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, int src, int dest) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		int nTargets = targets.size(); 
		TargetNode start = new TargetNode(targetmaps.get(src));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==dest) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//System.out.println();
				//printPath(node);
				pathcounter++;
				//System.out.println();
				goals.add(node);

				/*if(pathcounter>5000)
					break;*/
			}
			if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<TargetNode> succs = ExpandTarget(node, dmax, targetmaps, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}


	private static ArrayList<TargetNode> getPathsByBFS(double dmax,
			ArrayList<ArrayList<TargetNode>> paths, int[][] gamedata, ArrayList<TargetNode> targets) throws Exception 
			{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		int nTargets = targets.size(); 
		TargetNode start = new TargetNode(targets.get(0));
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==start.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//printPath(node, node);
				pathcounter++;
				//System.out.println();
				goals.add(node);

				/*if(pathcounter>5000)
					break;*/
			}
			else if(node.distancecoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<TargetNode> succs = Expand(node, dmax, targets);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}


		System.out.println("Total path "+ pathcounter);
		if(pathcounter==0)
		{
			throw new Exception("No route to select");
		}
		int Pmat[][] = new int[nTargets][pathcounter];
		int pathindex = 0;
		for(TargetNode n: goals)
		{
			TargetNode tmpgoal = n;
			TargetNode tmpstart = n.parent;
			Pmat[ map.get(tmpgoal.getTargetid())][pathindex]  = 1;
			while(true)
			{
				Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					break;
				}

			}
			pathindex++;
		}*/
		return goals;

			}

	public static void printPath(TargetNode node) {

		if(node.parent == null)
			return;
		printPath(node.parent);
		System.out.print(  node.getTargetid() +"->");

	}

	private static ArrayList<TargetNode> ExpandTarget(TargetNode node, double dmax, HashMap<Integer,TargetNode> targetmaps, ArrayList<TargetNode> targets) 
	{
		ArrayList<TargetNode> successors = new ArrayList<TargetNode>();

		/**
		 * find the index for  node
		 */
		TargetNode tmpnode = targetmaps.get(node.getTargetid());
		
		for(TargetNode nei: tmpnode.getNeighbors())
		{
			if(targets.contains(nei))
			{
				TargetNode newnei = new TargetNode(nei);
				newnei.distancecoveredyet = node.distancecoveredyet + tmpnode.getDistance(nei);
				newnei.parent = node;
				if(newnei.distancecoveredyet<=dmax)
				{
					successors.add(newnei);
				}
			}
		}
		return successors;

	}
	
	
	private static ArrayList<TargetNode> Expand(TargetNode node, double dmax, ArrayList<TargetNode> targets) 
	{
		ArrayList<TargetNode> successors = new ArrayList<TargetNode>();

		/**
		 * find the index for  node
		 */
		TargetNode tmpnode = new TargetNode(); 
		for(TargetNode t: targets)
		{
			if(t.getTargetid()==node.getTargetid())
			{
				tmpnode = t;
				break;
			}
		}
		for(TargetNode nei: tmpnode.getNeighbors())
		{
			TargetNode newnei = new TargetNode(nei);
			newnei.distancecoveredyet = node.distancecoveredyet + tmpnode.getDistance(nei);
			newnei.parent = node;
			if(newnei.distancecoveredyet<=dmax)
			{
				successors.add(newnei);
			}
		}
		return successors;

	}


	private static ArrayList<TargetNode> GreedyExpand(TargetNode node, 
			double dmax, ArrayList<TargetNode> targets, ArrayList<Integer> donetargets) {
		ArrayList<TargetNode> successors = new ArrayList<TargetNode>();

		/**
		 * find the index for  node
		 */
		TargetNode tmpnode = new TargetNode(); 
		for(TargetNode t: targets)
		{
			if(t.getTargetid()==node.getTargetid())
			{
				tmpnode = t;
				break;
			}
		}
		for(TargetNode nei: tmpnode.getNeighbors())
		{


			//TODO
			if(!donetargets.contains(nei.getTargetid()))
			{

				TargetNode newnei = new TargetNode(nei);
				newnei.distancecoveredyet = node.distancecoveredyet + tmpnode.getDistance(nei);
				newnei.parent = node;
				if(newnei.distancecoveredyet<=dmax)
				{
					successors.add(newnei);
				}
			}
		}
		return successors;

	}
	
	
	

	public void buildGraph(int numRow, int numCol, int[][] gamedata) 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, gamedata[target][0]);
			node.defenderreward = gamedata[target][0];
			node.defenderpenalty = gamedata[target][1];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = gamedata[target][3];

			targets.add(node);
			if(target==0)
			{
				graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;

		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				targets.get(targetid).setRowCol(row, col);
				targets.get(targetid).setCoinvalue(gamedata[targetid][0]);
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						targets.get(targetid).addNeighbor(targets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						targets.get(targetid).setPath(targets.get(neighborid), pathnodes);
						targets.get(targetid).setPathUtility(targets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						Double distance = 1.0;//rand.nextDouble()*10+5;
						targets.get(targetid).addDistance(targets.get(neighborid), Math.floor(distance));
						targets.get(neighborid).addDistance(targets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}

	}
	
	
	public static void buildGraph(int numRow, int numCol, int[][] gamedata, ArrayList<TargetNode> targets) 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, gamedata[target][0]);
			node.defenderreward = 0;
			node.defenderpenalty = -gamedata[target][2];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = 0;

			targets.add(node);
			if(target==0)
			{
				graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;

		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				targets.get(targetid).setRowCol(row, col);
				targets.get(targetid).setCoinvalue(gamedata[targetid][0]);
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						targets.get(targetid).addNeighbor(targets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						targets.get(targetid).setPath(targets.get(neighborid), pathnodes);
						targets.get(targetid).setPathUtility(targets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						Double distance = 1.0;//rand.nextDouble()*10+5;
						targets.get(targetid).addDistance(targets.get(neighborid), Math.floor(distance));
						targets.get(neighborid).addDistance(targets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}

	}
	



	public void buildGraph(int numRow, int numCol, double[][] u, double[][] e, int[][] gamedata) 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, 0);
			/*node.defenderreward = gamedata[target][0];
			node.defenderpenalty = gamedata[target][1];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = gamedata[target][3];*/

			targets.add(node);
			if(target==0)
			{
				graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;

		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				TargetNode tmp = targets.get(targetid);
				
				
				tmp.setRowCol(row, col);
				
				
				tmp.setCoinvalue(gamedata[targetid][0]);
				tmp.defenderreward = 0;
				tmp.defenderpenalty = -u[row][col];
				tmp.attackerreward = u[row][col];
				tmp.attackerpenalty = 0;
				tmp.setAnimaldensity(u[row][col]);
				
				
				
				
				
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						targets.get(targetid).addNeighbor(targets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						targets.get(targetid).setPath(targets.get(neighborid), pathnodes);
						targets.get(targetid).setPathUtility(targets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						//Double distance = 1.0;//rand.nextDouble()*10+5;
						
						double d1 = Math.sqrt((50*50) + e[neighborrow][neighborcol]*e[neighborrow][neighborcol]);
						
						/*System.out.println("target "+ targetid + ", neirow "+ neighborrow
								+ ", neicol "+ neighborcol + ", ele "+e[neighborrow][neighborcol]+", d = "+ d1);
						
						*/
						targets.get(targetid).addDistance(targets.get(neighborid), Math.floor(d1));
						//targets.get(neighborid).addDistance(targets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}

	}
	
	
	public static void buildcsvGraph(int numRow, int numCol, double[][] u, double[][] e, ArrayList<TargetNode> targets) 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, 0);
			/*node.defenderreward = gamedata[target][0];
			node.defenderpenalty = gamedata[target][1];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = gamedata[target][3];*/

			targets.add(node);
			if(target==0)
			{
				//graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;
		
		
		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				TargetNode tmp = targets.get(targetid);
				
				
				tmp.setRowCol(row, col);
				
				
				tmp.setCoinvalue(u[row][col]);
				tmp.defenderreward = 0;
				tmp.defenderpenalty = -u[row][col];
				tmp.attackerreward = u[row][col];
				tmp.attackerpenalty = 0;
				tmp.setAnimaldensity(u[row][col]);
				
				try {
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata3.csv"),true));
					
					pw.append(tmp.getTargetid()+","+u[row][col]+ ","+(row*50) + ","+(col*50)+ "\n");
					pw.close();
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

				
				
				
				
				
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						targets.get(targetid).addNeighbor(targets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						targets.get(targetid).setPath(targets.get(neighborid), pathnodes);
						targets.get(targetid).setPathUtility(targets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						//Double distance = 1.0;//rand.nextDouble()*10+5;
						
						double d1 = Math.sqrt((50*50) + e[neighborrow][neighborcol]*e[neighborrow][neighborcol]);
						
						/*System.out.println("target "+ targetid + ", neirow "+ neighborrow
								+ ", neicol "+ neighborcol + ", ele "+e[neighborrow][neighborcol]+", d = "+ d1);
						
						*/
						targets.get(targetid).addDistance(targets.get(neighborid), Math.floor(d1));
						//targets.get(neighborid).addDistance(targets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}

	}
	
	public static void buildcsvGraphExp(int numRow, int numCol, double[][] u,  ArrayList<TargetNode> targets, int iter) throws Exception 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);
		
		
		
		try {
			
			
			 File f = new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata");
			 
			 if(f.exists())
			 {
				 f.delete();
				 f.createNewFile();
			 }
			
			
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata"+iter+".csv"),true));
			
			pw.append("Id,U,X,Y"+"\n");
			pw.close();
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		for(int target=0; target< (numRow*numCol); target++)
		{
			//System.out.println("target "+ target + " animal density "+ gamedata[target][0]);

			TargetNode node = new TargetNode(target, 0);
			/*node.defenderreward = gamedata[target][0];
			node.defenderpenalty = gamedata[target][1];
			node.attackerreward = gamedata[target][2];
			node.attackerpenalty = gamedata[target][3];*/

			targets.add(node);
			if(target==0)
			{
				//graph = node;
				node.setStart(true);

			}
			if(target==((numRow*numCol)-1))
			{
				node.setGoal(true);
			}
		}

		//setDummyUtility();

		/**
		 * build the connections and graph
		 */
		int targetid = 0;
		
		
		for(int row=0; row<numRow; row++)
		{
			for(int col=0; col<numCol; col++)
			{
				/**
				 * add the neighbors and distances
				 */

				TargetNode tmp = targets.get(targetid);
				
				
				tmp.setRowCol(row, col);
				
				
				tmp.setCoinvalue(u[iter][targetid]);
				tmp.defenderreward = 0;
				tmp.defenderpenalty = -u[iter][targetid];
				tmp.attackerreward = u[iter][targetid];
				tmp.attackerpenalty = 0;
				tmp.setAnimaldensity(u[iter][targetid]);
				
				try {
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata"+iter+".csv"),true));
					
					pw.append(tmp.getTargetid()+","+u[iter][targetid]+ ","+(row*1) + ","+(col*1)+ "\n");
					pw.close();
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

				
				
				
				
				
				for(int neighborindex=0; neighborindex<8; neighborindex++)
				{
					int neighborrow = -1;
					int neighborcol = -1;
					if(neighborindex==0)
					{
						neighborrow = row-1;
						neighborcol = col-1;
					}
					else if(neighborindex==1)
					{
						neighborrow = row-1;
						neighborcol = col;
					}
					else if(neighborindex==2)
					{
						neighborrow = row-1;
						neighborcol = col+1;
					}
					else if(neighborindex==3)
					{
						neighborrow = row;
						neighborcol = col-1;
					}
					else if(neighborindex==4)
					{
						neighborrow = row;
						neighborcol = col+1;
					}
					else if(neighborindex==5)
					{
						neighborrow = row+1;
						neighborcol = col-1;
					}
					else if(neighborindex==6)
					{
						neighborrow = row+1;
						neighborcol = col;
					}
					else if(neighborindex==7)
					{
						neighborrow = row+1;
						neighborcol = col+1;
					}


					if(neighborrow >=0 && neighborrow <numRow && neighborcol >=0 && neighborcol < numCol)
					{
						//int targetid = targets.get(targetindex).getTargetid();
						int neighborid = (neighborrow* (numCol)) + neighborcol;
						targets.get(targetid).addNeighbor(targets.get(neighborid));
						ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
						//pathnodes.add(targets.get(targetid));
						//pathnodes.add(targets.get(neighborid));
						targets.get(targetid).setPath(targets.get(neighborid), pathnodes);
						targets.get(targetid).setPathUtility(targets.get(neighborid), 0.0);


						if(targetid==neighborid)
						{
							System.out.println("what !!!");
						}
						//System.out.println(" target "+ targetid + ", adding neighbor "+ neighborid);
						//Double distance = 1.0;//rand.nextDouble()*10+5;
						
						double d1 = 1;
						
						/*System.out.println("target "+ targetid + ", neirow "+ neighborrow
								+ ", neicol "+ neighborcol + ", ele "+e[neighborrow][neighborcol]+", d = "+ d1);
						
						*/
						targets.get(targetid).addDistance(targets.get(neighborid), Math.floor(d1));
						//targets.get(neighborid).addDistance(targets.get(targetid), Math.floor(distance));
					}

				}
				targetid++;
			}
		}

	}


	private void setDummyUtility() {

		targets.get(0).setAnimaldensity(10);
		targets.get(1).setAnimaldensity(0);
		targets.get(2).setAnimaldensity(0);
		targets.get(3).setAnimaldensity(10);
		targets.get(4).setAnimaldensity(0);
		targets.get(4).setAnimaldensity(0);
		targets.get(5).setAnimaldensity(0);
		targets.get(6).setAnimaldensity(0);
		targets.get(7).setAnimaldensity(0);
		targets.get(8).setAnimaldensity(0);
		targets.get(9).setAnimaldensity(0);
		targets.get(10).setAnimaldensity(0);
		targets.get(11).setAnimaldensity(0);
		targets.get(12).setAnimaldensity(0);
		targets.get(13).setAnimaldensity(10);
		targets.get(14).setAnimaldensity(0);
		targets.get(15).setAnimaldensity(0);
		targets.get(16).setAnimaldensity(0);
		targets.get(17).setAnimaldensity(0);
		targets.get(18).setAnimaldensity(0);
		targets.get(19).setAnimaldensity(10);



	}

	/**
	 * 
	 * @param numberoftargets number of unaccessible targets
	 */
	/*public void makeGeoConstraints(int numberoftargets, int nUnaccesstargets)
	{
		ArrayList<Integer> tmpunaccessibltargets = new ArrayList<Integer>();
		for(int i=0; i<numberoftargets; i++)
		{
			if((i!=0) && i != (numberoftargets-1))
			tmpunaccessibltargets.add(i);
		}
		Random rand = new Random(nUnaccesstargets);
		//SecurityGameContraction.unimportanttargets.clear();
		for(int i=0; i<nUnaccesstargets; i++)
		{
			int index = rand.nextInt((tmpunaccessibltargets.size()-1 - 0) + 1) + 0;
			int target = tmpunaccessibltargets.get(index);
			this.unaccessibltargets.put(target, target+1);
			SecurityGameContraction.unimportanttargets.add(target);
			//System.out.println(i+"th Unaccessible target "+ target);
			tmpunaccessibltargets.remove(index);

		}
	}
	 */
	public static ArrayList<TargetNode> chooseContractedNodes(int ndominatedtargets, int[][] gamedata)
	{






		ArrayList<TargetNode> dominatednodes = new ArrayList<TargetNode>(); // dominated
		if(ndominatedtargets==0)
			return dominatednodes;
		ArrayList<Integer> added = new ArrayList<Integer>();
		int dominatedcount = 0;
		Random rand = new Random();

		// nextInt is normally exclusive of the top value,
		// so add 1 to make it inclusive


		added.add(3);
		added.add(4);
		added.add(5);
		added.add(6);
		added.add(8);
		added.add(9);
		added.add(10);
		added.add(11);


		int l=0;
		while(true)
		{
			int randomNum = added.get(l++);//rand.nextInt((targets.size() - 1-1) + 1) + 1;
			//if(!added.contains(randomNum))
			//{
			added.add(randomNum);
			targets.get(randomNum).defenderreward=0;
			targets.get(randomNum).setAnimaldensity(0);

			gamedata[randomNum][0] = 0;
			dominatedcount++;
			//}
			if(dominatedcount==ndominatedtargets)
				break;
		}
		for(Integer x: added)
		{
			//dominatednodes.add(targets.get(x));
		}



		dominatedcount=0;


		for(int i=0; i<SecurityGameContraction.targets.size(); i++)
		{

			for(int j =0; j<SecurityGameContraction.targets.size(); j++)
			{
				if((i != j) && ((i!=0) && (i!=SecurityGameContraction.targets.size()-1)))
				{
					TargetNode inode = targets.get(i);
					TargetNode jnode = targets.get(j);
					double di = TargetNode.getAvgDistance(inode);
					double ui = inode.getAnimaldensity();
					double dj = TargetNode.getAvgDistance(jnode);
					double uj = jnode.getAnimaldensity();
					if((di>=dj) && (ui<uj))
					{
						dominatednodes.add(inode);
						//System.out.print(inode.getTargetid() + ", ");
						dominatedcount++;
						break;
					}
				}

			}
			if(dominatedcount==ndominatedtargets)
			{
				break;
			}
		}

		//System.out.println();
		return dominatednodes;

	}

	public static ArrayList<TargetNode> testContraction(int numberoftargets, 
			int nUnaccesstargets, int[][] gamedata, double dmax, int row, int col)
			{
		/**
		 * 
		 */
		SecurityGameContraction sgc = new SecurityGameContraction(row,col , gamedata);

		//sgc.makeGeoConstraints(numberoftargets, nUnaccesstargets);
		/*SecurityGameContraction.unimportanttargets.add(2); 
		SecurityGameContraction.unimportanttargets.add(3);
		 SecurityGameContraction.unimportanttargets.add(18);
		SecurityGameContraction.unimportanttargets.add(17);
		SecurityGameContraction.unimportanttargets.add(10);
		SecurityGameContraction.unimportanttargets.add(1);*/

		//sgc.removeUnaccessibleNodes();
		ArrayList<TargetNode> dominatednodes = chooseContractedNodes(nUnaccesstargets, gamedata);
		/*if(nUnaccesstargets==1)
		{
			dominatednodes.clear();
			dominatednodes.add(targets.get(1));
		}*/
		//dominatednodes.add(targets.get(2));
		//dominatednodes.add(targets.get(4));

		ArrayList<TargetNode> contractednodes = sgc.contractGraph(dominatednodes, targets, dmax);
		return contractednodes;
		//System.out.print( " HI " );
			}



	public static void setAnimalDensities(double[] animaldensity)
	{
		if(SecurityGameContraction.targets.size() != animaldensity.length)
			return;
		for(int i=0; i<animaldensity.length; i++)
		{
			SecurityGameContraction.targets.get(i).setAnimaldensity(animaldensity[i]);
		}


	}


	private ArrayList<TargetNode> contractGraph(ArrayList<TargetNode> dominatednodes, ArrayList<TargetNode> targets, double dmax) 
	{
		//double[]  animaldensity = {3, 0,0,0, 0, 4, 0, 0, 0,0, 5, 0,0 , 2, 0, 0, 0, 6, 7, 2};

		//setAnimalDensities(animaldensity);

		if(dominatednodes.size()==0)
		{
			return dominatednodes;
		}


		//ArrayList<TargetNode> dominatednodes = chooseContractedNodes();
		//System.out.println("# of dominated nodes " + dominatednodes.size());
		/*for(int i=0; i<SecurityGameContraction.targets.size(); i++)
		{
			System.out.println(" node " + targets.get(i).getTargetid() + " utility: " + targets.get(i).getAnimaldensity());
		}*/

		ArrayList<TargetNode> donedominatednodes = new ArrayList<TargetNode>();

		int c=0;
		for(TargetNode middle : dominatednodes)
		{
			cleanNeighbors(donedominatednodes, middle);
			/*if(donedominatednodes.size()>0)
			{
				System.out.println("Contraction result after node " + donedominatednodes.get(c).getTargetid());
				c++;
				cleanNeighbors(donedominatednodes, middle);
				//SecurityGameContraction.removePathsToDominatedNodes(donedominatednodes);
				SecurityGameContraction.printNodesWithNeighborsAndPath(donedominatednodes);
			}
			if(donedominatednodes.size()==0)
			{
				SecurityGameContraction.printNodesWithNeighborsAndPath(donedominatednodes);
			}*/
			ArrayList<Integer[]> donepair = new ArrayList<Integer[]>();

			//System.out.println("Contracting node " + middle.getTargetid());

			//System.out.println("Neighbors before contraction for node "+ middle.getTargetid());
			//printNeighbors(middle);

			ArrayList<TargetNode> srcnodes = middle.getNeighbors();
			ArrayList<TargetNode> destnodes = middle.getNeighbors();

			/*if(srcnodes.size()==0)
			{
				System.out.print("empty neighbor list");
			}*/



			for(int  srcindex=0; srcindex<srcnodes.size(); srcindex++)
			{
				TargetNode src = srcnodes.get(srcindex);
				for(int  destindex=0; destindex<destnodes.size(); destindex++)
				{
					TargetNode dest = destnodes.get(destindex);
					if((src.getTargetid()!= dest.getTargetid()) 
							&& (src.getTargetid() != middle.getTargetid()) 
							&& (dest.getTargetid() != middle.getTargetid()))
					{
						Integer[] pair = {src.getTargetid(), dest.getTargetid()}; 

						if(!donePair(pair, donepair))
						{


							//
							donepair.add(pair);

							/*System.out.println("middle node "+middle.getTargetid()+ ", evaluating edge "+src.getTargetid() +
									" -> "+middle.getTargetid() + " -> "+ dest.getTargetid());
							 */
							// calculate c(src, middle) + c(middle,dest)
							double distviamiddlenode = src.computeDistance(src,middle) + middle.computeDistance(middle,dest);

							/*System.out.println("c("+src.getTargetid() + ","+
									middle.getTargetid()+") + c("+ middle.getTargetid()+","+ 
									dest.getTargetid()+ ") = "+ distviamiddlenode);*/
							// calculate c(src,dest)
							double directdist = src.computeDistance(src,dest);

							//double shortestdistance  = isShortest(src,middle,dest,distviamiddlenode); // -1 if uvw is the shortest

							/*System.out.println("c("+src.getTargetid() + ","+
									dest.getTargetid()+") = " + shortestdistance);*/
							//if(shortestdistance==-1 && distviamiddlenode<dmax) // if middle node is in the shortest path

							if(distviamiddlenode<directdist && distviamiddlenode<dmax) // if middle node is in the shortest path
							{
								/**
								 * add the dest to neighbor list
								 * add middle node in the path to dest
								 */
								/**
								 * if dest is already in the neighbor list
								 * remove the dest from neighbor of src
								 */
								//printNeighbors(src);
								//printDistances(src);
								src.addNeighbor(dest);
								src.addDistance(dest, distviamiddlenode);
								//printNeighbors(src);
								//printDistances(src);
								/**
								 * build path from src -> path to middle -> dest -> path to middle in
								 * rev order
								 * add the path to src 
								 * and also to dest
								 */
								ArrayList<TargetNode> pathtodest = buildPath(src, middle, dest); 
								src.setPath(dest, pathtodest);
								/**
								 * compute path utility
								 */
								double pathutility = computePathUtility(src, pathtodest, dest);
								src.setPathUtility(dest, pathutility);

								//printNeighbors(dest);
								//printDistances(dest);
								dest.addNeighbor(src);
								dest.addDistance(src, distviamiddlenode);
								dest.setPathUtility(src, pathutility);
								//printNeighbors(dest);
								//printDistances(dest);

								/**
								 * rev the path to dest, add to dest
								 */
								ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
								for(int i=0; i<pathtodest.size(); i++)
								{
									pathtosrc.add(pathtodest.get(pathtodest.size()-i-1));
								}
								dest.setPath(src, pathtosrc);


								/*System.out.println("Adding node "+ middle.getTargetid() + " to hop node "
										+ "to reach node " + dest.getTargetid() + " from node "+ src.getTargetid() +" path utility : "+pathutility + "\n\n");
								 */	}
							else  // middle node is not in the shortest path
							{
								/*System.out.println("middle node " + middle.getTargetid() + " is not needed "
										+ "to reach node " + dest.getTargetid() + " from node "+ src.getTargetid() +"\n\n");*/

								/*src.removeNeighbor(middle);
							middle.removeNeighbor(src);
							middle.removeNeighbor(dest);
							dest.removeNeighbor(middle);

							src.removeDistance(middle);
							middle.removeDistance(src);
							middle.removeDistance(dest);
							dest.removeDistance(middle);

							src.removePath(middle);
							middle.removePath(src);
							middle.removePath(dest);
							dest.removePath(middle);*/
							}



						}
					}
				}
			}
			//System.out.println("Neighbors after contraction for node "+ middle.getTargetid()+", targets size "+ targets.size());
			//printNeighbors(middle);

			//cleanNeighbors(donedominatednodes, src);
			//cleanNeighbors(donedominatednodes, dest);	
			donedominatednodes.add(getTargetNode(middle.getTargetid(), targets));
		}
		return dominatednodes;
	}





	private ArrayList<TargetNode> contractGraphWithExtreamPruning(ArrayList<TargetNode> dominatednodes, ArrayList<TargetNode> targets, double dmax) 
	{
		//double[]  animaldensity = {3, 0,0,0, 0, 4, 0, 0, 0,0, 5, 0,0 , 2, 0, 0, 0, 6, 7, 2};

		//setAnimalDensities(animaldensity);

		if(dominatednodes.size()==0)
		{
			return dominatednodes;
		}

		TargetNode basenode = getTargetNode(0, targets);


		//ArrayList<TargetNode> dominatednodes = chooseContractedNodes();
		//System.out.println("# of dominated nodes " + dominatednodes.size());
		/*for(int i=0; i<SecurityGameContraction.targets.size(); i++)
		{
			System.out.println(" node " + targets.get(i).getTargetid() + " utility: " + targets.get(i).getAnimaldensity());
		}*/

		ArrayList<TargetNode> donedominatednodes = new ArrayList<TargetNode>();

		int c=0;
		for(TargetNode middle : dominatednodes)
		{
			cleanNeighbors(donedominatednodes, middle);
			/*if(donedominatednodes.size()>0)
			{
				System.out.println("Contraction result after node " + donedominatednodes.get(c).getTargetid());
				c++;
				cleanNeighbors(donedominatednodes, middle);
				//SecurityGameContraction.removePathsToDominatedNodes(donedominatednodes);
				SecurityGameContraction.printNodesWithNeighborsAndPath(donedominatednodes);
			}
			if(donedominatednodes.size()==0)
			{
				SecurityGameContraction.printNodesWithNeighborsAndPath(donedominatednodes);
			}*/
			ArrayList<Integer[]> donepair = new ArrayList<Integer[]>();

			System.out.println("Contracting node " + middle.getTargetid());

			//System.out.println("Neighbors before contraction for node "+ middle.getTargetid());
			//printNeighbors(middle);

			ArrayList<TargetNode> srcnodes = middle.getNeighbors();
			ArrayList<TargetNode> destnodes = middle.getNeighbors();

			/*if(srcnodes.size()==0)
			{
				System.out.print("empty neighbor list");
			}*/



			for(int  srcindex=0; srcindex<srcnodes.size(); srcindex++)
			{
				TargetNode src = srcnodes.get(srcindex);
				for(int  destindex=0; destindex<destnodes.size(); destindex++)
				{
					TargetNode dest = destnodes.get(destindex);
					if((src.getTargetid()!= dest.getTargetid()) 
							&& (src.getTargetid() != middle.getTargetid()) 
							&& (dest.getTargetid() != middle.getTargetid()))
					{
						Integer[] pair = {src.getTargetid(), dest.getTargetid()}; 

						if(!donePair(pair, donepair))
						{


							//
							donepair.add(pair);

							/*System.out.println("middle node "+middle.getTargetid()+ ", evaluating edge "+src.getTargetid() +
									" -> "+middle.getTargetid() + " -> "+ dest.getTargetid());
							 */
							// calculate c(src, middle) + c(middle,dest)
							double distviamiddlenode = src.computeDistance(src,middle) + middle.computeDistance(middle,dest);

							/*System.out.println("c("+src.getTargetid() + ","+
									middle.getTargetid()+") + c("+ middle.getTargetid()+","+ 
									dest.getTargetid()+ ") = "+ distviamiddlenode);*/
							// calculate c(src,dest)
							double directdist = src.computeDistance(src,dest);


							/**
							 * do the extream pruning 
							 * 
							 * find shortest dist from base to src and dest
							 * then if di+dj+distviamiddlenode< dmax add the edge
							 * 
							 */


							//ArrayList<Integer> pathnodes = new ArrayList<Integer>();
							double di= src.getDistfrombase();//findShortestPath(basenode, src, targets, dominatednodes, pathnodes);


							//ArrayList<Integer> pathnodes = new ArrayList<Integer>();
							//pathnodes.clear();
							double dj= dest.getDistfrombase();//findShortestPath(basenode, dest, targets, dominatednodes, pathnodes);



							//double shortestdistance  = isShortest(src,middle,dest,distviamiddlenode); // -1 if uvw is the shortest

							/*System.out.println("c("+src.getTargetid() + ","+
									dest.getTargetid()+") = " + shortestdistance);*/
							//if(shortestdistance==-1 && distviamiddlenode<dmax) // if middle node is in the shortest path

							if(distviamiddlenode<directdist && ((distviamiddlenode+di+dj) <=dmax)) // if middle node is in the shortest path
							{
								/**
								 * add the dest to neighbor list
								 * add middle node in the path to dest
								 */
								/**
								 * if dest is already in the neighbor list
								 * remove the dest from neighbor of src
								 */
								//printNeighbors(src);
								//printDistances(src);
								src.addNeighbor(dest);
								src.addDistance(dest, distviamiddlenode);
								//printNeighbors(src);
								//printDistances(src);
								/**
								 * build path from src -> path to middle -> dest -> path to middle in
								 * rev order
								 * add the path to src 
								 * and also to dest
								 */
								ArrayList<TargetNode> pathtodest = buildPath(src, middle, dest); 
								src.setPath(dest, pathtodest);
								/**
								 * compute path utility
								 */
								double pathutility = computePathUtility(src, pathtodest, dest);
								src.setPathUtility(dest, pathutility);

								//printNeighbors(dest);
								//printDistances(dest);
								dest.addNeighbor(src);
								dest.addDistance(src, distviamiddlenode);
								dest.setPathUtility(src, pathutility);
								//printNeighbors(dest);
								//printDistances(dest);

								/**
								 * rev the path to dest, add to dest
								 */
								ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
								for(int i=0; i<pathtodest.size(); i++)
								{
									pathtosrc.add(pathtodest.get(pathtodest.size()-i-1));
								}
								dest.setPath(src, pathtosrc);


								/*System.out.println("Adding node "+ middle.getTargetid() + " to hop node "
										+ "to reach node " + dest.getTargetid() + " from node "+ src.getTargetid() +" path utility : "+pathutility + "\n\n");
								 */	}
							else  // middle node is not in the shortest path
							{
								/*System.out.println("middle node " + middle.getTargetid() + " is not needed "
										+ "to reach node " + dest.getTargetid() + " from node "+ src.getTargetid() +"\n\n");*/

								/*src.removeNeighbor(middle);
							middle.removeNeighbor(src);
							middle.removeNeighbor(dest);
							dest.removeNeighbor(middle);

							src.removeDistance(middle);
							middle.removeDistance(src);
							middle.removeDistance(dest);
							dest.removeDistance(middle);

							src.removePath(middle);
							middle.removePath(src);
							middle.removePath(dest);
							dest.removePath(middle);*/
							}



						}
					}
				}
			}
			System.out.println("Neighbors after contraction for node "+ middle.getTargetid()+", targets size "+ targets.size());
			//printNeighbors(middle);

			//cleanNeighbors(donedominatednodes, src);
			//cleanNeighbors(donedominatednodes, dest);	
			donedominatednodes.add(getTargetNode(middle.getTargetid(), targets));
		}
		return dominatednodes;
	}





	private ArrayList<TargetNode> contractGraphV2(ArrayList<TargetNode> dominatednodes, ArrayList<TargetNode> targets) 
	{
		//double[]  animaldensity = {3, 0,0,0, 0, 4, 0, 0, 0,0, 5, 0,0 , 2, 0, 0, 0, 6, 7, 2};

		//setAnimalDensities(animaldensity);

		if(dominatednodes.size()==0)
		{
			return dominatednodes;
		}


		//ArrayList<TargetNode> dominatednodes = chooseContractedNodes();
		//System.out.println("# of dominated nodes " + dominatednodes.size());
		/*for(int i=0; i<SecurityGameContraction.targets.size(); i++)
		{
			System.out.println(" node " + targets.get(i).getTargetid() + " utility: " + targets.get(i).getAnimaldensity());
		}*/

		ArrayList<TargetNode> donedominatednodes = new ArrayList<TargetNode>();

		for(TargetNode middle : dominatednodes)
		{
			ArrayList<Integer[]> donepair = new ArrayList<Integer[]>();

			//System.out.println("Contracting node " + middle.getTargetid());
			cleanNeighbors(donedominatednodes, middle);
			//System.out.println("Neighbors before contraction for node "+ middle.getTargetid());
			//printNeighbors(middle);

			ArrayList<TargetNode> srcnodes = middle.getNeighbors();
			ArrayList<TargetNode> destnodes = middle.getNeighbors();

			/*if(srcnodes.size()==0)
			{
				System.out.print("empty neighbor list");
			}*/



			for(int  srcindex=0; srcindex<srcnodes.size(); srcindex++)
			{
				TargetNode src = srcnodes.get(srcindex);
				for(int  destindex=0; destindex<destnodes.size(); destindex++)
				{
					TargetNode dest = destnodes.get(destindex);
					if((src.getTargetid()!= dest.getTargetid()) 
							&& (src.getTargetid() != middle.getTargetid()) 
							&& (dest.getTargetid() != middle.getTargetid()))
					{
						Integer[] pair = {src.getTargetid(), dest.getTargetid()}; 

						if(!donePair(pair, donepair))
						{


							donedominatednodes.add(targets.get(middle.getTargetid()));
							donepair.add(pair);

							/*System.out.println("middle node "+middle.getTargetid()+ ", evaluating edge "+src.getTargetid() +
									" -> "+middle.getTargetid() + " -> "+ dest.getTargetid());*/

							// calculate c(src, middle) + c(middle,dest)
							double distviamiddlenode = src.computeDistance(src,middle) + middle.computeDistance(middle,dest);

							/*System.out.println("c("+src.getTargetid() + ","+
									middle.getTargetid()+") + c("+ middle.getTargetid()+","+ 
									dest.getTargetid()+ ") = "+ distviamiddlenode);*/
							// calculate c(src,dest)
							double directdist = src.computeDistance(src,dest);

							//double shortestdistance  = isShortest(src,middle,dest,distviamiddlenode); // -1 if uvw is the shortest

							/*System.out.println("c("+src.getTargetid() + ","+
									dest.getTargetid()+") = " + shortestdistance);*/
							if(distviamiddlenode<directdist) // if middle node is in the shortest path
							{
								/**
								 * add the dest to neighbor list
								 * add middle node in the path to dest
								 */
								/**
								 * if dest is already in the neighbor list
								 * remove the dest from neighbor of src
								 */
								//printNeighbors(src);
								//printDistances(src);
								src.addNeighbor(dest);
								src.addDistance(dest, distviamiddlenode);
								//printNeighbors(src);
								//printDistances(src);
								/**
								 * build path from src -> path to middle -> dest -> path to middle in
								 * rev order
								 * add the path to src 
								 * and also to dest
								 */
								ArrayList<TargetNode> pathtodest = buildPath(src, middle, dest); 
								src.setPath(dest, pathtodest);
								/**
								 * compute path utility
								 */
								double pathutility = computePathUtility(src, pathtodest, dest);
								src.setPathUtility(dest, pathutility);

								//printNeighbors(dest);
								//printDistances(dest);
								dest.addNeighbor(src);
								dest.addDistance(src, distviamiddlenode);
								dest.setPathUtility(src, pathutility);
								//printNeighbors(dest);
								//printDistances(dest);

								/**
								 * rev the path to dest, add to dest
								 */
								ArrayList<TargetNode> pathtosrc = new ArrayList<TargetNode>();
								for(int i=0; i<pathtodest.size(); i++)
								{
									pathtosrc.add(pathtodest.get(pathtodest.size()-i-1));
								}
								dest.setPath(src, pathtosrc);


								/*System.out.println("Adding node "+ middle.getTargetid() + " to hop node "
										+ "to reach node " + dest.getTargetid() + " from node "+ src.getTargetid() +" path utility : "+pathutility + "\n\n");
								 */	}
							else  // middle node is not in the shortest path
							{
								/*System.out.println("middle node " + middle.getTargetid() + " is not needed "
										+ "to reach node " + dest.getTargetid() + " from node "+ src.getTargetid() +"\n\n");*/

								/*src.removeNeighbor(middle);
							middle.removeNeighbor(src);
							middle.removeNeighbor(dest);
							dest.removeNeighbor(middle);

							src.removeDistance(middle);
							middle.removeDistance(src);
							middle.removeDistance(dest);
							dest.removeDistance(middle);

							src.removePath(middle);
							middle.removePath(src);
							middle.removePath(dest);
							dest.removePath(middle);*/
							}



						}
					}
				}
			}
			//System.out.println("Neighbors after contraction for node "+ middle.getTargetid());
			//printNeighbors(middle);

			//cleanNeighbors(donedominatednodes, src);
			//cleanNeighbors(donedominatednodes, dest);	
		}
		return dominatednodes;
	}






	public static void printDistances(TargetNode node) {

		System.out.println("Distances");
		for(TargetNode x: node.getDistances().keySet())
		{
			System.out.println(x.getTargetid()+" : "+ node.getDistance(x));
		}

	}

	public static void printNeighbors(TargetNode node) {

		System.out.print("neighbors of "+ node.getTargetid()+" : ");
		for(TargetNode x: node.getNeighbors())
		{
			System.out.print(x.getTargetid()+", ");
		}
		System.out.println();


	}

	private void cleanNeighbors(ArrayList<TargetNode> donedominatednodes,
			TargetNode middle) {

		for(TargetNode x: donedominatednodes)
		{
			middle.removeDistance(x);
			middle.removeNeighbor(x);
			middle.removePath(x);
			middle.removePathUtility(x);
		}

	}

	private boolean donePair(Integer[] pair, ArrayList<Integer[]> donepair) {

		if(donepair.size()==0)
			return false;
		for(Integer[] x: donepair)
		{
			if((x[0] == pair[0] && x[1] == pair[1]) ||
					(x[0] == pair[1] && x[1] == pair[0]))
			{
				return true;
			}
		}


		return false;
	}

	private static double computePathUtility(TargetNode src,
			ArrayList<TargetNode> pathtodest, TargetNode dest) {

		double sum =0;//src.getAnimaldensity();
		for(TargetNode n : pathtodest)
		{
			sum += n.getAnimaldensity();
		}
		//sum += dest.getAnimaldensity();





		return sum;
	}

	/**
	 * only distance constraints have been applied
	 * @param u
	 * @param v
	 * @param w
	 * @param distviamiddlenode
	 * @param targets2 
	 * @return -1 if u v w is the shortest path, otherwise returns the shortest distance
	 */
	private double isShortest(TargetNode u, TargetNode v,
			TargetNode w, double distviamiddlenode) 
	{
		/**
		 * we have to build a new tree
		 */

		TargetNode start = new TargetNode(u);
		double shortestdistancecoveredyet = 0.0;
		/**
		 * only tree nodes will be in the lists
		 */
		ArrayList<TargetNode> potentialnodes = new ArrayList<TargetNode>(); 
		ArrayList<TargetNode> closednodes = new ArrayList<TargetNode>();
		potentialnodes.add(start);
		//potentialnodes.add(start);
		//System.out.println("Start node "+ u.getTargetid());
		//System.out.println("u->v->w dist "+ distviamiddlenode);
		while((shortestdistancecoveredyet<distviamiddlenode) || potentialnodes.isEmpty())
		{


			if(potentialnodes.isEmpty())
			{
				// uvw is the shortest
				//System.out.println("uvw is the shortest, returning -1 ");
				return -1.0;
			}
			//get the node with shortest neighbor distance
			TargetNode tmpnodetoexpand = getNearestNode(potentialnodes);


			//use this to generate neighbors
			TargetNode nodetoexpand = getTargetNode(tmpnodetoexpand.getTargetid(), SecurityGameContraction.targets);
			shortestdistancecoveredyet = tmpnodetoexpand.getDistfromstart();

			potentialnodes.remove(tmpnodetoexpand);
			//removeNode(potentialnodes,nodetoexpand);
			closednodes.add(tmpnodetoexpand);
			//System.out.println("popping node "+ tmpnodetoexpand.getTargetid());
			//System.out.println("shortest dist from Start node yet "+ shortestdistancecoveredyet);


			if(shortestdistancecoveredyet>= distviamiddlenode)
			{
				// uvw is the shortest
				//System.out.println("uvw is the shortest, returning -1 ");
				return -1.0;
			}
			if(tmpnodetoexpand.getTargetid() == w.getTargetid())
			{

				//System.out.println("uvw is not the shortest edge dist, shortest is " +shortestdistancecoveredyet );

				return shortestdistancecoveredyet;
			}
			/*
			 * expand the nodes
			 * only generate successors other than v
			 */
			/**
			 * use the nodetoexpand to get the neighbors, tmpnodetoexpand does not has any connections other than parent
			 */
			for(TargetNode successor: nodetoexpand.getNeighbors())
			{
				TargetNode suc = new TargetNode(successor);
				if(suc.getTargetid() != v.getTargetid() && (!ifVisited(closednodes,suc)))
				{
					//System.out.println("adding sucessor node "+ suc.getTargetid() + " to "+ nodetoexpand.getTargetid());

					suc.setDistfromstart(tmpnodetoexpand.getDistfromstart() + nodetoexpand.getDistance(successor));
					//System.out.println("setting disst from start " + suc.getDistfromstart());
					suc.setParent(tmpnodetoexpand);
					//System.out.println("setting parent as " + suc.getParent().getTargetid() + " \n\n");
					potentialnodes.add(suc);
				}
			}


		}


		return -1.0;
	}

	private static boolean ifVisited(ArrayList<TargetNode> closednodes, TargetNode suc) {


		for(int i=0; i<closednodes.size(); i++)
		{
			if(closednodes.get(i).getTargetid()== suc.getTargetid())
				return true;
		}
		return false;
	}

	private void removeNode(ArrayList<TargetNode> potentialnodes,
			TargetNode nodetoexpand) {

		int size = potentialnodes.size();

		for(int i=0; i<size; i++)
		{
			if(potentialnodes.get(i).getTargetid() == nodetoexpand.getTargetid())
			{
				potentialnodes.remove(i);
				break;
			}

		}

	}

	private static TargetNode getNearestNode(ArrayList<TargetNode> potentialnodes) {

		if(potentialnodes.size()==0)
		{
			return null;
		}

		TargetNode closestneighbor = new TargetNode();
		double closestdist = Double.POSITIVE_INFINITY;
		for(TargetNode n : potentialnodes)
		{
			if(n.getDistfromstart()<closestdist)
			{
				closestdist = n.getDistfromstart();
				closestneighbor = n;
			}
		}
		//cn = closestneighbor;
		return closestneighbor; 


	}

	private ArrayList<TargetNode> buildPath(TargetNode src, TargetNode middle,
			TargetNode dest) {

		ArrayList<TargetNode> path = new ArrayList<TargetNode>();
		/**
		 * add path nodes from src to middle
		 */
		for(TargetNode tmpsrctomid : src.getPath(middle))
		{
			path.add(tmpsrctomid);

		}
		/**
		 * add middle node
		 */
		path.add(middle);
		/**
		 * i rev order add nodes from dest to middle
		 */

		/*for(TargetNode neighbor: dest.getNeighbors())
		{
			System.out.println("---Neighbor : "+ neighbor.getTargetid());
		 *//**
		 * print path
		 *//*
			ArrayList<TargetNode> path1 = dest.getPath(neighbor);
			System.out.print("Path : "+ dest.getTargetid()+ " --> ");
			for(TargetNode pathnode : path1)
			{
				System.out.print(pathnode.getTargetid()+" --> ");
			}

			System.out.print(neighbor.getTargetid()+ "\n");
			System.out.println("Distance : " + dest.getDistance(neighbor));
			System.out.print("Path utility : "+ dest.getPathUtility(neighbor)+"\n\n");

		}*/




		//printNeighbors(dest);
		//printDistances(dest);
		/*if(dest.getTargetid()==15 && middle.getTargetid()==7)
		{
			System.out.println();
		}*/
		int size = 0;
		if( !middle.getPath(dest).isEmpty())
		{
			size = middle.getPath(dest).size();
		}
		for(int i=0; i<size; i++)
		{
			path.add(middle.getPath(dest).get(i));
		}
		return path;



	}

	/*private void removeUnaccessibleNodes()
	{

		for(Integer fromnode: this.unaccessibltargets.keySet())
		{
			System.out.println("from "+fromnode.intValue() +", "+ this.unaccessibltargets.get(fromnode) + " is unaccessible");
			for(TargetNode neighbor: targets.get(fromnode.intValue()).getNeighbors())
			{
				System.out.print( neighbor.getTargetid() + ", ");

			}

			System.out.println();
			if(targets.get(fromnode.intValue()).getNeighbors().contains(targets.get(this.unaccessibltargets.get(fromnode))))
			{
				targets.get(fromnode.intValue()).removeNeighbor(targets.get(this.unaccessibltargets.get(fromnode)));
				System.out.print( " removed " );
				//break;
			}
			else 
			{
				System.out.print( " doesn't exist " );
			}
			for(TargetNode neighbor: targets.get(fromnode.intValue()).getNeighbors())
			{
				System.out.print( neighbor.getTargetid() + ", ");

			}
			System.out.println();

		}
	 *//**
	 * now try to remove those connections if exists as neighbor
	 *//*



	}
	  */

	public static void printNodesAsNeighbors( ArrayList<TargetNode> domindatednodes) 
	{
		System.out.println("Dominated intermediate nodes which are used in path");
		for(TargetNode node : targets)
		{
			//System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");
			if(!domindatednodes.contains(node))
			{

				for(TargetNode neighbor: node.getNeighbors())
				{
					//System.out.println("---Neighbor : "+ neighbor.getTargetid());
					/**
					 * print path
					 */
					ArrayList<TargetNode> path = node.getPath(neighbor);
					//System.out.print("Path : "+ node.getTargetid()+ " --> ");
					for(TargetNode pathnode : path)
					{
						System.out.print(pathnode.getTargetid()+",");
					}

					//System.out.print(neighbor.getTargetid()+ "\n");
					//System.out.println("Distance : " + node.getDistance(neighbor));
					//System.out.print("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");

				}
			}
		}
		System.out.println();


	}


	public static void printNodesWithNeighborsAndPath( ArrayList<TargetNode> domindatednodes,  ArrayList<TargetNode> targets) 
	{
		for(TargetNode node : targets)
		{
			System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");
			Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			if(!domindatednodes.contains(node))
			{

				for(TargetNode neighbor: node.getNeighbors())
				{
					System.out.println("---Neighbor : "+ neighbor.getTargetid());
					Logger.logit("---Neighbor : "+ neighbor.getTargetid()+"\n");
					/**
					 * print path
					 */
					ArrayList<TargetNode> path = node.getPath(neighbor);
					System.out.print("Path : "+ node.getTargetid()+ " --> ");
					Logger.logit("Path : "+ node.getTargetid()+ " --> ");
					for(TargetNode pathnode : path)
					{
						System.out.print(pathnode.getTargetid()+" --> ");
						Logger.logit(pathnode.getTargetid()+" --> ");
					}

					System.out.print(neighbor.getTargetid()+ "\n");
					System.out.println("Distance : " + node.getDistance(neighbor));
					Logger.logit(neighbor.getTargetid()+ "\n");
					Logger.logit("Distance : " + node.getDistance(neighbor));

					System.out.print("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");
					Logger.logit("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");

				}
			}
		}

	}


	public static void printEdges(ArrayList<TargetNode> targets) 
	{
		//int i=2;
		for(TargetNode node : targets)
		{
			//System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");


			for(TargetNode neighbor: node.getNeighbors())
			{
				//System.out.println("---Neighbor : "+ neighbor.getTargetid());
				/**
				 * print path
				 */
				//ArrayList<TargetNode> path = node.getPath(neighbor);
				System.out.print((node.getTargetid()+1)+" "+ (neighbor.getTargetid()+1) + " "+node.getDistance(neighbor)+"\n");
				//i=i+1;

			}

		}

	}



	public static void transformToDirectedGraph(ArrayList<TargetNode> targets) throws Exception {

		/**
		 * Generate all the duplicate nodes. 
		 * For each edge/neighbor make 2 nodes.
		 * one for out one for in
		 */
		duplicatetargets.clear();
		int transformedcounter = 0;
		for(int targetid = 0; targetid<targets.size(); targetid++)
		{

			int neighborsize = targets.get(targetid).getNeighbors().size();
			for(int neighbor = 0; neighbor<neighborsize; neighbor++)
			{
				/**
				 * create two nodes
				 * one for in 
				 * one for out
				 */
				boolean in = false;
				for(int i = 0; i<2; i++)
				{
					TargetNode node = new TargetNode(targets.get(targetid).getTargetid(), targets.get(targetid).getAnimaldensity());
					node.tranformedtargetid = transformedcounter++;
					node.in= in;
					node.out= (in==false)?true:false;
					in = true;
					node.edgeid = neighbor;
					node.usedforconnection = false;
					duplicatetargets.add(node);
					//System.out.println("adding node "+node.getTargetid() + " to duplicate");

				}
			}
		}

		/**
		 * connect all on and out nodes
		 */
		connectInOutNodes(targets);

		/**
		 * now connect the edges
		 */

		connectAllEdges(targets);

		//System.out.println("hi");




	}

	private static void connectAllEdges(ArrayList<TargetNode> targets) throws Exception {

		ArrayList<Integer[]> doneedges = new ArrayList<Integer[]>();

		for(int targetid = 0; targetid<targets.size(); targetid++) // i
		{
			for(TargetNode neighbor: targets.get(targetid).getNeighbors()) // j
			{
				/**
				 * if the edge(i, j) hasn't been processed yet
				 */
				if(!isEdgeProcessed(targets.get(targetid).getTargetid(), neighbor.getTargetid(), doneedges))
				{

					//System.out.println("done edge :"+ targets.get(targetid).getTargetid() +"-->"+ neighbor.getTargetid());
					doneedges.add(new Integer[]{targets.get(targetid).getTargetid(),neighbor.getTargetid() });

					/**
					 * establish the out connection
					 * 
					 */
					/**
					 * get an unused out node for i
					 */
					TargetNode outi = getOutINode(targets.get(targetid).getTargetid());
					if(null==outi)
					{
						throw new Exception("iout nnnol");
					}
					/**
					 * get an unused in node for j
					 */
					TargetNode inj = getInINode(neighbor.getTargetid());
					if(null==inj)
					{
						throw new Exception("jin nnnol");
					}
					/**
					 * add jin to i's neihbor
					 */
					outi.addNeighbor(inj);
					double dist = targets.get(targetid).getDistance(neighbor);
					outi.addDistance(inj, dist);
					double pathutility = targets.get(targetid).getPathUtility(neighbor);
					outi.setPathUtility(inj, pathutility);
					outi.usedforconnection  = true;
					inj.usedforconnection = true;

					/**
					 * establish the in connection
					 */

					/**
					 * get an unused out node for j
					 */
					TargetNode outj = getOutINode(neighbor.getTargetid());
					if(null==outj)
					{
						throw new Exception("jout nnnol");
					}
					/**
					 * get an unused in node for i
					 */
					TargetNode ini = getInINode(targets.get(targetid).getTargetid());
					if(null==ini)
					{
						throw new Exception("ini nnnol");
					}
					/**
					 * add ini to outj's neighbor
					 */
					outj.addNeighbor(ini);
					dist = neighbor.getDistance(targets.get(targetid));
					outj.addDistance(ini, dist);
					pathutility = neighbor.getPathUtility(targets.get(targetid));
					outj.setPathUtility(ini, pathutility);
					outj.usedforconnection = true;
					ini.usedforconnection = true;


				}
			}
		}

	}

	private static TargetNode getInINode(int targetid) {

		for(TargetNode x: duplicatetargets)
		{
			if((x.getTargetid()==targetid) && (x.in==true) && (x.usedforconnection==false))
				return x;
		}

		return null;
	}

	private static TargetNode getOutINode(int targetid) {

		for(TargetNode x: duplicatetargets)
		{
			if((x.getTargetid()==targetid) && (x.out==true) && (x.usedforconnection==false))
				return x;
		}

		return null;
	}

	private static boolean isEdgeProcessed(int targetid, int targetid2,
			ArrayList<Integer[]> doneedges) {

		for(Integer[] x: doneedges)
		{
			if((x[0]==targetid && x[1]==targetid2) || (x[1]==targetid && x[0]==targetid2))
			{
				return true;
			}
		}
		return false;
	}

	private static void connectInOutNodes(ArrayList<TargetNode> targets) {

		for(int targetid = 0; targetid<targets.size(); targetid++)
		{
			/**
			 * get all in nodes
			 */
			ArrayList<TargetNode> innodes = getAllInOutNodes(targets.get(targetid).getTargetid(), true);

			/**
			 * out nodes
			 */
			ArrayList<TargetNode> outnodes = getAllInOutNodes(targets.get(targetid).getTargetid(), false);

			/**
			 * connect in  to out nodes with 0 distance
			 */
			for( TargetNode in : innodes)
			{
				for(TargetNode out : outnodes)
				{
					in.addNeighbor(out);
					in.addDistance(out, 0.0);
					in.setPath(out, new ArrayList<TargetNode>());
					in.setPathUtility(out, 0.0);

				}
			}

		}


	}



	private static ArrayList<TargetNode> getAllInOutNodes(int targetid, boolean in) {

		ArrayList<TargetNode> nodes = new ArrayList<TargetNode>();

		for(TargetNode x: duplicatetargets)
		{

			if( in &&  (x.getTargetid()==targetid) && (x.in==true))
			{
				nodes.add(x);
			}
			else if(!in &&  (x.getTargetid()==targetid) && (x.out==true))
			{
				nodes.add(x);
			}
		}


		return nodes;
	}

	public static void addVirtualBase(int basenodeid, int vid, ArrayList<TargetNode> targets) {

		TargetNode basenode = getTargetNode(basenodeid, targets);

		TargetNode virtualbase = new TargetNode(vid, 0);

		ArrayList<TargetNode> pathtobase = new ArrayList<TargetNode>();



		virtualbase.setPath(basenode, pathtobase);
		virtualbase.setPathUtility(basenode, 0.0);
		virtualbase.addNeighbor(basenode);
		virtualbase.addDistance(basenode, 0.0);
		targets.add(virtualbase);

		basenode.addNeighbor(virtualbase);
		basenode.addDistance(virtualbase, 0.0);
		basenode.setPath(virtualbase, pathtobase);
		basenode.setPathUtility(virtualbase, 0.0);



	}

	public static int calculateEdgesWithCoin() {
		/**
		 * remove the dominated targets from targets arraylist
		 */
		ArrayList<Integer[]> donepairs = new ArrayList<Integer[]>();
		int qcounter = 0;

		for(int i=0; i<targets.size(); i++)
		{
			for(TargetNode neighbor : targets.get(i).getNeighbors())
			{

				if((targets.get(i).getPathUtility(neighbor)>0) && 
						!isEdgeDone(targets.get(i).getTargetid(), neighbor.getTargetid(), donepairs))
				{
					donepairs.add(new Integer[]{targets.get(i).getTargetid(), neighbor.getTargetid()});
					qcounter++;
				}

			}
		}
		return qcounter;
	}

	private static boolean isEdgeDone(int targetid, int targetid2,
			ArrayList<Integer[]> donepairs) {

		for(Integer[] x: donepairs)
		{
			if((x[0]==targetid && x[1]== targetid2) || (x[1]==targetid && x[0]== targetid2))
			{
				return true;
			}

		}

		return false;
	}

	public static HashMap<Integer, Integer[]> fillKToEdge(int nodewithcoins) {
		ArrayList<Integer[]> donepairs = new ArrayList<Integer[]>();
		HashMap<Integer, Integer[]> map = new HashMap<Integer, Integer[]>();// [i, j, utility(i,j)]
		int counter = nodewithcoins; // edges in v starts from targets.size
		for(int i=0; i<targets.size(); i++)
		{
			for(TargetNode neighbor : targets.get(i).getNeighbors())
			{
				if((targets.get(i).getPathUtility(neighbor)>0) && 
						!isEdgeDone(targets.get(i).getTargetid(), neighbor.getTargetid(), donepairs))
				{
					donepairs.add(new Integer[]{targets.get(i).getTargetid(), neighbor.getTargetid()});
					map.put(counter++, new Integer[]{targets.get(i).getTargetid(),neighbor.getTargetid(), targets.get(i).getPathUtility(neighbor).intValue() });
				}

			}
		}
		return map;



	}

	public static void removePathsToDominatedNodes(
			ArrayList<TargetNode> dominatednodes, ArrayList<TargetNode> targets) {

		for(int i=0; i<targets.size(); i++)
		{
			if(!dominatednodes.contains(targets.get(i)))
			{
				TargetNode tmp = targets.get(i);
				ArrayList<TargetNode> neibors = tmp.getNeighbors();
				for(int j=0; j<dominatednodes.size(); j++)
				{
					if(neibors.contains(dominatednodes.get(j)))
					{
						tmp.removeDistance(dominatednodes.get(j));
						tmp.removeNeighbor(dominatednodes.get(j));
						tmp.removePath(dominatednodes.get(j));

					}
				}
			}
		}

	}

	public static void removeDominatedTargets(ArrayList<TargetNode> dominatednodes, ArrayList<TargetNode> targets) {

		for(TargetNode x: dominatednodes)
		{
			if(targets.contains(x))
			{
				targets.remove(x);
			}
		}

	}

	public static void printDuplicateTargets() {

		for(int i=0; i<duplicatetargets.size(); i++)
		{
			TargetNode inode = duplicatetargets.get(i);
			System.out.println("id:"+inode.tranformedtargetid+ ", orig id: "+inode.getTargetid()+", in:"+inode.isIn()+",out:"+inode.isOut() );
			System.out.println("Neighbors: ");
			for(TargetNode nei: inode.getNeighbors())
			{
				System.out.println("id:"+nei.tranformedtargetid+ ", orig id: "+nei.getTargetid()+", in:"+nei.isIn()+",out:"+nei.isOut() );

			}
			System.out.println("##");
		}

	}

	public static HashMap<Integer, Integer> calculateNodesWithCoin(ArrayList<TargetNode> targets) {

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		int counter=0;
		for(TargetNode n : targets)
		{
			if(n.getAnimaldensity()>0)
			{
				System.out.println("Adding node "+ n.getTargetid()+ ", u ="+n.getAnimaldensity()+", to valuable entity");
				map.put(counter++, n.getTargetid());
			}
		}

		return map;




	}

	public static TargetNode getNodeWithCoin(Integer targetid) {

		for(int i=0; i<targets.size(); i++)
		{
			if(targets.get(i).getTargetid()==targetid)
			{
				return targets.get(i);
			}
		}
		return null;


	}

	public static void MPI() {

		int gridsize=24;
		int procgridsize=3;


		int[] displs = new int[9];
		int disp=0;


		for (int i=0; i<procgridsize; i++) {
			for (int j=0; j<procgridsize; j++) {
				displs[i*procgridsize+j] = disp;
				System.out.println("displs["+(i*procgridsize+j)+"]="+disp);
				disp += 1;
			}
			disp += ((gridsize/procgridsize)-1)*procgridsize;
		}


	}

	public static void basicInstantAbstractionAlgorithmWithExtrmPruningTest(double[][] density, int ITER, int nrow, 
			int ncol, int[] percentages, double dmax) throws Exception 
	{

		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				//long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					ArrayList<TargetNode> domindatednodes= new ArrayList<TargetNode>();
					preComputeShortestPaths(0, targets, domindatednodes);
					Date start = new Date();
					long l1 = start.getTime();

					ArrayList<Integer> greedypath = greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}



					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicInstantAbstractionWithExtreamPruning(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= targets.size();
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;
					solvingtime+=diff;
					revmaptime+=diff;



				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile((int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, Long.parseLong(x));
			}

		}
		//writeInFile(result,contractionsizes);


	}


	public static void basicInstantAbstractionAlgorithmTest(double[][] density, int ITER, int nrow, 
			int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{

		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();

					int base=0;
					ArrayList<Integer> greedypath = buildGreedyCover(targets, dmax, nTargets, base); //greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}
					 */


					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicInstantAbstraction(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));



				}
				//double avgtime = sumtime/(ITER);
				//double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("6",(int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}





	public static void basicGCSingleInstantAllPathsLPTest(double[][] density, int ITER, int nrow, 
			int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{

		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();

					int base=0;
					ArrayList<Integer> greedypath = buildGreedyCover(targets, dmax, nTargets, base); //greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}
					 */


					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicGCSingleInstantAllPathsLP(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));



				}
				//double avgtime = sumtime/(ITER);
				//double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("11",(int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}
	
	public static void basicAPSPGCMultiInstantAllPathsLPTest(double[][] density, int ITER, int nrow, 
			int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{

		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();

					int base=0;
					ArrayList<Integer> greedypath = buildGreedyCoverMultRes(targets, dmax, nTargets, base, nRes); //greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}
					 */


					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicAPSPGCMultiInstantAllPathsLP(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));



				}
				//double avgtime = sumtime/(ITER);
				//double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("12",(int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}



	

	public static void basicGCMultiInstantAllPathsLPTest(double[][] density, int ITER, int nrow, 
			int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{

		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();

					int base=0;
					ArrayList<Integer> greedypath = buildGreedyCoverMultRes(targets, dmax, nTargets, base, nRes); //greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}
					 */


					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicGCMultiInstantAllPathsLP(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));



				}
				//double avgtime = sumtime/(ITER);
				//double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("12",(int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}

	



	public static void basicGCMultiInstantGP3LPTest(double[][] density, int ITER, int nrow, 
			int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{

		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();

					int base=0;
					ArrayList<Integer> greedypath = buildGreedyCoverMultRes(targets, dmax, nTargets, base, nRes); //greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}
					 */


					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicGCMultiInstantGP3LP(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));



				}
				//double avgtime = sumtime/(ITER);
				//double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("13",(int)targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}




	public static void basicSeqAbstractionWithLPGreedyCover1AlgorithmTest(double[][] density, int ITER, int nrow, int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{



		//int[] contractint ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		int base=0;

		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//	makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ (iter+1));
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();


					/**
					 * we can change to a faster heuristic
					 */


					// 1. compute all pair shortest path

					ArrayList<Integer> greedypath = buildGreedyCover(targets, dmax, nTargets, base);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					System.out.println("\nInitial THreshold: "+ threshold);
					System.out.println("initial  tcur: ");
					for(Integer x: greedypath)
					{
						System.out.print(x+" ");
					}
					System.out.println();

					//printtargets(targets);



					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicSeqAbstractionLP(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));


				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("4",(int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}



	public static void basicSeqAbstractionWithLPGreedyCoverMultResAlgorithmTest(double[][] density, int ITER, int nrow, int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{



		//int[] contractint ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		int base=0;

		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//	makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ (iter+1));
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();


					/**
					 * we can change to a faster heuristic
					 */


					// 1. compute all pair shortest path

					ArrayList<Integer> greedypath = buildGreedyCoverMultRes(targets, dmax, nTargets, base, nRes);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					System.out.println("\nInitial THreshold: "+ threshold);
					System.out.println("initial  tcur: ");
					for(Integer x: greedypath)
					{
						System.out.print(x+" ");
					}
					System.out.println();

					//printtargets(targets);



					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicSeqAbstractionLP(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));


				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("5",(int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}






	public static void basicSeqAbstractionAlgorithmTest(double[][] density, int ITER, int nrow, int ncol, int[] percentages, double dmax, int nRes) throws Exception 
	{



		//int[] contractint ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		int base=0;

		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4]; //SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//	makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ (iter+1));
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();


					/**
					 * we can change to a faster heuristic
					 */


					// 1. compute all pair shortest path

					ArrayList<Integer> greedypath = buildGreedyCover(targets, dmax, nTargets, base);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					System.out.println("\nInitial THreshold: "+ threshold);
					System.out.println("initial  tcur: ");
					for(Integer x: greedypath)
					{
						System.out.print(x+" ");
					}
					System.out.println();

					//printtargets(targets);



					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.basicSeqAbstraction(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));


				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("3",(int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}






	public static void noContractionNoColumnGenerationTest(double[][] density, int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps) throws Exception 
	{



		//int[] contractint ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		int nTargets = nrow*ncol;
		int base=0;

		//for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{

					ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
					HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
					
					
					
					
					//printNodesWithNeighborsAndPath(targetmaps);

					int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					
					gamedata = constructGameData(targets);
					double [] result = new double[2];

					//int[][] gamedata = new int[nTargets][4]; //SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//	makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ (iter+1));
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					//SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					//assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();


					/**
					 * we can change to a faster heuristic
					 */


					// 1. compute all pair shortest path

					//ArrayList<Integer> greedypath = buildGreedyCover(targets, dmax, nTargets, base);


					threshold = 0;//initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					System.out.println("\nInitial THreshold: "+ threshold);
					System.out.println("initial  tcur: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+" ");
					}
					System.out.println();*/

					//printtargets(targets);



					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.noContractionNoColumnGeneration(targets, threshold, dmax, gamedata, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));


				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				//writeInFile("baseline",(int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
				
				SecurityGameContraction.writeInFile("baseline",(int)targetsize/ITER, avgsol, 0,solvingtime/ITER, 0, sumtime/ITER, nTargets);
				
				
			}

		}
		//writeInFile(result,contractionsizes);


	}
	
	


	public static void baselineForGroupTest( int ITER, int nTargets, double dmax, int nRes, 
			HashMap<Integer,ArrayList<TargetNode>> alltargets, HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps) throws Exception 
	{



		//int[] contractint ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};
		double threshold =0;

		//int nUnaccesstargets = percentage;
		//int nRes = 2;
		//int nTargets = nrow*ncol;
		int base=0;

		//for(double percentage: percentages)
		{
			//for(double threshold: thresholds)
			{

				double sumsol = 0;
				long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					//targets.clear();
					double [] result = new double[2];

					int[][] gamedata = new int[nTargets][4]; //SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					
					//ArrayList<Integer>[] clus = allclus.get(iter);
					ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
					HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
					
					
					//	makeZeroSum(gamedata,nTargets);
					gamedata = GroupingTargets.buildGamedata(targetmaps);

					System.out.println("\n Iter "+ (iter+1));
					//System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					//SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					
					
					//assignRandomDensityZeroSum(density, gamedata, targets, iter);

					Date start = new Date();
					long l1 = start.getTime();


					/**
					 * we can change to a faster heuristic
					 */


					// 1. compute all pair shortest path

					//ArrayList<Integer> greedypath = buildGreedyCover(targets, dmax, nTargets, base);


					threshold = 0;//initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					System.out.println("\nInitial THreshold: "+ threshold);
					System.out.println("initial  tcur: ");
					/*for(Integer x: greedypath)
					{
						System.out.print(x+" ");
					}
					System.out.println();*/

					//printtargets(targets);



					//sgc.contractGraph(domindatednodes, targets);

					result =	SecurityGameContraction.baseline(targets, threshold, dmax, gamedata, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= result[5];
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += result[3];
					solvingtime+= result[4];
					sumtime+=diff;

					//writeInFile(percentage,(int)sumthreshold/(iter+1), sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));
					//writeInFile(percentage,(int)result[1], sumiter/(iter+1),(int)targetsize/(iter+1), result[0], contractiontime/(iter+1), solvingtime/(iter+1), solvingtime/(iter+1));


				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile("baseline",(int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, sumtime/ITER);
			}

		}
		//writeInFile(result,contractionsizes);


	}





	private static ArrayList<Integer> buildGreedyCover(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);



		ArrayList<Integer> greedypath = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);

		greedypath =  greedyCover(base, targets, dmax, targetssorted, apspmat, map,mapback);

		return greedypath;
	}
	
	public static ArrayList<Integer> buildOneGreedyPathWithSrcDest(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int dest, int nRes, int[] disallocation) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		purifyAPSPMatrixZero(apspmat, targets, nTargets, map, mapback);
		



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);
		
		
		int s = map.get(base);
		int d = map.get(dest);
		dmax = apspmat[s][d];
		disallocation[0] = (int)dmax;

		tcur =  greedyOnePathWithSrcDest(base, dest, targets, dmax, targetssorted, apspmat, map,mapback, nRes);

		return tcur;
	}
	
	
	
	public static ArrayList<Integer> buildOneGreedyPathWithSrcDest(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int dest, int nRes) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		purifyAPSPMatrixZero(apspmat, targets, nTargets, map, mapback);
		



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);
		
		
		int s = map.get(base);
		int d = map.get(dest);
		

		tcur =  greedyOnePathWithSrcDest(base, dest, targets, dmax, targetssorted, apspmat, map,mapback, nRes);

		return tcur;
	}
	
	

	private static ArrayList<Integer> buildGreedyCoverWithSrcDest(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int dest, int nRes) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		purifyAPSPMatrixZero(apspmat, targets, nTargets, map, mapback);
		



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);

		tcur =  greedyCoverWithSrcDest(base, dest, targets, dmax, targetssorted, apspmat, map,mapback, nRes);

		return tcur;
	}





	public static ArrayList<Integer> buildGreedyCoverMultRes(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int nRes) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		purifyAPSPMatrixZeroGT(apspmat, targets, nTargets, map, mapback);
		



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);

		tcur =  greedyCoverMultRes(base, targets, dmax, targetssorted, apspmat, map,mapback, nRes);

		return tcur;
	}
	
	
	





	private static ArrayList<ArrayList<Integer>> buildGreedyPathMultRes(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int nRes) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);




		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


	//	System.out.println("93y6903486 "+ adjacencymatrix[map.get(34)][map.get(46)]);
		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);

		purifyAPSPMatrixWithNeighborAndZero(apspmat, targets, nTargets, map, mapback);

		//System.out.println("93y6903486 "+ apspmat[map.get(34)][map.get(46)]);

		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>(); //greedyFirstRoute(dmax,gamedata, targets);
		int[][] targetssorted = sortTargets(targets);

		Random r = new Random();
		for(int i=0; i<100; i++)
		{

			shuffle(targetssorted,r);
			ArrayList<Integer> path =  greedyPathMultRes(base, targets, dmax, targetssorted, apspmat, map,mapback, nRes);

			if(path.size()>=3)
			{
				paths.add(path);
			}
		}
		return paths;
	}







	private static void shuffle(int[][] targetssorted, Random r) {


		int i=1;

		while(i<targetssorted.length)
		{

			int j = r.nextInt(targetssorted.length);
			if(j>0 && j< targetssorted.length)
			{
				swap(targetssorted, i, j);
				i++;
			}
		}



	}

	private static void swap(int[][] targetssorted, int i, int j) {




		int [] tmp = {0,0};

		tmp = targetssorted[i];
		targetssorted[i] = targetssorted[j];
		targetssorted[j] = tmp;

	}

	private static ArrayList<ArrayList<Integer>> buildGreedyCoverMultRes2(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, 
			int nRes) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);

		purifyAPSPMatrixWithNeighborAndZero(apspmat, targets, nTargets, map, mapback);





		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);



		return greedyCoverMultRes2(base, targets, dmax, targetssorted, apspmat, map,mapback, nRes);
	}




	private static ArrayList<ArrayList<Integer>> buildGreedyPathMultRes2(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, 
			int nRes) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);

		purifyAPSPMatrixWithNeighborAndZero(apspmat, targets, nTargets, map, mapback);





		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);



		return greedyCoverMultRes2(base, targets, dmax, targetssorted, apspmat, map,mapback, nRes);
	}


	public static ArrayList<ArrayList<Integer>> buildSuperGreedyCoverMultRes2(
			HashMap<Integer, TargetNode> targetmaps, double dmax, int nTargets, int base, 
			int nRes, HashMap<Integer,Double> attackerstrategy, HashMap<Integer, SuperTarget> sts, 
			HashMap<Integer,Double> dstravel) {



		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
		
		
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		//ArrayList<Integer> graphint = new ArrayList<Integer>();
		
		
		int i = 1;
		for(Integer stid: sts.keySet())
		{
			map.put(stid, i);
			mapback.put(i, stid);
			//graphint.add(map.get(targets.get(i-1).getTargetid()));
			i++;
		}
		
		
		
		int[][] adjacencymatrix = new int[sts.size()+1][sts.size()+1];
		makeAdjacencyMatrixST(adjacencymatrix , sts, sts.size(), map, mapback);
		
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(sts.size());
		int[][] apsp = allPairShortestPath.allPairShortestPath(adjacencymatrix);
        //ArrayList<Integer> path = AllPairShortestPath.getPath(1, 2, allPairShortestPath.next);

		purifyAPSPSuperMatrixWithNeighborAndZero(apsp, targetmaps, sts.size(), map, mapback, sts);
		
		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortSuperTargets(sts);



		return greedySuperCoverMultRes2(base, targetmaps, dmax, targetssorted, apsp, map,mapback, nRes, attackerstrategy, sts, dstravel);
	}
	
	
	public static ArrayList<ArrayList<Integer>> buildGreedyCoverMultRes2(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, 
			int nRes, HashMap<Integer,Double> attackerstrategy) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);

		purifyAPSPMatrixWithNeighborAndZero(apspmat, targets, nTargets, map, mapback);





		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);



		return greedyCoverMultRes2(base, targets, dmax, targetssorted, apspmat, map,mapback, nRes, attackerstrategy);
	}


	private static ArrayList<ArrayList<Integer>> buildGreedyCoverMultResSrcDest(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int dest,
			int nRes, HashMap<Integer,Double> attackerstrategy) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		int icount =1;
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;
		}
		//makeAdjacencyMatrix(adjacencymatrix,targets,nTargets,map, mapback);

		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);


		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);

		purifyAPSPMatrixWithNeighborAndZero(apspmat, targets, nTargets, map, mapback);





		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = sortTargets(targets);



		return greedyCoverMultResSrcDest(base, dest, targets, dmax, targetssorted, apspmat, map,mapback, nRes, attackerstrategy);
	}
	
	
	



	private static ArrayList<Integer> greedyCover(int base, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback) 
			{


		int i=0; 

		ArrayList<Integer> tcur = new ArrayList<Integer>();
		ArrayList<Integer> greedypath = new ArrayList<Integer>();

		if(targetssorted[0][0]==base)
			i=1;
		else 
			i=0;


		if(apspmat[map.get(base)][map.get(targetssorted[i][0])]>dmax)
		{
			tcur.add(base);
			tcur.add(targetssorted[i][0]);
			greedypath.add(base);
			System.out.println();
			for(Integer n: greedypath)
			{
				System.out.print(n+" ");
			}
			System.out.println();
			return tcur;

		}

		greedypath.add(base);
		greedypath.add(targetssorted[i][0]);
		greedypath.add(base);


		int totaldist = apspmat[map.get(base)][map.get(targetssorted[i][0])] + 
				apspmat[map.get(targetssorted[i][0])][map.get(base)];

		i=i+1;

		while(i<=targets.size())
		{
			int besttotaldisttemp = -1;
			int bestj = -1;

			for(int j=1; j< greedypath.size(); j++ )
			{
				int s = map.get(greedypath.get(j-1));
				int d = map.get(greedypath.get(j));
				int totaldisttemp = totaldist - apspmat[s][d];

				s = map.get(greedypath.get(j-1));
				d = map.get(targetssorted[i][0]);
				totaldisttemp +=  apspmat[s][d];


				s = map.get(targetssorted[i][0]);
				d = map.get(greedypath.get(j));
				totaldisttemp +=  apspmat[s][d];

				if ((totaldisttemp<dmax) && 
						((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) )
				{
					bestj = j;
				}
			}
			if(bestj==-1)
			{
				int ulimit = targetssorted[i][1];

				for(int k=0; k<targetssorted.length; k++)
				{
					if(targetssorted[k][1]>=ulimit)
					{
						tcur.add(targetssorted[k][0]);
					}
				}
				System.out.println();
				for(Integer n: greedypath)
				{
					System.out.print(n+" ");
				}
				System.out.println();
				return tcur;
			}
			insertTsrt(greedypath, 0, bestj-1, targetssorted[i][0], bestj, greedypath.size()-1);

			int s = map.get(greedypath.get(bestj-1));
			int d = map.get(greedypath.get(bestj));
			totaldist = totaldist - apspmat[s][d];


			s = map.get(greedypath.get(bestj-1));
			d = map.get(targetssorted[i][0]);
			totaldist += apspmat[s][d];


			s = map.get(targetssorted[i][0]); 
			d = map.get(greedypath.get(bestj));
			totaldist += apspmat[s][d];

			i=i+1;
		} // end of while loop


		System.out.println();
		for(Integer n: greedypath)
		{
			System.out.print(n+" ");
		}
		System.out.println();



		return null;
			}




	public static ArrayList<Integer> greedyCoverMultRes(int base, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, int nRes) 
			{


		int i=0; 

		ArrayList<Integer> tcur = new ArrayList<Integer>();
		ArrayList<Integer> greedypath = new ArrayList<Integer>();

		if(targetssorted[0][0]==base)
			i=1;
		else 
			i=0;



		for(int res = 0; res<nRes; res++)
		{
			/*System.out.println("res = "+res);
			System.out.println("Adding base to greedy path : ");*/
			greedypath.clear();
			greedypath.add(base);
			greedypath.add(base);
			//printGreedyPath(greedypath);

			int totaldist = 0;
			boolean continueFlag= false;
			while (i<targetssorted.length && continueFlag==false)
			{
				int besttotaldisttemp = -1;
				int bestj = -1;
				/*System.out.println("i = "+i);
				System.out.println("continueFlag : "+continueFlag);
				System.out.println("besttotaldisttemp = "+besttotaldisttemp);
				System.out.println("bestj : "+bestj);*/

				for(int j=1; j< greedypath.size(); j++ )
				{
					int s = map.get(greedypath.get(j-1));
					int d = map.get(greedypath.get(j));
					int totaldisttemp = totaldist - apspmat[s][d];

					s = map.get(greedypath.get(j-1));
					d = map.get(targetssorted[i][0]);
					totaldisttemp +=  apspmat[s][d];


					s = map.get(targetssorted[i][0]);
					d = map.get(greedypath.get(j));
					totaldisttemp +=  apspmat[s][d];

					/*System.out.println("totaldisttemp = "+totaldisttemp);
					System.out.println("bestj : "+bestj);*/

					if ((totaldisttemp<dmax) && 
							((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) )
					{
						bestj = j;
						besttotaldisttemp = totaldisttemp;
						/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
						System.out.println("updating bestj : "+bestj);*/
					}
				}
				if(bestj==-1)
				{
					//System.out.println("could not insert, bestj : "+bestj);
					continueFlag = true;
					break;
				}
				//System.out.println("inserting to targetssorted : "+targetssorted[i][0]+" , greedy path : ");
				greedypath = insertTsrt(greedypath, 0, bestj-1, targetssorted[i][0], bestj, greedypath.size()-1);
				//printGreedyPath(greedypath);
				totaldist = besttotaldisttemp;
				i=i+1;
			}//while loop

			if(i==targetssorted.length)
			{
				/*System.out.println("i= "+i+"=targetssorted.length");
				System.out.println("All targets can be added to tcur: ");*/
				for(int k=0; k<targetssorted.length; k++)
				{
					tcur.add(targetssorted[k][0]);
				}
				//printGreedyPath(tcur);
				return tcur;
			}
			//System.out.println("res= "+res+", greedypath : ");
			//printGreedyPath(greedypath);
		}//end of for loop



		int ulimit = targetssorted[i][1];
		/*System.out.println("ulimit "+ulimit);
		System.out.println("tcur : ");*/

		for(int k=0; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				tcur.add(targetssorted[k][0]);
			}
		}
		//printGreedyPath(tcur);
		return tcur;


			}
	
	
	
	
	private static ArrayList<Integer> greedyCoverWithSrcDest(int base, int dest, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, int nRes) 
			{
		
		
		
		
		
		
		


		int i=0; 

		ArrayList<Integer> tcur = new ArrayList<Integer>();
		ArrayList<Integer> greedypath = new ArrayList<Integer>();

		if(targetssorted[0][0]==base)
			i=1;
		else 
			i=0;



		for(int res = 0; res<nRes; res++)
		{
			/*System.out.println("res = "+res);
			System.out.println("Adding base to greedy path : ");*/
			greedypath.clear();
			greedypath.add(base);
			greedypath.add(dest);
			//printGreedyPath(greedypath);

			int totaldist = 0;
			boolean continueFlag= false;
			while ( i<targetssorted.length && continueFlag==false)
			{
				
				int besttotaldisttemp = -1;
				int bestj = -1;
				/*System.out.println("i = "+i);
				System.out.println("continueFlag : "+continueFlag);
				System.out.println("besttotaldisttemp = "+besttotaldisttemp);
				System.out.println("bestj : "+bestj);*/

				for(int j=1; j< greedypath.size(); j++ )
				{
					int s = map.get(greedypath.get(j-1));
					int d = map.get(greedypath.get(j));
					int totaldisttemp = totaldist - apspmat[s][d];

					s = map.get(greedypath.get(j-1));
					d = map.get(targetssorted[i][0]);
					totaldisttemp +=  apspmat[s][d];


					s = map.get(targetssorted[i][0]);
					d = map.get(greedypath.get(j));
					totaldisttemp +=  apspmat[s][d];

					/*System.out.println("totaldisttemp = "+totaldisttemp);
					System.out.println("bestj : "+bestj);*/

					if (totaldisttemp>0 &&  (totaldisttemp<dmax) && 
							((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) )
					{
						bestj = j;
						besttotaldisttemp = totaldisttemp;
						/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
						System.out.println("updating bestj : "+bestj);*/
					}
				}
				if(bestj==-1)
				{
					//System.out.println("could not insert, bestj : "+bestj);
					continueFlag = true;
					break;
				}
				//System.out.println("inserting to targetssorted : "+targetssorted[i][0]+" , greedy path : ");
				greedypath = insertTsrt(greedypath, 0, bestj-1, targetssorted[i][0], bestj, greedypath.size()-1);
				//printGreedyPath(greedypath);
				totaldist = besttotaldisttemp;
				i=i+1;
			}//while loop

			if(i==targetssorted.length)
			{
				/*System.out.println("i= "+i+"=targetssorted.length");
				System.out.println("All targets can be added to tcur: ");*/
				for(int k=0; k<targetssorted.length; k++)
				{
					tcur.add(targetssorted[k][0]);
				}
				//printGreedyPath(tcur);
				return tcur;
			}
			//System.out.println("res= "+res+", greedypath : ");
			//printGreedyPath(greedypath);
		}//end of for loop


		
		
		

		int ulimit = targetssorted[i][1];
		/*System.out.println("ulimit "+ulimit);
		System.out.println("tcur : ");*/

		
		//if(greedypath.size()==2)
		
			//fidn the min
			// set the ulimit
			
		
			double min = getTargetNode(base, targets).attackerreward;
			for(int p=0; p<greedypath.size(); p++)
			{
			
			if(getTargetNode(greedypath.get(p), targets).attackerreward<min)
				min = (int)getTargetNode(greedypath.get(p), targets).attackerreward;
			
		
			}
		ulimit=(int)min;
		
		for(int k=0; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				tcur.add(targetssorted[k][0]);
			}
		}
		//printGreedyPath(tcur);
		return tcur;


			}

	
	
	public static ArrayList<Integer> greedyOnePathWithSrcDest(int base, int dest, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, int nRes) 
			{
		
		
		
		
		
		
		


		int i=0; 

		ArrayList<Integer> tcur = new ArrayList<Integer>();
		ArrayList<Integer> greedypath = new ArrayList<Integer>();

		if(targetssorted[0][0]==base)
			i=1;
		else 
			i=0;



		for(int res = 0; res<nRes; res++)
		{
			/*System.out.println("res = "+res);
			System.out.println("Adding base to greedy path : ");*/
			greedypath.clear();
			greedypath.add(base);
			greedypath.add(dest);
			//printGreedyPath(greedypath);

			int totaldist = 0;
			if(base!=dest)
			{
				totaldist = apspmat[map.get(base)][map.get(dest)];
			}
			boolean continueFlag= false;
			while ( i<targetssorted.length && continueFlag==false)
			{
				
				int besttotaldisttemp = -1;
				int bestj = -1;
				/*System.out.println("i = "+i);
				System.out.println("continueFlag : "+continueFlag);
				System.out.println("besttotaldisttemp = "+besttotaldisttemp);
				System.out.println("bestj : "+bestj);*/
				
				System.out.println("trying to insert  : "+targetssorted[i][0]);

				for(int j=1; j< greedypath.size(); j++ )
				{
					int s = map.get(greedypath.get(j-1));
					int d = map.get(greedypath.get(j));
					int totaldisttemp = totaldist - apspmat[s][d];

					s = map.get(greedypath.get(j-1));
					d = map.get(targetssorted[i][0]);
					totaldisttemp +=  apspmat[s][d];


					s = map.get(targetssorted[i][0]);
					d = map.get(greedypath.get(j));
					totaldisttemp +=  apspmat[s][d];

					/*System.out.println("totaldisttemp = "+totaldisttemp);
					System.out.println("bestj : "+bestj);*/

					if (totaldisttemp>0 &&  (totaldisttemp<=dmax) && 
							((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) )
					{
						bestj = j;
						besttotaldisttemp = totaldisttemp;
						/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
						System.out.println("updating bestj : "+bestj);*/
					}
				}
				if(bestj!=-1)
				{
					//System.out.println("could not insert, bestj : "+bestj);
					//continueFlag = true;
					//break;
					System.out.println("inserting : "+targetssorted[i][0]+" , greedy path : ");
					greedypath = insertTsrt(greedypath, 0, bestj-1, targetssorted[i][0], bestj, greedypath.size()-1);
					//printGreedyPath(greedypath);
					totaldist = besttotaldisttemp;
				}
				
				i=i+1;
			}//while loop
			
			
			
			
			

			/*if(i==targetssorted.length)
			{
				System.out.println("i= "+i+"=targetssorted.length");
				System.out.println("All targets can be added to tcur: ");
				for(int k=0; k<targetssorted.length; k++)
				{
					tcur.add(targetssorted[k][0]);
				}
				//printGreedyPath(tcur);
				return tcur;
			}*/
			//System.out.println("res= "+res+", greedypath : ");
			//printGreedyPath(greedypath);
		}//end of for loop


		return greedypath;
		
		

		/*int ulimit = targetssorted[i][1];
		System.out.println("ulimit "+ulimit);
		System.out.println("tcur : ");

		
		//if(greedypath.size()==2)
		
			//fidn the min
			// set the ulimit
			
		
			double min = getTargetNode(base, targets).attackerreward;
			for(int p=0; p<greedypath.size(); p++)
			{
			
			if(getTargetNode(greedypath.get(p), targets).attackerreward<min)
				min = (int)getTargetNode(greedypath.get(p), targets).attackerreward;
			
		
			}
		ulimit=(int)min;
		
		for(int k=0; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				tcur.add(targetssorted[k][0]);
			}
		}
		//printGreedyPath(tcur);
		return tcur;*/


			}







	private static ArrayList<Integer> greedyPathMultRes(int base, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, int nRes) 
			{


		int i=0; 

		ArrayList<Integer> tcur = new ArrayList<Integer>();
		ArrayList<Integer> greedypath = new ArrayList<Integer>();

		if(targetssorted[0][0]==base)
			i=1;
		else 
			i=0;



		for(int res = 0; res<nRes; res++)
		{
			//System.out.println("res = "+res);
			//System.out.println("Adding base to greedy path : ");
			greedypath.clear();
			greedypath.add(base);
			greedypath.add(base);
			//printGreedyPath(greedypath);

			int totaldist = 0;
			boolean continueFlag= false;
			while (i<targetssorted.length && continueFlag==false)
			{
				int besttotaldisttemp = -1;
				int bestj = -1;
				//System.out.println("i = "+i);
				//System.out.println("continueFlag : "+continueFlag);
				//System.out.println("besttotaldisttemp = "+besttotaldisttemp);
				//System.out.println("bestj : "+bestj);

				for(int j=1; j< greedypath.size(); j++ )
				{
					int s = map.get(greedypath.get(j-1));
					int d = map.get(greedypath.get(j));
					int totaldisttemp = totaldist - apspmat[s][d];

					s = map.get(greedypath.get(j-1));
					d = map.get(targetssorted[i][0]);
					totaldisttemp +=  apspmat[s][d];


					s = map.get(targetssorted[i][0]);
					d = map.get(greedypath.get(j));
					totaldisttemp +=  apspmat[s][d];

					//System.out.println("totaldisttemp = "+totaldisttemp);
					//System.out.println("bestj : "+bestj);

					if ((totaldisttemp<dmax) && 
							((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) && totaldisttemp>0 )
					{
						bestj = j;
						besttotaldisttemp = totaldisttemp;
						System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
						System.out.println("updating bestj : "+bestj);
					}
				}
				if(bestj==-1)
				{
					//System.out.println("could not insert, bestj : "+bestj);
					//continueFlag = true;
					//break;
				}
				else
				{
					System.out.println("inserting to targetssorted : "+targetssorted[i][0]+" , greedy path : ");
					//printGreedyPath(greedypath);
					greedypath = insertTsrt(greedypath, 0, bestj-1, targetssorted[i][0], bestj, greedypath.size()-1);
					//printGreedyPath(greedypath);
					totaldist = besttotaldisttemp;
				}
				i=i+1;
			}//while loop

			if(i==targetssorted.length)
			{
				//System.out.println("i= "+i+"=targetssorted.length");
				//System.out.println("All targets can be added to tcur: ");
				for(int k=0; k<targetssorted.length; k++)
				{
					tcur.add(targetssorted[k][0]);
				}
				//printGreedyPath(tcur);
				return greedypath;
			}
			//System.out.println("res= "+res+", greedypath : ");
			//printGreedyPath(greedypath);
		}//end of for loop



		int ulimit = targetssorted[i][1];
		//System.out.println("ulimit "+ulimit);
		//System.out.println("tcur : ");

		for(int k=0; k<targetssorted.length; k++)
		{
			if(targetssorted[k][1]>=ulimit)
			{
				tcur.add(targetssorted[k][0]);
			}
		}
		//printGreedyPath(tcur);
		return greedypath;


			}






	static void shuffleArray(int[] ar)
	{
		// If running on Java 6 or older, use `new Random()` on RHS here
		Random rnd = ThreadLocalRandom.current();
		for (int i = ar.length - 1; i > 0; i--)
		{
			int index = rnd.nextInt(i + 1);
			// Simple swap
			int a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}
	}
	
	
	
	private static ArrayList<ArrayList<Integer>> greedyCoverMultRes2(int base, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int nRes, HashMap<Integer,Double> attackerstrategy) 
			{


		ArrayList<Integer> tsrt = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> bestjointgr = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> jointgr = new ArrayList<ArrayList<Integer>>();
		double besttotalcoin = -1;
		//int[] coin = new int[targetssorted.length];

		HashMap<Integer, Double> coin = new HashMap<Integer, Double>();

		for(int k=0; k<targetssorted.length; k++)
		{
			//coin[targetssorted[k][0]]=;
			coin.put(targetssorted[k][0], targetssorted[k][1]*attackerstrategy.get(targetssorted[k][0]));
		}
		for(int k=1; k<targetssorted.length; k++)
		{
			tsrt.add(targetssorted[k][0]);
		}
		/*System.out.println("Tsrt : ");
		printGreedyPath(tsrt);*/

		for(int iter = 0; iter<10; iter++)
		{
			ArrayList<Integer> tmptsrt = new ArrayList<Integer>();
			if(iter==0)
			{
				tmptsrt=tsrt;
			}
			else
			{
				//shuffleArray(tmptsrt);
				tmptsrt = new ArrayList<Integer>(tsrt);
				Collections.shuffle(tmptsrt, new Random());
			}
			int i=0; 


			double totalcoin = 0;
			HashMap<Integer, Double> remaincoin = coin;
			ArrayList<Integer> greedypath = new ArrayList<Integer>();
			for(int res = 0; res<nRes; res++)
			{
				//System.out.println("res = "+res);
				//System.out.println("Adding base to greedy path : ");
				greedypath.clear();
				greedypath.add(base);
				greedypath.add(base);
				//printGreedyPath(greedypath);
				int totaldist = 0;
				totalcoin = totalcoin + remaincoin.get(base);
				remaincoin.remove(base);
				remaincoin.put(base, 0.0);
				while(i<tsrt.size())
				{
					int besttotaldisttemp = -1;
					int bestj = -1;
					for(int j=1; j<greedypath.size(); j++)
					{
						int s = map.get(greedypath.get(j-1));
						int d = map.get(greedypath.get(j));
						int totaldisttemp = totaldist - apspmat[s][d];

						s = map.get(greedypath.get(j-1));
						//d = map.get(targetssorted[i][0]);
						d = map.get(tmptsrt.get(i));
						totaldisttemp +=  apspmat[s][d];


						s = map.get(tmptsrt.get(i));
						d = map.get(greedypath.get(j));
						totaldisttemp +=  apspmat[s][d];

						/*System.out.println("totaldisttemp = "+totaldisttemp);
						System.out.println("bestj : "+bestj);*/

						if ((totaldisttemp<dmax) && 
								((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) && totaldisttemp>0)
						{
							bestj = j; // trying to inseert as much target as possible by taking target which will result in minimum distance 
							/**
							 * a b c d ....we will try to insert x between ab or bc or cd
							 * we will choose which will result in the minimum total distance
							 */
							besttotaldisttemp = totaldisttemp;
							/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
							System.out.println("updating bestj : "+bestj);*/
						}



					}
					if(bestj>-1)
					{
						greedypath = insertTsrt(greedypath, 0, bestj-1, tmptsrt.get(i) , bestj, greedypath.size()-1);
						if(greedypath.size()>dmax)
						{
							printGreedyPath(greedypath);
						}
						//printGreedyPath(greedypath);
						totaldist = besttotaldisttemp;

						totalcoin += remaincoin.get(tmptsrt.get(i));
						//remaincoin[tmptsrt.get(i)] = 0;
						remaincoin.remove(tmptsrt.get(i));
						remaincoin.put(tmptsrt.get(i), 0.0);
					}
					i++;
				} // while loop
				if(greedypath.size()==2 && greedypath.get(0)==base && greedypath.get(0)==base)
				{
					greedypath.clear();
					greedypath.add(base);
				}
				if(greedypath.size()>=3)
				{
					ArrayList<Integer> tmp = new ArrayList<Integer>(greedypath);
					if(greedypath.size()>dmax)
					{
						printGreedyPath(greedypath);
					}
					jointgr.add(tmp);
					//printPaths(jointgr);
				}
			}
			if(totalcoin>besttotalcoin)
			{
				bestjointgr = jointgr;
				besttotalcoin = totalcoin;
			}
			//printPaths(jointgr);

		}
		return bestjointgr;


			}
	
	
	private static ArrayList<ArrayList<Integer>> greedySuperCoverMultRes2(int base, HashMap<Integer, TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int nRes, HashMap<Integer,Double> attackerstrategy, HashMap<Integer, SuperTarget> sts, HashMap<Integer,Double> dstravel) 
			{


		ArrayList<Integer> tsrt = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> bestjointgr = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> jointgr = new ArrayList<ArrayList<Integer>>();
		double besttotalcoin = -1;
		//int[] coin = new int[targetssorted.length];

		HashMap<Integer, Double> coin = new HashMap<Integer, Double>();

		for(int k=0; k<targetssorted.length; k++)
		{
			//coin[targetssorted[k][0]]=;
			coin.put(targetssorted[k][0], targetssorted[k][1]*attackerstrategy.get(targetssorted[k][0]));
		}
		for(int k=1; k<targetssorted.length; k++)
		{
			tsrt.add(targetssorted[k][0]);
		}
		/*System.out.println("Tsrt : ");
		printGreedyPath(tsrt);*/

		for(int iter = 0; iter<10; iter++)
		{
			ArrayList<Integer> tmptsrt = new ArrayList<Integer>();
			if(iter==0)
			{
				tmptsrt=tsrt;
			}
			else
			{
				//shuffleArray(tmptsrt);
				tmptsrt = new ArrayList<Integer>(tsrt);
				Collections.shuffle(tmptsrt, new Random());
			}
			int i=0; 


			double totalcoin = 0;
			HashMap<Integer, Double> remaincoin = coin;
			ArrayList<Integer> greedypath = new ArrayList<Integer>();
			for(int res = 0; res<nRes; res++)
			{
				//System.out.println("res = "+res);
				//System.out.println("Adding base to greedy path : ");
				greedypath.clear();
				greedypath.add(base);
				greedypath.add(base);
				//printGreedyPath(greedypath);
				int totaldist = 0;
				totalcoin = totalcoin + remaincoin.get(base);
				remaincoin.remove(base);
				remaincoin.put(base, 0.0);
				while(i<tsrt.size())
				{
					int besttotaldisttemp = -1;
					int bestj = -1;
					for(int j=1; j<greedypath.size(); j++)
					{
						
						
						//TODO consider distance, intra cluster
						int s = map.get(greedypath.get(j-1));
						int d = map.get(greedypath.get(j));
						//System.out.println(" s : "+ greedypath.get(j-1) + " , d : "+greedypath.get(j) + " j-1 "+ (j-1) + " j "+ j +  " gs "+ greedypath.size());
						int totaldisttemp = totaldist - apspmat[s][d] - dstravel.get(greedypath.get(j-1)).intValue() - dstravel.get(greedypath.get(j)).intValue();

						s = map.get(greedypath.get(j-1));
						//d = map.get(targetssorted[i][0]);
						d = map.get(tmptsrt.get(i));
						//System.out.println(" s : "+ greedypath.get(j-1) + " , d : "+tmptsrt.get(i) + " j-1 "+ (j-1) + " i "+ i +  " gs "+ greedypath.size());
						totaldisttemp +=  apspmat[s][d];
						
						
						totaldisttemp +=  dstravel.get(greedypath.get(j-1)).intValue() ;
						
						totaldisttemp +=  dstravel.get(tmptsrt.get(i)).intValue();


						s = map.get(tmptsrt.get(i));
						d = map.get(greedypath.get(j));
						//System.out.println(" s : "+ greedypath.get(i) + " , d : "+greedypath.get(j) + " i "+ i + " j "+ j +  " gs "+ greedypath.size());
						totaldisttemp +=  apspmat[s][d] + dstravel.get(tmptsrt.get(i)).intValue() + dstravel.get(greedypath.get(j)).intValue();;

						/*System.out.println("totaldisttemp = "+totaldisttemp);
						System.out.println("bestj : "+bestj);*/

						if ((totaldisttemp<dmax) && 
								((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) && totaldisttemp>0)
						{
							bestj = j; // trying to inseert as much target as possible by taking target which will result in minimum distance 
							/**
							 * a b c d ....we will try to insert x between ab or bc or cd
							 * we will choose which will result in the minimum total distance
							 */
							besttotaldisttemp = totaldisttemp;
							/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
							System.out.println("updating bestj : "+bestj);*/
						}



					}
					if(bestj>-1)
					{
						greedypath = insertTsrt(greedypath, 0, bestj-1, tmptsrt.get(i) , bestj, greedypath.size()-1);
						if(greedypath.size()>dmax)
						{
							//printGreedyPath(greedypath);
						}
						//printGreedyPath(greedypath);
						totaldist = besttotaldisttemp;

						totalcoin += remaincoin.get(tmptsrt.get(i));
						//remaincoin[tmptsrt.get(i)] = 0;
						remaincoin.remove(tmptsrt.get(i));
						remaincoin.put(tmptsrt.get(i), 0.0);
					}
					i++;
				} // while loop
				if(greedypath.size()==2 && greedypath.get(0)==base && greedypath.get(0)==base)
				{
					greedypath.clear();
					greedypath.add(base);
				}
				if(greedypath.size()>=3)
				{
					ArrayList<Integer> tmp = new ArrayList<Integer>(greedypath);
					if(greedypath.size()>dmax)
					{
						printGreedyPath(greedypath);
					}
					jointgr.add(tmp);
					//printPaths(jointgr);
				}
			}
			if(totalcoin>besttotalcoin)
			{
				bestjointgr = jointgr;
				besttotalcoin = totalcoin;
			}
			//printPaths(jointgr);

		}
		return bestjointgr;


			}



	private static ArrayList<ArrayList<Integer>> greedyCoverMultResSrcDest(int base, int dest, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int nRes, HashMap<Integer,Double> attackerstrategy) 
			{


		ArrayList<Integer> tsrt = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> bestjointgr = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> jointgr = new ArrayList<ArrayList<Integer>>();
		double besttotalcoin = -1;
		//int[] coin = new int[targetssorted.length];

		HashMap<Integer, Double> coin = new HashMap<Integer, Double>();

		for(int k=0; k<targetssorted.length; k++)
		{
			//coin[targetssorted[k][0]]=;
			coin.put(targetssorted[k][0], targetssorted[k][1]*attackerstrategy.get(targetssorted[k][0]));
		}
		for(int k=1; k<targetssorted.length; k++)
		{
			tsrt.add(targetssorted[k][0]);
		}
		/*System.out.println("Tsrt : ");
		printGreedyPath(tsrt);*/

		for(int iter = 0; iter<10; iter++)
		{
			ArrayList<Integer> tmptsrt = new ArrayList<Integer>();
			if(iter==0)
			{
				tmptsrt=tsrt;
			}
			else
			{
				//shuffleArray(tmptsrt);
				tmptsrt = new ArrayList<Integer>(tsrt);
				Collections.shuffle(tmptsrt, new Random());
			}
			int i=0; 


			double totalcoin = 0;
			HashMap<Integer, Double> remaincoin = coin;
			ArrayList<Integer> greedypath = new ArrayList<Integer>();
			for(int res = 0; res<nRes; res++)
			{
				//System.out.println("res = "+res);
				//System.out.println("Adding base to greedy path : ");
				greedypath.clear();
				greedypath.add(base);
				greedypath.add(dest);
				//printGreedyPath(greedypath);
				int totaldist = 0;
				totalcoin = totalcoin + remaincoin.get(base);
				remaincoin.remove(base);
				remaincoin.put(base, 0.0);
				while(i<tsrt.size())
				{
					int besttotaldisttemp = -1;
					int bestj = -1;
					for(int j=1; j<greedypath.size(); j++)
					{
						int s = map.get(greedypath.get(j-1));
						int d = map.get(greedypath.get(j));
						int totaldisttemp = totaldist - apspmat[s][d];

						s = map.get(greedypath.get(j-1));
						//d = map.get(targetssorted[i][0]);
						d = map.get(tmptsrt.get(i));
						totaldisttemp +=  apspmat[s][d];


						s = map.get(tmptsrt.get(i));
						d = map.get(greedypath.get(j));
						totaldisttemp +=  apspmat[s][d];

						/*System.out.println("totaldisttemp = "+totaldisttemp);
						System.out.println("bestj : "+bestj);*/

						if ((totaldisttemp<dmax) && 
								((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) && totaldisttemp>0)
						{
							bestj = j; // trying to inseert as much target as possible by taking target which will result in minimum distance 
							/**
							 * a b c d ....we will try to insert x between ab or bc or cd
							 * we will choose which will result in the minimum total distance
							 */
							besttotaldisttemp = totaldisttemp;
							/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
							System.out.println("updating bestj : "+bestj);*/
						}



					}
					if(bestj>-1)
					{
						greedypath = insertTsrt(greedypath, 0, bestj-1, tmptsrt.get(i) , bestj, greedypath.size()-1);
						if(greedypath.size()>dmax)
						{
							printGreedyPath(greedypath);
						}
						//printGreedyPath(greedypath);
						totaldist = besttotaldisttemp;

						totalcoin += remaincoin.get(tmptsrt.get(i));
						//remaincoin[tmptsrt.get(i)] = 0;
						remaincoin.remove(tmptsrt.get(i));
						remaincoin.put(tmptsrt.get(i), 0.0);
					}
					i++;
				} // while loop
				if(greedypath.size()==2 && greedypath.get(0)==base && greedypath.get(1)==base)
				{
					greedypath.clear();
					greedypath.add(base);
				}
				if(greedypath.size()>=3)
				{
					ArrayList<Integer> tmp = new ArrayList<Integer>(greedypath);
					if(greedypath.size()>dmax)
					{
						printGreedyPath(greedypath);
					}
					jointgr.add(tmp);
					//printPaths(jointgr);
				}
			}
			if(totalcoin>besttotalcoin)
			{
				bestjointgr = jointgr;
				besttotalcoin = totalcoin;
			}
			//printPaths(jointgr);

		}
		return bestjointgr;


			}




	private static ArrayList<ArrayList<Integer>> greedyCoverMultRes2(int base, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int nRes) 
			{


		ArrayList<Integer> tsrt = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> bestjointgr = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> jointgr = new ArrayList<ArrayList<Integer>>();
		double besttotalcoin = -1;
		//int[] coin = new int[targetssorted.length];

		HashMap<Integer, Double> coin = new HashMap<Integer, Double>();

		for(int k=0; k<targetssorted.length; k++)
		{
			//coin[targetssorted[k][0]]=;
			coin.put(targetssorted[k][0], Math.floor(targetssorted[k][1]));
		}
		for(int k=1; k<targetssorted.length; k++)
		{
			tsrt.add(targetssorted[k][0]);
		}
		/*System.out.println("Tsrt : ");
		printGreedyPath(tsrt);*/

		for(int iter = 0; iter<100; iter++)
		{
			ArrayList<Integer> tmptsrt = new ArrayList<Integer>();
			if(iter==0)
			{
				tmptsrt=tsrt;
			}
			else
			{
				//shuffleArray(tmptsrt);
				tmptsrt = new ArrayList<Integer>(tsrt);
				Collections.shuffle(tmptsrt, new Random());
			}
			int i=0; 


			double totalcoin = 0;
			HashMap<Integer, Double> remaincoin = coin;
			ArrayList<Integer> greedypath = new ArrayList<Integer>();
			for(int res = 0; res<nRes; res++)
			{
				//System.out.println("res = "+res);
				//System.out.println("Adding base to greedy path : ");
				greedypath.clear();
				greedypath.add(base);
				greedypath.add(base);
				//printGreedyPath(greedypath);
				int totaldist = 0;
				totalcoin = totalcoin + remaincoin.get(base);
				remaincoin.remove(base);
				remaincoin.put(base, 0.0);
				while(i<tsrt.size())
				{
					int besttotaldisttemp = -1;
					int bestj = -1;
					for(int j=1; j<greedypath.size(); j++)
					{
						int s = map.get(greedypath.get(j-1));
						int d = map.get(greedypath.get(j));
						int totaldisttemp = totaldist - apspmat[s][d];

						s = map.get(greedypath.get(j-1));
						//d = map.get(targetssorted[i][0]);
						d = map.get(tmptsrt.get(i));
						totaldisttemp +=  apspmat[s][d];


						s = map.get(tmptsrt.get(i));
						d = map.get(greedypath.get(j));
						totaldisttemp +=  apspmat[s][d];

						//System.out.println("totaldisttemp = "+totaldisttemp);
						//System.out.println("bestj : "+bestj);

						if ((totaldisttemp<=dmax) && 
								((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) && totaldisttemp>0)
						{
							bestj = j;
							besttotaldisttemp = totaldisttemp;
							//System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
							//System.out.println("updating bestj : "+bestj);
						}



					}
					if(bestj>-1)
					{
						greedypath = insertTsrt(greedypath, 0, bestj-1, tmptsrt.get(i) , bestj, greedypath.size()-1);
						if(greedypath.size()>dmax)
						{
							printGreedyPath(greedypath);
						}
						printGreedyPath(greedypath);
						totaldist = besttotaldisttemp;

						totalcoin += remaincoin.get(tmptsrt.get(i));
						//remaincoin[tmptsrt.get(i)] = 0;
						remaincoin.remove(tmptsrt.get(i));
						remaincoin.put(tmptsrt.get(i), 0.0);
					}
					i++;
				} // while loop
				if(greedypath.size()==2 && greedypath.get(0)==base && greedypath.get(0)==base)
				{
					greedypath.clear();
					greedypath.add(base);
				}
				if(greedypath.size()>=3)
				{
					ArrayList<Integer> tmp = new ArrayList<Integer>(greedypath);
					if(greedypath.size()>dmax)
					{
						printGreedyPath(greedypath);
					}
					jointgr.add(tmp);
					//printPaths(jointgr);
				}
			}
			if(totalcoin>besttotalcoin)
			{
				bestjointgr = jointgr;
				besttotalcoin = totalcoin;
			}
		//	printPaths(jointgr);

		}
		return bestjointgr;


			}



	private static ArrayList<ArrayList<Integer>> greedyPathMultRes2(int base, ArrayList<TargetNode> targets,
			double dmax, int[][] targetssorted, int[][] apspmat,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int nRes) 
			{


		ArrayList<Integer> tsrt = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> bestjointgr = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> jointgr = new ArrayList<ArrayList<Integer>>();
		double besttotalcoin = -1;
		//int[] coin = new int[targetssorted.length];

		HashMap<Integer, Double> coin = new HashMap<Integer, Double>();

		for(int k=0; k<targetssorted.length; k++)
		{
			//coin[targetssorted[k][0]]=;
			coin.put(targetssorted[k][0], Math.floor(targetssorted[k][1]));
		}
		for(int k=1; k<targetssorted.length; k++)
		{
			tsrt.add(targetssorted[k][0]);
		}
		/*System.out.println("Tsrt : ");
		printGreedyPath(tsrt);*/

		for(int iter = 0; iter<100; iter++)
		{
			ArrayList<Integer> tmptsrt = new ArrayList<Integer>();
			if(iter==0)
			{
				tmptsrt=tsrt;
			}
			else
			{
				//shuffleArray(tmptsrt);
				tmptsrt = new ArrayList<Integer>(tsrt);
				Collections.shuffle(tmptsrt, new Random());
			}
			int i=0; 


			double totalcoin = 0;
			HashMap<Integer, Double> remaincoin = coin;
			ArrayList<Integer> greedypath = new ArrayList<Integer>();
			for(int res = 0; res<nRes; res++)
			{
				//System.out.println("res = "+res);
				//System.out.println("Adding base to greedy path : ");
				greedypath.clear();
				greedypath.add(base);
				greedypath.add(base);
				//printGreedyPath(greedypath);
				int totaldist = 0;
				totalcoin = totalcoin + remaincoin.get(base);
				remaincoin.remove(base);
				remaincoin.put(base, 0.0);
				while(i<tsrt.size())
				{
					int besttotaldisttemp = -1;
					int bestj = -1;
					for(int j=1; j<greedypath.size(); j++)
					{
						int s = map.get(greedypath.get(j-1));
						int d = map.get(greedypath.get(j));
						int totaldisttemp = totaldist - apspmat[s][d];

						s = map.get(greedypath.get(j-1));
						//d = map.get(targetssorted[i][0]);
						d = map.get(tmptsrt.get(i));
						totaldisttemp +=  apspmat[s][d];


						s = map.get(tmptsrt.get(i));
						d = map.get(greedypath.get(j));
						totaldisttemp +=  apspmat[s][d];

						/*System.out.println("totaldisttemp = "+totaldisttemp);
						System.out.println("bestj : "+bestj);*/

						if ((totaldisttemp<dmax) && 
								((besttotaldisttemp==-1) || (totaldisttemp<besttotaldisttemp)) && totaldisttemp>0)
						{
							bestj = j;
							besttotaldisttemp = totaldisttemp;
							/*System.out.println("updating besttotaldisttemp = "+besttotaldisttemp);
							System.out.println("updating bestj : "+bestj);*/
						}



					}
					if(bestj>-1)
					{
						greedypath = insertTsrt(greedypath, 0, bestj-1, tmptsrt.get(i) , bestj, greedypath.size()-1);
						if(greedypath.size()>dmax)
						{
							printGreedyPath(greedypath);
						}
						//printGreedyPath(greedypath);
						totaldist = besttotaldisttemp;

						totalcoin += remaincoin.get(tmptsrt.get(i));
						//remaincoin[tmptsrt.get(i)] = 0;
						remaincoin.remove(tmptsrt.get(i));
						remaincoin.put(tmptsrt.get(i), 0.0);
					}
					i++;
				} // while loop
				if(greedypath.size()==2 && greedypath.get(0)==base && greedypath.get(0)==base)
				{
					greedypath.clear();
					greedypath.add(base);
				}
				if(greedypath.size()>=3)
				{
					ArrayList<Integer> tmp = new ArrayList<Integer>(greedypath);
					if(greedypath.size()>dmax)
					{
						printGreedyPath(greedypath);
					}
					jointgr.add(tmp);
					//printPaths(jointgr);
				}
			}
			if(totalcoin>besttotalcoin)
			{
				bestjointgr = jointgr;
				besttotalcoin = totalcoin;
			}
			//printPaths(jointgr);

		}
		return bestjointgr;


			}

	private static void printGreedyPath(ArrayList<Integer> greedypath) {

		System.out.println();
		for(Integer n: greedypath)
		{
			System.out.print(n+" ");
		}
		System.out.println();


	}

	public static ArrayList<Integer> insertTsrt(ArrayList<Integer> greedypath, int i, int j,
			int k, int l, int m) {


		
		
		
		for(int p=0; p<greedypath.size(); p++)
		{
			if(k==greedypath.get(p))
				return greedypath;
		}
		

		ArrayList<Integer> tmp = new ArrayList<Integer>();

		for(int p = i; p<=j; p++)
		{
			tmp.add(greedypath.get(p));
		}

		tmp.add(k);


		for(int p = l; p<=m; p++)
		{
			tmp.add(greedypath.get(p));
		}
		return tmp;


	}
	
	

	private static void printtargets(ArrayList<TargetNode> targets2) {


		for(TargetNode n: targets2)
		{
			System.out.println("taregt "+ n.getTargetid()+", u: "+ n.getAnimaldensity());
		}

	}

	public static void basicSeqAbstractionAlgorithmWithExtreamPrunningTest(double[][] density, int ITER, int nrow, int ncol, int[] percentages, double dmax) throws Exception 
	{



		//int[] contractionsizes = {0,2,5,8,10};
		//int ITER = 20;
		//double[] result = new double[contractionsizes.length];
		int rindex=0;
		//double percentages[] = {20};
		double thresholds[] = {0};
		/*int nrow= 5;
		int ncol = 5;
		int dmax = 8;*/
		//int nUnaccesstargets = percentage;
		int nRes = 2;
		int nTargets = nrow*ncol;
		for(double percentage: percentages)
		{
			for(double threshold: thresholds)
			{

				double sumsol = 0;
				//long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				long targetsize=0;
				long sumthreshold=0;
				long sumiter =0;


				for(int iter=0; iter<ITER; iter++)
				{
					targets.clear();
					double [] result = new double[2];

					int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
					//	makeZeroSum(gamedata,nTargets);

					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);
					//System.out.println("Unnecessary targets "+ nUnaccesstargets);
					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);

					//chooseDummyNodes(nUnaccesstargets);

					//makeStarGraph(gamedata, nTargets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					preComputeShortestPaths(0, targets, domindatednodes );

					Date start = new Date();
					long l1 = start.getTime();

					ArrayList<Integer> greedypath = greedyFirstRoute(dmax,gamedata, targets);


					threshold = initializeThreshold(greedypath);
					//sumthreshold+= threshold;

					//System.out.println("\nInitial THreshold: "+ threshold);
					//	System.out.println("Greedy path: ");
					for(Integer x: greedypath)
					{
						System.out.print(x+"->");
					}



					//sgc.contractGraph(domindatednodes, targets);



					result =	SecurityGameContraction.basicSeqAbstractionExtreamPruning(targets, threshold, dmax, gamedata, sgc, nRes);
					sumsol += result[0];
					sumthreshold += result[1];
					sumiter += result[2];
					targetsize+= targets.size();
					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;
					solvingtime+=diff;
					revmaptime+=diff;



				}
				//double avgtime = sumtime/(ITER);
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				String x = df.format(revtime);
				writeInFile((int)targetsize/ITER, avgsol, contractiontime/ITER, solvingtime/ITER, Long.parseLong(x));
			}

		}
		//writeInFile(result,contractionsizes);


	}

	public static void spanningTreeTest(double[][] density, int ITER) {




		int nrow= 10;
		int ncol = 10;
		int dmax = 18;
		int nRes = 2;
		//double perc = .20;
		int nTargets = nrow*ncol;
		double thresholds[]={2,3,4,5};
		//int ITER=4;

		double percentages[] = {20};
		double[] result = new double[percentages.length];
		int rindex=0;
		double defexp=0;

		//int lstart=1,  lend=4, hstart=8, hend=10;
		//double[][] density=generateRandomDensity( percentages[0], ITER, lstart, lend,  hstart, hend, nTargets);


		for(double perc: percentages)
		{

			for(double threshold: thresholds)
			{
				double sumsol = 0;
				//long sumtime = 0;
				long contractiontime=0;
				long solvingtime=0;
				long revmaptime=0;
				int targetsize=0;


				for(int iter=0; iter<ITER; iter++)
				{

					targets.clear();

					int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

					//makeZeroSum(gamedata,nTargets);
					System.out.println("\n Iter "+ iter);
					System.out.println("Number of targets "+ nrow*ncol);
					System.out.println("dmax "+ dmax);

					System.out.println("nRes "+ nRes);
					SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
					//printEdges(targets);
					assignRandomDensityZeroSum(density, gamedata, targets, iter);

					ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
					selectDominatedTargets(targets, domindatednodes, threshold);
					//printNodesWithNeighborsAndPath(domindatednodes, targets);

					Date start = new Date();
					long l1 = start.getTime();


					instantContraction(domindatednodes, targets, dmax);




					Date stop = new Date();
					long l2 = stop.getTime();
					long diff = l2 - l1;
					contractiontime += diff;

					SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, targets);
					SecurityGameContraction.removeDominatedTargets(domindatednodes, targets);
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
					targetsize+= targets.size();




					/**
					 * make spanning tree
					 * 
					 */

					ArrayList<TargetNode> spanningtree = makeSpanningTree(targets, targets.size(), domindatednodes);
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, spanningtree);
					domindatednodes = SecurityGameContraction.removeDominatedTargetsUsingID(domindatednodes, spanningtree);
					SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, spanningtree);

					System.out.println("Resulting dspanning tree size "+ spanningtree.size());




					int[][] p = new int[spanningtree.size()][];
					try 
					{
						start = new Date();
						l1 = start.getTime();
						ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, spanningtree);
						ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
						/**
						 * map has present id
						 * mapback gives the original ids
						 */
						HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
						HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
						makePathSeq(pathseq, goals, goals.size(), spanningtree.size(), map, mapback, spanningtree );
						System.out.println("Total path with duplicates "+pathseq.size());
						printPaths(pathseq);
						pathseq = removeDuplicatePathSimple(pathseq);
						printPaths(pathseq);
						System.out.println("Total path without duplicates "+pathseq.size());

						//int k = 2;
						Integer[] input = new Integer[pathseq.size()];
						int[] branch = new int[nRes];//{0,0};//new char[k];


						for(int i=0; i<input.length; i++)
						{
							input[i] = i;
						}
						HashSet jSet=new HashSet();
						if(pathseq.size()==0)
						{
							//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
							//choose the worst payoff for defender

							Double mAxpayoff = Double.MIN_VALUE;
							Double defpayoff = 0.0;
							for(int i=0; i<domindatednodes.size(); i++)
							{
								spanningtree.add(domindatednodes.get(i));
							}
							for(TargetNode x: spanningtree)
							{
								if(x.attackerreward>mAxpayoff)
								{
									mAxpayoff= x.attackerreward;
									defpayoff = x.defenderpenalty;
								}
							}

							sumsol += defpayoff;
							System.out.println("Defender expected payoff "+ defpayoff);
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;


						}
						else
						{
							//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
							if(pathseq.size()<nRes)
							{

								branch = new int[pathseq.size()];
								jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
							}
							else
							{
								jSet=combine(input, nRes, 0, branch, 0, jSet);
							}

							List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
							/**
							 * columns will be combination of paths for each resources. 
							 */
							/**
							 * pmat, where columns will be combination of paths. 
							 * rows are targets. 
							 * each entry will say whether the target is in the joint schedule
							 */
							//jSet.

							//printJointSchedule(jset);

							//printNodesAsNeighbors(dominatednodes);

							p = makePmat(pathseq, jset, mapback, spanningtree);
							//printPathMat(p);

							/**
							 * remove duplicates from p
							 */
							//removeDuplicatesFromP(p);
							//System.out.println();
							//printPathMat(p);

							//System.out.println("Number of targets after contraction "+ targets.size());
							//System.out.println("mip in... ");
							System.out.println("Iter "+iter+", mip in... ");
							double[] probdistribution = MIPSolver4.solve(p, gamedata, spanningtree, nRes);
							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							solvingtime += diff;
							if(probdistribution.equals(null))
							{
								throw new Exception("Prob null...");
							}
							System.out.println("Iter "+iter+", mip out... ");
							//printPathWithPositiveCoverage(p, coverage, jset, pathseq, map);


							start = new Date();
							l1 = start.getTime();
							int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, spanningtree);
							//removeDuplicatesFromP(origpmat);
							//printPathMat(origpmat);
							//System.out.println("\n after mapping back");
							//printPmatWithPositiveCoverage(origpmat, coverage, jset, pathseq,map);

							stop = new Date();
							l2 = stop.getTime();
							diff = l2 - l1;
							if(threshold>0)
							{
								revmaptime += diff;
							}


							/*for(int i=0; i<coverage.length; i++)
	{
		for(int j=0; j<coverage[i].length; j++)
		{
			if(coverage[i][j]>0)
			{
				System.out.println("selected path : " + j);

				for(int k=0; k<origpmat.length; k++)
				{
					System.out.print(origpmat[k][j] + " ");
				}
				System.out.println();

			}
		}
	}*/			

							int maxtargetforattacker = findAttackTarget(origpmat, probdistribution, gamedata);

							double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);
							defexp = defexpectedpayoff;
							System.out.println("Attacked target is "+ maxtargetforattacker);
							System.out.println("Defender expected payoff "+ defexpectedpayoff);


							Logger.logit("instant contraction :  \n");
							for(TargetNode n: domindatednodes)
							{
								Logger.logit(" dominated node "+ n.getTargetid()+"\n");
							}
							Logger.logit("expected payoff "+ defexpectedpayoff+ "\n" + "threshold "+ threshold+"\n");



							sumsol += defexp;
							/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defexpectedpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/

						}

					} 
					catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

					//	writeInFile(iter,(int) threshold, targetsize/(iter+1), defexp, contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));

					//writeInFile(iter,(int) threshold, targetsize/(iter+1), sumsol/(iter+1), contractiontime/(iter+1), solvingtime/(iter+1), revmaptime/(iter+1));
					//writeInFile(iter,(int) threshold, targetsize/(ITER), sumsol/(ITER), contractiontime/(ITER), solvingtime/(ITER), revmaptime/(ITER));


				}
				double avgsol = sumsol/ITER;
				//result[rindex++] = avgsol;
				DecimalFormat df = new DecimalFormat("#.#######");
				double revtime = revmaptime/ITER;
				//String x = df.format(revtime);
				writeInFile(targetsize/ITER, sumsol/ITER, contractiontime/ITER, solvingtime/ITER, (long)revtime);
			}
		}



	}

	private static ArrayList<TargetNode> removeDominatedTargetsUsingID(
			ArrayList<TargetNode> domindatednodes,
			ArrayList<TargetNode> spanningtree) {

		ArrayList<TargetNode> dom = new ArrayList<TargetNode>();


		for(int i=0; i<domindatednodes.size(); i++)
		{
			for(int j=0;j<spanningtree.size(); j++)
			{
				if(domindatednodes.get(i).getTargetid()==spanningtree.get(j).getTargetid())
				{
					dom.add(spanningtree.get(j));
					spanningtree.remove(spanningtree.get(j));
					break;
				}

			}
		}

		return dom;

	}

	public static void contractionWithSingleOracleGreedyTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = contractionWithGreedySingleOracle(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER);

		}

		System.out.println("\nDef exp utility : "+ sumsol/ITER);

		writeInFile("9" ,(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}


	public static void contractionWithSingleOracleTest1(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = contractionWithSingleOracle(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef exp utility : "+ sumsol/ITER);

		writeInFile("10", (int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}


	
	



	public static void contractionWithDoubelOracleGreedyPath3TestLP(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = contractionWithGreedyPath3DoubleOracleLP(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("14",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}








	public static void contractionWithDoubelOracleTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = contractionWithGreedyDoubleOracle(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("9", (int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}


	
	


	public static void doubleOracleGCSingleGPMultiLPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCSingleGPMultiLPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("1",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}


	
	


	public static void doubleOracleGCSingleGPMultiLPSamplePathTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCSingleGPMultiLPSamplePath(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("5",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}


	
	
	


	public static void doubleOracleGCSingleGP3LPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCSingleGP3LPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("3",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	
	


	public static void doubleOracleGCSingleGP3LPSamplePathTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCSingleGP3LPSamplePath(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("7",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	
	
	

	public static void doubleOracleGCMultiGPMultiLPSamplePathTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGPMultiLPSamplePath(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("6",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}



	
	

	public static void doubleOracleGCMultiGPMultiLPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGPMultiLPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("2",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	
	

	public static void doubleOracleGCMultiGP3LPGCMultiRealWorldDataTest(
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{
			
			double u [][] = new double[nrow][ncol];
			double e [][] = new double[nrow][ncol];
			double utility [][] = new double[560][560];
			double elevation [][] = new double[560][560];
			
			
			
			ReadData.readData(560, 560, utility, elevation);
			ReadData.getChunk( utility, elevation, 0, 0, nrow, ncol, u, e);
			int[][] gamedata = constructGameData(u);//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			
			
			
			
			
			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LPGCMultiRealWorldData(nTargets, nRes, dmax, iter, nrow, ncol, u, e, gamedata);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile(String.valueOf(nTargets),(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	
	

	public static int[][] constructGameData(double[][] u) {
		
		
		int[][] gd = new int[u.length*u[0].length][u[0].length];
		
		int index =0;
		
		for(int i=0; i<u.length; i++)
		{
			for(int j=0; j<u[0].length; j++)
			{
				gd[index][0] = 0;
				gd[index][1] = (int)Math.floor(-u[i][j]);
				gd[index][2] = -gd[index][1];
				gd[index][3] = 0;
				index++;
			}
		}
		
		
		
		return gd;
	}
	
	
public static int[][] constructGameData(ArrayList<TargetNode> u) {
		
		
		int[][] gd = new int[u.size()][u.size()];
		
		int index =0;
		
		for(int i=0; i<u.size(); i++)
		{
			//for(int j=0; j<u.size(); j++)
			{
				gd[index][0] = 0;
				gd[index][1] = (int)Math.floor(-u.get(i).attackerreward);
				gd[index][2] = -gd[index][1];
				gd[index][3] = 0;
				index++;
			}
		}
		
		
		
		return gd;
	}
	
	
	
	public static void GroupingTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;
		int srcid=1;
		int destid=3;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = GroupingWithDO(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, srcid, destid);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("4",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER,totaltime/ITER, nTargets);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	
	
	public static void doubleOracleGCMultiGP3LPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("4",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER,totaltime/ITER, nTargets);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	public static void DOTest(double[][] density,
			int ITER, int nrow, int ncol,
			double dmax, int nRes, HashMap<Integer,ArrayList<TargetNode>> alltargets, HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{

			
			
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			//printNodesWithNeighborsAndPath(targetmaps);

			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
			
			gamedata = constructGameData(targets);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = DO(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol, targets, targetmaps);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			
			//SecurityGameContraction.writeRes("DO", iter, (int)sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);

			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("DO",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER,totaltime/ITER, nTargets);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	
	public static void doubleOracleTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracle(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];
			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("DO",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER,totaltime/ITER, nTargets);
		//writeInFile("4",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10);
		//(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10

	}
	
	
	
	public static void doubleOracleGCMultiExactLPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;
		long sumslavetime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiExactLPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime += res[5];

			//writeInFile(Integer.toString(iter),  (int)res[3], res[0], sumcontractiontime/iter, sumsolvtime/iter, sumslavetime/10, totaltime/10);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("Exact",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER, totaltime/ITER, nTargets );
		//writeInFile("Exact",(int)sumfinaltargetsize/10, sumsol/10, sumcontractiontime/10, sumsolvtime/10, sumslavetime/10, totaltime/10 );



	}
	
	
	
	
	public static void doubleOracleGCMultiGP3LP_OPTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime=0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LP_OP(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime+= res[5];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("OP",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER, totaltime/ITER, nTargets);



	}
	
	
	public static void doubleOracleGCMultiGP3LP_lexicoOPTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime=0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LP_lexicoOP(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime+= res[5];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("lexicoOP",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER,sumslavetime/ITER,  totaltime/ITER, nTargets);



	}
	
	
	public static void doubleOracleGCMultiGP3LP_modifiedOPTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime=0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LP_modifiedOP(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime+= res[5];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("modifiedOP",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER,sumslavetime/ITER,  totaltime/ITER, nTargets);



	}
	
	
	public static void doubleOracleGCMultiGP3LPGCMultiTest(HashMap<Integer, ArrayList<TargetNode>> alltargets,
			HashMap<Integer, HashMap<Integer, TargetNode>> alltargetmaps,
			int ITER, int nTargets,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime=0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{
			ArrayList<TargetNode> targets = alltargets.get(iter);
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); 
			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LPGCMulti(gamedata, nTargets, nRes, dmax, iter, targetmaps, targets);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime+= res[5];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("DO",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER,sumslavetime/ITER,  totaltime/ITER, nTargets);



	}
	
	
	
	public static void doubleOracleGCMultiGP3LP_TOPTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long sumslavetime=0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LP_TOP(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];
			sumslavetime+= res[5];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("TOP",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, sumslavetime/ITER,  totaltime/ITER, nTargets);



	}
	
	
	
	
	
	

	public static void doubleOracleGCMultiGP3LPOPTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			
			//long t1= System.nanoTime();

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("4",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	

	public static void doubleOracleGCMultiGP3APSPLPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3APSPLPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("4-GP3-APSP",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	




	public static void doubleOracleAPSPGCMultiGP3LPGCMultiTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleAPSPGCMultiGP3LPGCMulti(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("APSP",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	public static void doubleOracleGCMultiGP3LPSamplePathTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = doubleOracleGCMultiGP3LPSamplePath(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("8",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}






	public static void contractionWithDoubelOracleGreedyCover3LPTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = contractionWithGreedyDoubleOracleWGreedyCover3LP(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("13",(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}




	public static void noContractionWithColumnGenerationTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = noContractionWithColumnGeneration(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("15", (int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}
	
	
	public static void baselineTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = noContractionWithColumnGeneration(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("15", (int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}







	public static void contractionWithDoubelOracleGreedyCoverTest(double[][] density,
			int ITER, int nrow, int ncol, int[] percentages, int[] thresholds,
			double dmax, int nRes) throws Exception {
		// TODO Auto-generated method stub

		int nTargets = nrow*ncol;
		double sumsol=0;
		long sumcontractiontime = 0;
		long sumsolvtime =0;
		long sumfinaltargetsize = 0;
		long sumthreshold = 0;
		long totaltime = 0;


		for(int iter=0; iter<ITER; iter++)
		{


			int[][] gamedata = new int[nTargets][4];//SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);

			Date start = new Date();
			long l1 = start.getTime();
			double[] res = contractionWithGreedyDoubleOracleWGreedyCover(gamedata, nTargets, nRes, density, dmax, iter, nrow, ncol);
			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			totaltime += diff;


			System.out.println("\nDef exp utility : "+ res[0]);
			sumsol += res[0];
			sumcontractiontime += res[1];
			sumsolvtime += res[2];
			sumfinaltargetsize += res[3];
			sumthreshold += res[4];

			//writeInFile(sumthreshold/ITER, sumfinaltargetsize/ITER, res[0], sumcontractiontime/ITER, sumsolvtime/ITER, 0);

		}

		System.out.println("\nDef avg exp utility : "+ sumsol/ITER);

		writeInFile("10" ,(int)sumfinaltargetsize/ITER, sumsol/ITER, sumcontractiontime/ITER, sumsolvtime/ITER, totaltime/ITER);



	}



	private static void writeInFile(String expno, int finalsize, double avgsol, long contracttime,
			long solvingtime, long totaltime) 
	{


		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("grp-result.csv"),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(expno+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," + totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	
	


	public static void writeRes(String expno, int itr, int finalsize, double avgsol, long contracttime,
			long solvingtime, long totaltime) 
	{


		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("itr-result.csv"),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(expno+","+itr+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," + totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	
	
	
	public static void writeInFile(String expno, int finalsize, double avgsol, long contracttime,
			long solvingtime, long slavetime, long totaltime, int nTargets ) 
	{


		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("grp-result.csv"),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(expno+","+nTargets+","+finalsize+ ","+ avgsol+ ","+contracttime+"," + solvingtime+"," +slavetime+","+ totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}





	}
	
	

	private static void writeInFile(String filename, int iter, int outerstage, double  outertime) 
	{


		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File(filename),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(iter+","+outerstage+ ","+ outertime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	
	
	private static void writeInFile(String filename, int iter, int outerstage, int innerstage, double  innertime) 
	{


		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File(filename),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(iter+","+outerstage+","+innerstage+ ","+ innertime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}




	private static double[] contractionWithSingleOracle(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = new ArrayList<Integer>();
		currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);



		int currentPlace = 1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;


		while(true)
		{

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);



			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);



			p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);


				start = new Date();
				l1 = start.getTime();

				probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);



				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvingtime += diff;


				attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);

				attackedtarget = mapback.get(attackedtarget);

				System.out.println("attack target "+ attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("attacker u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}

				if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
				{
					System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

					break;
				}
				else 
				{
					int prevcur = currentPlace;
					currentPlace += 1;
					if(currentPlace>targetssorted.length)
					{
						currentPlace = targetssorted.length;
					}
					//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

					for(int k= prevcur+1; k<=currentPlace; k++ )
					{
						currenttargets.add(targetssorted[k][0]);
					}
				}



				System.out.println();


			}


		}

		System.out.println("Final target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}

/*
	private static double[] doubleOracleGCMultiGPMultiLPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		*//**
		 * 1. sort the targets
		 *//*
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0); //  new ArrayList<Integer>();
		
		
		currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;
		double attackerv;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat;






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			
			pathseq = buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			
			//pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			*//**
			 * keep only nRes*3 paths from the end
			 *//*

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);



			while(true)
			{


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					*//**
					 * columns will be combination of paths for each resources. 
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 *//*
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();
					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,attackerstrategy );



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					*//**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 *//*
					System.out.println("attacked target after rev map "+ attackedtarget);

					//ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);

					
					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					
					
					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					int newsize = pathseq.size();
					System.out.println("haa ");


					if(oldsize==newsize)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);

				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<3 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 3-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	

*/




	private static double[] contractionWithGreedyPath3DoubleOracleLP(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;
		double attackerv;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat;






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			
			//pathseq = buildGreedyPathMultRes(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			
			pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);



			while(true)
			{


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();
					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes,attackerstrategy );



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();

					ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);

					makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					int newsize = pathseq.size();
					System.out.println("haa ");


					if(oldsize==newsize)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);

				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			/*if(addcount<3 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 3-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}*/






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	









	private static double[] contractionWithGreedyDoubleOracle(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;
		double attackerv;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat;






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			ArrayList<Integer> pathnodes = new ArrayList<Integer>();
			pathseq = generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			
			//double sd = findShortestPathThrougGraph(tmpgraph.get(0), tmpgraph.get(2), tmpgraph, pathnodes );
			//findShortestPath(tmpgraph.get(0), tmpgraph.get(1), tmpgraph, domindatednodes, pathnodes);
			
			
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);



			while(true)
			{


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();
					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					//attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					//attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();

					ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);

					makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					int newsize = pathseq.size();
					System.out.println("haa ");


					if(oldsize==newsize)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
						break;

				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			/*if(addcount<3 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 3-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}*/






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	


	
	

	private static double[] doubleOracleGCSingleGP3LPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}*/


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			


			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}
	
	
	
	
	
	

	private static double[] doubleOracleGCSingleGP3LPSamplePath(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}*/


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();
					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}
	
	
	
	private static double[] doubleOracleGCSingleGPMultiLPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		
		int stagecounter = 0;
		while(true)
		{
			
			long outertime = 0;
			long innertime = 0;
			
			Date outerstart = new Date();
			long outerl1 = outerstart.getTime();

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			/*System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}*/


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			//pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}
*/


			int innerstagecounter=0;
			while(true)
			{
				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/
				
				Date innerstart = new Date();
				long innerl1 = innerstart.getTime();
				
				
				
				
				innertime=0;

				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						Date innerend = new Date();
						long innerl2 = innerend.getTime();
						
						innertime = innerl2- innerl1;
						
						writeInFile("inner.csv", iter, stagecounter, innerstagecounter, innertime);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						Date innerend = new Date();
						long innerl2 = innerend.getTime();
						
						innertime = innerl2- innerl1;
						
						writeInFile("inner.csv", iter, stagecounter, innerstagecounter, innertime);
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (innerstagecounter>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						Date innerend = new Date();
						long innerl2 = innerend.getTime();
						
						innertime = innerl2- innerl1;
						
						writeInFile("inner.csv", iter, stagecounter, innerstagecounter, innertime);
						break;
					}

					//printPaths(pathseq);


				} // end if else
				//System.out.println("iter"+ innerstagecounter++);
				
				Date innerend = new Date();
				long innerl2 = innerend.getTime();
				
				innertime = innerl2- innerl1;
				
				writeInFile("inner.csv", iter, stagecounter, innerstagecounter, innertime);
				
				
				
				
				innerstagecounter++;
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			
			Date outerend = new Date();
			long outerl2 = outerend.getTime();
			
			diff = outerl2 - outerl1;

			outertime = diff; 
			
			
			writeInFile("outer.csv", iter, stagecounter, outertime);
			
			
			
			
			stagecounter++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	


	
	
	
	private static double[] doubleOracleGCSingleGPMultiLPSamplePath(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		
		int stagecounter = 0;
		while(true)
		{
			
			long outertime = 0;
			long innertime = 0;
			
			Date outerstart = new Date();
			long outerl1 = outerstart.getTime();

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}
*/

			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			//pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			

			int innerstagecounter=0;
			while(true)
			{
				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/
				
				Date innerstart = new Date();
				long innerl1 = innerstart.getTime();
				
				
				
				
				innertime=0;

				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();
					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (innerstagecounter>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				//System.out.println("iter"+ innerstagecounter++);
				
				Date innerend = new Date();
				long innerl2 = innerend.getTime();
				
				innertime = innerl2- innerl1;
				
				writeInFile("inner.csv", iter, stagecounter, innerstagecounter, innertime);
				
				
				
				
				innerstagecounter++;
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}





			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			
			Date outerend = new Date();
			long outerl2 = outerend.getTime();
			
			diff = outerl2 - outerl1;

			outertime = diff; 
			
			
			writeInFile("outer.csv", iter, stagecounter, outertime);
			
			
			
			
			stagecounter++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}
	
	
	
	


	private static double[] doubleOracleGCMultiGPMultiLPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}*/


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			//pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	




	
	


	private static double[] doubleOracleGCMultiGPMultiLPSamplePath(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}
*/

			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();
*/
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			//pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();
					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}
	
	
	


	
	private static double[] doubleOracleGCMultiGP3LPGCMultiRealWorldData(
			int nTargets, int nRes, double
			dmax, int iter, int nrow, int ncol, double[][] u, double[][] e, int[][] gamedata) throws Exception {


		
		
		//int base = 1842;

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, u, e, gamedata);
		//assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();
		
		
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, targets);
		

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					
					 if( (currentPlace < targetssorted.length-1) && (attackeru<targetssorted[currentPlace+1][1]))
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>=targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);
				
				if(prevcur<targetssorted.length)
				{

					for(int k= prevcur+1; k<=currentPlace; k++ )
					{

						//System.out.println("adding target  "+ targetssorted[k][0]);
						currenttargets.add(targetssorted[k][0]);
					}
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}
	
	
	
	
	private static double[] doubleOracleGCMultiGP3LP_OP(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		long slavetime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//test
					
					start = new Date();
					l1 = start.getTime();
					
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					slavetime += diff;
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);


					
					/** test
					 * 
					 */
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					
					

					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					/*System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	

	private static double[] doubleOracleGCMultiGP3LPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double
			dmax, int iter, HashMap<Integer,TargetNode> targetmaps,
			ArrayList<TargetNode> targets) throws Exception {


		

		//targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nTargets, gamedata);
		assignRandomDensityZeroSum(targetmaps, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		long slavetime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//test
					
					start = new Date();
					l1 = start.getTime();
					
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.modifiedTOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					slavetime += diff;
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);


					
					/** test
					 * 
					 */
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				  //  newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					/*if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					*/
					

					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	private static double[] doubleOracleGCMultiGP3LP_modifiedOP(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		long slavetime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//test
					
					start = new Date();
					l1 = start.getTime();
					
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.modifiedTOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					slavetime += diff;
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);


					
					/** test
					 * 
					 */
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				  //  newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					/*if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					*/
					

					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	
	
	
	
	private static double[] doubleOracleGCMultiGP3LP_lexicoOP(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		long slavetime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//test
					
					start = new Date();
					l1 = start.getTime();
					
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.lexicographicOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					slavetime += diff;
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);


					
					/** test
					 * 
					 */
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				  //  newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					/*if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					*/
					

					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	
	
	private static double[] doubleOracleGCMultiGP3LP_TOP(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		long slavetime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//test
					
					start = new Date();
					l1 = start.getTime();
					
					ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalTOP(nRes, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					slavetime += diff;
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);


					
					/** test
					 * 
					 */
					
					/*System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					*/
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					
					

					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	

	
	private static double[] doubleOracle(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		//targets.clear();
		//SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		
		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		SecurityGameContraction.buildGraph(nrow, ncol, gamedata, targets);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;

		


		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/**test
					 * 
					 */
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					
					

					//test
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					/*pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/

					printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	private static double[] doubleOracleGCMultiGP3LPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		
		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;

		


		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/**test
					 * 
					 */
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					
					

					//test
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					/*pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/

					printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	
	private static double[] DO(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps) throws Exception {


		

		/*targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		
		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		
		assignRandomDensityZeroSum(density, gamedata, targets, iter);
*/

		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;

		


		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					
					attackedtarget = mapback.get(attackedtarget);
					int attackedtargetrestrgraph = attackedtarget;
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);
					
					
					ArrayList<Integer> attackpath = pathForAT(dmax, tmpgraph, currenttargets, nRes, attackedtargetrestrgraph);
					
					newpathseq.add(attackpath);
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/**test
					 * 
					 */
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					
					
					

					//test
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					/*pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	
	
	
	

	
	
	private static double[] GroupingWithDO(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol, int srcid, int destid) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		
		//targets = dummygraph1(gamedata);
		targets = dummygraph2(gamedata);
		targets= dummygraph1(gamedata);
		
		
		//assignRandomDensityZeroSum(density, gamedata, targets, iter);


		printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverWithSrcDest(targets, dmax, nTargets, srcid, destid, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			ArrayList<TargetNode> goals = generatePathsGreedy2SrcDest(dmax, gamedata, tmpgraph, currenttargets, nRes, srcid, destid);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			//pathseq =  generatePathsGreedy3WithSrcDest(dmax, gamedata, tmpgraph, currenttargets, nRes, srcid, destid);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			/*int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}*/
			makePathSeqSrcDest(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff && x.getTargetid()!=srcid)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
							attackedtarget=x.getTargetid();
						}
					}
					canaddpath=false;
					break;
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultResSrcDest(tmpgraph, dmax, tmpgraph.size(), srcid, destid, nRes, attackerstrategy);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/**test
					 * 
					 */
					/*System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					*/
					
					

					//test
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ############### \n final paths: ");
						printPaths(pathseq);
						break;
					}

					


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	

	
	

	
	
	public static double[] subGraphDOSolver(ArrayList<TargetNode> subgraph, int[][] gamedata,
			int nTargets, int nRes, double dmax,  int srcid, int destid) throws Exception 
	{


		if(dmax==0)
		{
			double max = Double.MIN_VALUE;
			for(TargetNode t: subgraph)
			{
				if(t.getTargetid()!=srcid && max<t.attackerreward)
				{
					max = t.attackerreward;
				}
			}
			return (new double[]{-max, max});
		}

		//targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction();
		
		//targets = dummygraph1(gamedata);
		//targets = dummygraph2(gamedata);
		//targets= dummygraph1(gamedata);
		
		targets = subgraph;
		
		
		//assignRandomDensityZeroSum(density, gamedata, targets, iter);


		printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverWithSrcDest(targets, dmax, nTargets, srcid, destid, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/
		
		if(currenttargets.size()==1)
		{
			double max = Double.MIN_VALUE;
			for(TargetNode t: subgraph)
			{
				if(t.getTargetid()!=currenttargets.get(0) && max<t.attackerreward)
				{
					max = t.attackerreward;
				}
			}
			return (new double[]{-max, max});
		}



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;
		
		
		HashMap<Integer, Integer> originalmap = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> originalmapback = new HashMap<Integer, Integer>();
		
		makeOriginalMapping(originalmap, originalmapback, targets);



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;
		
		int masteritr=0;

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			HashMap<Integer, TargetNode> tmpgraphmap = new HashMap<Integer, TargetNode>();
			
			for(TargetNode n: tmpgraph)
			{
				tmpgraphmap.put(n.getTargetid(), n);
			}
			
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodesWMap(targetssorted, currentPlace+1, tmpgraph,tmpgraphmap);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			ArrayList<TargetNode> goals = generatePathsGreedy2SrcDest(dmax, gamedata, tmpgraph, currenttargets, nRes, srcid, destid);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			//pathseq =  generatePathsGreedy3WithSrcDest(dmax, gamedata, tmpgraph, currenttargets, nRes, srcid, destid);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			/*int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}*/
			makePathSeqSrcDest(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				System.out.println("Entered inner loop...slave");
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff && x.getTargetid()!=srcid)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
							attackedtarget=x.getTargetid();
							attackeru = x.attackerreward;
						}
					}
					canaddpath=false;
					break;
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMat(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					attackedtarget = mapback.get(attackedtarget);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				//	attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					attackerv = expectedAttackerPayoff(attackedtarget, origpmat, probdistribution, gamedata, originalmap);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(currentPlace<targetssorted.length-1 && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);
					
					
					start = new Date();
					l1 = start.getTime();


					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultResSrcDest(tmpgraph, dmax, tmpgraph.size(), srcid, destid, nRes, attackerstrategy);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/**test
					 * 
					 */
					/*System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");
					*/
					
					

					//test
					//ArrayList<ArrayList<Integer>> newpathseq = MIPSolver4.originalOP(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ############### \n final paths: ");
						printPaths(pathseq);
						break;
					}

					


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}



			

			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff,  attackeru};
		return res;
	}
	
	



	
	private static ArrayList<TargetNode> dummygraph2(int[][] gamedata) {
		
		
		ArrayList<TargetNode> graph = new ArrayList<TargetNode>();
		
		TargetNode d = new TargetNode(0, 5);
		TargetNode e = new TargetNode(1, 4);
		TargetNode f = new TargetNode(2, 5);
		TargetNode g = new TargetNode(3, 4);
		
		
		//d -> e
		ArrayList<TargetNode> path = new ArrayList<TargetNode>();
		
		d.addNeighbor(e);
		d.addDistance(e, 1.0);
		d.setPath(e, path);
		
		e.addNeighbor(d);
		e.addDistance(d, 1.0);
		e.setPath(d, path);
		
		//d->g
		
		d.addNeighbor(g);
		d.addDistance(g, 1.0);
		d.setPath(g, path);
		
		g.addNeighbor(d);
		g.addDistance(d, 1.0);
		g.setPath(d, path);
		
		
		// d->f
		d.addNeighbor(f);
		d.addDistance(f, 2.0);
		d.setPath(f, path);
		
		f.addNeighbor(d);
		f.addDistance(d, 2.0);
		f.setPath(d, path);
		
		
		// g->f
				g.addNeighbor(f);
				g.addDistance(f, 1.0);
				g.setPath(f, path);
				
				f.addNeighbor(g);
				f.addDistance(g, 1.0);
				f.setPath(g, path);
		
		
				// e->f
				e.addNeighbor(f);
				e.addDistance(f, 1.0);
				e.setPath(f, path);
				
				f.addNeighbor(e);
				f.addDistance(e, 1.0);
				f.setPath(e, path);
		
		gamedata[0][0] = 0;//(int)density[iter][i] ;
		gamedata[0][1] = -5;
		gamedata[0][2] = (int)5;  // uncovered
		gamedata[0][3] = 0; //covered
		
		
		
		
		gamedata[1][0] = 0;//(int)density[iter][i] ;
		gamedata[1][1] = -4;
		gamedata[1][2] = (int)4;  // uncovered
		gamedata[1][3] = 0; //covered
		
		
		gamedata[2][0] = 0;//(int)density[iter][i] ;
		gamedata[2][1] = -5;
		gamedata[2][2] = (int)5;  // uncovered
		gamedata[2][3] = 0; //covered
		
		gamedata[3][0] = 0;//(int)density[iter][i] ;
		gamedata[3][1] = -4;
		gamedata[3][2] = (int)4;  // uncovered
		gamedata[3][3] = 0; //covered
		
		
		
		d.defenderreward=0;
		d.defenderpenalty= -5;
		d.attackerreward = 5;
		d.attackerpenalty= 0;
		
		e.defenderreward=0;
		e.defenderpenalty= -4;
		e.attackerreward = 4;
		e.attackerpenalty= 0;
		
		f.defenderreward=0;
		f.defenderpenalty= -5;
		f.attackerreward = 5;
		f.attackerpenalty= 0;
		
		g.defenderreward=0;
		g.defenderpenalty= -4;
		g.attackerreward = 4;
		g.attackerpenalty= 0;
		
		
		
		graph.add(d);
		graph.add(e);
		graph.add(f);
		graph.add(g);
		
		
		return graph;
	}
	

	
	private static ArrayList<TargetNode> dummygraph3(int[][] gamedata) {
		
		
		ArrayList<TargetNode> graph = new ArrayList<TargetNode>();
		
		TargetNode a = new TargetNode(0, 10); //base
		TargetNode b = new TargetNode(1, 10); // sp1
		TargetNode c = new TargetNode(2, 5); // sp2
		
		
		
		//a -> b
		ArrayList<TargetNode> path = new ArrayList<TargetNode>();
		
		a.addNeighbor(b);
		a.addDistance(b, 1.0);
		a.setPath(b, path);
		
		b.addNeighbor(a);
		b.addDistance(a, 1.0);
		b.setPath(a, path);
		
		//a->c
		
		a.addNeighbor(c);
		a.addDistance(c, 2.0);
		a.setPath(c, path);
		
		c.addNeighbor(a);
		c.addDistance(a, 2.0);
		c.setPath(a, path);
		
		
		// c->b
		c.addNeighbor(b);
		c.addDistance(b, 1.0);
		c.setPath(b, path);
		
		b.addNeighbor(c);
		b.addDistance(c, 1.0);
		b.setPath(c, path);
		
		
		
		gamedata[0][0] = 0;//(int)density[iter][i] ;
		gamedata[0][1] = -10;
		gamedata[0][2] = (int)10;  // uncovered
		gamedata[0][3] = 0; //covered
		
		
		a.defenderreward=0;
		a.defenderpenalty= -10;
		a.attackerreward = 10;
		a.attackerpenalty= 0;
		
		b.defenderreward=0;
		b.defenderpenalty= -9;
		b.attackerreward = 9;
		b.attackerpenalty= 0;
		
		c.defenderreward=0;
		c.defenderpenalty= -9;
		c.attackerreward = 9;
		c.attackerpenalty= 0;
		
		
		
		
		gamedata[1][0] = 0;//(int)density[iter][i] ;
		gamedata[1][1] = -9;
		gamedata[1][2] = (int)9;  // uncovered
		gamedata[1][3] = 0; //covered
		
		
		gamedata[2][0] = 0;//(int)density[iter][i] ;
		gamedata[2][1] = -9;
		gamedata[2][2] = (int)9;  // uncovered
		gamedata[2][3] = 0; //covered
		
		
		
		graph.add(a);
		graph.add(b);
		graph.add(c);
		
		
		return graph;
	}
	
	
	
	private static ArrayList<TargetNode> dummygraph1(int[][] gamedata) {
		
		
		ArrayList<TargetNode> graph = new ArrayList<TargetNode>();
		
		TargetNode a = new TargetNode(0, 10);
		TargetNode b = new TargetNode(1, 9);
		TargetNode c = new TargetNode(2, 9);
		
		
		
		//a -> b
		ArrayList<TargetNode> path = new ArrayList<TargetNode>();
		
		a.addNeighbor(b);
		a.addDistance(b, 1.0);
		a.setPath(b, path);
		
		b.addNeighbor(a);
		b.addDistance(a, 1.0);
		b.setPath(a, path);
		
		//a->c
		
		a.addNeighbor(c);
		a.addDistance(c, 2.0);
		a.setPath(c, path);
		
		c.addNeighbor(a);
		c.addDistance(a, 2.0);
		c.setPath(a, path);
		
		
		// c->b
		c.addNeighbor(b);
		c.addDistance(b, 1.0);
		c.setPath(b, path);
		
		b.addNeighbor(c);
		b.addDistance(c, 1.0);
		b.setPath(c, path);
		
		
		
		gamedata[0][0] = 0;//(int)density[iter][i] ;
		gamedata[0][1] = -10;
		gamedata[0][2] = (int)10;  // uncovered
		gamedata[0][3] = 0; //covered
		
		
		a.defenderreward=0;
		a.defenderpenalty= -10;
		a.attackerreward = 10;
		a.attackerpenalty= 0;
		
		b.defenderreward=0;
		b.defenderpenalty= -9;
		b.attackerreward = 9;
		b.attackerpenalty= 0;
		
		c.defenderreward=0;
		c.defenderpenalty= -9;
		c.attackerreward = 9;
		c.attackerpenalty= 0;
		
		
		
		
		gamedata[1][0] = 0;//(int)density[iter][i] ;
		gamedata[1][1] = -9;
		gamedata[1][2] = (int)9;  // uncovered
		gamedata[1][3] = 0; //covered
		
		
		gamedata[2][0] = 0;//(int)density[iter][i] ;
		gamedata[2][1] = -9;
		gamedata[2][2] = (int)9;  // uncovered
		gamedata[2][3] = 0; //covered
		
		
		
		graph.add(a);
		graph.add(b);
		graph.add(c);
		
		
		return graph;
	}

	private static double[] doubleOracleGCMultiExactLPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {


		

		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];
		long slavetime = 0;






		boolean canaddpath = true;
		
		int masteritr=0;

		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			
			SecurityGameContraction.origtargets.clear();
			copyInOrigTargets(tmpgraph);
			
			
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);
			
			
			///// originalgraph 
			
			

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				
				
				itr++;
				System.out.println("Entered inner loop...slave itr: "+itr);

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if((currentPlace<targetssorted.length-1) && attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<v : "+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					//System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();//buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//test
					//ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>(); //MIPSolver4.slavePath(1, gamedata, tmpgraph, nRes, nTargets, dmax);
					
					
					
					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);


					
					//copyInOrigTargets();
					
					
					
					///slave
					start = new Date();
					l1 = start.getTime();


					
					
					
					
					
					SecurityGameContraction.addVirtualBase(0,nTargets, tmpgraph);

					//System.out.println("Added virtual base , tmpgraph size "+ tmpgraph.size());

					//int q = SecurityGameContraction.calculateEdgesWithCoin();
					HashMap<Integer, Integer> nodewithcoins= SecurityGameContraction.calculateNodesWithCoin(tmpgraph);

					SecurityGameContraction.transformToDirectedGraph(tmpgraph);
					
					//System.out.println("after transformed to directed , tmpgraph size "+ tmpgraph.size());
					//System.out.println("after transformed to directed , duplicate size "+ SecurityGameContraction.duplicatetargets.size());


					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, SecurityGameContraction.duplicatetargets);
					ArrayList<Integer> pathtoconsider = MIPSolver3.solve(SecurityGameContraction.duplicatetargets, tmpgraph,
							nodewithcoins, domindatednodes, dmax, nTargets);
					removeVirtualBase(tmpgraph, nTargets);
					//System.out.println("After removing virtual base tmpgrap size "+ tmpgraph.size());
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					

					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					
					
					
					
					newpathseq.add(pathtoconsider);
					
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq by slave : ");
					printPaths(newpathseq);
					
					/**
					 * need to determine which path seq can be added depending on probability distribution and origpmat
					 */
					
					System.out.println("newpathseq size before purify : "+newpathseq.size());
				    newpathseq = determineNewPaths(newpathseq, origpmat, probdistribution);
					System.out.println("newpathseq size after purify : "+newpathseq.size());
					
					
					if((newpathseq.size()==0) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("New whole path seq ");


					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					/*pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");
					
					if((newsize==oldsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}*/
					


					
					printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v&&cantaddpath/alltargetsadded"+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/
			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru, slavetime};
		return res;
	}
	
	
	public static ArrayList<ArrayList<Integer>> determineNewPaths(ArrayList<ArrayList<Integer>> newpathseq, int[][] origpmat,
			double[] probdistribution) {
		
		ArrayList<ArrayList<Integer>> purifiedpaths = new ArrayList<ArrayList<Integer>>();
		
		
		/*for(int i=0; i<origpmat.length; i++)
		{
			for(int j=0; j<origpmat[i].length; j++)
			{
				System.out.print(origpmat[i][j]+"  ");
			}
			System.out.print("\n");
		}*/
		
		//compare path not joint schedule
		
		for(ArrayList<Integer> p: newpathseq)
		{
			boolean f=false;
			
			
			
			
			for(int j=0; j<origpmat[0].length; j++)
			{
				f=false;
				for(Integer t: p)
				{
					if(origpmat[t][j]!=1)
					{
						f=true;
						break;
					}
				}
				if(!f)
				{
					break;
				}
				
				
			}
			
			if(f)
			{
				purifiedpaths.add(newpathseq.get(newpathseq.indexOf(p)));
			}
		}
		
		
		return purifiedpaths;
	}

	private static boolean isCovered(Integer t, int[][] origpmat, double[] probdistribution) {
		
		
		for(int j=0; j<origpmat[t].length; j++)
		{
			if(origpmat[t][j]==1 && probdistribution[j]>0)
				return true;
		}
		
		return false;
	}

	private static double[] doubleOracleGCMultiGP3APSPLPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}*/


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			
			
			/*ArrayList<Integer> domint = new ArrayList<Integer>();
			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
				domint.add(s.getTargetid());
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			//System.out.println("tmpgraph size "+ tmpgraph.size());
			//System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			//System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			//System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					System.out.println("iter "+ iter);
					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					//origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					//System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					//System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					//System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					//System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}
	
	
	
	
	private static double[] doubleOracleAPSPGCMultiGP3LPGCMulti(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}
*/

			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			/*if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}*/



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	






	private static double[] doubleOracleGCMultiGP3LPSamplePath(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		//printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		//ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			/*System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}*/


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			/*System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();*/

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContractionWithAPSP(domindatednodes, tmpgraph, dmax);
			//sgc.contractGraph(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3WithAPSP(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = new ArrayList<ArrayList<Integer>>();
					//ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<5 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 5-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}





			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	







	private static double[] contractionWithGreedyDoubleOracleWGreedyCover3LP(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution= new double[jset.size()];
		double attackeru= -999;
		double attackerv = -999;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat = new int[nTargets][];






		boolean canaddpath = true;

		while(true)
		{

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq =  generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			if(pathseq.size()==0)
			{
				System.out.println("pathseq 0, iter..nqqqqq"+ iter);
			}



			int itr=0;
			while(true)
			{
				
				itr++;

				

				/*if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ iter);
				}*/


				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter"+ iter);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, tmpgraph, nRes, attackerstrategy);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes, attackerstrategy);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("haa ");


					if((oldsize==newsize) || (itr>=10))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
				System.out.println("iter"+ itr++);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			/*if(addcount<3 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 3-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}*/






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	










	private static double[] noContractionWithColumnGeneration(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = new ArrayList<Integer>();//buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/
		for(int i=0; i<targetssorted.length; i++)
		{
			currenttargets.add(targetssorted[i][0]);
		}



		//int currentPlace = currenttargets.size()-1;


		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;
		double attackerv;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat;






		boolean canaddpath = true;


		pathseq = new ArrayList<ArrayList<Integer>>();

		//System.out.println("\nCurrent place : "+ currentPlace);

		System.out.print("Current target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}


		//tmpgraph = getDuplicateGraph(targets);


		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


		Date start = new Date();
		long l1 = start.getTime();


		//instantContraction(domindatednodes, tmpgraph, dmax);


		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;

		contractiontime += diff;
		//SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
		//SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

		System.out.println("tmpgraph size "+ targets.size());
		//System.out.println("dom size "+ domindatednodes.size());
		//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
		p = new int[targets.size()][]; // p matrix

		//apply greedy approach
		//TODO generate paths where there will be at least one target
		ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, targets);

		makePathSeq(pathseq, goals, goals.size(), targets.size(), map, mapback, targets);
		//printPaths(pathseq);
		System.out.println("Total path with duplicates "+pathseq.size());
		pathseq = removeDuplicatePathSimple(pathseq);
		System.out.println("Total path without duplicates "+pathseq.size()+"\n");
		//printPaths(pathseq);

		/**
		 * keep only nRes*3 paths from the end
		 */

		//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
		//System.out.println("Initial number of paths "+ pathseq.size());
		//printPaths(pathseq);

		int internaliter = 0;


		while(true)
		{
			internaliter++;

			canaddpath = true;

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;

				for(TargetNode x: targets)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}
			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				p = makePmat(pathseq, jset, mapback, targets);
				//printPathMat(p);

				start = new Date();
				l1 = start.getTime();

				HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

				probdistribution = MIPSolver4.solveForAttackerLP(p, gamedata, targets, nRes, attackerstrategy);



				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvingtime += diff;

				attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
				//attackedtarget = mapback.get(attackedtarget);
				System.out.println("attack target before rev map "+ attackedtarget);
				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("attacker u= "+attackeru);

				//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);





				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}



				/**
				 * apply greedy slave
				 * 
				 * find the attack target and find a path that includes that target

				 */
				System.out.println("attacked target after rev map "+ attackedtarget);

				ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(targets,dmax, targets.size(), 0, nRes);

				if(newpathseq.size()==0)
				{
					canaddpath = false;
					System.out.println("Slave can't add any new path ###############");
					break;
				}
				System.out.println("tcur: ");
				//printGreedyPath(currenttargets);
				System.out.println("newpathseq: ");
				//printPaths(newpathseq);

				System.out.println("Old path seq size "+ pathseq.size());

				int oldsize = pathseq.size();
				for(ArrayList<Integer> q: newpathseq)
				{
					pathseq.add(q);
				}

				System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

				pathseq = removeDuplicatePathSimple(pathseq);
				System.out.println("New path seq size "+ pathseq.size());
				//printPaths(pathseq);
				int newsize = pathseq.size();
				System.out.println("haa ");


				if((oldsize==newsize) || (internaliter>=10))
				{
					canaddpath = false;
					System.out.println("Slave can't add any new path ###############");
					break;
				}




			} // end if else
			
		} // inner while loop 





		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, p, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	









	private static double[] contractionWithGreedyDoubleOracleWGreedyCover(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		//ArrayList<Integer> currenttargets = buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes); //  new ArrayList<Integer>();
		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;
		double attackerv;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat;






		boolean canaddpath = true;
		
		

		while(true)
		{
			

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;
			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[targets.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			//pathseq = buildGreedyPathMultRes(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq=buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);

			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<tmpgraph.size(); i++)
			{

				map.put(tmpgraph.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, tmpgraph.get(i).getTargetid());
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			//ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);
			//System.out.println("Initial number of paths "+ pathseq.size());
			//printPaths(pathseq);


			int itr = 0;

			while(true)
			{
				itr++;

				canaddpath = true;

				Integer[] input = new Integer[pathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(pathseq.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}
				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(pathseq.size()<nRes)
					{

						branch = new int[pathseq.size()];
						jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(pathseq, jset, mapback, tmpgraph);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					//HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					//attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					//attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(MIPSolver4.attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target after rev map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					/*if(attackeru>= targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.$$$$$$$$$$$$$$$$$..attacker u>=v="+attackeru);
						break;
					}*/
					if(attackeru<targetssorted[currentPlace+1][1])
					{
						System.out.println("inner loop ....breaking.%%%%%%%%%%..attacker u<=v="+attackeru);
						break;
					}
					/**
					 * apply greedy slave
					 * 
					 * find the attack target and find a path that includes that target

					 */
					System.out.println("attacked target after rev map "+ attackedtarget);

					ArrayList<ArrayList<Integer>> newpathseq = buildGreedyCoverMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);

					//ArrayList<TargetNode> goal = generatePathsSlave(dmax, gamedata, tmpgraph, attackedtarget, nRes, currenttargets);



					//makeSlavePathSeq(newpathseq, goal);
					//removeDuplicatePathSimple(newpathseq);
					if(newpathseq.size()==0)
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}
					System.out.println("tcur: ");
					//printGreedyPath(currenttargets);
					System.out.println("newpathseq: ");
					//printPaths(newpathseq);

					System.out.println("Old path seq size "+ pathseq.size());

					int oldsize = pathseq.size();
					for(ArrayList<Integer> q: newpathseq)
					{
						pathseq.add(q);
					}

					System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

					pathseq = removeDuplicatePathSimple(pathseq);
					System.out.println("New path seq size "+ pathseq.size());
					//printPaths(pathseq);
					int newsize = pathseq.size();
					System.out.println("itr "+itr);


					if((oldsize==newsize) || (itr>=15))
					{
						canaddpath = false;
						System.out.println("Slave can't add any new path ###############");
						break;
					}

					//printPaths(pathseq);


				} // end if else
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					
					addcount++;
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);

					if(addcount>=5)
					{
						break;
					}
				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			/*if(addcount<3 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 3-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}*/






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}	









	public static double[] contractionWithDoubleOracle(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double
			dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0); //  new ArrayList<Integer>();
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;
		double attackerv;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		int [][] origpmat;








		while(true)
		{

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.origtargets.clear();
			copyInOrigTargets(tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);

			System.out.println("tmpgraph size "+ tmpgraph.size());
			System.out.println("dom size "+ domindatednodes.size());




			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);



			p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			/**
			 * keep only nRes*3 paths from the end
			 */

			ArrayList<ArrayList<Integer>> initpaths =	filterPaths(pathseq, 3*nRes, currenttargets);

			System.out.println("Initial number of paths "+ initpaths.size());

			printPaths(initpaths);


			while(true)
			{


				Integer[] input = new Integer[initpaths.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				jSet=new HashSet();
				if(initpaths.size()==0)
				{
					//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
					//choose the worst payoff for defender

					Double mAxpayoff = Double.MIN_VALUE;
					Double defpayoff = 0.0;
					for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}
					for(TargetNode x: tmpgraph)
					{
						if(x.attackerreward>mAxpayoff)
						{
							mAxpayoff= x.attackerreward;
							defpayoff = x.defenderpenalty;
						}
					}



				}
				else
				{
					//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
					if(initpaths.size()<nRes)
					{

						branch = new int[initpaths.size()];
						jSet=combine(input, initpaths.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}

					jset = new ArrayList<ArrayList<Integer>>(jSet);
					/**
					 * columns will be combination of paths for each resources. 
					 */
					/**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = makePmat(initpaths, jset, mapback, tmpgraph);


					start = new Date();
					l1 = start.getTime();

					probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);



					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;

					attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target wo map "+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
					System.out.println("attacker u= "+attackeru);



					origpmat = makeOrigPMatWOMap(p, initpaths, jset, nTargets, domindatednodes, map, mapback, tmpgraph);
					attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);
					System.out.println("attack target w map"+ attackedtarget);
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);
					System.out.println("attacker v= "+attackerv);



					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}*/

					if(attackeru>= attackerv)
					{
						System.out.println("inner loop ....breaking...attacker u>=v="+attackeru);
						break;
					}


					/**
					 * solve the coin collection problem to add new path to s
					 */

					//copyInOrigTargets();
					SecurityGameContraction.addVirtualBase(0,nTargets, tmpgraph);

					System.out.println("Added virtual base , tmpgraph size "+ tmpgraph.size());

					//int q = SecurityGameContraction.calculateEdgesWithCoin();
					HashMap<Integer, Integer> nodewithcoins= SecurityGameContraction.calculateNodesWithCoin(tmpgraph);

					SecurityGameContraction.transformToDirectedGraph(tmpgraph);


					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, SecurityGameContraction.duplicatetargets);
					ArrayList<Integer> pathtoconsider = MIPSolver3.solve(SecurityGameContraction.duplicatetargets, tmpgraph,
							nodewithcoins, domindatednodes, dmax, nTargets);
					removeVirtualBase(tmpgraph, nTargets);
					System.out.println("After removing virtual base tmpgrap size "+ tmpgraph.size());

					System.out.print("\npath to consider : ");
					for(int i=0; i<pathtoconsider.size(); i++)
					{
						System.out.print(pathtoconsider.get(i) + " ");
					}
					System.out.println();

					boolean yes = canBeAdded(origpmat, pathtoconsider, probdistribution);


					if(yes)
					{
						System.out.println("new path added : ");
						initpaths.add(pathtoconsider);

						printPaths(initpaths);
					}
					else
					{
						System.out.println("new path was not added : ");
						initpaths.add(pathtoconsider);

						printPaths(initpaths);
						break;

					}

				} // end if else

				SecurityGameContraction.duplicatetargets.clear();





			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if(currentPlace==targetssorted.length-1 || (attackeru>= attackerv))
			{
				System.out.println("outer loop ....breaking...attacker u>=v="+attackeru);
				break;
			}




			double ulimit = getTargetNode(attackedtarget, targets).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					currenttargets.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);

				}
			}

			System.out.println("addcount : "+ addcount);

			currentPlace = currenttargets.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			/*if(addcount<3 || addcount==0)
			{
				System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += 3-addcount;

				System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length)
				{
					currentPlace = targetssorted.length;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					currenttargets.add(targetssorted[k][0]);
				}
			}*/






			/*int prevcur = currentPlace;
			currentPlace += 3;
			if(currentPlace>targetssorted.length)
			{
				currentPlace = targetssorted.length;
			}
			System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

			for(int k= prevcur+1; k<=currentPlace; k++ )
			{
				currenttargets.add(targetssorted[k][0]);
			}


			break;*/


		} // outer while loop

		System.out.println("Final target list size : "+ currenttargets.size());

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}


	private static boolean canBeAdded(int[][] pmat,
			ArrayList<Integer> pathtoconsider, double[] probdistribution) {


		boolean covered = false;

		for(Integer nodeid: pathtoconsider)
		{
			covered = false;
			// check if nodeid is covered if not return true
			for(int i=0; i<probdistribution.length; i++)
			{
				if(pmat[nodeid][i]==1)
				{
					covered=true;
					break;
				}
			}
			if(!covered)
			{
				return true;
			}

		}
		return false;
	}



	private static void removeVirtualBase(ArrayList<TargetNode> tmpgraph,
			int virtualbaseid) {

		TargetNode v = getTargetNode(virtualbaseid, tmpgraph);
		TargetNode b = getTargetNode(0, tmpgraph);

		b.removeDistance(v);
		b.removeNeighbor(v);
		tmpgraph.remove(v);



	}

	public static void copyInOrigTargets(ArrayList<TargetNode> targets) {

		for(TargetNode x:  targets)
		{
			SecurityGameContraction.origtargets.add(x);
		}

	}



	private static ArrayList<ArrayList<Integer>> filterPaths(ArrayList<ArrayList<Integer>> pathseq, int n, ArrayList<Integer> currenttargets) {


		ArrayList<ArrayList<Integer>> tmp = new ArrayList<ArrayList<Integer>>();


		int count=0;


		for(int i = 0; i< pathseq.size(); i++)
		{

			//check if atleast 30 percent of targets should be in the list

			int shouldcount = (int)Math.ceil(0.4*currenttargets.size());
			boolean ok = checkIfCountIsOk(shouldcount, currenttargets, pathseq.get(i));
			if(ok)
			{
				tmp.add(pathseq.get(i));
			}
		}


		return tmp;


	}

	private static boolean checkIfCountIsOk(int shouldcount,
			ArrayList<Integer> currenttargets, ArrayList<Integer> path) {



		int count = 0;

		for(Integer x: currenttargets)
		{
			if(path.contains(x))
			{
				count++;
			}
		}

		if(count>=shouldcount)
			return true;
		return false;
	}

	private static double[] contractionWithSingleOracleRevMap(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = new ArrayList<Integer>();
		currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);



		int currentPlace = 1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;

		int[][] origpmat;


		while(true)
		{

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);



			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);



			p = new int[targets.size()][]; // p matrix
			ArrayList<TargetNode> goals = generatePaths(dmax, gamedata, tmpgraph);
			pathseq = new ArrayList<ArrayList<Integer>>();
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);


				start = new Date();
				l1 = start.getTime();

				probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);



				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvingtime += diff;


				/**
				 * do rev map
				 */

				origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

				attackedtarget = findAttackTarget(origpmat, probdistribution, gamedata);





				/*double defexpectedpayoff = expectedPayoffDef(maxtargetforattacker, origpmat, gamedata, probdistribution);
				defexp = defexpectedpayoff;





				attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);

				attackedtarget = mapback.get(attackedtarget);*/

				System.out.println("attack target "+ attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedPayoffAtt(attackedtarget, origpmat, gamedata, probdistribution);  //expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("attacker u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}

				if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
				{
					System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

					break;
				}
				else 
				{
					int prevcur = currentPlace;
					currentPlace += 3;
					if(currentPlace>targetssorted.length)
					{
						currentPlace = targetssorted.length;
					}
					System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

					for(int k= prevcur+1; k<=currentPlace; k++ )
					{
						currenttargets.add(targetssorted[k][0]);
					}
				}



				System.out.println();


			}


		}

		System.out.println("Final target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		//double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);
		double defpayoff = expectedPayoffDef(attackedtarget, origpmat, gamedata, probdistribution);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}





	private static double[] contractionWithGreedySingleOracle(int[][] gamedata,
			int nTargets, int nRes, double[][] density, double dmax, int iter, int nrow, int ncol) throws Exception {



		targets.clear();
		SecurityGameContraction sgc = new SecurityGameContraction(nrow, ncol, gamedata);
		assignRandomDensityZeroSum(density, gamedata, targets, iter);


		//printtargets(targets);

		/**
		 * 1. sort the targets
		 */
		int[][] targetssorted = sortTargets(targets);
		printSortedTargets(targetssorted);

		ArrayList<Integer> currenttargets = buildGreedyCover(targets, dmax, nTargets, 0);
		/*currenttargets.add(targetssorted[0][0]);
		currenttargets.add(targetssorted[1][0]);*/



		int currentPlace = currenttargets.size()-1;


		ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		int attackedtarget=-1;
		int[][] p;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		HashSet jSet=new HashSet();
		ArrayList<ArrayList<Integer>> pathseq = new ArrayList<ArrayList<Integer>>();
		List<ArrayList<Integer>> jset = new ArrayList<ArrayList<Integer>>(jSet);
		double[] probdistribution;
		double attackeru;



		long contractiontime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;


		while(true)
		{

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<currenttargets.size(); i++)
			{
				System.out.print(currenttargets.get(i)+",");
			}


			tmpgraph = getDuplicateGraph(targets);
			if(currentPlace<targetssorted.length-1)
				domindatednodes = selectDominatedNodes(targetssorted, currentPlace+1, tmpgraph);
			else
			{
				domindatednodes.clear();
			}

			System.out.print("\nDom targets : ");
			for(TargetNode s: domindatednodes)
			{
				System.out.print(s.getTargetid()+" ");
			}
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			instantContraction(domindatednodes, tmpgraph, dmax);


			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			contractiontime += diff;



			SecurityGameContraction.removePathsToDominatedNodes(domindatednodes, tmpgraph);
			SecurityGameContraction.removeDominatedTargets(domindatednodes, tmpgraph);



			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);



			p = new int[targets.size()][]; // p matrix
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			pathseq = new ArrayList<ArrayList<Integer>>();
			pathseq = generatePathsGreedy3(dmax, gamedata, tmpgraph, currenttargets, nRes);
			
			
			/**
			 * map has present id
			 * mapback gives the original ids
			 */
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount =0;
			for(int i=0; i<targets.size(); i++)
			{
				map.put(targets.get(i).getTargetid(), icount);
				//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;
			}

			
			
		//	makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");
			//printPaths(pathseq);

			Integer[] input = new Integer[pathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			jSet=new HashSet();
			if(pathseq.size()==0)
			{
				//System.out.println("pathseq 0, iter"+ iter+", contrac "+ contractionsize);
				//choose the worst payoff for defender

				Double mAxpayoff = Double.MIN_VALUE;
				Double defpayoff = 0.0;
				for(int i=0; i<domindatednodes.size(); i++)
				{
					tmpgraph.add(domindatednodes.get(i));
				}
				for(TargetNode x: tmpgraph)
				{
					if(x.attackerreward>mAxpayoff)
					{
						mAxpayoff= x.attackerreward;
						defpayoff = x.defenderpenalty;
					}
				}


				//System.out.println("Defender expected payoff "+ defpayoff);
				/*try
				{
					PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+contractionsize+".csv"),true));
					pw.append(iter+ "," + defpayoff+"\n");
					pw.close();

				}
				catch(Exception e)
				{

				}*/



			}
			else
			{
				//System.out.println("pathseq "+pathseq.size()+", iter"+ iter+", contrac "+ contractionsize);
				if(pathseq.size()<nRes)
				{

					branch = new int[pathseq.size()];
					jSet=combine(input, pathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}

				jset = new ArrayList<ArrayList<Integer>>(jSet);
				/**
				 * columns will be combination of paths for each resources. 
				 */
				/**
				 * pmat, where columns will be combination of paths. 
				 * rows are targets. 
				 * each entry will say whether the target is in the joint schedule
				 */
				//jSet.

				//printJointSchedule(jset);

				//printNodesAsNeighbors(dominatednodes);

				p = makePmat(pathseq, jset, mapback, tmpgraph);
				//printPathMat(p);


				start = new Date();
				l1 = start.getTime();

				probdistribution = MIPSolver4.solveForAttacker(p, gamedata, tmpgraph, nRes);



				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;

				solvingtime += diff;


				attackedtarget = findAttackTargetWMapping(p, probdistribution, gamedata, map, mapback);

				attackedtarget = mapback.get(attackedtarget);

				System.out.println("attack target "+ attackedtarget);

				//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
				attackeru = expectedAttackerPayoff(attackedtarget, p, probdistribution, gamedata, map);
				System.out.println("attacker u= "+attackeru);

				if(probdistribution.equals(null))
				{
					throw new Exception("Prob null...");
				}

				if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
				{
					System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

					break;
				}
				else 
				{
					int prevcur = currentPlace;
					currentPlace += 1;
					if(currentPlace>targetssorted.length)
					{
						currentPlace = targetssorted.length;
					}
					System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

					for(int k= prevcur+1; k<=currentPlace; k++ )
					{
						currenttargets.add(targetssorted[k][0]);
					}
				}



				System.out.println();


			}


		}

		System.out.println("Final target list : ");

		for(int i=0; i<currenttargets.size(); i++)
		{
			System.out.print(currenttargets.get(i)+",");
		}

		double defpayoff = expectedDefenderPayoff(attackedtarget, p, probdistribution, gamedata, map);




		//int[][] origpmat = makeOrigPMatWOMap(p, pathseq, jset, nTargets, domindatednodes, map, mapback, targets);

		double[] res = {defpayoff, contractiontime, solvingtime, currenttargets.size(), attackeru};
		return res;
	}

	public static ArrayList<TargetNode> selectDominatedNodes(
			int[][] targetssorted, int start, ArrayList<TargetNode> tmpgraph) {


		ArrayList<TargetNode> dom = new ArrayList<TargetNode>();


		for(int i=start; i<targetssorted.length; i++)
		{
			dom.add(tmpgraph.get(targetssorted[i][0]));
		}


		return dom;
	}
	
	
	private static ArrayList<TargetNode> selectDominatedNodesWMap(
			int[][] targetssorted, int start, ArrayList<TargetNode> tmpgraph,  HashMap<Integer, TargetNode> tmpgraphmap) {


		ArrayList<TargetNode> dom = new ArrayList<TargetNode>();


		for(int i=start; i<targetssorted.length; i++)
		{
			dom.add(tmpgraphmap.get(targetssorted[i][0]));
		}


		return dom;
	}

	public static void printSortedTargets(int[][] targetssorted) {

		System.out.println(" Sorted Targets : ");

		for(int i=0; i<targetssorted.length; i++)
		{

			System.out.println("Target "+targetssorted[i][0] + ", u = "+ targetssorted[i][1]);

		}


	}

	public static int[][] sortTargets(ArrayList<TargetNode> targets) {

		int[][] srted = new int[targets.size()][2];
		for(int i=0; i<srted.length; i++)
		{
			srted[i][0] = targets.get(i).getTargetid();
			srted[i][1] = (int)targets.get(i).attackerreward;
		}
		int[] swap = {0,0};

		for (int k = 0; k < srted.length; k++) 
		{
			for (int d = 1; d < srted.length-k; d++) 
			{
				if (srted[d-1][1] < srted[d][1])    // ascending order
				{
					swap = srted[d];
					srted[d]  = srted[d-1];
					srted[d-1] = swap;
				}
			}
		}
		return srted;
	}
	
	
	public static int[][] sortSuperTargets(HashMap<Integer, SuperTarget> sts) {

		int[][] srted = new int[sts.size()][2];
		int i = 0;
		for(SuperTarget st: sts.values())
		{
			srted[i][0] = st.stid;
			srted[i][1] = (int)st.attackerreward;
			i++;
		}
		int[] swap = {0,0};

		for (int k = 0; k < srted.length; k++) 
		{
			for (int d = 1; d < srted.length-k; d++) 
			{
				if (srted[d-1][1] < srted[d][1])    // ascending order
				{
					swap = srted[d];
					srted[d]  = srted[d-1];
					srted[d-1] = swap;
				}
			}
		}
		return srted;
	}






}
