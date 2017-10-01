package groupingtargets;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;

import cs.Interval.ILP.MIPSolver4;
import cs.Interval.contraction.Logger;
import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;
import cs.com.allpair.AllPairShortestPath;
import cs.com.realworld.ReadData;
import weka.clusterers.EM;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;


public class GroupingTargets {

	public static Random rand1 = new Random(100);
	public static final int INFINITY = 99999;
	
	
	
	public static void testWeka() throws Exception {
		
		/*FileInputStream fstream = new FileInputStream("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/realdata2.csv");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		Instances instances = new Instances(br);
		 */
		 
		 
		 // Print header and instances.
		// System.out.println("\nDataset:\n");
		 //System.out.println(instances.toSummaryString());
		 
		/* SimpleKMeans model = new SimpleKMeans();
		 model.setNumClusters(10);
		 model.buildClusterer(instances);
		 model.setDistanceFunction(new weka.core.ManhattanDistance());
		 
		 System.out.println(model);*/
		 
		/* MakeDensityBasedClusterer dc = new MakeDensityBasedClusterer();
		 dc.setNumClusters(10);
		
		 dc.buildClusterer(instances);
		 System.out.println(dc);
		 
		 
		 
		 for(int i=0; i<instances.size(); i++)
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(instances.get(i)));
		 
			 
		 }*/
		 
		 
		 
		 
		 
		/* EM em = new EM();
		 em.setNumClusters(10);
		 em.buildClusterer(instances);
		 System.out.println(em);*/
		 
		 
		 
		 CSVLoader csvload = new CSVLoader();
		 csvload.setSource(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/realdata3.csv"));
		 Instances data = csvload.getDataSet();
		 
		 
		 ArffSaver arf = new ArffSaver();
		 arf.setInstances(data);
		 arf.setFile(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/newdata3.arff"));
		 arf.writeBatch();
		
		
	}
	
	
	public static void createCSVTestData()
	{
		int nrow = 10;
		int ncol = 10;
		double utility [][] = new double[560][560];
		double elevation [][] = new double[560][560];
		double u [][] = new double[nrow][ncol];
		double e [][] = new double[nrow][ncol];
		
		
		
		
		ReadData.readData(560, 560, utility, elevation);
		//ReadData.getChunk(560, 560, utility, elevation, 0, 0, nrow, ncol, u, e);
		
		
		//ReadData.createCSVData(560, 560, utility, elevation);
		//int[][] gamedata = SecurityGameContraction.constructGameData(utility);
		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		SecurityGameContraction.buildcsvGraph(nrow,ncol,utility, elevation,targets );
	}
	
	
	
	public static HashMap<Integer, SuperTarget> clusterTargets(ArrayList<Integer> targetstocluster, 
			ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps, double dmax, int k, int radius,
			HashMap<Integer, Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths)
	{
		
		
		/**
		 * 1. start with every node in its own cluster
		 */
		
		
		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[targetmaps.size()];
		for(int t=0; t<targetmaps.size(); t++)
		{
			clusters[t] = new ArrayList<Integer>();
			clusters[t].add(t);
		}
		HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,targetmaps);
		// printSuperTargets(sts);
		 int stssize = 1; // to keep track if any targets were clustered
		 while(true)
		 {
			 if(stssize == sts.size())
				 break;
			 
			 stssize = sts.size();
			 /**
				 * 2. For every pair of cluster: minimize the distance d : d1 + d(a1,a2) + d2 
				 * 3. Find a1 and a2 which will minimize the distance
				 */
			 
			 /**
			  *  Make a temporary merged supertarget
			  *  Then for every pair of access points
			  *  try to find the minimum one
			  *  save it the config for a particular d
			  *  
			  */
			 
			 System.out.println("Computing dmin");
			 double mindi = Double.MAX_VALUE;
			 int stid1 = -1;
			 int stid2 = -1;
			 int aid1 = -1;
			 int aid2 = -1;
			 double sda1a2 = -1;
			 ArrayList<Integer> spath = new ArrayList<Integer>();
			 SuperTarget dminst = new SuperTarget();
			 for(SuperTarget st1 : sts.values())
			 {
				 
				if(canBeMerged(st1, targetstocluster))
				{
				 
				 for(SuperTarget st2: sts.values())
				 {
					 
					 if(canBeMerged(st2, targetstocluster))
					 {
					 
					 //if they are not same
					 // check if st1 and st2 are neighbprs
					 // have nodes in n eighbor
					 boolean isneighbor = isNeighbor(st1, st2);
					 // min dist between two st needs to be within threshold to be considered as potential merge
					 double mindist = 500000;
					 if(isneighbor)
						 mindist = minimumDist(st1, st2);
					if( (st1.stid != st2.stid) && (isneighbor == true) 
							 && st1.stid != 0 && st2.stid != 0 && (mindist <= radius ))
					 {
						 //System.out.println("SuperTarget "+ st1.stid+ ", and "+st2.stid + " will be merged");
						 //merge them and create a temp supertarget
						 SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
						 tempst.stid = 200 + st1.stid + st2.stid;
						// printSuperTarget(tempst);
						/**
						  * for every pair of access points
						  */
						 for(TargetNode a1: tempst.ap.values())
						 {
							 /**
							  * 	for every pair of other supertargets which are not st1 and st2
							  * */
							 for(TargetNode a2: tempst.ap.values())
							 {
								 /**
								  * 		find the min d = d1 + dist(a1,a2) + d2
								  */
								 // should the ap be same for entry and exit ?
								 if(a1.getTargetid() != a2.getTargetid())
								 {
									 // for every pair of supertargets which are not st1 or st2
									 for(SuperTarget s1: sts.values())
									 {
										 for(SuperTarget s2: sts.values())
										 {
											 if(s1.stid != st1.stid && s1.stid != st2.stid && 
													 s2.stid != st1.stid && s2.stid != st2.stid)
											 {
												 // measure the distance between a1->s1 and a2->s2
												 
												 //a1->s1 distance between a target and a supertarget
//												 System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
//														 "\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
//														 "\n tempst "+ tempst.stid);
												 
												 
												 double d1 = shortestdist(a1,s1);
												// System.out.println("shortestdist(a1,s1) "+ d1);
												 
												 double d2 = shortestdist(a2,s2);
												// System.out.println("shortestdist(a2,s2) "+ d2);
												 
												 // next measure the intra cluster shortest traveling path using a1 and a2
												 ArrayList<Integer> tmpspath = new ArrayList<Integer>();
												 double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpspath);
												// System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);
												 
												 // if any of the dist is <0 we know that it's not possible to have a path
												 
												 if(d1 > 0 && d2 > 0 && dista1a2 > 0)
												 {
													 double totaldi = d1 + dista1a2 + d2;
													 if(totaldi < mindi)
													 {
														 stid1 = st1.stid;
														 stid2 = st2.stid;
														 dminst = tempst;
														 aid1 = a1.getTargetid();
														 aid2 = a2.getTargetid();
														 mindi = totaldi;
														 sda1a2 = dista1a2;
														 spath.clear();
														 for(Integer in: tmpspath)
														 {
															 spath.add(in);
														 }
														 
														 System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
																 "\n a1 "+ aid1 + ", a2 "+ aid2);
														 
														 
													 }
												 }
												 
											 }
										 }
									 }
								 }
							 }
						 }

					 }
				 }
					  
				 }
			 	}
			 }
			 
			 
			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
				//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			 /**
			  * Merge two cluster which have min d
			  * For the new id of supertargets concat the strings with comma
			  * remove the old two clusters
			  * add the new one. 
			  */
			 
			 if(stid1 != -1)
			 {

				 dstravel.put(200+stid1+stid2, sda1a2);
				 stpaths.put(200+stid1+stid2, spath);
				 SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
				// printSuperTarget(newst);
				 sts.remove(stid1);
				 sts.remove(stid2);
				 updateNeighbors(newst, sts, stid1, stid2);
				 sts.put(newst.stid, newst);
				 //update the neighbors of ST
				// System.out.println("\n\n After merging # supertargets : "+ sts.size());
				// System.out.println("\n After merging new supertargets : ");
				// printSuperTargets(sts);
			 }
			 System.out.println("hi");
			 
		 }
		 return sts;
	}
	
	
	

	
	

	public static HashMap<Integer, SuperTarget> clusterTargetsWekaRW(ArrayList<Integer> targetstocluster, 
			ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps, double dmax, int k, int radius,
			HashMap<Integer, Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, EM dc, Instances instances) throws Exception
	{
		
		
		
		
		
		/**
		 * make new instances with targets to cluster
		 */
		
		Instances newinstance = new Instances(instances);
		
		
		
		
		
		/**
		 * remove unwanted instance
		 * index is the target id
		 */
		
		for(Integer tid: targetmaps.keySet())
		{
			if(!targetstocluster.contains(tid))
			{
				
				
				
				
				int index = findIndex(newinstance, tid);
				
				if(index != -1)
				{
					newinstance.remove(index);
				}
				
			}
		}
		
		
		System.out.println("New instance size  "+ newinstance.size());
		
		/**
		 * cluster
		 */
		
		int totalcluster = k + (instances.size()-newinstance.size());
		
		ArrayList<Integer>[] clusters = clusterWithWekaRW(k, newinstance, totalcluster, targetstocluster,targetmaps, dc);
		
		
		//ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[targetmaps.size()];
		
		/*for(int t=0; t<targetmaps.size(); t++)
		{
			clusters[t] = new ArrayList<Integer>();
			clusters[t].add(t);
		}
		*/
		
		
		
		
		
		
		HashMap<Integer, SuperTarget> sts = SuperTarget.buildSuperTargets(clusters,targetmaps);
		
		/**
		 * for every super target determine 2 AP
		 */
		
		



		 printSuperTargets(sts);
		int stssize = 1; // to keep track if any targets were clustered
		
		/**
		 * 2. For every pair of cluster: minimize the distance d : d1 + d(a1,a2) + d2 
		 * 3. Find a1 and a2 which will minimize the distance
		 */

		/**
		 *  Make a temporary merged supertarget
		 *  Then for every pair of access points
		 *  try to find the minimum one
		 *  save it the config for a particular d
		 *  
		 */

		System.out.println("Computing dmin");
		
		
		
		
		

		for(SuperTarget tempst: sts.values())
		{
			//SuperTarget tempst = SuperTarget.mergeSuperTargets(st1,st2);
			//tempst.stid = 200 + st1.stid + st2.stid;
			// printSuperTarget(tempst);
			/**
			 * for every pair of access points
			 */
			
			System.out.println("Computing ap for ST "+ tempst.stid);
			
			
			double mindi = Double.MAX_VALUE;
			int stid1 = -1;
			//int stid2 = -1;
			int aid1 = -1;
			int aid2 = -1;
			int s1id = -1;
			int s2id = -2;
			double sda1a2 = -1;
			ArrayList<Integer> spath = new ArrayList<Integer>();
			SuperTarget dminst = new SuperTarget();
			
			double mindist = 500000;
			
			
			
			

			if(tempst.nodes.size()>1)
			{
				
				if(tempst.stid==14)
				{
					System.out.println("shortestdist(a1,s1) ");
				}

				for(TargetNode a1: tempst.ap.values())
				{
					/**
					 * 	for every pair of other supertargets which are not st1 and st2
					 * */
					for(TargetNode a2: tempst.ap.values())
					{
						/**
						 * 		find the min d = d1 + dist(a1,a2) + d2
						 */
						// should the ap be same for entry and exit ?
						if(a1.getTargetid() != a2.getTargetid())
						{
							// for every pair of supertargets which are not st1 or st2
							for(SuperTarget s1: tempst.neighbors.values()) // need neighbor cluster of a1
							{
								for(SuperTarget s2: tempst.neighbors.values()) // need neighbor cluster  of a2?
								{
									
									if(s1.stid != s2.stid)
									{


										boolean arebothnei = areBothNei(s1,s2,tempst);

										if(s1.stid != tempst.stid && 
												s2.stid != tempst.stid  && arebothnei)
										{
											// measure the distance between a1->s1 and a2->s2

											//a1->s1 distance between a target and a supertarget
											System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);


											double d1 = shortestdist(a1,s1); // should i use bfs for longer path other than near neighbor?
											//TODO
											System.out.println("shortestdist(a1,s1) "+ d1);

											double d2 = shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2
											
											
											
											
											ArrayList<Integer> tmpspath = new ArrayList<Integer>();
											double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpspath);

											if(dista1a2 ==0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												System.out.println("No path found to compute AP for st "+ tempst.stid);
											}

											System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1 > 0 && d2 > 0 && dista1a2 > 0)
											{
												double totaldi = d1 + dista1a2 + d2;
												if(totaldi < mindi)
												{
													//stid1 = .stid;
													//stid2 = st2.stid;
													dminst = tempst;
													aid1 = a1.getTargetid();
													aid2 = a2.getTargetid();
													mindi = totaldi;
													sda1a2 = dista1a2;
													s1id = s1.stid;
													s2id = s2.stid;
													//spath.add(tmpspath.get(0));
													spath.clear();
													for(Integer in: tmpspath)
													{
														spath.add(in);
													}

													System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);


												}
											}

										}
									}
								}
							}
						}
					}
				}
			}

			//}
			//						}
			//
			//					}
			//				}
			//			}


			// System.out.println("Current mindi "+ mindi + " \n ST1 "+ stid1 + " ST2 "+ stid2 + 
			//	 "\n a1 "+ aid1 + ", a2 "+ aid2);
			/**
			 * Merge two cluster which have min d
			 * For the new id of supertargets concat the strings with comma
			 * remove the old two clusters
			 * add the new one. 
			 */
			
			
			//TODO what to do when there is no path from a1 to a2?????
			
			
			
			if(tempst.stid==14 && sda1a2==-1)
			{
				System.out.println("shortestdist(a1,s1) ");
				printSuperTarget(tempst);
			}

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				System.out.println("AP done for st "+ tempst.stid);
				
				dstravel.put(tempst.stid, sda1a2);
				stpaths.put(tempst.stid, spath);
				//SuperTarget newst = SuperTarget.mergeSuperTargets(sts.get(stid1), sts.get(stid2), aid1, aid2, targetmaps);
				// printSuperTarget(newst);
				/*sts.remove(stid1);
				sts.remove(stid2);*/
				//printSuperTargets(sts);
				
				updateAP(tempst, sts, aid1, aid2);
				//sts.put(newst.stid, newst);
				//update the neighbors of ST
				// System.out.println("\n\n After merging # supertargets : "+ sts.size());
				// System.out.println("\n After merging new supertargets : ");
				 
			}
			System.out.println("hi");
		}

		//}
		updateNeighbors(sts);
		printSuperTargets(sts);
		 return sts;
	}
	
	
	
	


	private static boolean areBothNei(SuperTarget s1, SuperTarget s2, SuperTarget tempst) {
		
		
		if(isNeighbor(s1, tempst) && isNeighbor(s2, tempst))
			return true;
		
		return false;
	}


	private static int findIndex(Instances newinstance, Integer tid) {
		// TODO Auto-generated method stub
		
		int index = 0;
		
		Iterator<Instance> ins = newinstance.iterator();
		
		while(ins.hasNext())
		{
			Instance newins = ins.next();
			if(newins.value(0) == tid.intValue())
				return index;
			index++;
		}
		
		
		return -1;
	}


	
	
	
	private static ArrayList<Integer>[] clusterWithWekaRW(int k, Instances newinstance, int totalcluster,
			ArrayList<Integer> targetstocluster, HashMap<Integer, TargetNode> targetmaps, EM dc) throws Exception {
		
		
		
		ArrayList<Integer>[] clusters = (ArrayList<Integer>[])new ArrayList[totalcluster];
		
		for(int i=0; i<totalcluster; i++)
		{
			clusters[i] = new ArrayList<Integer>();
		}
		
		
	/*	FileInputStream fstream = new FileInputStream("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/realdata2.csv");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		Instances instances = new Instances(br);*/
		 
		 
		 
		 // Print header and instances.
		// System.out.println("\nDataset:\n");
		 //System.out.println(instances.toSummaryString());
		 
		/* SimpleKMeans model = new SimpleKMeans();
		 model.setNumClusters(10);
		 model.buildClusterer(instances);
		 model.setDistanceFunction(new weka.core.ManhattanDistance());
		 
		 System.out.println(model);*/
		
		
		//EM dc = new EM();
		
		 
		//SimpleKMeans dc = new SimpleKMeans();
		
		 
		 if(newinstance.get(0).value(0) != 0)
		 {
			 throw new Exception("0 is not the base");
		 }
		 
		 newinstance.remove(0); // remove base
		 
		 dc.setNumClusters(k-1); // 0 for base
		
		 dc.buildClusterer(newinstance);
		 System.out.println(dc);
		 
		 
		 
		 for(int i=0; i<newinstance.size(); i++)
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(newinstance.get(i))+1);
			 int clusterid = dc.clusterInstance(newinstance.get(i));
			 int tid = (int)newinstance.get(i).value(0);
			 if(tid!=0)
			 {
				 clusters[clusterid+1].add(tid);
			 }
		 
			 
		 }
		 
		 clusters[0].add(0); //base
		 
		 printClusters(clusters);
		 
		 int j = k;
		 for(Integer t: targetmaps.keySet())
		 {
			 if(!targetstocluster.contains(t))
			 {
				 clusters[j++].add(t.intValue());
			 }
		 }
		 

		 printClusters(clusters);
		
		return clusters;
	}


	private static boolean canBeMerged(SuperTarget st1, ArrayList<Integer> targetstocluster) {
		
		//if(st1.nodes.size() ==1)
		
		
		if(st1.nodes.size()>1)
			return true;



		for(TargetNode t: st1.nodes.values())
		{
			//boolean f = false;
			for(Integer n: targetstocluster)
			{
				if(n.equals(t.getTargetid()))
				{
					//f = true;
					return true;
				}
			}
		}



		return false;
	}




	public static double minimumDist(SuperTarget st1, SuperTarget st2) {
		
		double mind = Double.MAX_VALUE;
		
		// min dist using ap
		
		for(TargetNode t1: st1.ap.values())
		{
			for(TargetNode t2: st2.ap.values())
			{
				if(t1.getNeighbors().contains(t2))
				{
					if(mind > t1.getDistance(t2))
					{
						mind = t1.getDistance(t2);
					}
				}
			}
		}
		
		
		
		
		return mind;
	}




	private static void updateNeighbors(SuperTarget newst, HashMap<Integer, SuperTarget> sts, int stid1, int stid2) {
		
		// remove old sts as neighbors
		
		for(SuperTarget st: sts.values())
		{
			if(st.neighbors.keySet().contains(stid1))
			{
				st.neighbors.remove(stid1);
			}
			if(st.neighbors.keySet().contains(stid2))
			{
				st.neighbors.remove(stid2);
			}
		}
		
		newst.neighbors.clear();
		
		//update new neighbor for every st
		for(SuperTarget st: sts.values())
		{
			if(isPotentialNeighbor(newst, st))
			{
				st.neighbors.put(newst.stid, newst);
				newst.neighbors.put(st.stid, st);
				
				
			}
		}
		
	}
	
	
private static void updateAP(SuperTarget curst, HashMap<Integer, SuperTarget> sts, int aid1, int aid2 ) {
		
		// remove old sts as neighbors
		
		
		
		//curst.neighbors.clear();
		curst.ap.clear();
		
		curst.ap.put(aid1, curst.nodes.get(aid1));
		curst.ap.put(aid2, curst.nodes.get(aid2));
		
		
		//update new neighbor for every st and curst
		for(SuperTarget st: sts.values())
		{
			if(curst.stid != st.stid)
			{
				if(!isPotentialNeighbor(curst, st))
				{
					st.neighbors.remove(curst.stid);
					st.distances.remove(curst);
					st.path.remove(curst);
					
					
					curst.neighbors.remove(st.stid);
					curst.distances.remove(st);
					curst.path.remove(st);
					
				}
			}
		}
		
	}


private static void updateNeighbors(HashMap<Integer, SuperTarget> sts) {
	
	// remove old sts as neighbors
	
	
	
	
	//update new neighbor for every st and curst
	for(SuperTarget st: sts.values())
	{


		for(SuperTarget st2: sts.values())
		{

			if(st2.stid != st.stid)
			{
				if(!isPotentialNeighbor(st2, st))
				{

					//should enter neighbor again?
					//st.neighbors.put(curst.stid, curst);
					//st2.neighbors.put(st.stid, st);
					
					st.neighbors.remove(st2.stid);
					st2.neighbors.remove(st.stid);


				}
				else
				{
					
				}
			}
		}
	}
	
}


	private static boolean isPotentialNeighbor(SuperTarget newst, SuperTarget st) {
		
		
		if(newst.stid == st.stid)
			return false;
		
		
		for(TargetNode nei: newst.ap.values())
		{
			for(TargetNode nei2 : st.ap.values())
			{
				if(nei.getNeighbors().contains(nei2))
					return true;
			}
		}
		
		return false;
	}


	private static boolean isNeighbor(SuperTarget st1, SuperTarget st2) {
		
		
		// if they have nodes in neighbor
		if(st1.neighbors.values().contains(st2))
			return true;
		
		return false;
	}


	private static double shortestdist(TargetNode a1, TargetNode a2, SuperTarget tempst, double dmax, ArrayList<Integer> spath) {
		
		
		
		int nTargets = tempst.nodes.size(); 
		TargetNode start = new TargetNode(a1);
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		ArrayList<TargetNode> closed = new ArrayList<TargetNode>();
		
		fringequeue.add(start);
		int pathcounter = 0;
		int nodestocover = tempst.nodes.size();
		while(fringequeue.size()>0)
		{
			System.out.println("Polling from queue ");
			System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==a2.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//System.out.println();
				//printPath(node);
				//pathcounter++;
				//System.out.println();
				
				goals.add(node);
				SecurityGameContraction.makeClusterPathSeq(goals, spath);
				return node.distancecoveredyet;
				//break;

				/*if(pathcounter>5000)
					break;*/
			}
			else if(node.distancecoveredyet<=dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				closed.add(node);
				ArrayList<TargetNode> succs = ExpandTarget(node, tempst.nodes, dmax);
				/**
				 * add nodes to queue
				 */
				for(TargetNode suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					if(!isInclosed(closed,node))
					{
						fringequeue.add(suc);
					}

				}
				
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}

		
		
		return 0;
	}
	
	
	
	public static boolean isInclosed(ArrayList<TargetNode> closed, TargetNode node) {

		
		for(TargetNode t: closed)
		{
			if(t.getTargetid() == node.getTargetid())
				return true;
		}
		return false;
		
	}


	private static ArrayList<TargetNode> ExpandTarget(TargetNode node, HashMap<Integer,TargetNode> nodes, double dmax ) 
	{
		ArrayList<TargetNode> successors = new ArrayList<TargetNode>();

		/**
		 * find the index for  node
		 */
		TargetNode tmpnode = nodes.get(node.getTargetid());
		
		for(TargetNode nei: tmpnode.getNeighbors())
		{
			if(nodes.values().contains(nei))
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


	
	
	
	
	
	
	private static double shortestdist(TargetNode a1, SuperTarget s1) {
		// TODO Auto-generated method stub



		double dmin = Double.MAX_VALUE;
		boolean isnei = false;
		for(TargetNode t: s1.ap.values())
		{
			// are each other's neighbor
			for(TargetNode nei : t.getNeighbors())
			{
				
					if(nei.getTargetid() == a1.getTargetid())
					{
						double tmpd = a1.getDistance(t);
						if(tmpd<dmin)
						{
							dmin = a1.getDistance(t);
						}
						isnei = true;
					}
				
			}
		}
		
		if(isnei)
			return dmin;
		else 
			return -1;
	}


	public static void testClustering()
	{
		
		//testing merge conflict



		// We will work tomorrow on the main alorithm

		int base = 0;
		int dest = 0;
		int k = 5;
		int radius = 3;
		int dmax = 100;
		int nRes=1;
		int dlim = 10;
		int nTargets = 20;
		int ITER = 1;
		int ap = 4; // should be </= than cluster size
		int utiliy_l=0;
		int utility_h=10;
		






		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		ArrayList<Integer>[] clus = GroupingTargets.makeGraph(k, radius, dlim , nTargets, utiliy_l, utility_h, ap, targets, targetmaps);

	//	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
	//	HashMap<Integer,TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		
		
		//ArrayList<Integer>[] clus = GroupingTargets.createGraph3(targets, targetmaps);
		printClusters(clus);
		printNodesWithNeighborsAndPath(targets);
		
		ArrayList<Integer> targetstocluster = new ArrayList<Integer>();
		
		/*for(int i=0; i<11; i++)
		{
			targetstocluster.add(i);
		}*/
		
		targetstocluster.add(0);
		
		targetstocluster.add(1);
		targetstocluster.add(4);
		targetstocluster.add(5);

		targetstocluster.add(6);
		targetstocluster.add(9);
		targetstocluster.add(10);
		
		HashMap<Integer, Double> d = new HashMap<>();
		HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
		
		HashMap<Integer, SuperTarget> sts = GroupingTargets.clusterTargets(targetstocluster, targets, targetmaps, dmax, k+1, radius, d, stpaths);
		printSuperTargets(sts);
		System.out.println("hii");
		
	}


	
	
	
	
	
	


	public static ArrayList<Integer>[] makeDummyCluster(int k)
	{
		ArrayList<Integer>[] cluster = (ArrayList<Integer>[])new ArrayList[k];
		for(int i=0; i<k; i++)
		{
			cluster[i] = new ArrayList<Integer>();
		}

		cluster[0].add(0);
		cluster[1].add(1);
		cluster[1].add(2);
		cluster[1].add(3);
		cluster[2].add(4);
		cluster[2].add(5);
		cluster[2].add(6);





		return cluster;
	}
	
	
	public static ArrayList<Integer>[] makeGraph(int k, int radius, int dlim, int nTargets, int utiliy_l,
			int utility_h, int ap, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps)
	{
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		// create targets


		targetmaps.clear();
		targets.clear();

		for(int i=0; i<nTargets; i++)
		{
			
			
			// select a random number for density
			int ru = randInt(utiliy_l, utility_h);
			if(i==0)
				ru=utility_h;

			TargetNode t = new TargetNode(i, ru);
			if(t.getTargetid()==0)
			{
				t.attackerreward = utility_h;
				t.attackerpenalty = 0;
				t.defenderpenalty = -utility_h;
				t.defenderreward = 0;
				t.setAnimaldensity(ru);
			}
			else
			{
				t.attackerreward = ru;
				t.attackerpenalty = 0;
				t.defenderpenalty = -ru;
				t.defenderreward = 0;
				t.setAnimaldensity(ru);
			}
			targetmaps.put(i, t);
			targets.add(t);
		}

		// select 0 as base, that means set it more than radius distance away and don't include it in grouping
		//select how many connection you want with base, a random number dependent on nTarget

		int nedges = randInt(nTargets/2, nTargets-1);

		// 

		int[] nodes = new int[nTargets-1]; 

		for(int i=0; i<nodes.length; i++)
		{
			nodes[i] = i+1;
		}

		shuffleArray(nodes);

		TargetNode base = targetmaps.get(0);
		for(int eid=0; eid<nedges; eid++)
		{
			// pick a target
			TargetNode t = targetmaps.get(nodes[eid]);

			ArrayList<TargetNode> path = new ArrayList<TargetNode>();

			double dist = randInt(radius+1, dlim);  // more distance than radius

			base.addNeighbor(t);
			base.addDistance(t, dist);
			base.setPath(t, path);

			t.addNeighbor(base);
			t.addDistance(base, dist);
			t.setPath(base, path);


		}

		ArrayList<Integer>[] cluster = (ArrayList<Integer>[])new ArrayList[k];
		for(int i=0; i<k; i++)
		{
			cluster[i] = new ArrayList<Integer>();
		}
		cluster[0].add(0); // base

		int curindex = 0;
		int tleft = nTargets-1;

		nodes = new int[nTargets-1];
		for(int i=0; i<nodes.length; i++)
		{
			nodes[i] = i+1;
		}
		shuffleArray(nodes);

		int nodepercluster = (nTargets-1)/(k);
		int remainder = (nTargets-1)%(k);
		if(remainder>0)
		{
			nodepercluster+=1;
		}

		for(int cid=1; cid<k; cid++)
		{
			// pick how many targets to include in cluster cid

			int count = 0;
			if(cid==k-1 && remainder>0)
				nodepercluster += (nTargets - ((k-1)*nodepercluster + 1));
			for(int j=curindex; j<(nTargets-1); j++)
			{
				cluster[cid].add(nodes[j]);
				count++;
				if(count==nodepercluster)
					break;
			}
			curindex += nodepercluster;


		}

		printClusters(cluster);
		for(int i=0; i<nTargets; i++)
		{
			
			for(ArrayList<Integer> clus: cluster)
			{
				int c = randInt(0, 1);
				int low = 0;
				int high = 0;
				if(c==0)
				{
					low = utiliy_l;
					high = utility_h/2;
				}
				else
				{
					low = utility_h/2;
					high = utility_h;
				}
				for(Integer n: clus)
				{
					
					if(n!=0)
					{
						int ru = randInt(low, high);
						TargetNode t = targetmaps.get(n);
						t.attackerreward = ru;
						t.attackerpenalty = 0;
						t.defenderpenalty = -ru;
						t.defenderreward = 0;
						t.setAnimaldensity(ru);
					}

				}
			}
			
			/*// select a random number for density
			int ru = randInt(utiliy_l, utility_h);
			if(i==0)
				ru=utility_h;

			TargetNode t = new TargetNode(i, ru);
			t.attackerreward = ru;
			t.attackerpenalty = 0;
			t.defenderpenalty = -ru;
			t.defenderreward = 0;
			t.setAnimaldensity(ru);
			targetmaps.put(i, t);
			targets.add(t);*/
		}
		
		
		
		
		
		
		// for each of k clusters
		//     select some random nodes to be within radius

		for(int clusterid = 1; clusterid<k; clusterid++)
		{
			int centerid = cluster[clusterid].get(0);
			TargetNode center = targetmaps.get(centerid);
			for(Integer nei: cluster[clusterid])
			{
				TargetNode neinode = targetmaps.get(nei);
				if(nei!=centerid)
				{
					ArrayList<TargetNode> path = new ArrayList<TargetNode>();

					double dist = randInt(1, radius);  // more distance than radius

					center.addNeighbor(neinode);
					center.addDistance(neinode, dist);
					center.setPath(neinode, path);

					neinode.addNeighbor(center);
					neinode.addDistance(center, dist);
					neinode.setPath(center, path);
				}
			}
		}




		// for each k group select some random nodes and connect them with more than radius and withim dlim




		for(int cid=1; cid<k-1; cid++)
		{
			for(int cid2=cid+1; cid2<k; cid2++)
			{
				for(int con=0; con<ap; con++)
				{
					//System.out.println("cid "+cid+ " cid2 "+ cid2+ " con "+con) ;
					
					TargetNode t1 = targetmaps.get(cluster[cid].get(con));
					TargetNode t2 = targetmaps.get(cluster[cid2].get(con));


					ArrayList<TargetNode> path = new ArrayList<TargetNode>();

					double dist = randInt(radius+1, dlim);  // more distance than radius

					t1.addNeighbor(t2);
					t1.addDistance(t2, dist);
					t1.setPath(t2, path);

					t2.addNeighbor(t1);
					t2.addDistance(t1, dist);
					t2.setPath(t1, path);
				}
			}
		}


		return cluster;

	}



	public static ArrayList<Integer>[] makeCluster(int k, int nTargets, int utility_l, 
			int utility_h, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps, double[][] density, int iter)
	{
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		// create targets


		//targetmaps.clear();
		//targets.clear();

		

		ArrayList<Integer>[] cluster = (ArrayList<Integer>[])new ArrayList[k];
		for(int i=0; i<k; i++)
		{
			cluster[i] = new ArrayList<Integer>();
		}
		cluster[0].add(0); // base

		int curindex = 0;
		int tleft = nTargets-1;

		int[] nodes = new int[nTargets-1];
		for(int i=0; i<nodes.length; i++)
		{
			nodes[i] = i+1;
		}
		shuffleArray(nodes);

		int nodepercluster = (nTargets-1)/(k);
		int remainder = (nTargets-1)%(k);
		if(remainder>0)
		{
			nodepercluster+=1;
		}

		
		ArrayList<Integer> done = new ArrayList<Integer>();
		done.add(0);
		
		
		int cid = 1;
		while(true)
		{
			// pick how many targets to include in cluster cid

			int count = 0;
			if(cid==k-1 && remainder>0)
				nodepercluster += (nTargets - ((k-1)*nodepercluster + 1));
			
			if(done.size()==targetmaps.size())
			{
				break;
			}
			while(true)
			{
				//pick a node to add in cluster
				
				
				
				int nodeid = pickNode(targetmaps, done, nodepercluster, cluster[cid]);
				if((nodeid == -1) || (done.size()==targetmaps.size()))
				{
					break;
				}
				cluster[cid].add(nodeid);
				if(nodeid==0)
				{
					System.out.println("fffff count  "+ count + ", nodeperclus "+ nodepercluster);

				}
				done.add(nodeid);
				
				count++;
				if(count == nodepercluster)
				{
					System.out.println("count  "+ count + ", nodeperclus "+ nodepercluster);
					break;
				}
			}
			//curindex += nodepercluster;
			if(done.size()<targetmaps.size())
			{
				cid ++;
				if(cid == k)
				{
					cid =1;
				}
			}

		}

		printClusters(cluster);
		
		//Random rand = new Random();
		
		for(int i=0; i<nTargets; i++)
		{
			
			for(ArrayList<Integer> clus: cluster)
			{
				if(clus.get(0) == 0)
				{
					density[iter][0] = utility_h;
					
				}
				else
				{
					int utility = randInt(utility_l, utility_h);
					for(Integer n: clus)
					{
						density[iter][n] = Math.abs(utility - randInt(0,2));
					}
				}
				
			}

		}

		return cluster;

	}
	
	
	
	public static ArrayList<Integer>[] makeClusterWithRange(int k, int nTargets, int utility_l, 
			int utility_h, ArrayList<TargetNode> targets, HashMap<Integer,TargetNode> targetmaps,
			double[][] density, int iter, int[][] ranges, int[] percforranges)
	{
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		// create targets


		//targetmaps.clear();
		//targets.clear();

		

		ArrayList<Integer>[] cluster = (ArrayList<Integer>[])new ArrayList[k];
		for(int i=0; i<k; i++)
		{
			cluster[i] = new ArrayList<Integer>();
		}
		cluster[0].add(0); // base

		int curindex = 0;
		int tleft = nTargets-1;

		int[] nodes = new int[nTargets-1];
		for(int i=0; i<nodes.length; i++)
		{
			nodes[i] = i+1;
		}
		shuffleArray(nodes);

		int nodepercluster = (nTargets-1)/(k);
		int remainder = (nTargets-1)%(k);
		if(remainder>0)
		{
			nodepercluster+=1;
		}

		
		ArrayList<Integer> done = new ArrayList<Integer>();
		done.add(0);
		
		
		int cid = 1;
		while(true)
		{
			// pick how many targets to include in cluster cid

			int count = 0;
			if(cid==k-1 && remainder>0)
				nodepercluster += (nTargets - ((k-1)*nodepercluster + 1));
			
			if(done.size()==targetmaps.size())
			{
				break;
			}
			while(true)
			{
				//pick a node to add in cluster
				
				
				
				int nodeid = pickNode(targetmaps, done, nodepercluster, cluster[cid]);
				if((nodeid == -1) || (done.size()==targetmaps.size()))
				{
					break;
				}
				cluster[cid].add(nodeid);
				if(nodeid==0)
				{
					System.out.println("fffff count  "+ count + ", nodeperclus "+ nodepercluster);

				}
				done.add(nodeid);
				
				count++;
				if(count == nodepercluster)
				{
					System.out.println("count  "+ count + ", nodeperclus "+ nodepercluster);
					break;
				}
			}
			//curindex += nodepercluster;
			if(done.size()<targetmaps.size())
			{
				cid ++;
				if(cid == k)
				{
					cid =1;
				}
			}

		}

		printClusters(cluster);
		
		//Random rand = new Random();
		
		
		int nlclus = (int)Math.floor(nTargets*(percforranges[0]/100.0))-1;
		int nmidclus = (int)Math.ceil(nTargets*(percforranges[1]/100.0));
		int nhighclus = nTargets - nlclus - nmidclus;
		int limit[] = {nlclus, nmidclus, nhighclus};
		
		System.out.println("low clus: "+ nlclus + ", mid clus "+ nmidclus + ", high clus "+ nhighclus);
		
		int[] counter = new int[3];
		
		for(int clusid=0; clusid<k; clusid++)
		{
			
			if(clusid==0)
			{
				density[iter][0] = ranges[2][1];
			}
			else
			{
				int cls = 1;
				
				while(true)
				{
					cls = randInt(0, 2);
					if(counter[cls]<limit[cls])
					{
						break;
					}
				}
				
				for(Integer n: cluster[clusid])
				{
					
					int utility = randInt(ranges[cls][0], ranges[cls][1]);
					density[iter][n] = utility;
					counter[cls]++;
					
				}
				
			}
			
			
			/*if(counter[2]<nhighclus)
			{
				counter[2]++;
				for(Integer n: cluster[clusid])
				{
					
					int utility = randInt(ranges[2][0], ranges[2][1]);
					density[iter][n] = utility;
					if(n==0)
					{
						density[iter][n] = ranges[2][1];
					}
					
					
				}
			}
			else if(counter[1]<nmidclus)
			{
				counter[1]++;
				for(Integer n: cluster[clusid])
				{
					int utility = randInt(ranges[1][0], ranges[1][1]);
					density[iter][n] = utility;
					if(n==0)
					{
						density[iter][n] = ranges[2][1];
					}
					
				}
			}
			else
			{
				counter[0]++;
				for(Integer n: cluster[clusid])
				{
					int utility = randInt(ranges[0][0], ranges[0][1]);
					density[iter][n] = utility;
					if(n==0)
					{
						density[iter][n] = ranges[2][1];
					}
					
				}
			}*/
		}

		
		
		
		
		
		
		/*for(int i=0; i<nTargets; i++)
		{
			
			for(ArrayList<Integer> clus: cluster)
			{
				if(clus.get(0) == 0)
				{
					density[iter][0] = utility_h;
					
				}
				else
				{
					int utility = randInt(utility_l, utility_h);
					for(Integer n: clus)
					{
						density[iter][n] = Math.abs(utility - randInt(0,2));
					}
				}
				
			}

		}*/

		return cluster;

	}
	
	
	

	private static int pickNode(HashMap<Integer, TargetNode> targetmaps, ArrayList<Integer> done, int nodepercluster,
			ArrayList<Integer> cluster) {
		
		
		int nodeid = -1;
		
		int[] nodes = new int[targetmaps.size()-done.size()];
		
		int index = 0;
		for(Integer i: targetmaps.keySet())
		{
			if(!done.contains((i)))
			{
				nodes[index++] = i;
				//System.out.println("Adding node "+ i);
				
			}
		}
		shuffleArray(nodes);

		
		
		
		// if cluster empty, choose any node
		if(cluster.size()==0)
		{
			int nodeidindex = randInt(0, nodes.length-1);
			nodeid = nodes[nodeidindex];
		}
		else
		{
			// pick a node which is not in cluster, and done
			//which is neighbor of any node in cluster
			//nodes = new int[targetmaps.size()-1-done.size()];
			
			ArrayList<Integer> potentialnode = new ArrayList<Integer>(); // neighbbor of cluster and not in cluster and done
			
			for(Integer clusnode: cluster)
			{
				for(TargetNode potnode: targetmaps.get(clusnode).getNeighbors())
				{
					// see if it's in done 
					if(!done.contains(potnode.getTargetid()) && (potnode.getTargetid() != 0))
					{
						potentialnode.add(potnode.getTargetid());
					}
				}
				
				
				
			}
			
			
			nodes = new int[potentialnode.size()];
			int ind = 0;
			for(int n: potentialnode)
			{
				nodes[ind++] = n;
			}
			
			shuffleArray(nodes);
			if(nodes.length==0)
			{
				System.out.println("ouch! found no node to add");
				return -1;
			}
			int nodeidindex = randInt(0, nodes.length-1);
			nodeid = nodes[nodeidindex];
			
			if(nodeid ==0)
			{
				System.out.println("ouch! found no node to add");
			}

		}
		
		
		
		return nodeid;
		
		
		
	}


	static void shuffleArray(int[] ar)
	{
		// If running on Java 6 or older, use `new Random()` on RHS here
		//Random rnd = ThreadLocalRandom.current();
		for (int i = ar.length - 1; i > 0; i--)
		{
			int index = randInt(0, i);
			// Simple swap
			int a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}
	}
	
	
	
	public static double[] groupingBaseline2NewMILP(int base, int dest, int k, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, 
			ArrayList<Integer>[] clusters, int dmaxsuper, int dminsuper) throws Exception
	{
		double[] res = new double[5];
		
		long revmaptime = 0;
		long solvingtime = 0;


		/*
		int base = 0;
		int dest = 0;
		int k = 5;
		int radius = 2;
		int dmax = 25;
		int nRes=1;
		int dlim = 10;
		int nTargets = 20;*/

		//HashMap<Integer, TargetNode>  ta = makeGraph(5, 5, 20, 20, 2, 10, 2);
		//printNodesWithNeighborsAndPath(ta);

		// create graph
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		/*for(TargetNode t: targets)
		{
			targetmaps.put(t.getTargetid(), t);
		}
		 */

		//makeGraph(k, radius, dlim , nTargets, 2, 10, 2, targets, targetmaps);

		printNodesWithNeighborsAndPath(targets);


		HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>(); // map for apsp
		HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>(); //mapback for apsp
		
		
		 //MultiValuedMap<String, String> map = new MultiValuedHashMap<String, String>();
		

		//fix it
		HashMap<String,Integer> subgamemap = new HashMap<String, Integer>(); // map for all the subgames
		HashMap<Integer,String> subgamemapback = new HashMap<Integer, String>(); // mapback for all the subgames
		HashMap<String,HashMap<ArrayList<Integer>,Double >> subgamepaths = new HashMap<String,HashMap<ArrayList<Integer>,Double>>();

		ArrayList<ArrayList<Double>> origpaths = new ArrayList<ArrayList<Double>>();
		
		HashMap<Integer, Double> superuncovattckr =new HashMap<Integer, Double>(); 
		
		
		// which subgames belongs to which super targets
		HashMap<String, Integer> subgamecontainer = new HashMap<String, Integer>();



		int[][] apspmat = buildAPSP(targets.size(), apspmap, apspmapback, targets, targetmaps);
		//ArrayList<Integer>[] clusters = groupTargets(base, dest, targetmaps, k, radius, targets, apspmat, apspmap);
		//printClusters(clusters);
		HashMap<Integer, SuperTarget> sts =  SuperTarget.buildSuperTargets2(clusters,targetmaps);  // sts are supre targets
		printSuperTargets(sts);

		int[][] gamedata = buildGamedata(targetmaps);

		// build tables for each of the supertargets

		double[][][][] sttable = new double[k][][][];



		// do i keep track of which targets were covered for each entry?
		// supertable contains all the subgame info's
		
		Date start = new Date();
		long l1 = start.getTime();

		//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
		
		//fix it duplicate key can happen
		HashMap<String, double[]>  supertable = buildSuperTable2(sttable, targets, targetmaps, sts, dmaxsuper, gamedata,
				apspmat, apspmap, subgamemap, subgamemapback, subgamepaths, dminsuper,subgamecontainer);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;
		solvingtime += diff;
		
		computeSuperAttckrUncovPayoffs(sts, superuncovattckr);
		
		

		printSubGamePathsPayoffs(subgamepaths, supertable);

		// generate pathsxxxxxxx
		try {
			ArrayList<ArrayList<String>> superpathseq = new ArrayList<ArrayList<String>>();
			ArrayList<SuperTarget> supergoals = generatePaths(dmax, sts, base, dest, supertable, dmaxsuper);
			makeSuperPathSeqSrcDest(superpathseq, supergoals);
			if(superpathseq.size()==0)
			{
				throw new Exception("No superpath");
			}
			printSuperPaths(superpathseq);


			int [][] p = new int[supertable.size()][]; // p matrix for subgame
			Integer[] input = new Integer[superpathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(superpathseq.size()==0)
			{

			}
			else
			{
				// make joint schedule
				// solve using an LP

				if(superpathseq.size()<nRes)
				{

					branch = new int[superpathseq.size()];
					jSet=combine(input, superpathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}
				List<ArrayList<Integer>> superjset = new ArrayList<ArrayList<Integer>>(jSet);
				//printJointSchedule(superjset);
				p = makeSuperPmat(superpathseq, superjset, subgamemap, supertable);

				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				
				start = new Date();
				l1 = start.getTime();

				
				double[] probdistribution = MIPSolver4.solveForSuperAttackerNewMILP(p, sts, supertable, 
						nRes, astrategy, subgamemap, subgamemapback, subgamecontainer, superuncovattckr );
				
				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;
				solvingtime += diff;
				
				
				
				start = new Date();
				l1 = start.getTime();

				
				

				String attackedsupertarget = findPlayedSubGameWMapping(p, probdistribution, subgamemap, subgamemapback, supertable);
				

				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;
				revmaptime += diff;

				System.out.println("attacked subgame  "+ attackedsupertarget);

				//printSubGamePath(subgamepaths,attackedsupertarget);

				//String sattackedsupertarget = String.valueOf(attackedsupertarget);// attackedsupertarget

				/*if(attackedsupertarget==0)
				{
					sattackedsupertarget += "0"+"0"; 
				}
				 */

				int attakedtarget = (int)supertable.get(attackedsupertarget)[2];

				System.out.println("attacked target  "+ attakedtarget);


				// don't know how will work
				//int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps);


				if(attakedtarget!=-1)
				{
					// Integer is the path index
					// Double is the path prob
					
					start = new Date();
					l1 = start.getTime();

					
					
					HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
					buildOrigPaths(probdistribution, superjset, superpathseq, subgamepaths, targetmaps, finalpaths);

					//TODO check the probability distribution
					printStrategyPaths(probdistribution, superjset, finalpaths);
					int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps, finalpaths);
					//double sprob = getSubGamePathProb(superjset.get(0), superpathseq, finalpaths);
					HashMap<Integer, Double>  attckrexppayoffs = computeAttackerExpectedPayoffs(targetmaps, probdistribution, origPmat,
							superjset, superpathseq, finalpaths);
					// find the [attacked target, def payoffs, attckr payoffs]
					double[] attackedtarget = findAttackedTarget(attckrexppayoffs, targetmaps, origPmat, probdistribution, superjset, superpathseq, finalpaths);


					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					revmaptime += diff;
					
					
					System.out.println("\nAtatcked target: "+ (attackedtarget[0])+ "\n Def Exp: "+ (attackedtarget[1]) + "\n Attckr Exp: "+ (attackedtarget[2]));

					res[0] = attackedtarget[1];
					res[1] = attackedtarget[2];
					//System.out.println("hiii ");

				}


				//System.out.println("hiii ");

			}




		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}




		res[3] = revmaptime;
		res[4] = solvingtime;
		return res;
	}
	
	
	public static double[] groupingBaseline2(int base, int dest, int k, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, 
			ArrayList<Integer>[] clusters, int dmaxsuper, int dminsuper) throws Exception
	{
		double[] res = new double[5];
		
		long revmaptime = 0;
		long solvingtime = 0;


		/*
		int base = 0;
		int dest = 0;
		int k = 5;
		int radius = 2;
		int dmax = 25;
		int nRes=1;
		int dlim = 10;
		int nTargets = 20;*/

		//HashMap<Integer, TargetNode>  ta = makeGraph(5, 5, 20, 20, 2, 10, 2);
		//printNodesWithNeighborsAndPath(ta);

		// create graph
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		/*for(TargetNode t: targets)
		{
			targetmaps.put(t.getTargetid(), t);
		}
		 */

		//makeGraph(k, radius, dlim , nTargets, 2, 10, 2, targets, targetmaps);

		printNodesWithNeighborsAndPath(targets);


		HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>(); // map for apsp
		HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>(); //mapback for apsp
		
		
		 //MultiValuedMap<String, String> map = new MultiValuedHashMap<String, String>();
		

		//fix it
		HashMap<String,Integer> subgamemap = new HashMap<String, Integer>(); // map for all the subgames
		HashMap<Integer,String> subgamemapback = new HashMap<Integer, String>(); // mapback for all the subgames
		HashMap<String,HashMap<ArrayList<Integer>,Double >> subgamepaths = new HashMap<String,HashMap<ArrayList<Integer>,Double>>();

		ArrayList<ArrayList<Double>> origpaths = new ArrayList<ArrayList<Double>>();
		
		HashMap<String, Integer> subgamecontainer = new HashMap<String, Integer>();


		int[][] apspmat = buildAPSP(targets.size(), apspmap, apspmapback, targets, targetmaps);
		//ArrayList<Integer>[] clusters = groupTargets(base, dest, targetmaps, k, radius, targets, apspmat, apspmap);
		//printClusters(clusters);
		HashMap<Integer, SuperTarget> sts =  SuperTarget.buildSuperTargets2(clusters,targetmaps);  // sts are supre targets
		printSuperTargets(sts);

		int[][] gamedata = buildGamedata(targetmaps);

		// build tables for each of the supertargets

		double[][][][] sttable = new double[k][][][];



		// do i keep track of which targets were covered for each entry?
		// supertable contains all the subgame info's
		
		Date start = new Date();
		long l1 = start.getTime();

		//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
		
		//fix it duplicate key can happen
		HashMap<String, double[]>  supertable = buildSuperTable2(sttable, targets, targetmaps, sts, dmaxsuper, gamedata,
				apspmat, apspmap, subgamemap, subgamemapback, subgamepaths, dminsuper, subgamecontainer);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;
		solvingtime += diff;
		
		

		printSubGamePathsPayoffs(subgamepaths, supertable);

		// generate pathsxxxxxxx
		try {
			ArrayList<ArrayList<String>> superpathseq = new ArrayList<ArrayList<String>>();
			ArrayList<SuperTarget> supergoals = generatePaths(dmax, sts, base, dest, supertable, dmaxsuper);
			makeSuperPathSeqSrcDest(superpathseq, supergoals);
			if(superpathseq.size()==0)
			{
				throw new Exception("No superpath");
			}
			printSuperPaths(superpathseq);


			int [][] p = new int[supertable.size()][]; // p matrix for subgame
			Integer[] input = new Integer[superpathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(superpathseq.size()==0)
			{

			}
			else
			{
				// make joint schedule
				// solve using an LP

				if(superpathseq.size()<nRes)
				{

					branch = new int[superpathseq.size()];
					jSet=combine(input, superpathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}
				List<ArrayList<Integer>> superjset = new ArrayList<ArrayList<Integer>>(jSet);
				//printJointSchedule(superjset);
				p = makeSuperPmat(superpathseq, superjset, subgamemap, supertable);

				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				
				start = new Date();
				l1 = start.getTime();

				
				double[] probdistribution = MIPSolver4.solveForSuperAttacker(p, sts, supertable, nRes, astrategy, subgamemap, subgamemapback );
				
				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;
				solvingtime += diff;
				
				
				
				start = new Date();
				l1 = start.getTime();

				
				

				String attackedsupertarget = findPlayedSubGameWMapping(p, probdistribution, subgamemap, subgamemapback, supertable);
				

				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;
				revmaptime += diff;

				System.out.println("attacked subgame  "+ attackedsupertarget);

				//printSubGamePath(subgamepaths,attackedsupertarget);

				//String sattackedsupertarget = String.valueOf(attackedsupertarget);// attackedsupertarget

				/*if(attackedsupertarget==0)
				{
					sattackedsupertarget += "0"+"0"; 
				}
				 */

				int attakedtarget = (int)supertable.get(attackedsupertarget)[2];

				System.out.println("attacked target  "+ attakedtarget);


				// don't know how will work
				//int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps);


				if(attakedtarget!=-1)
				{
					// Integer is the path index
					// Double is the path prob
					
					start = new Date();
					l1 = start.getTime();

					
					
					HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
					buildOrigPaths(probdistribution, superjset, superpathseq, subgamepaths, targetmaps, finalpaths);

					//TODO check the probability distribution
					printStrategyPaths(probdistribution, superjset, finalpaths);
					int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps, finalpaths);
					//double sprob = getSubGamePathProb(superjset.get(0), superpathseq, finalpaths);
					HashMap<Integer, Double>  attckrexppayoffs = computeAttackerExpectedPayoffs(targetmaps, probdistribution, origPmat,
							superjset, superpathseq, finalpaths);
					// find the [attacked target, def payoffs, attckr payoffs]
					double[] attackedtarget = findAttackedTarget(attckrexppayoffs, targetmaps, origPmat, probdistribution, superjset, superpathseq, finalpaths);


					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					revmaptime += diff;
					
					
					System.out.println("\nAtatcked target: "+ (attackedtarget[0])+ "\n Def Exp: "+ (attackedtarget[1]) + "\n Attckr Exp: "+ (attackedtarget[2]));

					res[0] = attackedtarget[1];
					res[1] = attackedtarget[2];
					//System.out.println("hiii ");

				}


				//System.out.println("hiii ");

			}




		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}




		res[3] = revmaptime;
		res[4] = solvingtime;
		return res;
	}


	
	
	public static double[] groupingBaseline3(int base, int dest, int k, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, 
			ArrayList<Integer>[] clusters, int dmaxsuper, int dminsuper) throws Exception
	{
		double[] res = new double[5];
		
		long revmaptime = 0;
		long solvingtime = 0;


		/*
		int base = 0;
		int dest = 0;
		int k = 5;
		int radius = 2;
		int dmax = 25;
		int nRes=1;
		int dlim = 10;
		int nTargets = 20;*/

		//HashMap<Integer, TargetNode>  ta = makeGraph(5, 5, 20, 20, 2, 10, 2);
		//printNodesWithNeighborsAndPath(targetmaps);

		// create graph
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		/*for(TargetNode t: targets)
		{
			targetmaps.put(t.getTargetid(), t);
		}
		 */

		//makeGraph(k, radius, dlim , nTargets, 2, 10, 2, targets, targetmaps);

		printNodesWithNeighborsAndPath(targets);


		HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>(); // map for apsp
		HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>(); //mapback for apsp
		
		
		ArrayList<String> doneimproveingsubgames = new ArrayList<String>();
		
		 //MultiValuedMap<String, String> map = new MultiValuedHashMap<String, String>();
		

		HashMap<String,Double> maxdistallocation = new HashMap<String, Double>(); // map for all the subgames
		
		
		//fix it
		HashMap<String,Integer> subgamemap = new HashMap<String, Integer>(); // map for all the subgames
		HashMap<Integer,String> subgamemapback = new HashMap<Integer, String>(); // mapback for all the subgames
		HashMap<String,HashMap<ArrayList<Integer>,Double >> subgamepaths = new HashMap<String,HashMap<ArrayList<Integer>,Double>>();
		
		
		HashMap<String, Integer> subgamecontainer = new HashMap<String, Integer>();
		
		HashMap<Integer, Double> superuncovattckr =new HashMap<Integer, Double>(); 

		ArrayList<ArrayList<Double>> origpaths = new ArrayList<ArrayList<Double>>();
		
		
		HashMap<Integer, HashMap<Integer, Integer>> stcountoutnei = new HashMap<Integer, HashMap<Integer, Integer>>();


		int[][] apspmat = buildAPSP(targets.size(), apspmap, apspmapback, targets, targetmaps);
		//ArrayList<Integer>[] clusters = groupTargets(base, dest, targetmaps, k, radius, targets, apspmat, apspmap);
		printClusters(clusters);
		HashMap<Integer, SuperTarget> sts =  SuperTarget.buildSuperTargets3(clusters,targetmaps, stcountoutnei);  // sts are supre targets
		printSuperTargets(sts);

		int[][] gamedata = buildGamedata(targetmaps);

		// build tables for each of the supertargets

		double[][][][] sttable = new double[k][][][];



		// do i keep track of which targets were covered for each entry?
		// supertable contains all the subgame info's
		
		Date start = new Date();
		long l1 = start.getTime();

		//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
		
		//fix it duplicate key can happen
		HashMap<String, double[]>  supertable = buildSuperTable3(sttable, targets, targetmaps, sts, dmaxsuper, gamedata,
				apspmat, apspmap, subgamemap, subgamemapback, subgamepaths, dminsuper, maxdistallocation, subgamecontainer);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;
		solvingtime += diff;
		
		computeSuperAttckrUncovPayoffs(sts, superuncovattckr);
		
		
		
		/**
		 * start a loop
		 * loop untill no distance to increase
		 */
		
		int loopiter = 0;
		
		while(true)
		{

			loopiter++;
			if(loopiter==10)
				break;


			printSubGamePathsPayoffs(subgamepaths, supertable);

			// generate pathsxxxxxxx
			try 
			{
				ArrayList<ArrayList<String>> superpathseq = new ArrayList<ArrayList<String>>();
				ArrayList<SuperTarget> supergoals = generatePaths(dmax, sts, base, dest, supertable, dmaxsuper);
				makeSuperPathSeqSrcDest(superpathseq, supergoals);
				if(superpathseq.size()==0)
				{
					throw new Exception("No superpath");
				}
				printSuperPaths(superpathseq);


				int [][] p = new int[supertable.size()][]; // p matrix for subgame
				Integer[] input = new Integer[superpathseq.size()];
				int[] branch = new int[nRes];//{0,0};//new char[k];

				for(int i=0; i<input.length; i++)
				{
					input[i] = i;
				}
				HashSet jSet=new HashSet();


				if(superpathseq.size()==0)
				{

				}
				else
				{
					// make joint schedule
					// solve using an LP

					if(superpathseq.size()<nRes)
					{

						branch = new int[superpathseq.size()];
						jSet=combine(input, superpathseq.size(), 0, branch, 0, jSet);
					}
					else
					{
						jSet=combine(input, nRes, 0, branch, 0, jSet);
					}
					List<ArrayList<Integer>> superjset = new ArrayList<ArrayList<Integer>>(jSet);
					ArrayList<Integer> first = superjset.get(0);
					ArrayList<Integer> last = superjset.get(superjset.size()-1);
					
					
					/*if(superjset.size()==2)
					{
						superjset.remove(0);
						superjset.remove(0);
						superjset.add(0, last);
						superjset.add(superjset.size(), first);
					}
					else if(superjset.size()>2)
					{
						//superjset.remove(0);
						superjset.remove(superjset.size()-1);
						superjset.add(0, last);
						//superjset.add(superjset.size(), first);
					}
					*/
					
					
					printJointSchedule(superjset);
					p = makeSuperPmat(superpathseq, superjset, subgamemap, supertable);

					HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();

					start = new Date();
					l1 = start.getTime();


					double[] probdistribution = MIPSolver4.solveForSuperAttackerNewMILP(p, sts, supertable, nRes, astrategy, 
							subgamemap, subgamemapback, subgamecontainer, superuncovattckr);
					
					/*if(loopiter==2)
					{
						probdistribution[3] = 1;
						probdistribution[0] = 0;

					}*/

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					solvingtime += diff;



					start = new Date();
					l1 = start.getTime();



					// need to find attacked supertarget
					Integer attackedsupertarget = findAttackedSuperTargetWMapping(p, probdistribution, subgamemap,
							subgamemapback, supertable, subgamecontainer, sts);


					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					revmaptime += diff;

					System.out.println("attacked supertarget  "+ attackedsupertarget);

					//printSubGamePath(subgamepaths,attackedsupertarget);

					//String sattackedsupertarget = String.valueOf(attackedsupertarget);// attackedsupertarget

					/*if(attackedsupertarget==0)
				{
					sattackedsupertarget += "0"+"0"; 
				}
					 */

					//int attakedsupertarget = (int)supertable.get(attackedsupertarget)[2];

					//System.out.println("attacked target  "+ attakedsupertarget);


					// don't know how will work
					//int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps);


					if(attackedsupertarget!=-1)
					{
						// Integer is the path index
						// Double is the path prob

						start = new Date();
						l1 = start.getTime();



						HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
						buildOrigPaths(probdistribution, superjset, superpathseq, subgamepaths, targetmaps, finalpaths);

						//TODO check the probability distribution
						printStrategyPaths(probdistribution, superjset, finalpaths);
						int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps, finalpaths);
						//double sprob = getSubGamePathProb(superjset.get(0), superpathseq, finalpaths);
						HashMap<Integer, Double>  attckrexppayoffs = computeAttackerExpectedPayoffs(targetmaps, probdistribution, origPmat,
								superjset, superpathseq, finalpaths);
						// find the [attacked target, def payoffs, attckr payoffs]
						double[] result = findAttackedTarget(attckrexppayoffs, targetmaps, origPmat, probdistribution, superjset, superpathseq, finalpaths);


						stop = new Date();
						l2 = stop.getTime();
						diff = l2 - l1;
						revmaptime += diff;


						System.out.println("\nAtatcked target: "+ (result[0])+ "\n Def Exp: "+ (result[1]) + "\n Attckr Exp: "+ (result[2]));

						res[0] = result[1];
						res[1] = result[2];
						//System.out.println("hiii ");
						
						/**
						 * get the maximum distance allocation for the attacked supertarget
						 * 
						 */
						// choose a subgame and then improve the supertarget based on the last subgame distance allocation
						// how to make a better choice ???? 
						String chosensubgame = findAttackedSubgame(subgamecontainer, attackedsupertarget, dmaxsuper, maxdistallocation, doneimproveingsubgames);
						if(!chosensubgame.equals(""))
						{
							System.out.println("Choosing subgame "+ chosensubgame + " to improve");
							
							String[] keys = chosensubgame.split(",");
							
							String key = keys[0] + "," + keys[1];
							
							int attackedentry = Integer.parseInt(keys[0]);
							int attackedexit = Integer.parseInt(keys[1]);
							double attackedtargetdistallctn = maxdistallocation.get(key);
							int stid = getStId(sts, attackedentry, attackedexit);
							HashMap<ArrayList<Integer>, Double> attackedtargetpaths = subgamepaths.get(chosensubgame);
							
							
							if(attackedtargetdistallctn==dmaxsuper)
							{
								break;
								// we can add entry points later
								// or exit points
							}
							
							for(int dist = (int)(attackedtargetdistallctn + 1); dist<= dmaxsuper; dist++)
							{
								
								

								ArrayList<TargetNode> subgraph = buildSubGraph(sts.get(stid).nodes);
								HashMap<Integer, TargetNode> subtargetmaps = buildSubTargetMaps(subgraph);
							
								ArrayList<Integer> newpath = getNewPath(subgraph,subtargetmaps ,dist, attackedentry, attackedexit);
								
								// check if the path adds any new target
								
								boolean addsnewtarget = doesAddNewTarget(attackedtargetpaths, newpath);
								// if so add the path to subgamepath and add the subgame to supertable
								// otherwise increase distance allocation
								
								if(addsnewtarget)
								{
									// subgamemap mapback
									
									maxdistallocation.put(key, (double)dist);
									
									int count = Collections.max(subgamemap.values());
									
									String newkey = attackedentry +","+ attackedexit+ ","+ dist;
									
									subgamemap.put(newkey, (count+1));
									subgamemapback.put(count+1, newkey);
									
									// subgamepaths
									
									HashMap<ArrayList<Integer>, Double> newp = new HashMap<ArrayList<Integer>, Double>();
									newp.put(newpath, 1.0);
									subgamepaths.put(newkey, newp);
									
									// supertable 
									// Should solve the game using a solver....
									
									
									int maxtarget = -1;
									double maxattackerpayoff = Double.NEGATIVE_INFINITY;
									for(TargetNode t: subtargetmaps.values())
									{
										// if  covered then 
										double tmppayoff = -1;
										if(newpath.contains(t.getTargetid()))
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
									

									
									supertable.put(newkey, new double[]{-maxattackerpayoff, maxattackerpayoff, maxtarget}) ;
									subgamecontainer.put(newkey, stid);
									System.out.println("Adding new subgame "+ newkey);
									break;
									
									
								}
								else if(dist==dmaxsuper && !addsnewtarget)
								{
									doneimproveingsubgames.add(key);
									System.out.println("Done improving subgame "+ key);
								}
							}
							
							//throw new Exception("No subgame to improve ...");
						}
						else
						{
							if(attackedsupertarget!=0)
							{
								printPossibleAP(stcountoutnei.get(attackedsupertarget));
							}
							System.out.println("hi");
							//TODO
							//add entry and exit points 
							// allocate the shortest distance  , NOTE: if entry and exit are same then allocate distance so that another target can be included. 
							// make another entry in the supertable
							// add the entry and/or exit point in the sts
							// the improve the subgame incrementally
						}

						
						
						
						
						
						
						
						
						/**
						 * 
						 * with that entry and exit
						 * try to increase distance allocation for same entry and exit point
						 * 
						 * if adds the same path then increase distance more. 
						 */
						
						
						

					}


					//System.out.println("hiii ");

				}




			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}


		res[3] = revmaptime;
		res[4] = solvingtime;
		return res;
	}


	
	
	
	public static double[] groupingWithDO(int base, int dest, int ncluster, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, 
			ArrayList<Integer>[] clusters) throws Exception
	{
		
		
		
		
		
		
		
		
		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		//SecurityGameContraction.printSortedTargets(targetssorted);
		
		
		//Get the list of initial targets using GCR from Tsrt, Tcur = GreedyCoverR()
		
		ArrayList<Integer> targetstocluster = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes);
		
		int currentPlace = targetstocluster.size()-1;
		
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		//HashMap<Integer, TargetNode> tmptargetmaps = new HashMap<Integer, TargetNode>();
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

		


		long clusteringtime=0;
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

			for(int i=0; i<targetstocluster.size(); i++)
			{
				System.out.print(targetstocluster.get(i)+",");
			}


			/*tmpgraph = SecurityGameContraction.getDuplicateGraph(targets);
			for(TargetNode t: tmpgraph)
			{
				tmptargetmaps.put(t.getTargetid(), t);
			}*/
			
			
			
			
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			//Build an abstract graph Gt through target clustering given Tcur, Gt
			
			
			HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
			HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
			
			HashMap<Integer, SuperTarget> currentst = GroupingTargets.clusterTargets(targetstocluster, targets, 
					targetmaps, dmax, ncluster+1, radius, dstravel, stpaths);
			targetsize= currentst.size();
			

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			clusteringtime += diff;
			
			
			//TODO save distance traveled for each cluster
			// remove the unncessary ones. 
			
			ArrayList<Integer> notin = new ArrayList<Integer>();
			
			for(SuperTarget st: currentst.values())
			{
				if(st.nodes.size()==1)
				{
					dstravel.put(st.stid, 0.0);
					
				}
				
			}
			
			
			int ind = 0;
			for(Integer t: dstravel.keySet())
			{
				if(!currentst.keySet().contains(t))
				{
					notin.add(t);
				}
				ind++;
			}
			
			for(Integer x: notin)
			{
				dstravel.remove(x);
				stpaths.remove(x);
			}
			
			
			//HashMap<Integer, Double> stvalue
			assignSTValues(currentst, targetmaps);
			
			//System.out.println("olaa ");
			//printSuperTargets(currentst);
			
			
			//Generate initial set of paths using GreedyPathR, scur = GPR(Gt)
			
			

			

			System.out.println("current st size "+ currentst.size());
			
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[currentst.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, currentst, targetmaps, nRes, dstravel);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: currentst.values())
			{

				map.put(st.stid, icount);
				System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			SecurityGameContraction.printPaths(pathseq);
			
			System.out.println("hi");

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

				

				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ masteritr);
				}


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
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: currentst.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
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
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, currentst, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, currentst, targetmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, currentst, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, currentst, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = SecurityGameContraction.makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							targetmaps, currentst, stpaths);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, targetmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, targetmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
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

					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSuperGreedyCoverMultRes2(targetmaps, 
							dmax, currentst.size(), 0, nRes, attackerstrategy, currentst, dstravel);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/*System.out.println("newpathseq size before purify : "+newpathseq.size());
					    //newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
						System.out.println("newpathseq size after purify : "+newpathseq.size());*/
						
						
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
						SecurityGameContraction.printPaths(newpathseq);

						System.out.println("Old path seq size "+ pathseq.size());

						int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

						pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						int newsize = pathseq.size();
						//System.out.println("haa ");


						if((oldsize==newsize) || (itr>=10))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ###############");
							break;
						}

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				break;
			}
			
			double ulimit = targetmaps.get(attackedtarget).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;
			int ADD_C = 5;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					targetstocluster.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=ADD_C)
					{
						break;
					}
				}
			}
			
			System.out.println("addcount : "+ addcount);

			currentPlace = targetstocluster.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<ADD_C || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += ADD_C-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					targetstocluster.add(targetssorted[k][0]);
				}
			}

			masteritr++;


		} // outer while loop

		System.out.println("Final target list size : "+ targetstocluster.size());

		for(int i=0; i<targetstocluster.size(); i++)
		{
			System.out.print(targetstocluster.get(i)+",");
		}

		
		double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, targetmaps, probdistribution);




		

		double[] res1 = {defpayoff, clusteringtime, solvingtime, targetsize, attackeru, slavetime, revmaptime};
		return res1;
	}
	
	

	
	public static double[] wekaClusteringWithDORW(int base, int dest, int ncluster, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps) throws Exception
	{
		
		
		
		
		
		
		
		
		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		//SecurityGameContraction.printSortedTargets(targetssorted);
		
		
		//Get the list of initial targets using GCR from Tsrt, Tcur = GreedyCoverR()
		
		ArrayList<Integer> targetstocluster = SecurityGameContraction.buildGreedyCoverMultRes(targets, dmax, nTargets, 0, nRes);
		
		int currentPlace = targetstocluster.size()-1;
		
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		//ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		//HashMap<Integer, TargetNode> tmptargetmaps = new HashMap<Integer, TargetNode>();
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

		


		long clusteringtime=0;
		long solvingtime=0;
		long revmaptime=0;
		int targetsize=0;
		double slavetime = 0;
		int [][] origpmat = new int[nTargets][];
		boolean canaddpath = true;
		
		int masteritr=0;
		
		
		
		
		/**
		 * make instances for clustering with weka
		 */
		
		
		FileInputStream fstream = new FileInputStream("C:\\Users\\IASRLusers\\eclipse-workspace\\IntervalSGAbstraction\\newdata3.arff");
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		// Read all the instances in the file (ARFF, CSV, XRFF, ...)
		 //DataSource source = new DataSource(br);
		 Instances instances = new Instances(br);
		 // Print header and instances.
		// System.out.println("\nDataset:\n");
		 System.out.println(instances.toSummaryString());
		 
		/* SimpleKMeans model = new SimpleKMeans();
		 model.setNumClusters(10);
		 model.buildClusterer(instances);
		 model.setDistanceFunction(new weka.core.ManhattanDistance());
		 System.out.println(model);*/
		 //MakeDensityBasedClusterer dc = new MakeDensityBasedClusterer();
		 EM dc = new EM();
		// instances.remove(0); // remove the base
		 dc.setNumClusters(ncluster-1);
		 dc.buildClusterer(instances);
		 System.out.println(dc);
		 
		 for(int i=0; i<instances.size(); i++) //without base
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(instances.get(i)));
		  
		 }
		
		

		
		while(true)
		{
			
			System.out.println("Outer loop...Master");

			pathseq = new ArrayList<ArrayList<Integer>>();

			System.out.println("\nCurrent place : "+ currentPlace);

			System.out.print("Current target list : ");

			for(int i=0; i<targetstocluster.size(); i++)
			{
				System.out.print(targetstocluster.get(i)+",");
			}


			/*tmpgraph = SecurityGameContraction.getDuplicateGraph(targets);
			for(TargetNode t: tmpgraph)
			{
				tmptargetmaps.put(t.getTargetid(), t);
			}*/
			
			
			
			
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			//Build an abstract graph Gt through target clustering given Tcur, Gt
			
			
			HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
			HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
			
			HashMap<Integer, SuperTarget> currentst = GroupingTargets.clusterTargetsWekaRW(targetstocluster, targets, 
					targetmaps, dmax, ncluster, radius, dstravel, stpaths, dc, instances);
			targetsize= currentst.size();
			
			
			printSuperTargets(currentst);

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;

			clusteringtime += diff;
			
			
			//TODO save distance traveled for each cluster
			// remove the unncessary ones. 
			
			ArrayList<Integer> notin = new ArrayList<Integer>();
			
			for(SuperTarget st: currentst.values())
			{
				if(st.nodes.size()==1)
				{
					dstravel.put(st.stid, 0.0);
					
				}
				
			}
			
			
			int ind = 0;
			for(Integer t: dstravel.keySet())
			{
				if(!currentst.keySet().contains(t))
				{
					notin.add(t);
				}
				ind++;
			}
			
			for(Integer x: notin)
			{
				dstravel.remove(x);
				stpaths.remove(x);
			}
			
			
			//HashMap<Integer, Double> stvalue
			assignSTValues(currentst, targetmaps);
			
			//System.out.println("olaa ");
			//printSuperTargets(currentst);
			
			
			//Generate initial set of paths using GreedyPathR, scur = GPR(Gt)
			
			

			

			System.out.println("current st size "+ currentst.size());
			
			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);
			p = new int[currentst.size()][]; // p matrix

			//apply greedy approach
			//TODO generate paths where there will be at least one target
			//ArrayList<TargetNode> goals = generatePathsGreedy2(dmax, gamedata, tmpgraph, currenttargets, nRes);
			//pathseq =  buildGreedyPathMultRes2(tmpgraph, dmax, tmpgraph.size(), 0, nRes);
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, currentst, targetmaps, nRes, dstravel);
			map = new HashMap<Integer, Integer>();
			mapback = new HashMap<Integer, Integer>();
			int icount = 0;
			for(SuperTarget st: currentst.values())
			{

				map.put(st.stid, icount);
				System.out.println("SuperTarget "+ st.stid +" --> "+icount);
				mapback.put(icount, st.stid);
				icount++;

			}
			//makePathSeq(pathseq, goals, goals.size(), tmpgraph.size(), map, mapback, tmpgraph);
			//printPaths(pathseq);
			System.out.println("Total path with duplicates "+pathseq.size());
			pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
			System.out.println("Total path without duplicates "+pathseq.size()+"\n");




			SecurityGameContraction.printPaths(pathseq);
			
			System.out.println("hi");

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

				

				if(pathseq.size()==0)
				{
					System.out.println("pathseq 0, iter.ohhh"+ masteritr);
				}


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
					/*for(int i=0; i<domindatednodes.size(); i++)
					{
						tmpgraph.add(domindatednodes.get(i));
					}*/
					
					for(SuperTarget st: currentst.values())
					{
						for(TargetNode x: st.nodes.values())
						{
							if(x.attackerreward>mAxpayoff)
							{
								mAxpayoff= x.attackerreward;
								defpayoff = x.defenderpenalty;
							}
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
					 *//*
					*//**
					 * pmat, where columns will be combination of paths. 
					 * rows are targets. 
					 * each entry will say whether the target is in the joint schedule
					 */
					//jSet.

					//printJointSchedule(jset);

					p = SecurityGameContraction.makeSuperPmat(pathseq, jset, mapback, currentst, map);
					//printPathMat(p);

					start = new Date();
					l1 = start.getTime();

					HashMap<Integer, Double> attackerstrategy = new HashMap<Integer, Double>();

					System.out.println("Solving LP");
					probdistribution = MIPSolver4.solveForAttackerLPST(p, currentst, targetmaps, nRes, attackerstrategy);

					

					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					solvingtime += diff;
					
					
					start = new Date();
					l1 = start.getTime();


					attackedtarget = SecurityGameContraction.findAttackSuperTargetWMapping(p, probdistribution, currentst, map, mapback);
					attackedtarget = mapback.get(attackedtarget);
					System.out.println("attack target before rev map "+ attackedtarget);
					
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackeru = SecurityGameContraction.expectedAttackerSTPayoff(attackedtarget, p, probdistribution, currentst, map);
					//System.out.println("attacker u= "+attackeru);

					//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);

					origpmat = SecurityGameContraction.makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							targetmaps, currentst, stpaths);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, targetmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, targetmaps, probdistribution);
					//System.out.println("attacker v= "+attackerv);
					
					System.out.println("master "+masteritr+", slave "+itr+", u= "+attackeru+", v= "+attackerv);
					System.out.println("attack target after rev map"+ attackedtarget);
					
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					revmaptime += diff;

					

					if(probdistribution.equals(null))
					{
						throw new Exception("Prob null...");
					}

					/*if(attackeru>=targetssorted[currentPlace+1][1] || currentPlace==targetssorted.length)
					{
						System.out.println("attacker u "+ attackeru +" is greater than u("+targetssorted[currentPlace+1][0]+")="+targetssorted[currentPlace+1][1]);

						break;
					}

					if(attackeru>= targetssorted[currentPlace+1][1])
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

					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSuperGreedyCoverMultRes2(targetmaps, 
							dmax, currentst.size(), 0, nRes, attackerstrategy, currentst, dstravel);
					
					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;

					slavetime += diff;
					
					
					/*System.out.println("newpathseq size before purify : "+newpathseq.size());
					    //newpathseq = SecurityGameContraction.determineNewPaths(newpathseq, p, probdistribution);
						System.out.println("newpathseq size after purify : "+newpathseq.size());*/
						
						
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
						SecurityGameContraction.printPaths(newpathseq);

						System.out.println("Old path seq size "+ pathseq.size());

						int oldsize = pathseq.size();
						for(ArrayList<Integer> q: newpathseq)
						{
							pathseq.add(q);
						}

						System.out.println("new paths added by slave *************, attacked target "+ attackedtarget);

						pathseq = SecurityGameContraction.removeDuplicatePathSimple(pathseq);
						System.out.println("New path seq size "+ pathseq.size());
						//printPaths(pathseq);
						int newsize = pathseq.size();
						//System.out.println("haa ");


						if((oldsize==newsize) || (itr>=20))
						{
							canaddpath = false;
							System.out.println("Slave can't add any new path ############### or iteration>20");
							printSuperTargets(currentst);
							SecurityGameContraction.printPaths(pathseq);
							break;
						}

						//SecurityGameContraction.printPaths(pathseq);

				} // end if else
				System.out.println("iter"+ itr);
				
			} // inner while loop 




			// add all targets all targets with utility >= U(a')


			if((currentPlace==targetssorted.length-1 || (attackeru>= attackerv)) && !canaddpath)
			{
				System.out.println("outer loop ....breaking.@@@@@@@@@@@@@@@..attacker u>=v="+attackeru);
				printSuperTargets(currentst);
				break;
			}
			
			double ulimit = targetmaps.get(attackedtarget).attackerreward;

			System.out.println("attacked target "+ attackedtarget+", adding all target w u >= "+ ulimit);


			int addcount=0;
			int ADD_C = 5;

			for(int k=currentPlace+1; k<targetssorted.length; k++)
			{
				if(targetssorted[k][1]>=ulimit)
				{
					addcount++;
					
					targetstocluster.add(targetssorted[k][0]);
					System.out.println("adding target "+targetssorted[k][0] +", u = "+ targetssorted[k][1]);
					if(addcount>=ADD_C)
					{
						break;
					}
				}
			}
			
			System.out.println("addcount : "+ addcount);

			currentPlace = targetstocluster.size()-1;

			System.out.println("currentplace  : "+ currentPlace);

			if(addcount<ADD_C || addcount==0)
			{
				//System.out.println("adding more ");

				int prevcur = currentPlace;
				currentPlace += ADD_C-addcount;

				//System.out.println("currentplace  : "+ currentPlace);
				if(currentPlace>targetssorted.length-1)
				{
					currentPlace = targetssorted.length-1;
				}
				//System.out.println("attacker u "+ attackeru +" is less than u("+targetssorted[currentPlace][0]+")="+targetssorted[currentPlace][1]);

				for(int k= prevcur+1; k<=currentPlace; k++ )
				{

					System.out.println("adding target  "+ targetssorted[k][0]);
					targetstocluster.add(targetssorted[k][0]);
				}
			}

			masteritr++;
			
		


		} // outer while loop

		System.out.println("Final target list size : "+ targetstocluster.size());

		for(int i=0; i<targetstocluster.size(); i++)
		{
			System.out.print(targetstocluster.get(i)+",");
		}

		
		double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, targetmaps, probdistribution);




		

		double[] res1 = {defpayoff, clusteringtime, solvingtime, targetsize, attackeru, slavetime, revmaptime};
		return res1;
	}
	
	
	public static ArrayList<Integer> buildGreedyCoverMultResWeka(
			ArrayList<TargetNode> targets, double dmax, int nTargets, int base, int nRes,
			HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			int[][] apspmat, AllPairShortestPath apsp ) {



		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];

		/**
		 * make mapping
		 */

		// map = new HashMap<Integer, Integer>();
		// mapback = new HashMap<Integer, Integer>();
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


		
		 int[][] apspmat1 =  apsp.allPairShortestPath(adjacencymatrix);
		
		
		SecurityGameContraction.purifyAPSPMatrixZeroGT(apspmat1, targets, nTargets, map, mapback);
		
		//apspmat = new int[apspmat1.length][apspmat1[0].length];
		
		
		
		for(int i=0; i<apspmat1.length; i++)
		{
			for(int j=0; j<apspmat1[0].length; j++)
			{
				apspmat[i][j] = apspmat1[i][j];
			}
		}



		ArrayList<Integer> tcur = new ArrayList<Integer>(); //greedyFirstRoute(dmax,gamedata, targets);

		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);

		tcur =  SecurityGameContraction.greedyCoverMultRes(base, targets, dmax, targetssorted, apspmat1, map,mapback, nRes);

		return tcur;
	}

	
	


	private static void assignSTValues(HashMap<Integer, SuperTarget> currentst,
			HashMap<Integer, TargetNode> targetmaps) {

		
		
		for(SuperTarget st: currentst.values())
		{
			double maxval = Double.MIN_VALUE;
			for(TargetNode t: st.nodes.values())
			{
				if(t.attackerreward>maxval)
				{
					maxval = t.attackerreward;
					st.attackerreward = t.attackerreward;
					st.attackerpenalty = 0;
					
					st.defenderreward = 0;
					st.defenderpenalty = -t.attackerreward;
				}
			}
		}
		
	}




	private static void printPossibleAP(HashMap<Integer, Integer> hashMap) {
		// TODO Auto-generated method stub
		for(Integer n: hashMap.keySet())
		{
			System.out.println("Node " + n + ", #connection "+ hashMap.get(n));
		}
		
	}


	private static String findAttackedSubgame(HashMap<String, Integer> subgamecontainer, Integer attackedsupertarget,
			int dmaxsuper, HashMap<String,Double> maxdistallocation, ArrayList<String> doneimproveingsubgames) {
		
		
		Integer maxd = -1;
		String sub = "";
		for(String subgame: subgamecontainer.keySet())
		{
			if(subgamecontainer.get(subgame)==attackedsupertarget)
			{
				String [] keys = subgame.split(",");
				String k = keys[0]+","+keys[1];
				if((maxd< Integer.parseInt(keys[2]))  && (Integer.parseInt(keys[2])<dmaxsuper) && !doneimproveingsubgames.contains(k))
				{
					maxd = Integer.parseInt(keys[2]);
					sub = subgame;
				}
				
			}
		}
		return sub;
		
	}


	private static void computeSuperAttckrUncovPayoffs(HashMap<Integer, SuperTarget> sts,
			HashMap<Integer, Double> superuncovattckr) {
		// TODO Auto-generated method stub
		for(SuperTarget st: sts.values())
		{
			int maxt = -1;
			Double maxu = Double.NEGATIVE_INFINITY;
			for(TargetNode n: st.nodes.values())
			{
				if(n.attackerreward>maxu)
				{
					maxu = n.attackerreward;
				}
			}
			superuncovattckr.put(st.stid, maxu);
		}
		
	}


	private static boolean doesAddNewTarget(HashMap<ArrayList<Integer>, Double> attackedtargetpaths,
			ArrayList<Integer> newpath)
	{


		for(int n: newpath)
		{

			boolean issame = true;

			for(ArrayList<Integer> p: attackedtargetpaths.keySet())
			{
				issame = true;
				if(!p.contains(n))
				{
					issame = false;
				}

			}
			if(!issame)
			{
				return true;
			}
		}





		return false;
	}


	private static ArrayList<Integer> getNewPath(ArrayList<TargetNode> subgraph,HashMap<Integer, TargetNode> subtargetmaps ,int dist, int attackedentry,
			int attackedexit) {
		
		
		
		
		ArrayList<Integer> newpath = SecurityGameContraction.buildOneGreedyPathWithSrcDest(subgraph, dist, subgraph.size(), 
				attackedentry, attackedexit, 1);
		
		
		
		
		
		
		
		return newpath;
	}


	private static int getStId(HashMap<Integer, SuperTarget> sts, int attackedentry, int attackedexit) {
		
		
		for(SuperTarget st: sts.values())
		{
			if(st.ap.keySet().contains(attackedentry) && st.ap.keySet().contains(attackedexit)) 
			{
				return st.stid;
			}
		}
		
		
		return -1;
	}


	public static double[] groupingBaseline(int base, int dest, int k, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, 
			ArrayList<Integer>[] clusters, int dmaxsuper, int dminsuper) throws Exception
	{
		double[] res = new double[5];
		
		long revmaptime = 0;
		long solvingtime = 0;


		/*
		int base = 0;
		int dest = 0;
		int k = 5;
		int radius = 2;
		int dmax = 25;
		int nRes=1;
		int dlim = 10;
		int nTargets = 20;*/

		//HashMap<Integer, TargetNode>  ta = makeGraph(5, 5, 20, 20, 2, 10, 2);
		//printNodesWithNeighborsAndPath(ta);

		// create graph
		//ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		//HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		/*for(TargetNode t: targets)
		{
			targetmaps.put(t.getTargetid(), t);
		}
		 */

		//makeGraph(k, radius, dlim , nTargets, 2, 10, 2, targets, targetmaps);

		//printNodesWithNeighborsAndPath(targets);


		HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>(); // map for apsp
		HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>(); //mapback for apsp
		
		
		 //MultiValuedMap<String, String> map = new MultiValuedHashMap<String, String>();
		

		//fix it
		HashMap<String,Integer> subgamemap = new HashMap<String, Integer>(); // map for all the subgames
		HashMap<Integer,String> subgamemapback = new HashMap<Integer, String>(); // mapback for all the subgames
		HashMap<String,HashMap<ArrayList<Integer>,Double >> subgamepaths = new HashMap<String,HashMap<ArrayList<Integer>,Double>>();

		ArrayList<ArrayList<Double>> origpaths = new ArrayList<ArrayList<Double>>();


		int[][] apspmat = buildAPSP(targets.size(), apspmap, apspmapback, targets, targetmaps);
		//ArrayList<Integer>[] clusters = groupTargets(base, dest, targetmaps, k, radius, targets, apspmat, apspmap);
		//printClusters(clusters);
		HashMap<Integer, SuperTarget> sts =  SuperTarget.buildSuperTargets(clusters,targetmaps);  // sts are supre targets
	//	printSuperTargets(sts);

		int[][] gamedata = buildGamedata(targetmaps);

		// build tables for each of the supertargets

		double[][][][] sttable = new double[k][][][];



		// do i keep track of which targets were covered for each entry?
		// supertable contains all the subgame info's
		
		Date start = new Date();
		long l1 = start.getTime();

		//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
		
		//fix it duplicate key can happen
		HashMap<String, double[]>  supertable = buildSuperTable(sttable, targets, targetmaps, sts, dmaxsuper, gamedata,
				apspmat, apspmap, subgamemap, subgamemapback, subgamepaths, dminsuper);

		Date stop = new Date();
		long l2 = stop.getTime();
		long diff = l2 - l1;
		solvingtime += diff;
		
		

		//printSubGamePathsPayoffs(subgamepaths, supertable);

		// generate pathsxxxxxxx
		try {
			ArrayList<ArrayList<String>> superpathseq = new ArrayList<ArrayList<String>>();
			ArrayList<SuperTarget> supergoals = generatePaths(dmax, sts, base, dest, supertable, dmaxsuper);
			makeSuperPathSeqSrcDest(superpathseq, supergoals);
			if(superpathseq.size()==0)
			{
				throw new Exception("No superpath");
			}
			//printSuperPaths(superpathseq);


			int [][] p = new int[supertable.size()][]; // p matrix for subgame
			Integer[] input = new Integer[superpathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(superpathseq.size()==0)
			{

			}
			else
			{
				// make joint schedule
				// solve using an LP

				if(superpathseq.size()<nRes)
				{

					branch = new int[superpathseq.size()];
					jSet=combine(input, superpathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}
				List<ArrayList<Integer>> superjset = new ArrayList<ArrayList<Integer>>(jSet);
				//printJointSchedule(superjset);
				p = makeSuperPmat(superpathseq, superjset, subgamemap, supertable);

				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				
				start = new Date();
				l1 = start.getTime();

				
				double[] probdistribution = MIPSolver4.solveForSuperAttacker(p, sts, supertable, nRes, astrategy, subgamemap, subgamemapback );
				
				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;
				solvingtime += diff;
				
				
				
				start = new Date();
				l1 = start.getTime();

				
				

				String attackedsupertarget = findPlayedSubGameWMapping(p, probdistribution, subgamemap, subgamemapback, supertable);
				

				stop = new Date();
				l2 = stop.getTime();
				diff = l2 - l1;
				revmaptime += diff;

				System.out.println("attacked subgame  "+ attackedsupertarget);

				//printSubGamePath(subgamepaths,attackedsupertarget);

				//String sattackedsupertarget = String.valueOf(attackedsupertarget);// attackedsupertarget

				/*if(attackedsupertarget==0)
				{
					sattackedsupertarget += "0"+"0"; 
				}
				 */

				int attakedtarget = (int)supertable.get(attackedsupertarget)[2];

				System.out.println("attacked target  "+ attakedtarget);


				// don't know how will work
				//int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps);


				if(attakedtarget!=-1)
				{
					// Integer is the path index
					// Double is the path prob
					
					start = new Date();
					l1 = start.getTime();

					
					
					HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
					buildOrigPaths(probdistribution, superjset, superpathseq, subgamepaths, targetmaps, finalpaths);

					//TODO check the probability distribution
					printStrategyPaths(probdistribution, superjset, finalpaths);
					int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps, finalpaths);
					//double sprob = getSubGamePathProb(superjset.get(0), superpathseq, finalpaths);
					HashMap<Integer, Double>  attckrexppayoffs = computeAttackerExpectedPayoffs(targetmaps, probdistribution, origPmat,
							superjset, superpathseq, finalpaths);
					// find the [attacked target, def payoffs, attckr payoffs]
					double[] attackedtarget = findAttackedTarget(attckrexppayoffs, targetmaps, origPmat, probdistribution, superjset, superpathseq, finalpaths);


					stop = new Date();
					l2 = stop.getTime();
					diff = l2 - l1;
					revmaptime += diff;
					
					
					System.out.println("\nAtatcked target: "+ (attackedtarget[0])+ "\n Def Exp: "+ (attackedtarget[1]) + "\n Attckr Exp: "+ (attackedtarget[2]));

					res[0] = attackedtarget[1];
					res[1] = attackedtarget[2];
					//System.out.println("hiii ");

				}


				//System.out.println("hiii ");

			}




		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}




		res[3] = revmaptime;
		res[4] = solvingtime;
		return res;
	}
	
	
	public static void groupingWithDOExp(int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap,HashMap<Integer,ArrayList<Integer>[]> allclus, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;
		long clusteringtime = 0;
		long slavetime = 0;
		int finalsize = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			
			
			
			
		/*	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
			
			
			ArrayList<Integer>[] clus = GroupingTargets.createGraph3(targets, targetmaps);*/



			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = groupingWithDO(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, clus);
			//double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime};
			
			
			sumdefexp += res[0];
			solvingtime+= res[2];
			revmaptime += res[6];
			clusteringtime += res[1];
			slavetime += res[5];
			finalsize += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		finalsize /= LIMIT;
		clusteringtime /= LIMIT;
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		writeInFileST("ClusteringWithDO",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime);
		

	}
	
	public static void wekaClusteringWithDOExpRW(int nrow, int ncol, int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap,HashMap<Integer,ArrayList<Integer>[]> allclus, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;
		long clusteringtime = 0;
		long slavetime = 0;
		int finalsize = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			/*ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			*/
			
			
			
			
			
			
			
		/*	ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
			
			
			ArrayList<Integer>[] clus = GroupingTargets.createGraph3(targets, targetmaps);*/
			
			
			
			
			/*int nrow = 50;
			int ncol = 50;*/
			double utility [][] = new double[560][560];
			double elevation [][] = new double[560][560];
			double u [][] = new double[nrow][ncol];
			double e [][] = new double[nrow][ncol];
			
			
			
			
			ReadData.readData(560, 560, utility, elevation);
			ReadData.getChunk(utility, elevation, 0, 0, nrow, ncol, u, e);
			
			
			ReadData.createCSVData(560, 560, utility, elevation);
			//int[][] gamedata = SecurityGameContraction.constructGameData(utility);
			ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
			SecurityGameContraction.buildcsvGraph(nrow,ncol,u, e,targets );
			
			HashMap<Integer,TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
			for(TargetNode t: targets)
			{
				targetmaps.put(t.getTargetid(), t);
			}
			
			
			printNodesWithNeighborsAndPath(targetmaps);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = wekaClusteringWithDORW(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps);
			//double[] res1 = {defpayoff, clusteringtime, solvingtime, targetstocluster.size(), attackeru, slavetime, revmaptime};
			
			
			sumdefexp += res[0];
			solvingtime+= res[2];
			revmaptime += res[6];
			clusteringtime += res[1];
			slavetime += res[5];
			finalsize += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		finalsize /= LIMIT;
		clusteringtime /= LIMIT;
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		writeInFileST("AutoClusteringWithDO",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime);
		

	}
	
	
	
	
	
	
	public static void groupTargetBaseline3Exp(int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap,HashMap<Integer,ArrayList<Integer>[]> allclus, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int dmaxsuper, int dminsuper ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			//ArrayList<Integer>[] clus = allclus.get(iter);
			//ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			//HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			
			ArrayList<TargetNode> targets = new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
			
			
			ArrayList<Integer>[] clus = GroupingTargets.createGraph3(targets, targetmaps);



			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = groupingBaseline3(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, clus, dmaxsuper, dminsuper);
			sumdefexp += res[0];
			solvingtime+= res[4];
			revmaptime += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		writeInFile(nTargets, k, sumdefexp, solvingtime, revmaptime, totaltime);
		

	}


	public static void groupTargetBaseline2Exp(int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap,HashMap<Integer,ArrayList<Integer>[]> allclus, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int dmaxsuper, int dminsuper ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = groupingBaseline2(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, clus, dmaxsuper, dminsuper);
			sumdefexp += res[0];
			solvingtime+= res[4];
			revmaptime += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		writeInFile(nTargets, k, sumdefexp, solvingtime, revmaptime, totaltime);
		

	}
	
	
	public static void groupTargetBaseline2ExpNewMILP(int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap,HashMap<Integer,ArrayList<Integer>[]> allclus, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int dmaxsuper, int dminsuper ) throws Exception
	{



		
		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = groupingBaseline2NewMILP(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, clus, dmaxsuper, dminsuper);
			sumdefexp += res[0];
			solvingtime+= res[4];
			revmaptime += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		writeInFile(nTargets, k, sumdefexp, solvingtime, revmaptime, totaltime);
		

	}
	
	
	

	public static void groupTargetBaselineExp(int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap,HashMap<Integer,ArrayList<Integer>[]> allclus, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
			HashMap<Integer,HashMap<Integer,TargetNode>> alltargetmaps, int dmaxsuper, int dminsuper ) throws Exception
	{



		/*int base = 0;
		int dest = 0;
		int k = 4;
		int radius = 1;
		int dmax = 10;
		int nRes=1;
		int dlim = 4;
		int nTargets = 12;
		int LIMIT = 1;
		int ap = 2;
		 */
		//HashMap<Integer, TargetNode>  ta = makeGraph(5, 5, 20, 20, 2, 10, 2);
		//printNodesWithNeighborsAndPath(ta);

		// create graph

		/*for(TargetNode t: targets)
		{
			targetmaps.put(t.getTargetid(), t);
		}
		 */


		double sumdefexp = 0.0;
		long totaltime = 0;
		long solvingtime = 0;
		long revmaptime = 0;

		for(int iter=0; iter<LIMIT; iter++)
		{
			ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = groupingBaseline(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, clus, dmaxsuper, dminsuper);
			sumdefexp += res[0];
			solvingtime+= res[4];
			revmaptime += res[3];

			Date stop = new Date();
			long l2 = stop.getTime();
			long diff = l2 - l1;
			totaltime += diff;

		}
		
		sumdefexp /= LIMIT;
		solvingtime /= LIMIT;
		revmaptime /= LIMIT;
		totaltime /= LIMIT; 
		
		//System.out.println("Defender exp "+ (double)sumdefexp/LIMIT + ", time : "+ (long)totaltime/LIMIT);
		writeInFile(nTargets, k, sumdefexp, solvingtime, revmaptime, totaltime);
		

	}
	
	
	private static void writeInFile(int nTargets, int k, double defexp, long solvingtime, long revmaptime, long totaltime) 
	{

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"grp-result.csv"),true));
			pw.append(nTargets+ "," + k+ ","+defexp+"," + solvingtime+ ","+revmaptime+"," + totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	
	
	private static void writeInFileST(String algo, int finalsize, double defexp, long solvingtime, long revmaptime, long clusteringtime, long totaltime) 
	{
		
		//ClusteringWithDO",finalsize,sumdefexp, solvingtime, revmaptime, clusteringtime ,totaltime

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"grp-result.csv"),true));
			pw.append(algo+ ","+finalsize+","+defexp+"," + clusteringtime+ ","+solvingtime+ ","+revmaptime+","+totaltime+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}


	private static void writeLinearityCheck(double coverage, double expectedpayoff) 
	{


		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("LCHK.csv"),true));
			//PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append( coverage+"," + expectedpayoff+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}



	public static void groupTargetTest() throws Exception
	{



		int base = 0;
		int dest = 0;
		int k = 5;
		int radius = 2;
		int dmax = 25;
		int nRes=1;
		int dlim = 10;
		int nTargets = 20;
		int dmaxsuper = 3;

		//HashMap<Integer, TargetNode>  ta = makeGraph(5, 5, 20, 20, 2, 10, 2);
		//printNodesWithNeighborsAndPath(ta);

		// create graph
		ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		/*for(TargetNode t: targets)
		{
			targetmaps.put(t.getTargetid(), t);
		}
		 */

		makeGraph(k, radius, dlim , nTargets, 2, 10, 2, targets, targetmaps);

		printNodesWithNeighborsAndPath(targets);


		HashMap<Integer, Integer> apspmap = new HashMap<Integer, Integer>(); // map for apsp
		HashMap<Integer, Integer> apspmapback = new HashMap<Integer, Integer>(); //mapback for apsp

		HashMap<String,Integer> subgamemap = new HashMap<String, Integer>(); // map for all the subgames
		HashMap<Integer,String> subgamemapback = new HashMap<Integer, String>(); // mapback for all the subgames

		HashMap<String,HashMap<ArrayList<Integer>,Double >> subgamepaths = new HashMap<String,HashMap<ArrayList<Integer>,Double>>();

		ArrayList<ArrayList<Double>> origpaths = new ArrayList<ArrayList<Double>>();


		int[][] apspmat = buildAPSP(targets.size(), apspmap, apspmapback, targets, targetmaps);
		ArrayList<Integer>[] clusters = groupTargets(base, dest, targetmaps, k, radius, targets, apspmat, apspmap);
		printClusters(clusters);
		HashMap<Integer, SuperTarget> sts =  SuperTarget.buildSuperTargets(clusters,targetmaps);  // sts are supre targets
		printSuperTargets(sts);

		int[][] gamedata = buildGamedata(targetmaps);

		// build tables for each of the supertargets

		double[][][][] sttable = new double[k][][][];



		// do i keep track of which targets were covered for each entry?
		// supertable contains all the subgame info's
		HashMap<String, double[]>  supertable = buildSuperTable(sttable, targets, targetmaps, sts, 4, gamedata,
				apspmat, apspmap, subgamemap, subgamemapback, subgamepaths, 14);


		printSubGamePathsPayoffs(subgamepaths, supertable);

		// generate paths
		try {
			ArrayList<ArrayList<String>> superpathseq = new ArrayList<ArrayList<String>>();
			ArrayList<SuperTarget> supergoals = generatePaths(dmax, sts, base, dest, supertable, dmaxsuper);
			makeSuperPathSeqSrcDest(superpathseq, supergoals);
			printSuperPaths(superpathseq);


			int [][] p = new int[supertable.size()][]; // p matrix for subgame
			Integer[] input = new Integer[superpathseq.size()];
			int[] branch = new int[nRes];//{0,0};//new char[k];

			for(int i=0; i<input.length; i++)
			{
				input[i] = i;
			}
			HashSet jSet=new HashSet();


			if(superpathseq.size()==0)
			{

			}
			else
			{
				// make joint schedule
				// solve using an LP

				if(superpathseq.size()<nRes)
				{

					branch = new int[superpathseq.size()];
					jSet=combine(input, superpathseq.size(), 0, branch, 0, jSet);
				}
				else
				{
					jSet=combine(input, nRes, 0, branch, 0, jSet);
				}
				List<ArrayList<Integer>> superjset = new ArrayList<ArrayList<Integer>>(jSet);
				//printJointSchedule(jset);
				p = makeSuperPmat(superpathseq, superjset, subgamemap, supertable);

				HashMap<Integer, Double> astrategy = new HashMap<Integer, Double>();
				double[] probdistribution = MIPSolver4.solveForSuperAttacker(p, sts, supertable, nRes, astrategy, subgamemap, subgamemapback );

				String attackedsupertarget = findPlayedSubGameWMapping(p, probdistribution, subgamemap, subgamemapback, supertable);

				System.out.println("attacked subgame  "+ attackedsupertarget);

				String sattackedsupertarget = String.valueOf(attackedsupertarget);// attackedsupertarget

				/*if(attackedsupertarget==0)
				{
					sattackedsupertarget += "0"+"0"; 
				}
				 */

				int attakedtarget = (int)supertable.get(sattackedsupertarget)[2];

				System.out.println("attacked supertarget  "+ attakedtarget);


				// don't know how will work
				//int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps);


				if(attakedtarget!=-1)
				{
					// Integer is the path index
					// Double is the path prob
					HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
					buildOrigPaths(probdistribution, superjset, superpathseq, subgamepaths, targetmaps, finalpaths);
					int[][] origPmat = buildOrigPmat(superjset, superpathseq, subgamepaths, targetmaps, finalpaths);
					//double sprob = getSubGamePathProb(superjset.get(0), superpathseq, finalpaths);
					HashMap<Integer, Double>  attckrexppayoffs = computeAttackerExpectedPayoffs(targetmaps, probdistribution, origPmat,
							superjset, superpathseq, finalpaths);
					// find the [attacked target, def payoffs, attckr payoffs]
					double[] attackedtarget = findAttackedTarget(attckrexppayoffs, targetmaps, origPmat, probdistribution, superjset, superpathseq, finalpaths);

					System.out.println("\nAtatcked target: "+ attackedtarget[0]+ "\n Def Exp: "+ attackedtarget[1] + "\n Attckr Exp: "+ attackedtarget[2]);

					System.out.println("hiii ");

				}


				System.out.println("hiii ");

			}




		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}

	private static double getSubGamePathProb(ArrayList<Integer> jointschedule, ArrayList<ArrayList<String>> superpathseq,
			HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths, int target) {


		double prob = 0.0;



		for(Integer pathindex: jointschedule) // there is only one joint patrolling path
		{
			// get the path and check if target is covered
			double maxprob = 0.0;
			for(ArrayList<Double> path:  finalpaths.get(pathindex))
			{
				//ArrayList<Double> path = finalpaths.get(pathindex).get(i);
				for(int j=1; j<path.size(); j++)
				{
					if(path.get(j)==target)
					{
						maxprob = Math.max(maxprob, path.get(0));
						break;
					}
				}
			}
			prob += maxprob;
		}
		//prob /= jointschedule.size(); // divide by number of resources


		return prob;
	}

	/**
	 * 
	 * @param attckrexppayoffs
	 * @param targetmaps
	 * @param origPmat
	 * @param probdistribution
	 * @param finalpaths 
	 * @param superpathseq 
	 * @param superjset 
	 * @return attacked target and defender exputility and attacker exp
	 */


	private static double[] findAttackedTarget(HashMap<Integer, Double> attckrexppayoffs,
			HashMap<Integer, TargetNode> targetmaps, int[][] origPmat, double[] probdistribution, List<ArrayList<Integer>> superjset,
			ArrayList<ArrayList<String>> superpathseq, HashMap<Integer,ArrayList<ArrayList<Double>>> finalpaths ) {


		double attackrmax = Double.NEGATIVE_INFINITY;
		int maxtarget = -1;
		
		

		for(Integer t: attckrexppayoffs.keySet())
		{
			
			System.out.println("Atckr payoff target "+ t + " = "+ attckrexppayoffs.get(t));
			
			//int target = Integer.parseInt(t);
			if(attackrmax<attckrexppayoffs.get(t))
			{
				attackrmax = attckrexppayoffs.get(t);
				maxtarget = t;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(Integer t: attckrexppayoffs.keySet())
		{
			//int target = Integer.parseInt(t);
			if((t!=maxtarget) && (attckrexppayoffs.get(t) == attackrmax))
			{
				System.out.println("tied target "+t +", attkr payoff "+ attckrexppayoffs.get(t));

				tiedtargets.add(t);
			}
		}

		int selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{

			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			HashMap<Integer, Double>  defenderexpectedpayoffs =  defenderExpectedPayoffs(tiedtargets, probdistribution, origPmat, 
					targetmaps, superjset, superpathseq, finalpaths);
			/*for(int target: defenderexpectedpayoffs.keySet())
			{
				System.out.println("target "+target +", defender payoff "+ defenderexpectedpayoffs.get(target));
			}
			 */
			double defmax = Double.NEGATIVE_INFINITY;
			int defmaxtarget = -1;

			for(Integer t: defenderexpectedpayoffs.keySet())
			{
				//int target = Integer.parseInt(t);
				System.out.println("target "+t +", defender payoff "+ defenderexpectedpayoffs.get(t));
				if(defmax<defenderexpectedpayoffs.get(t))
				{
					defmax = defenderexpectedpayoffs.get(t);
					defmaxtarget = t;
				}
			}
			selectedtarget = defmaxtarget;
			return new double[] {selectedtarget, defmax, attackrmax};
		}



		double defmax = defenderExpectedPayoff(selectedtarget, probdistribution, origPmat, targetmaps, superjset, superpathseq, finalpaths);

		return new double[] {selectedtarget, defmax, attackrmax};



	}



	private static HashMap<Integer, Double> defenderExpectedPayoffs(ArrayList<Integer> tiedtargets,
			double[] probdistribution, int[][] origPmat, HashMap<Integer,TargetNode> targetmaps,
			List<ArrayList<Integer>> superjset, ArrayList<ArrayList<String>> superpathseq, HashMap<Integer,ArrayList<ArrayList<Double>>> finalpaths ) {



		HashMap<Integer, Double> expectedpayoffs = new HashMap<Integer, Double>();

		double c = 0;
		for(Integer target: tiedtargets)
		{
			//System.out.println("\n\nTarget "+ target);



			c=0;
			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{
				double subgameprob = getSubGamePathProb(superjset.get(jointschedule), superpathseq, finalpaths, target);
				c += origPmat[target][jointschedule]* probdistribution[jointschedule]* subgameprob;

			}

			double x = c*0.0;
			double y = (1-c)*targetmaps.get(target).defenderpenalty;
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs.put(target, Math.round((x+y)*100)/100.00) ; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}


	private static Double defenderExpectedPayoff(int target,
			double[] probdistribution, int[][] origPmat, HashMap<Integer,TargetNode> targetmaps,
			List<ArrayList<Integer>> superjset, ArrayList<ArrayList<String>> superpathseq, HashMap<Integer,ArrayList<ArrayList<Double>>> finalpaths ) {



		double expectedpayoff = 0;

		double c = 0;
		//for(Integer target: tiedtargets)
		{
			//System.out.println("\n\nTarget "+ target);



			c=0;
			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{
				double subgameprob = getSubGamePathProb(superjset.get(jointschedule), superpathseq, finalpaths, target);

				c += origPmat[target][jointschedule]* probdistribution[jointschedule]* subgameprob;

			}

			double x = c*targetmaps.get(target).defenderreward;
			double y = (1-c)*targetmaps.get(target).defenderpenalty;
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoff = Math.round((x+y)*100)/100.00 ; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoff;
	}


	private static HashMap<Integer, Double> computeAttackerExpectedPayoffs(HashMap<Integer, TargetNode> targetmaps,
			double[] probdistribution, int[][] origPmat, List<ArrayList<Integer>> superjset, ArrayList<ArrayList<String>> superpathseq, 
			HashMap<Integer,ArrayList<ArrayList<Double>>> finalpaths) {


		HashMap<Integer, Double> expectedpayoffs = new HashMap<Integer, Double>();

		double c = 0;
		for(Integer target: targetmaps.keySet())
		{
			//System.out.println("\n\nTarget "+ target);
			c=0;
			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{

				// 
				double subgameprob = getSubGamePathProb(superjset.get(jointschedule), superpathseq, finalpaths, target);
				// need to use max prob obtained from bfs
				c += origPmat[target][jointschedule]* probdistribution[jointschedule]*subgameprob;

			}

			double x = c*0.0;
			double y = (1-c)*targetmaps.get(target).attackerreward;
			System.out.println("target "+ target+ ", coverage "+ c+ ", r:"+y+",p:"+x);
			expectedpayoffs.put(target, Math.round((x+y)*100)/100.00) ; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;




	}



	private static int[][] buildOrigPmat(List<ArrayList<Integer>> superjset, ArrayList<ArrayList<String>> superpathseq,
			HashMap<String, HashMap<ArrayList<Integer>, Double>> subgamepaths,
			HashMap<Integer, TargetNode> targetmaps, HashMap<Integer,ArrayList<ArrayList<Double>>> finalpaths) {


		int[][] pmat = new int[targetmaps.size()][superjset.size()];


		// for every target

		for(Integer tid: targetmaps.keySet())
		{
			// find out if it's covered by the joint schedule
			for(ArrayList<Integer> jointschedule: superjset)
			{


				// find out if it's covered by the paths
				for(Integer superpath: jointschedule)
				{
					boolean iscovered = isCoveredByPath(finalpaths.get(superpath), tid);
					if(iscovered)
					{
						pmat[tid][superjset.indexOf(jointschedule)] = 1;
						break;
					}
					else
					{
						pmat[tid][superjset.indexOf(jointschedule)] = 0;
					}
				}



			}
		}


		return pmat;
	}



	private static boolean isCoveredByPath(ArrayList<ArrayList<Double>> paths, Integer tid) {


		for(ArrayList<Double> p: paths)
		{
			for(int i=1; i<p.size(); i++)
			{
				if(tid == p.get(i).intValue())
					return true;
			}
		}
		return false;
	}



	private static void buildOrigPathsWithProb(ArrayList<ArrayList<Double>> origpaths, double[] probdistribution,
			List<ArrayList<Integer>> superjset, ArrayList<ArrayList<String>> superpathseq,
			HashMap<String,HashMap<ArrayList<Integer>,Double>> subgamepaths, HashMap<Integer, TargetNode> targetmaps, 
			HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths) throws Exception {



		/*for(int superprobindex=0; superprobindex<probdistribution.length; superprobindex++)
		{
			for(int superpathindex: superjset.get(superprobindex))
			{
				ArrayList<Integer> origpath = new ArrayList<Integer>();
				// get the path
				for(Integer supernode: superpathseq.get(superpathindex))
				{
					// for every super target reversemap the paths
		 *//**
		 * add every node 
		 *//*
					for(Double subprob: subgamepaths.get(supernode).keySet())
					{
						// wont work
						// use BFS
					}

				}
			}

		}*/

		buildOrigPaths( probdistribution,superjset, superpathseq, subgamepaths, targetmaps, finalpaths);


		//return goals;

	}



	private static void buildOrigPaths(
			double[] probdistribution, List<ArrayList<Integer>> superjset, ArrayList<ArrayList<String>> superpathseq,
			HashMap<String, HashMap<ArrayList<Integer>, Double>> subgamepaths,HashMap<Integer, TargetNode> targetmaps, 
			HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths) throws Exception {

		// run bfs for every path

		ArrayList<SuperTarget> goals = new ArrayList<SuperTarget>();

		//HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();

		for(ArrayList<String> spath : superpathseq)
		{

			/*System.out.println("Expanding path: ");
			for(Integer p: spath)
			{
				System.out.print(p+"->");
			}*/
			ArrayList<ArrayList<Double>> finalpath = new ArrayList<ArrayList<Double>>();
			ArrayList<SuperTarget> subgoals =  bfs(spath, subgamepaths, targetmaps, finalpath);
			finalpaths.put(superpathseq.indexOf(spath), finalpath);
			// for each of the subgoals add the path with probs in final paths
			//print paths with probs

		}




		//return goals;
	}


	private static void printStrategyPaths(
			double[] probdistribution, List<ArrayList<Integer>> superjset, 

			HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths) {

		// run bfs for every path

		//ArrayList<SuperTarget> goals = new ArrayList<SuperTarget>();

		//HashMap<Integer, ArrayList<ArrayList<Double>>> finalpaths = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		System.out.print("\nDEfender patrollling path :\n\n");
		for(ArrayList<Integer> jschedule : superjset)
		{

			//if(probdistribution[superjset.indexOf(jschedule)]>0)
			{
				for(Integer pinddex: jschedule)
				{
					System.out.println("Path index "+ pinddex);

					//finalpaths.put(superpathseq.indexOf(spath), finalpath);
					// for each of the subgoals add the path with probs in final paths
					//print paths with probs

					int f=1;
					for(ArrayList<Double> origpath : finalpaths.get(pinddex))
					{
						
						f=1;
						if(origpath.get(0)>0)
						{
							for(Double n: origpath)
							{
								if(origpath.indexOf(n)>0)
									System.out.print(n.intValue()+"->");
								else if(f==1)
								{
									f=0;
									System.out.print(n*probdistribution[superjset.indexOf(jschedule)]+"->");
								}
								else
									System.out.print(n+"->");

							}
							System.out.print("\n");
						}
					}
					System.out.print("\n");

				}
			}

		}
		//return goals;
	}






	private static ArrayList<SuperTarget> bfs(ArrayList<String> spath,
			HashMap<String, HashMap<ArrayList<Integer>, Double>> subgamepaths, HashMap<Integer, TargetNode> targtemaps, 
			ArrayList<ArrayList<Double>> finalpath) throws Exception {


		SuperTarget start = new SuperTarget();

		// get the targetnode of the base super node
		int basetnode = -1;

		String k = spath.get(0);
		String key = String.valueOf(k);

	/*	if(k==0)
		{
			key = "0,0,0";
		}
*/

		for(ArrayList<Integer> n: subgamepaths.get(key).keySet())
		{
			basetnode = n.get(0);
			break;
		}

		//start.nodes.put(basetnode, targtemaps.get(basetnode));
		start.pathnodesid.add(basetnode);
		if(spath.get(0).equals("0,0,0"))
		{
			start.stid = 0;
		}
		else
		{
			throw new Exception("nooooppppp"); // need to write code for handling when base is not 0
		}

		Queue<SuperTarget> fringequeue = new LinkedList<SuperTarget>();
		ArrayList<SuperTarget> goals = new ArrayList<SuperTarget>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);


		while(fringequeue.size()>0)
		{

			SuperTarget node = fringequeue.poll();


			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.stid==start.stid) && (node.distcoveredyet>0))
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//System.out.println();
				//SuperTarget.printPath(node);
				//System.out.println();
				ArrayList<Double> fpath = new ArrayList<Double>();
				fpath.add(node.currentprob);
				SuperTarget.printFinalPathWithProb(node, fpath);
				//System.out.println();
				finalpath.add(fpath);
				goals.add(node);

				/*if(pathcounter>5000)
					break;*/
			}
			else
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.stid);
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<SuperTarget> succs = Expand(node, start.stid, node.currentindex, node.currentprob, spath, subgamepaths, targtemaps);
				/**
				 * add nodes to queue
				 */
				for(SuperTarget suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}



		}





		return goals;
	}



	private static ArrayList<SuperTarget> Expand(SuperTarget node, int dest, int currentindex, double currentprob,
			ArrayList<String> spath, HashMap<String, HashMap<ArrayList<Integer>, Double>> subgamepaths, HashMap<Integer,TargetNode> targtemaps) {



		// get the neighbor
		if(currentindex==spath.size())
			return null;
		
		// get the paths 
		String neiid = spath.get(currentindex+1);
		
		
		// convert the neiid to an integer
		String sepneiid [] = neiid.split(",");
		String snid = "";
		for(String sn: sepneiid)
		{
			snid+= sn;
		}
		int intneiid = Integer.parseInt(snid);
		ArrayList<SuperTarget> sucs = new ArrayList<SuperTarget>();
		

		//System.out.println("hhhhheeeloo  "+ key);
		
		
		//TODO in all the paths compute which targets are common
		// compute the maximum prob
		// for common nodes use max prob
		// in that case in the same path two nodes may have different coverages. 

		for(ArrayList<Integer> path: subgamepaths.get(neiid).keySet())
		{
			// for every path create a super target
			// update the current prob
			// add the nodes of the path
			// set the parent
			
			
			SuperTarget tmp = new SuperTarget();
			tmp.currentindex = currentindex+1;
			tmp.stid= intneiid;//Integer.parseInt(neiid);
			tmp.currentprob = node.currentprob * subgamepaths.get(neiid).get(path);
			// need to have a data structure to save different probs for each target n in path if they are common node
			// need to have max prob if it's a  common nodes
			tmp.distcoveredyet++;
			//System.out.println("Adding nodes ");
			for(int n: path)
			{
				//System.out.print(n+" ");

				//tmp.nodes.put(n, targtemaps.get(n)); // lol problem if has same target id
				tmp.pathnodesid.add(n);
			}
			tmp.parent = node;
			//System.out.println("Parent node  "+ node.stid);

			sucs.add(tmp);

		}




		return sucs;
	}



	private static void printSubGamePathsPayoffs(HashMap<String, HashMap<ArrayList<Integer>,Double>> subgamepaths, HashMap<String,double[]> supertable) {


		for(String subgame: subgamepaths.keySet())
		{
			System.out.println("\nSubgame : "+ subgame);
			//int f=0;
			for(ArrayList<Integer> nodes: subgamepaths.get(subgame).keySet())
			{
				System.out.print(subgamepaths.get(subgame).get(nodes)+" : ");
				for(Integer node: nodes)
				{

					System.out.print(node+"->");
				}
				System.out.println();

			}
			//System.out.println();
			double[] res = supertable.get(subgame);
			System.out.print("def = "+ res[0]+ "  attkr = "+ res[1]+"\n");
		}

	}

	private static void printSubGamePath(HashMap<String, HashMap<ArrayList<Integer>,Double>> subgamepaths, String subgame) {


		//	for(String subgame: subgamepaths.keySet())
		{
			System.out.println("\nSubgame : "+ subgame);
			//int f=0;
			for(ArrayList<Integer> nodes: subgamepaths.get(subgame).keySet())
			{
				System.out.print(subgamepaths.get(subgame).get(nodes)+" : ");
				for(Integer node: nodes)
				{

					System.out.print(node+"->");
				}
				System.out.println();

			}
			System.out.println();
		}

	}



	private static Integer findAttackedSuperTargetWMapping(int[][] pmat, double[] probdistribution,
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback, HashMap<String,double[]> supertable,
			HashMap<String,Integer> subgamecontainer, HashMap<Integer,SuperTarget> sts) throws Exception {


		HashMap<Integer, Double> expectedpayoffs = expectedAttackerSuperPayoffsWithMapping(pmat, probdistribution, subgamemap, 
				subgamemapback,supertable, subgamecontainer, sts);
		System.out.println();
		for(Integer supertarget : expectedpayoffs.keySet())
		{
			System.out.println("super target "+supertarget+", attkr payoff "+ expectedpayoffs.get(supertarget));
		}
		double max = Double.NEGATIVE_INFINITY;
		Integer maxtarget = -1;

		for(Integer t: expectedpayoffs.keySet())
		{
			if(max<expectedpayoffs.get(t))
			{
				max = expectedpayoffs.get(t);
				maxtarget = t;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<Integer> tiedtargets = new ArrayList<Integer>();
		tiedtargets.add(maxtarget);
		for(Integer t: expectedpayoffs.keySet())
		{
			//int target = Integer.parseInt(t);
			if((!t.equals(maxtarget)) && (expectedpayoffs.get(t) == max))
			{
				System.out.println("tied supertarget "+t +", attkr payoff "+ expectedpayoffs.get(t));

				tiedtargets.add(t);
			}
		}

		Integer selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{

			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			HashMap<Integer, Double>  defenderexpectedpayoffs =  defenderSuperExpectedPayoffsWithMapping(tiedtargets, probdistribution, 
					pmat, subgamemap, subgamemapback,supertable,subgamecontainer, sts);
			
			
			for(Integer supertarget :defenderexpectedpayoffs.keySet())
			{
				System.out.println("supertarget "+supertarget +", defender payoff "+ defenderexpectedpayoffs.get(supertarget));
			}
			double defmax = Double.NEGATIVE_INFINITY;
			Integer defmaxtarget = -1;;

			for(Integer t: defenderexpectedpayoffs.keySet())
			{
				//int target = Integer.parseInt(t);
				if(defmax<defenderexpectedpayoffs.get(t))
				{
					defmax = defenderexpectedpayoffs.get(t);
					defmaxtarget = t;
				}
			}
			
			selectedtarget = defmaxtarget;
			return selectedtarget;
		}



		return selectedtarget;
	}
	
	
	private static String findPlayedSubGameWMapping(int[][] pmat, double[] probdistribution,
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback, HashMap<String,double[]> supertable) throws Exception {


		HashMap<String, Double> expectedpayoffs = expectedAttackerPayoffsWithMapping(pmat, probdistribution, subgamemap, subgamemapback,supertable);

		for(String subgame : expectedpayoffs.keySet())
		{
			System.out.println("target "+subgame+", attkr payoff "+ expectedpayoffs.get(subgame));
		}
		double max = Double.NEGATIVE_INFINITY;
		String maxtarget = "";

		for(String t: expectedpayoffs.keySet())
		{
			
			/*String [] x = t.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}
			*/
			
			if(max<expectedpayoffs.get(t))
			{
				max = expectedpayoffs.get(t);
				maxtarget = t;
			}
		}

		//System.out.println("max target "+maxtarget+", attkr payoff "+ expectedpayoffs[maxtarget]);


		/**
		 * check if there are ties
		 */

		ArrayList<String> tiedtargets = new ArrayList<String>();
		tiedtargets.add(maxtarget);
		for(String t: expectedpayoffs.keySet())
		{
			//int target = Integer.parseInt(t);
			if((!t.equals(maxtarget)) && (expectedpayoffs.get(t) == max))
			{
				System.out.println("tied target "+t +", attkr payoff "+ expectedpayoffs.get(t));

				tiedtargets.add(t);
			}
		}

		String selectedtarget = maxtarget;

		if(tiedtargets.size()>1)
		{

			/**
			 * there are tied targets
			 * choose in favor of defender
			 */
			HashMap<String, Double>  defenderexpectedpayoffs =  defenderExpectedPayoffsWithMapping(tiedtargets, probdistribution, pmat, subgamemap, subgamemapback,supertable);
			for(String subgame :defenderexpectedpayoffs.keySet())
			{
				System.out.println("target "+subgame +", defender payoff "+ defenderexpectedpayoffs.get(subgame));
			}
			double defmax = Double.NEGATIVE_INFINITY;
			String defmaxtarget = "";

			for(String t: defenderexpectedpayoffs.keySet())
			{
				//int target = Integer.parseInt(t);
				if(defmax<defenderexpectedpayoffs.get(t))
				{
					defmax = defenderexpectedpayoffs.get(t);
					defmaxtarget = t;
				}
			}
			/**
			 * select all the subgames with defmax
			 */
			ArrayList<String> deftiedtargets = new ArrayList<String>();
			for(String t: defenderexpectedpayoffs.keySet())
			{
				//int target = Integer.parseInt(t);
				if(defmax==defenderexpectedpayoffs.get(t))
				{
					deftiedtargets.add(t);
				}
			}
			/**
			 * choose subgame with max d
			 */
			
			int maxd = -1;
			for(String starget: deftiedtargets)
			{
				int subd = Integer.parseInt(starget.split(",")[2]);
				if(maxd<subd)
				{
					maxd=subd;
					defmaxtarget = starget;
				}
			}
			selectedtarget = defmaxtarget;
			return selectedtarget;
		}



		return selectedtarget;
	}

	private static HashMap<String, Double> defenderExpectedPayoffsWithMapping(
			ArrayList<String> tiedtargets, double[] coverage,
			int [][] pmat, 
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback, HashMap<String,double[]> supertable) throws Exception {

		HashMap<String, Double> expectedpayoff = new HashMap<String, Double>();

		double c = 0;



		for(String target: tiedtargets)
		{
			//int target = tiedtargets.get(t);
			c=0;
			for(int path=0; path<coverage.length; path++)
			{
				c += pmat[subgamemap.get(target)][path]*coverage[path];


				/*double defreward = origpmat[target][path]*coverage[path]*gamedata[target][0];
					double defpenalty = (1-origpmat[target][path])*(1-coverage[path])*gamedata[target][1];
					sum = sum + defpenalty + defreward;*/
			}

			//String starget = String.valueOf(target);

			//System.out.println(" subgame $$$$  "+ target);

			if(!supertable.containsKey(target))
			{
				throw new Exception("Subgame does not exists....");
			}

			double tmp = c*0.0 + (1-c)*(-supertable.get(target)[1]);
			expectedpayoff.put(target,  Math.round(tmp*100)/100.00);

		}


		return expectedpayoff;
	}
	
	
	private static HashMap<Integer, Double> defenderSuperExpectedPayoffsWithMapping(
			ArrayList<Integer> tiedtargets, double[] probdistribution,
			int [][] pmat, 
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback, 
			HashMap<String,double[]> supertable, HashMap<String,Integer> subgamecontainer, HashMap<Integer,SuperTarget> sts) throws Exception {

		
		HashMap<Integer, Double> expectedpayoffs = new HashMap<Integer, Double>();// for supertargets

		for(Integer stid: tiedtargets)
		{
			
			double exp = 0;
			String maxsubgame = "";
			double maxc = -1;
			for(String subgame: subgamecontainer.keySet())
			{
				double c = 0;
				if(subgamecontainer.get(subgame)==stid) // if that's the super target where subgame belongs
				{
					
					for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
					{
						c += pmat[subgamemap.get(subgame)][jointschedule]* probdistribution[jointschedule];

					}
					if(maxc<c)
					{
						maxc = c;
						maxsubgame = subgame;
					}

					
				}

				//}
			}
			double x = maxc*0.0/*(-supertable.get(starget)[0])*/;
			double y = (maxc)*(-supertable.get(maxsubgame)[1]);
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			//expectedpayoffs.put(starget, Math.round((x+y)*100)/100.00) ; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);
			exp += Math.round((x+y)*100)/100.00;
			expectedpayoffs.put(stid, exp);

		
		}
		return expectedpayoffs;
		
	}



	private static HashMap<String, Double> expectedAttackerPayoffsWithMapping(int[][] pmat,
			double[] probdistribution, HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback, HashMap<String,double[]> supertable) {


		HashMap<String, Double> expectedpayoffs = new HashMap<String, Double>();

		double c = 0;
		for(String starget: supertable.keySet())
		{
			//System.out.println("\n\nTarget "+ target);

			/*String [] z = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: z)
			{
				stringid.append(y);
			}*/
			

			c=0;
			for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
			{
				c += pmat[subgamemap.get(starget)][jointschedule]* probdistribution[jointschedule];

			}

			double x = c*0.0/*(-supertable.get(starget)[0])*/;
			double y = (1-c)*(supertable.get(starget)[1]);
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			expectedpayoffs.put(starget, Math.round((x+y)*100)/100.00) ; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);

		}


		return expectedpayoffs;
	}
	
	
	private static HashMap<Integer, Double> expectedAttackerSuperPayoffsWithMapping(int[][] pmat,
			double[] probdistribution, HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback, 
			HashMap<String,double[]> supertable, HashMap<String,Integer> subgamecontainer, HashMap<Integer,SuperTarget> sts) {


		HashMap<Integer, Double> expectedpayoffs = new HashMap<Integer, Double>();// for supertargets

		for(Integer stid: sts.keySet())
		{
			
			double exp = 0;
			System.out.print("\n\nST "+ stid);
			String maxsubgame = "";
			double maxc = -1;
			for(String subgame: subgamecontainer.keySet())
			{
				double c = 0;
				if(subgamecontainer.get(subgame)==stid) // if that's the super target where subgame belongs
				{
					System.out.print("\nSubgame "+ subgame + " , c:  ");
					for(int jointschedule=0; jointschedule<probdistribution.length; jointschedule++)
					{
						c += pmat[subgamemap.get(subgame)][jointschedule]* probdistribution[jointschedule];

					}
					System.out.print(c+"\n ");
					if(maxc<c)
					{
						maxc = c;
						maxsubgame = subgame;
					}
					
					
				}

				//}
			}
			double x = maxc*0.0/*(-supertable.get(starget)[0])*/;
			double y = (maxc)*(supertable.get(maxsubgame)[1]);
			//System.out.println(" attk: target "+ target+ ", coverage "+ c[target]+ ", r:"+y+",p:"+x);
			//expectedpayoffs.put(starget, Math.round((x+y)*100)/100.00) ; //c[target]*gamedata[target][3] + (1-c[target]*gamedata[target][2]);
			System.out.println("maxsubgame "+maxsubgame+", c : "+maxc+" , exp : "+ (Math.round((x+y)*100)/100.00));
			exp += Math.round((x+y)*100)/100.00;
			expectedpayoffs.put(stid, exp);

		
		}
		return expectedpayoffs;
	}



	// map for every entry to an id of subgames

	private static int[][] makeSuperPmat(ArrayList<ArrayList<String>> superpathseq,
			List<ArrayList<Integer>> jset, HashMap<String,Integer> subgamemap, HashMap<String, double[]>  supertable) {


		int[][] pmat = new int[supertable.size()][jset.size()];

		System.out.println("number of subgames "+ supertable.size());
		/**
		 * mapping happened
		 */
		for(String t: supertable.keySet())
		{
			for(int j=0; j<jset.size(); j++)
			{
				/**
				 * check if target t is in j schedule
				 */
				
				/*String [] x = t.split(",");
				StringBuilder stringid = new StringBuilder();
				for(String y: x)
				{
					stringid.append(y);
				}*/
				
				//int tmpid = Integer.parseInt(stringid.toString());
				int tindex = subgamemap.get(t);//targets.get(t).getTargetid();
				int isinj = isInJointSchedule(t, jset.get(j), superpathseq);
				//System.out.println("tindex "+ tindex);

				pmat[tindex][j] = isinj;

			}
		}


		return pmat;
	}


	private static int isInJointSchedule(String combid, ArrayList<Integer> jointschedule,
			ArrayList<ArrayList<String>> superpathseq) {


		for(Integer j: jointschedule)
		{
			ArrayList<String> path = superpathseq.get(j);
			for(String tt: path)
			{
				if(combid.equals(tt))
					return 1;
			}
		}


		return 0;
	}


	private static void printJointSchedule(List<ArrayList<Integer>> jset) {

		System.out.println();
		for(int i=0; i<jset.size(); i++)
		{
			System.out.print(i + " : [");
			for(int j=0; j<jset.get(i).size(); j++)
			{
				System.out.print(jset.get(i).get(j)+",");
			}
			System.out.println("]");
		}
		System.out.println();


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


	private static void printSuperPaths(ArrayList<ArrayList<String>> pathseq) {


		System.out.println("Paths ");
		//Logger.logit("paths : \n");

		int i=0;
		for(ArrayList<String> path: pathseq)
		{
			System.out.print(i+" : ");

			for(String x: path)
			{
				System.out.print(x+"->");
				//Logger.logit(x+"->");
			}
			System.out.println();
			//Logger.logit("\n");
			i++;

		}



	}


	public static int[][] buildGamedata(HashMap<Integer, TargetNode> targetmaps) {


		int[][] gamedata = new int[targetmaps.values().size()][4];

		for(TargetNode t: targetmaps.values())
		{
			gamedata[t.getTargetid()][0] = (int)t.defenderreward;
			gamedata[t.getTargetid()][1] = (int)t.defenderpenalty;
			gamedata[t.getTargetid()][2] = (int)t.attackerreward;
			gamedata[t.getTargetid()][3] = (int)t.attackerpenalty;

		}
		return gamedata;

	}
	
	
	private static HashMap<String, double[]>  buildSuperTable(double[][][][] sttable, ArrayList<TargetNode> targets,
			HashMap<Integer, TargetNode> targetmaps, HashMap<Integer, SuperTarget> sts, int dmax, int[][] gamedata,
			int[][] apspmat, HashMap<Integer,Integer> apspmap, 
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback,
			HashMap<String,HashMap<ArrayList<Integer>,Double>> subgamepaths, int dminsuper) throws Exception
	{



		HashMap<String, double[]> supertable = new  HashMap<String, double[]>();

		int subgamecount = 0;
		for(Integer stid: sts.keySet())
		{
			SuperTarget st = sts.get(stid);
			//if(st.nodes.size()>1)
			{
				//System.out.println("\nST "+stid);
				sttable[stid] = new double[st.ap.size()][st.ap.size()][dmax]; 
				ArrayList<TargetNode> subgraph = buildSubGraph(st.nodes);
				HashMap<Integer, TargetNode> subtargetmaps = buildSubTargetMaps(subgraph);
				printNodesWithNeighborsAndPath(subgraph);
				for(int d=0; d<=dmax; d++)
				{
					if(st.nodes.size()>1 && d<dminsuper)
					{
						continue;
					}
					//System.out.print("d = "+d+", ");
					for(TargetNode entry: st.ap.values())
					{
						for(TargetNode exit: st.ap.values())
						{
							//System.out.print("*************Doing DO for ("+entry.getTargetid()+","+exit.getTargetid()+ ", d = "+d+ ")\n");
							/*
							double res[] = SecurityGameContraction.subGraphDOSolver(subgraph, gamedata, st.nodes.size(), 
									1, d, entry.getTargetid(), exit.getTargetid());


							if((entry.getTargetid() != exit.getTargetid())  &&  (apspmat[entry.getTargetid()][exit.getTargetid()]>d) )
							{

								System.out.println("@@@@@@@@@@@Invalid d<dmax \n");
								sttable[stid][entry.getTargetid()][exit.getTargetid()][d] = -7777;
							}
							else
							{

							 */
							//if(!(d>0 && (entry.getTargetid() == exit.getTargetid())))
							{

								//d=2;
								
								//for(double r=0; r<=1; r+=0.01)
								
								//TEST 
								if( d==5)
								{
									System.out.println("hello...");
								}
								
								
								
								{
									
									double[] res = SecurityGameContraction.solvingSubgraphBaseline(subgraph, subtargetmaps, d, gamedata, 1, 
											entry.getTargetid(),  exit.getTargetid(), apspmat, apspmap, subgamepaths/* r*/);

									if(res[2] != -1)
									{

										/*res[0]+=10;
										res[1]+=10;
										*/
										//writeLinearityCheck(r, res[0]);

										System.out.println("("+entry.getTargetid()+","+exit.getTargetid()+")="+res[1]+"");

										//sttable[stid][entry.getTargetid()][exit.getTargetid()][d] = res[1];
										String key = entry.getTargetid()+","+exit.getTargetid()+","+d;
										supertable.put(key, res);
										System.out.println(key+","+ res[0] + ","+ res[1] + ","+ res[2]);
										String k = (entry.getTargetid()+","+exit.getTargetid()+","+d);
										
										//TEST
										/*if(subgamemap.containsKey(k))
										{
											*//**
											 * major problem, not unique id ........ Need to solve this
											 *//*
											
											break;
										}*/
										subgamemap.put(k, subgamecount);
										subgamemapback.put(subgamecount, k);
										subgamecount += 1;
										if(st.nodes.size()==1 && d==0)
											break;
									}
								}
								//System.out.println("done");
								
							}
							//}

							//System.out.println("Hi");


						}
						if(st.nodes.size()==1 && d==0)
							break;
					}
					if(st.nodes.size()==1 && d==0)
						break;
					//System.out.println();
				}

			}
		}

		//System.out.println();
		printSuperTable(supertable, sts, dmax);

		return supertable;



	}

	private static HashMap<String, double[]>  buildSuperTable2(double[][][][] sttable, ArrayList<TargetNode> targets,
			HashMap<Integer, TargetNode> targetmaps, HashMap<Integer, SuperTarget> sts, int dmax, int[][] gamedata,
			int[][] apspmat, HashMap<Integer,Integer> apspmap, 
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback,
			HashMap<String,HashMap<ArrayList<Integer>,Double>> subgamepaths, int dminsuper, 
			HashMap<String,Integer> subgamecontainer) throws Exception
	{



		HashMap<String, double[]> supertable = new  HashMap<String, double[]>();

		int subgamecount = 0;
		for(Integer stid: sts.keySet())
		{
			SuperTarget st = sts.get(stid);

			//System.out.println("\nST "+stid);
			sttable[stid] = new double[st.ap.size()][st.ap.size()][dmax]; 
			ArrayList<TargetNode> subgraph = buildSubGraph(st.nodes);
			HashMap<Integer, TargetNode> subtargetmaps = buildSubTargetMaps(subgraph);
			printNodesWithNeighborsAndPath(subgraph);

			for(TargetNode entry: st.ap.values())
			{
				for(TargetNode exit: st.ap.values())
				{

					
					
					
					int [] disallocation = new int[1];

					double[] res = SecurityGameContraction.solvingSubgraphBaseline2(subgraph, subtargetmaps, dmax, gamedata, 1, 
							entry.getTargetid(),  exit.getTargetid(), apspmat, apspmap, subgamepaths/* r*/, disallocation);

					if(res[2] != -1)
					{

						/*res[0]+=10;
										res[1]+=10;
						 */
						//writeLinearityCheck(r, res[0]);

						System.out.println("("+entry.getTargetid()+","+exit.getTargetid()+")="+res[1]+"");

						//sttable[stid][entry.getTargetid()][exit.getTargetid()][d] = res[1];
						String key = entry.getTargetid()+","+exit.getTargetid()+","+disallocation[0];
						supertable.put(key, res);
						System.out.println(key+","+ res[0] + ","+ res[1] + ","+ res[2]);
						String k = (entry.getTargetid()+","+exit.getTargetid()+","+disallocation[0]);


						subgamemap.put(k, subgamecount);
						subgamemapback.put(subgamecount, k);
						subgamecontainer.put(k, stid);
						subgamecount += 1;
						
						String key2 = entry.getTargetid()+","+exit.getTargetid();
						
						//maxdistallocation.put(key2, (double)disallocation[0]);
						
						
						

					}
					if(st.nodes.size()==1)
						break;

				}
				if(st.nodes.size()==1)
					break;


			}

		}

		//System.out.println();
		printSuperTable(supertable, sts, dmax);

		return supertable;



	}
	
	
	private static HashMap<String, double[]>  buildSuperTable3(double[][][][] sttable, ArrayList<TargetNode> targets,
			HashMap<Integer, TargetNode> targetmaps, HashMap<Integer, SuperTarget> sts, int dmax, int[][] gamedata,
			int[][] apspmat, HashMap<Integer,Integer> apspmap, 
			HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback,
			HashMap<String,HashMap<ArrayList<Integer>,Double>> subgamepaths, int dminsuper, HashMap<String,Double> maxdistallocation,
			HashMap<String,Integer> subgamecontainer) throws Exception
	{



		HashMap<String, double[]> supertable = new  HashMap<String, double[]>();

		int subgamecount = 0;
		for(Integer stid: sts.keySet())
		{
			SuperTarget st = sts.get(stid);

			//System.out.println("\nST "+stid);
			sttable[stid] = new double[st.ap.size()][st.ap.size()][dmax]; 
			ArrayList<TargetNode> subgraph = buildSubGraph(st.nodes);
			HashMap<Integer, TargetNode> subtargetmaps = buildSubTargetMaps(subgraph);
			printNodesWithNeighborsAndPath(subgraph);

			for(TargetNode entry: st.ap.values())
			{
				for(TargetNode exit: st.ap.values())
				{

					
					
					
					int [] disallocation = new int[1];

					
					// TODO if entry and exit points are same we need different type of distance allocation where there will 
					// be another target in the path.
					double[] res = SecurityGameContraction.solvingSubgraphBaseline2(subgraph, subtargetmaps, dmax, gamedata, 1, 
							entry.getTargetid(),  exit.getTargetid(), apspmat, apspmap, subgamepaths/* r*/, disallocation);

					if(res[2] != -1)
					{

						/*res[0]+=10;
										res[1]+=10;
						 */
						//writeLinearityCheck(r, res[0]);

						System.out.println("("+entry.getTargetid()+","+exit.getTargetid()+")="+res[1]+"");

						//sttable[stid][entry.getTargetid()][exit.getTargetid()][d] = res[1];
						String key = entry.getTargetid()+","+exit.getTargetid()+","+disallocation[0];
						supertable.put(key, res);
						System.out.println(key+","+ res[0] + ","+ res[1] + ","+ res[2]);
						String k = (entry.getTargetid()+","+exit.getTargetid()+","+disallocation[0]);


						subgamemap.put(k, subgamecount);
						subgamemapback.put(subgamecount, k);
						subgamecontainer.put(k, stid);
						subgamecount += 1;
						
						String key2 = entry.getTargetid()+","+exit.getTargetid();
						
						maxdistallocation.put(key2, (double)disallocation[0]);
						
						
						

					}
					if(st.nodes.size()==1)
						break;

				}
				if(st.nodes.size()==1)
					break;


			}

		}

		//System.out.println();
		printSuperTable(supertable, sts, dmax);

		return supertable;



	}


	private static HashMap<Integer, TargetNode> buildSubTargetMaps(ArrayList<TargetNode> subgraph) {
		
		
		HashMap<Integer, TargetNode> maps = new HashMap<Integer, TargetNode>();
		
		for(TargetNode n: subgraph)
		{
			maps.put(n.getTargetid(), n);
		}
		
		return maps;
	}


	private static void makeSuperPathSeqSrcDest(ArrayList<ArrayList<String>> pathseq
			,ArrayList<SuperTarget> goals) {



		//int Pmat[][] = new int[nTargets][pathcounter];
		//int pathindex = 0;
		for(SuperTarget n: goals)
		{

			ArrayList<String> tmppathseq = new ArrayList<String>();
			builSuperpathSeqq(n, tmppathseq);

			//pathindex++;
			pathseq.add(tmppathseq);
		}

		//return Pmat;





	}


	private static ArrayList<String> builSuperpathSeqq(SuperTarget n, ArrayList<String> tmppathseq) {

		if(n==null)
			return null;
		//System.out.println("parent  "+ n.parent.getTargetid());
		builSuperpathSeqq(n.parent, tmppathseq);
		//FIXIT combid is not unique
		tmppathseq.add(n.combid);
		//System.out.println("Adding "+ n.combid);



		return tmppathseq;
	}


	private static void printSuperTable(HashMap<String, double[]> supertable, HashMap<Integer,SuperTarget> sts, int dmax) {


		for(Integer stid: sts.keySet())
		{
			SuperTarget st = sts.get(stid);
			//if(st.nodes.size()>1)
			{
				System.out.println("\n\n%%%%%%%%%%%%%%  ST "+stid);


				for(int d=0; d<=dmax; d++)
				{
					//System.out.print("d = "+d+", ");
					for(TargetNode entry: st.ap.values())
					{
						for(TargetNode exit: st.ap.values())
						{

							//sttable[stid][entry.getTargetid()][exit.getTargetid()][d] = res[1];
							String key = entry.getTargetid()+","+exit.getTargetid()+","+d;
							//String key1 = entry.getTargetid()+"-"+exit.getTargetid()+"-"+d;

							if(supertable.containsKey(key))
							{
								System.out.println(key+" =("+supertable.get(key)[0]+","+supertable.get(key)[1]+") ") ;
							}



						}
					}
					//System.out.println();
				}

			}
		}


	}




	private static ArrayList<TargetNode> buildSubGraph(HashMap<Integer, TargetNode> src) {


		///
		ArrayList<TargetNode> duplicategraph = new ArrayList<TargetNode>();

		for(TargetNode s: src.values())
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
		for(TargetNode s: src.values())
		{
			TargetNode t = duplicategraph.get(tindex);
			// add all neighbors

			for(TargetNode nei: s.getNeighbors())
			{


				//if nei is inside the supertarget
				if(src.containsValue(nei))
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
			}

			tindex++;
		}








		return duplicategraph;
	}


	public static ArrayList<SuperTarget> generatePaths(double dmax, HashMap<Integer, SuperTarget> sts, int base, int dest,
			HashMap<String,double[]> supertable, int dmaxsuper) throws Exception
	{


		ArrayList<SuperTarget> path = getPathsByBFS(dmax,sts, base, dest, supertable, dmaxsuper);
		//printPathMat(pathsmat);

		return path;

	}

	private static ArrayList<SuperTarget> getPathsByBFS(double dmax, HashMap<Integer, SuperTarget> sts,
			int base, int dest, HashMap<String,double[]> supertable, int dmaxsuper) throws Exception 
	{
		//TargetNode start = SecurityGameContraction.targets.get(0);
		int nTargets = sts.size(); 
		SuperTarget start = new SuperTarget(sts.get(base), sts.get(base).nodes.get(base), sts.get(base).nodes.get(base)); // add all the start nodes possible to the queue,
		/**
		 * only one node is there
		 * so it will be both entry and exit
		 * and there is only one start node
		 * not multiple, unless there are more than one node in the base-st
		 */

		start.combid = "0,0,0";
		Queue<SuperTarget> fringequeue = new LinkedList<SuperTarget>();
		ArrayList<SuperTarget> goals = new ArrayList<SuperTarget>();
		//System.out.println("start node "+ start.getTargetid());
		fringequeue.add(start);
		int pathcounter = 0;
		while(fringequeue.size()>0)
		{
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			SuperTarget node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.stid==start.stid) &&
					(node.distcoveredyet>0) && node.distcoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
				//System.out.println();
				//	SuperTarget.printPath(node);
				pathcounter++;
				//System.out.println();
				goals.add(node);

				/*if(pathcounter>5000)
					break;*/
			}
			if(node.distcoveredyet<dmax)
			{
				/**
				 * expand the node
				 */
				//System.out.println("Expanding node "+ node.getTargetid());
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
				ArrayList<SuperTarget> succs = Expand(node, dmax, sts, base, dest, supertable, dmaxsuper);
				/**
				 * add nodes to queue
				 */
				for(SuperTarget suc: succs)
				{
					//System.out.println("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue");
					//Logger.logit("Adding node "+ suc.getTargetid() +", distance covered  "+suc.distancecoveredyet+", to queue"+"\n");
					fringequeue.add(suc);

				}
				//System.out.println("Queue size after adding "+ fringequeue.size());

			}

		}
		return goals;

	}

	private static ArrayList<SuperTarget> Expand(SuperTarget node, double dmax, HashMap<Integer, SuperTarget> sts, 
			int base, int dest, HashMap<String,double[]> supertable, int dmaxsuper) 
	{
		ArrayList<SuperTarget> successors = new ArrayList<SuperTarget>();


		SuperTarget tmpnode = sts.get(node.stid);


		for(SuperTarget nei: tmpnode.neighbors.values())
		{

			// for every neighbor 
			/**
			 * for every edge : node.exit -> (nei.entry, nei.exit) 
			 */
			//TargetNode exitpoint = null;
			for(TargetNode srcexit: node.exitpoints.values()) // should be only one exit
			{
				// find targetnode neibor in supertarget nei
				for(TargetNode neientry : nei.ap.values())
				{
					if(neientry.getNeighbors().contains(srcexit)) // if there is an edge src.exit -> nei.entry
					{
						// for every exit point in nei. create a super target with nei.entry
						for(TargetNode neiexit: nei.ap.values()) 
						{
							// for every distance create a supertarget
							for(int d=0; d<=dmaxsuper; d++)
							{
								

								if((node.distcoveredyet + d +srcexit.getDistance(neientry)) <= dmax)
								{
									/**
									 * create new supertarget 
									 * allocate distance
									 * assign nodes
									 * update entry exit point
									 * assign parent
									 * assign distancecveredyet
									 * 
									 */
									String key = (neientry.getTargetid()+","+neiexit.getTargetid()+","+d);
									//System.out.println("Key&&&&&& "+ key);
									if(supertable.containsKey(key))
									{

										if(!(supertable.get(key)[1]==777 || supertable.get(key)[1]==-777) )
										{

											SuperTarget newst = new SuperTarget(nei, neientry, neiexit);
											newst.stid = nei.stid;
											newst.allocateddistance = d;
											newst.parent = node;
											newst.distcoveredyet = node.distcoveredyet + d + srcexit.getDistance(neientry);
											newst.attackerreward = supertable.get(key)[1];
											newst.attackerpenalty = 0;
											newst.defenderreward = 0;
											newst.defenderpenalty = -newst.attackerreward;
											// can have duplicates combid, make it string with comma
											newst.combid = /*Integer.parseInt*/(neientry.getTargetid()+","+neiexit.getTargetid()+","+d);
											successors.add(newst);
											if(newst.stid==dest)
											{
												break;
											}
										}
									}

								}

							}
						}
					}

				}

			}


		}
		return successors;

	}



	
	private static void printSuperTarget(SuperTarget sts) {



		//Iterator<SuperTarget> itr = sts.values().iterator();

		//for(SuperTarget st : itr)
		//while(itr.hasNext())
		//{
			SuperTarget node = sts;
			System.out.println("\n\n******Super target node " + node.stid+"******");
			//Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				// print the nodes
				System.out.print("---Nodes : ");
				for(TargetNode n: node.nodes.values())
				{
					System.out.print( n.getTargetid()+ " ");

				}
				System.out.print("\n---Neighbors : ");
				for(SuperTarget neighbor: node.neighbors.values())
				{
					System.out.print(neighbor.stid+ " ");

				}
				System.out.print("\n---AP : ");
				for(TargetNode a: node.ap.values())
				{
					System.out.print(a.getTargetid()+ " ");

				}
			}
		//}

	}

	
	
	

	private static void printSuperTargets(HashMap<Integer, SuperTarget> sts) {



		Iterator<SuperTarget> itr = sts.values().iterator();

		//for(SuperTarget st : itr)
		while(itr.hasNext())
		{
			SuperTarget node = itr.next();
			System.out.println("\n\n******Super target node " + node.stid+"******");
			//Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				// print the nodes
				System.out.print("---Nodes : ");
				for(TargetNode n: node.nodes.values())
				{
					System.out.print( n.getTargetid()+ " ");

				}
				System.out.print("\n---Neighbors : ");
				for(SuperTarget neighbor: node.neighbors.values())
				{
					System.out.print(neighbor.stid+ " ");

				}
				ArrayList<TargetNode> already = new ArrayList<TargetNode>();
				System.out.print("\n---Neighbor Targets : ");
				for(TargetNode neighbor: node.ap.values())
				{
					for(TargetNode nei: neighbor.getNeighbors())
					{
						if(!already.contains(nei) && !node.nodes.values().contains(nei))
						{
							System.out.print(nei.getTargetid()+ " ");
							already.add(nei);
						}
					}

				}
				System.out.print("\n---AP : ");
				for(TargetNode a: node.ap.values())
				{
					System.out.print(a.getTargetid()+ " ");

				}
			}
		}

	}



	private static ArrayList<TargetNode> createGraph() {


		ArrayList<TargetNode> graph = new ArrayList<TargetNode>();


		TargetNode base = new TargetNode(0, 10);

		TargetNode a = new TargetNode(1, 10);
		TargetNode b = new TargetNode(2, 9);
		TargetNode c = new TargetNode(3, 9);







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



		// base ->a
		base.addNeighbor(a);
		base.addDistance(a, 5.0);
		base.setPath(a, path);

		a.addNeighbor(base);
		a.addDistance(base, 5.0);
		a.setPath(base, path);




		// base ->b
		base.addNeighbor(b);
		base.addDistance(b, 6.0);
		base.setPath(b, path);

		b.addNeighbor(base);
		b.addDistance(base, 6.0);
		b.setPath(base, path);









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


		/*

		gamedata[1][0] = 0;//(int)density[iter][i] ;
		gamedata[1][1] = -9;
		gamedata[1][2] = (int)9;  // uncovered
		gamedata[1][3] = 0; //covered


		gamedata[2][0] = 0;//(int)density[iter][i] ;
		gamedata[2][1] = -9;
		gamedata[2][2] = (int)9;  // uncovered
		gamedata[2][3] = 0; //covered
		 */


		graph.add(base);
		graph.add(a);
		graph.add(b);
		graph.add(c);



		TargetNode d = new TargetNode(4, 5);
		TargetNode e = new TargetNode(5, 4);
		TargetNode f = new TargetNode(6, 5);
		TargetNode g = new TargetNode(7, 4);


		//d -> e
		// path = new ArrayList<TargetNode>();

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




		// base ->d
		base.addNeighbor(d);
		base.addDistance(d, 7.0);
		base.setPath(d, path);

		d.addNeighbor(base);
		d.addDistance(base, 7.0);
		d.setPath(base, path);



		// base ->g
		base.addNeighbor(g);
		base.addDistance(g, 6.0);
		base.setPath(g, path);

		g.addNeighbor(base);
		g.addDistance(base, 6.0);
		g.setPath(base, path);



		// b ->d
		b.addNeighbor(d);
		b.addDistance(d, 3.0);
		b.setPath(d, path);

		d.addNeighbor(b);
		d.addDistance(b, 3.0);
		d.setPath(b, path);

		// c ->e
		c.addNeighbor(e);
		c.addDistance(e, 4.0);
		c.setPath(e, path);

		e.addNeighbor(c);
		e.addDistance(c, 4.0);
		e.setPath(c, path);




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
	
	
	public static ArrayList<Integer>[] createGraph3(ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps) 
	{


		//ArrayList<TargetNode> graph = new ArrayList<TargetNode>();


		TargetNode base = new TargetNode(0, 10);

		TargetNode a = new TargetNode(1, 7);
		TargetNode b = new TargetNode(2, 7);
		TargetNode c = new TargetNode(3, 7);
		TargetNode d = new TargetNode(4, 7);
		TargetNode e = new TargetNode(5, 7);
		
		
		
		TargetNode f = new TargetNode(6, 2);
		TargetNode g = new TargetNode(7, 1);
		TargetNode h = new TargetNode(8, 1);
		TargetNode i = new TargetNode(9, 0);
		TargetNode j = new TargetNode(10, 0);
		
		
		
		
		
		







		//a -> b
		ArrayList<TargetNode> path = new ArrayList<TargetNode>();

		a.addNeighbor(b);
		a.addDistance(b, 1.0);
		a.setPath(b, path);
		b.addNeighbor(a);
		b.addDistance(a, 1.0);
		b.setPath(a, path);
		
		e.addNeighbor(b);
		e.addDistance(b, 1.0);
		e.setPath(b, path);
		b.addNeighbor(e);
		b.addDistance(e, 1.0);
		b.setPath(e, path);

		//a->c
		a.addNeighbor(c);
		a.addDistance(c, 2.0);
		a.setPath(c, path);
		c.addNeighbor(a);
		c.addDistance(a, 2.0);
		c.setPath(a, path);
		// c->e
		c.addNeighbor(e);
		c.addDistance(e, 2.0);
		c.setPath(e, path);
		e.addNeighbor(c);
		e.addDistance(c, 2.0);
		e.setPath(c, path);
		
		
		
		//a->d
		a.addNeighbor(d);
		a.addDistance(d, 3.0);
		a.setPath(d, path);
		d.addNeighbor(a);
		d.addDistance(a, 3.0);
		d.setPath(a, path);
		
		
		
		
		
		
		// d->e
		d.addNeighbor(e);
		d.addDistance(e, 3.0);
		d.setPath(e, path);
		e.addNeighbor(d);
		e.addDistance(d, 3.0);
		e.setPath(d, path);




		// base ->a
		base.addNeighbor(a);
		base.addDistance(a, 5.0);
		base.setPath(a, path);
		a.addNeighbor(base);
		a.addDistance(base, 5.0);
		a.setPath(base, path);




		
		base.defenderreward=0;
		base.defenderpenalty= -10;
		base.attackerreward = 10;
		base.attackerpenalty= 0;



		a.defenderreward=0;
		a.defenderpenalty= -7;
		a.attackerreward = 7;
		a.attackerpenalty= 0;

		b.defenderreward=0;
		b.defenderpenalty= -7;
		b.attackerreward = 7;
		b.attackerpenalty= 0;

		c.defenderreward=0;
		c.defenderpenalty= -7;
		c.attackerreward = 7;
		c.attackerpenalty= 0;


		
		
		
		
		
		
		
		graph.add(base);
		graph.add(a);
		graph.add(b);
		graph.add(c);

		//d -> e
		// path = new ArrayList<TargetNode>();

		d.addNeighbor(e);
		d.addDistance(e, 1.0);
		d.setPath(e, path);

		e.addNeighbor(d);
		e.addDistance(d, 1.0);
		e.setPath(d, path);

		//d->g

		/*d.addNeighbor(g);
		d.addDistance(g, 1.0);
		d.setPath(g, path);

		g.addNeighbor(d);
		g.addDistance(d, 1.0);
		g.setPath(d, path);
		 */

		// i->f
		i.addNeighbor(f);
		i.addDistance(f, 1.0);
		i.setPath(f, path);
		f.addNeighbor(i);
		f.addDistance(i, 1.0);
		f.setPath(i, path);

		
		i.addNeighbor(j);
		i.addDistance(j, 1.0);
		i.setPath(j, path);
		j.addNeighbor(i);
		j.addDistance(i, 1.0);
		j.setPath(i, path);

		
		
		h.addNeighbor(f);
		h.addDistance(f, 2.0);
		h.setPath(f, path);
		f.addNeighbor(h);
		f.addDistance(h, 2.0);
		f.setPath(h, path);

		
		h.addNeighbor(j);
		h.addDistance(j, 2.0);
		h.setPath(j, path);
		j.addNeighbor(h);
		j.addDistance(h, 2.0);
		j.setPath(h, path);

		
		
		g.addNeighbor(f);
		g.addDistance(f, 3.0);
		g.setPath(f, path);
		f.addNeighbor(g);
		f.addDistance(g, 3.0);
		f.setPath(g, path);

		
		g.addNeighbor(j);
		g.addDistance(j, 3.0);
		g.setPath(j, path);
		j.addNeighbor(g);
		j.addDistance(g, 3.0);
		j.setPath(g, path);


		



		// base ->d
		base.addNeighbor(j);
		base.addDistance(j, 7.0);
		base.setPath(j, path);

		j.addNeighbor(base);
		j.addDistance(base, 7.0);
		j.setPath(base, path);




		// c ->e
		f.addNeighbor(e);
		f.addDistance(e, 4.0);
		f.setPath(e, path);

		e.addNeighbor(f);
		e.addDistance(f, 4.0);
		e.setPath(f, path);
		
		
		
		i.addNeighbor(d);
		i.addDistance(d, 6.0);
		i.setPath(d, path);
		d.addNeighbor(i);
		d.addDistance(i, 6.0);
		d.setPath(i, path);
		
		
		
		//a->j
		a.addNeighbor(j);
		a.addDistance(j, 8.0);
		a.setPath(j, path);
		j.addNeighbor(a);
		j.addDistance(a, 8.0);
		j.setPath(a, path);




		d.defenderreward=0;
		d.defenderpenalty= -1;
		d.attackerreward = 1;
		d.attackerpenalty= 0;

		e.defenderreward=0;
		e.defenderpenalty= -1;
		e.attackerreward = 1;
		e.attackerpenalty= 0;
		
		
		
		
		

		f.defenderreward=0;
		f.defenderpenalty= -5;
		f.attackerreward = 5;
		f.attackerpenalty= 0;

		g.defenderreward=0;
		g.defenderpenalty= -4;
		g.attackerreward = 4;
		g.attackerpenalty= 0;
		
		h.defenderreward=0;
		h.defenderpenalty= -4;
		h.attackerreward = 4;
		h.attackerpenalty= 0;

		i.defenderreward=0;
		i.defenderpenalty= -1;
		i.attackerreward = 1;
		i.attackerpenalty= 0;

		j.defenderreward=0;
		j.defenderpenalty= -1;
		j.attackerreward = 1;
		j.attackerpenalty= 0;
		 


		graph.add(d);
		graph.add(e);
		graph.add(f);
		graph.add(g);
		graph.add(h);
		graph.add(i);
		graph.add(j);
		
		//graph.add(g);


		for(TargetNode n: graph)
		{
			targetmaps.put(n.getTargetid(), n);
		}
		
		
		ArrayList<Integer>[] cluster = (ArrayList<Integer>[])new ArrayList[3];
		for(int t=0; t<3; t++)
		{
			cluster[t] = new ArrayList<Integer>();
		}
		
		cluster[0].add(0);
		
		cluster[1].add(1);
		cluster[1].add(2);
		cluster[1].add(3);
		cluster[1].add(4);
		cluster[1].add(5);
		
		
		cluster[2].add(6);
		cluster[2].add(7);
		cluster[2].add(8);
		cluster[2].add(9);
		cluster[2].add(10);
		


		return cluster;


	}






	public static ArrayList<TargetNode> createGraph2(ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps) {


		//ArrayList<TargetNode> graph = new ArrayList<TargetNode>();


		TargetNode base = new TargetNode(0, 10);

		TargetNode a = new TargetNode(1, 7);
		TargetNode b = new TargetNode(2, 7);
		TargetNode c = new TargetNode(3, 7);







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
		a.addDistance(c, 1.0);
		a.setPath(c, path);

		c.addNeighbor(a);
		c.addDistance(a, 1.0);
		c.setPath(a, path);


		// c->b
		c.addNeighbor(b);
		c.addDistance(b, 1.0);
		c.setPath(b, path);

		b.addNeighbor(c);
		b.addDistance(c, 1.0);
		b.setPath(c, path);



		// base ->a
		base.addNeighbor(a);
		base.addDistance(a, 5.0);
		base.setPath(a, path);

		a.addNeighbor(base);
		a.addDistance(base, 5.0);
		a.setPath(base, path);




		// base ->b
		/*base.addNeighbor(b);
		base.addDistance(b, 6.0);
		base.setPath(b, path);

		b.addNeighbor(base);
		b.addDistance(base, 6.0);
		b.setPath(base, path);
		 */








		a.defenderreward=0;
		a.defenderpenalty= -7;
		a.attackerreward = 7;
		a.attackerpenalty= 0;

		b.defenderreward=0;
		b.defenderpenalty= -4;
		b.attackerreward = 4;
		b.attackerpenalty= 0;

		c.defenderreward=0;
		c.defenderpenalty= -5;
		c.attackerreward = 5;
		c.attackerpenalty= 0;


		/*

		gamedata[1][0] = 0;//(int)density[iter][i] ;
		gamedata[1][1] = -9;
		gamedata[1][2] = (int)9;  // uncovered
		gamedata[1][3] = 0; //covered


		gamedata[2][0] = 0;//(int)density[iter][i] ;
		gamedata[2][1] = -9;
		gamedata[2][2] = (int)9;  // uncovered
		gamedata[2][3] = 0; //covered
		 */


		graph.add(base);
		graph.add(a);
		graph.add(b);
		graph.add(c);



		TargetNode d = new TargetNode(4, 5);
		TargetNode e = new TargetNode(5, 5);
		TargetNode f = new TargetNode(6, 5);
		//TargetNode g = new TargetNode(7, 5);


		//d -> e
		// path = new ArrayList<TargetNode>();

		d.addNeighbor(e);
		d.addDistance(e, 1.0);
		d.setPath(e, path);

		e.addNeighbor(d);
		e.addDistance(d, 1.0);
		e.setPath(d, path);

		//d->g

		/*d.addNeighbor(g);
		d.addDistance(g, 1.0);
		d.setPath(g, path);

		g.addNeighbor(d);
		g.addDistance(d, 1.0);
		g.setPath(d, path);
		 */

		// d->f
		d.addNeighbor(f);
		d.addDistance(f, 1.0);
		d.setPath(f, path);

		f.addNeighbor(d);
		f.addDistance(d, 1.0);
		f.setPath(d, path);


		// g->f
		/*g.addNeighbor(f);
		g.addDistance(f, 1.0);
		g.setPath(f, path);

		f.addNeighbor(g);
		f.addDistance(g, 1.0);
		f.setPath(g, path);
		 */

		// e->f
		e.addNeighbor(f);
		e.addDistance(f, 1.0);
		e.setPath(f, path);

		f.addNeighbor(e);
		f.addDistance(e, 1.0);
		f.setPath(e, path);




		// base ->d
		base.addNeighbor(d);
		base.addDistance(d, 7.0);
		base.setPath(d, path);

		d.addNeighbor(base);
		d.addDistance(base, 7.0);
		d.setPath(base, path);



		/*// base ->g
		base.addNeighbor(g);
		base.addDistance(g, 6.0);
		base.setPath(g, path);

		g.addNeighbor(base);
		g.addDistance(base, 6.0);
		g.setPath(base, path);
		 */


		// b ->d
		a.addNeighbor(d);
		a.addDistance(d, 3.0);
		a.setPath(d, path);

		d.addNeighbor(a);
		d.addDistance(a, 3.0);
		d.setPath(a, path);

		/*// c ->e
		c.addNeighbor(e);
		c.addDistance(e, 4.0);
		c.setPath(e, path);

		e.addNeighbor(c);
		e.addDistance(c, 4.0);
		e.setPath(c, path);
		 */



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

		/*g.defenderreward=0;
		g.defenderpenalty= -4;
		g.attackerreward = 4;
		g.attackerpenalty= 0;
		 */


		graph.add(d);
		graph.add(e);
		graph.add(f);
		//graph.add(g);


		for(TargetNode n: graph)
		{
			targetmaps.put(n.getTargetid(), n);
		}


		return graph;


	}






	public static ArrayList<Integer>[] groupTargets(int base, int dest, HashMap<Integer,TargetNode> targetmaps, 
			int k, int radius, ArrayList<TargetNode> targets, int[][] apspmat,  HashMap<Integer,Integer> map)
	{
		ArrayList<Integer>[] cluster = (ArrayList<Integer>[])new ArrayList[k];
		int centers[] = new int[k];
		int prevcenter[] = new int[k];
		for(int i=0; i<k; i++)
		{
			cluster[i] = new ArrayList<Integer>();
		}

		// pick center randomly which are not neighbor and more than 2*radius away





		centers = pickCenter(targetmaps, radius,k,base,dest, targets, apspmat, map);
		//printCenters(centers);

		//System.out.println();


		int iter = 0;






		while(true)
		{
			System.out.println("Iter "+ iter+ "\n");

			// save centers to prev center

			for(int i=0; i<k; i++)
			{
				prevcenter[i] = centers[i];
			}

			// clear the clusters
			for(int i=0; i<k; i++)
			{
				cluster[i].clear();
			}


			// for every target find the distance from every center


			for(TargetNode t: targets)
			{
				//if base or dest assign that to 

				if(t.getTargetid() == base )
				{
					// assign to cluster 0
					cluster[0].add(t.getTargetid());
				}
				if(t.getTargetid() == dest && base != dest)
				{
					//assign to cluster 1
					cluster[1].add(t.getTargetid());
				}

				if(t.getTargetid() != base || t.getTargetid() != dest)
				{
					double dist [][] = new double[k][2];
					int distindex = 0;
					//System.out.println("Target "+ t.getTargetid());
					for(Integer center: centers)
					{
						TargetNode centernode = targetmaps.get(center);
						if(centernode.getTargetid()==base || centernode.getTargetid()==dest)
						{
							dist[distindex][0] = Double.MAX_VALUE;
							dist[distindex][1] = centernode.getTargetid();
							distindex++;
						}
						else if(t.getTargetid()==centernode.getTargetid())
						{
							// dsitance between center and targetnode t
							dist[distindex][0] = 0;
							dist[distindex][1] = centernode.getTargetid();
							distindex++;
						}
						else
						{
							// dsitance between center and targetnode t
							dist[distindex][0] = apspmat[map.get(t.getTargetid())][map.get(center)];
							dist[distindex][1] = centernode.getTargetid();
							distindex++;
						}
					}
					// find the min distance and the center
					double mindist = Double.MAX_VALUE;
					int minindex = -1;
					for(int j=0; j<k; j++)
					{
						//System.out.println("dist("+ t.getTargetid()+","+dist[j][1]+") = "+ dist[j][0]);
						if(mindist>dist[j][0])
						{
							mindist = dist[j][0];
							minindex = j;
						}
					}

					//System.out.println("Min distance : "+ dist[minindex][0]+ "\nclosest center "+ dist[minindex][1] + ", cluster "+ minindex);
					// assign the target to the closest center
					cluster[minindex].add(t.getTargetid());
					//System.out.println("Target : "+ t.getTargetid()+ "assigned to cluster "+ minindex + "\n\n");


				}

			}
			//print clusters

			//printClusters(cluster);


			// compute the center

			computeCenter(cluster, centers, targetmaps, apspmat);
			//System.out.println("New centers\n");
			//	printCenters(centers);

			// now check with prev center. If they aren't same. then move to next itearation. otherwise stop

			boolean stop = checkIfStop(centers, prevcenter);

			if(stop)
			{
				//System.out.println("New centers are same as prevcenter....stopping....\n");
				break;
			}

			//System.out.println("New centers are different than prevcenter....continuing....\n");

			//System.out.println();
			iter++;
			if(iter==100)
				break;
		}


		return cluster;
	}






	private static boolean checkIfStop(int[] centers, int[] prevcenter) {



		for(int i=0; i<centers.length; i++)
		{
			// check if the center can be found
			boolean found = false;
			for(int j=0; j<prevcenter.length; j++)
			{
				if(centers[i]==prevcenter[j])
				{
					found=true;
					break;
				}
			}
			if(!found)
				return false;
		}


		return true;
	}




	private static void computeCenter(ArrayList<Integer>[] clusters, int[] centers,
			HashMap<Integer, TargetNode> targetmaps, int[][] apspmat) {

		// for every cluster
		int index = 0; // cluster index
		for(ArrayList<Integer> cluster: clusters)
		{

			if(cluster.size()>1)
			{

				//System.out.println("cluster "+ index);
				// compute farthest distance
				int farthestdist = Integer.MIN_VALUE;
				int fx=-1;
				int fy=-1;
				for(int x: cluster)
				{
					for(int y: cluster)
					{

						if(x!=y && (farthestdist<apspmat[x][y]))
						{
							farthestdist = apspmat[x][y];
							fx=x;
							fy=y;
						}
					}
				}
				//System.out.println("farthest dist "+ farthestdist);
				//System.out.println("farthest target x "+ fx);
				//System.out.println("farthest target y "+ fy);
				// now find the closest  target to farthestdist/2
				// find z which is not x and not y
				// where d(x,z)-d(x,y) smallest

				double halffd = farthestdist/2.0;

				int closestdist = Integer.MAX_VALUE;
				int closestz = -1;
				for(int z: cluster)
				{

					if(z!=fx  && z!= fy)
					{
						//System.out.println("target z = "+ z);
						int dzx = apspmat[z][fx];
						int dzy = apspmat[z][fy];

						//System.out.println("dist(z,fx) "+z +", "+ fx+" = "+ dzx);
						//System.out.println("dist(z,fy) "+z +", "+ fy+" = "+ dzy);

						if((dzx-dzy)<closestdist)
						{
							closestdist = Math.abs(dzx-dzy);
							closestz = z;
						}
					}

				}

				//	System.out.println("closest target z "+ closestz + ", dist diff = "+ closestdist);
				// assign the z as center
				centers[index] = closestz;

			}
			index++;
		}

	}




	private static void printClusters(ArrayList<Integer>[] cluster) {


		for(int i=0; i<cluster.length; i++)
		{
			System.out.print("cluster "+ i + ": ");
			for(Integer t: cluster[i])
			{
				System.out.print(t+ ", ");
			}
			System.out.println();
		}

	}




	private static int[][] buildAPSP(int nTargets, HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback,
			ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps)
	{


		int[][] adjacencymatrix = new int[nTargets+1][nTargets+1];



		int icount =1;
		for(int i=0; i<targetmaps.size(); i++)
		{
			map.put(targetmaps.get(i).getTargetid(), icount);
			//System.out.println("Target "+ targets.get(i).getTargetid() +" --> "+icount);
			mapback.put(icount, targetmaps.get(i).getTargetid());
			icount++;
		}


		makeAdjacencyMatrix(adjacencymatrix, targets, nTargets, map, mapback);
		// make the adjacency matrix

		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		int[][] apspmat =  apsp.allPairShortestPath(adjacencymatrix);
		purifyAPSPMatrixZero(apspmat, targets, nTargets, map, mapback);



		return apspmat;
	}




	public static void printCenters(int[] centers) {

		int i=0;
		for(Integer c: centers)
		{
			System.out.println(" Center "+ i + " : target " + c);
			i++;
		}

	}

	/**
	 * 
	 * @param targetmaps list of targets
	 * @param radius radius limit for cluster
	 * @param k how many cluster
	 * @param base base id
	 * @param dest destination id
	 * @return
	 */
	public static int[] pickCenter(HashMap<Integer,TargetNode> targetmaps, int radius, int k, int base,
			int dest, ArrayList<TargetNode> targets, int[][] apspmat, HashMap<Integer, Integer> map) {



		int centers[] = new int[k];
		int nTargets = targetmaps.size();
		int startingcenterindex = 0; // from which index to start assigning
		ArrayList<Integer> donecenters = new ArrayList<Integer>();


		// assign the base to first center
		centers[0] = base;
		donecenters.add(0);
		startingcenterindex++;

		//if base != dest assign dest to next center, increment the starting index. 
		if(base != dest)
		{
			centers[1] = dest;
			donecenters.add(1);
			startingcenterindex++;
		}

		//ArrayList<Integer> centeroptions = new ArrayList<Integer>();

		//	System.out.println("starting center index "+ startingcenterindex);


		for(TargetNode n: targets)
		{
			if(n.getTargetid() != base && n.getTargetid()!= dest)
			{

				int count = startingcenterindex;
				//	System.out.println("Considering target "+ n.getTargetid());

				while(true)
				{
					// pick a random index for center from centeroptions
					int centerindex = count;
					//System.out.println("Considering center index "+ centerindex);
					//System.out.println("random center picked "+ centeroptions.get(centerindex));

					// now check if it's more than 2*radius distance away from other centers
					boolean ok = checkIfOk(centers, apspmat, n.getTargetid(), donecenters, 
							targetmaps, radius, base, dest, map, targets);
					if(ok)
					{

						//System.out.println("Assigning target"+ n.getTargetid() + " in center "+ centerindex);
						donecenters.add(centerindex);
						centers[centerindex] = n.getTargetid();
						//centeroptions.remove(centerindex);
						startingcenterindex++;
						break;
					}
					count++;
					if(count==k)
						break;

				}
				if(startingcenterindex==k)
				{
					//System.out.println("Every center has been assigned \n startingcenterindex  "+ startingcenterindex);
					break;
				}

			}

		}


		return centers;
	}

	public static boolean checkIfOk(int[] centers, int[][] apspmat, Integer candidatecenter, ArrayList<Integer> 
	donecenters, HashMap<Integer, TargetNode> targetmaps,
	int radius, int base, int dest, HashMap<Integer,Integer> map, ArrayList<TargetNode> targets) {




		for(int c: donecenters)
		{
			if(c!=base && c!= dest)
			{


				/*TargetNode tmp = targetmaps.get(candidatecenter);
				TargetNode tmp2 = targetmaps.get(c);

				if(tmp2.getNeighbors().contains(tmp))
				{
					if(tmp2.getDistance(tmp)<=radius)
						return false;
				}*/

				//System.out.println("cendidate center "+ candidatecenter+ " is "+ apspmat[map.get(candidatecenter)][map.get(c)] + " away from center target "+ c);
				if(apspmat[map.get(candidatecenter)][map.get(c)]<=radius)
				{
					return false;

				}
			}
		}


		return true;
	}


	private static void purifyAPSPMatrixZero(int[][] adjacencymatrix,
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
					//adjacencymatrix[source][destination] = INFINITY;
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



	private static void makeAdjacencyMatrix(int[][] adjacencymatrix,
			ArrayList<TargetNode> targetmaps, int nTargets, HashMap<Integer,Integer> map, HashMap<Integer,Integer> mapback) {


		int i=1; 
		for(TargetNode n: targetmaps)
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




	public static void printNodesWithNeighborsAndPath(  ArrayList<TargetNode> targets) 
	{
		for(TargetNode node : targets)
		{
			System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");
			Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				for(TargetNode neighbor: node.getNeighbors())
				{
					System.out.println("---Neighbor : "+ neighbor.getTargetid());
					//Logger.logit("---Neighbor : "+ neighbor.getTargetid()+"\n");
					/**
					 * print path
					 */
					ArrayList<TargetNode> path = node.getPath(neighbor);
					System.out.print("Path : "+ node.getTargetid()+ " --> ");
					//Logger.logit("Path : "+ node.getTargetid()+ " --> ");
					for(TargetNode pathnode : path)
					{
						System.out.print(pathnode.getTargetid()+" --> ");
						Logger.logit(pathnode.getTargetid()+" --> ");
					}

					System.out.print(neighbor.getTargetid()+ "\n");
					System.out.println("Distance : " + node.getDistance(neighbor));
					//Logger.logit(neighbor.getTargetid()+ "\n");
					//Logger.logit("Distance : " + node.getDistance(neighbor));

					System.out.print("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");
					//Logger.logit("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");

				}
			}
		}

	}
	
	
	


	public static void printNodesWithNeighborsAndPath(  HashMap<Integer, TargetNode> targets) 
	{
		for(TargetNode node : targets.values())
		{
			System.out.println("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******");
			Logger.logit("\n\n****** target node " + node.getTargetid()+", utility : "+node.getAnimaldensity() +"******\n");
			//if(!domindatednodes.contains(node))
			{

				for(TargetNode neighbor: node.getNeighbors())
				{
					System.out.println("---Neighbor : "+ neighbor.getTargetid());
					//Logger.logit("---Neighbor : "+ neighbor.getTargetid()+"\n");
					/**
					 * print path
					 */
					ArrayList<TargetNode> path = node.getPath(neighbor);
					System.out.print("Path : "+ node.getTargetid()+ " --> ");
					//Logger.logit("Path : "+ node.getTargetid()+ " --> ");
					for(TargetNode pathnode : path)
					{
						System.out.print(pathnode.getTargetid()+" --> ");
						Logger.logit(pathnode.getTargetid()+" --> ");
					}

					System.out.print(neighbor.getTargetid()+ "\n");
					System.out.println("Distance : " + node.getDistance(neighbor));
					//Logger.logit(neighbor.getTargetid()+ "\n");
					//Logger.logit("Distance : " + node.getDistance(neighbor));

					System.out.print("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");
					//Logger.logit("Path utility : "+ node.getPathUtility(neighbor)+"\n\n");

				}
			}
		}

	}


	public static int randInt(int min, int max) {

		// Usually this should be a field rather than a method variable so
		// that it is not re-seeded every call.


		// nextInt is normally exclusive of the top value,
		// so add 1 to make it inclusive
		int randomNum = rand1.nextInt((max - min) + 1) + min;

		return randomNum;
	}



	



}
