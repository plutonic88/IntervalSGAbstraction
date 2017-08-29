package groupingtargets;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
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
import weka.clusterers.EM;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;

public class ClusterTargets {
	
	
	public static Random rand1 = new Random(100);
	public static final int INFINITY = 99999;
	
	
	
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
	
	
	public static void buildcsvGraphExp(int numRow, int numCol, double[][] u,  ArrayList<TargetNode> targets, int iter) throws Exception 
	{
		/**
		 * create the nodes and add to the target list
		 */
		Random rand = new Random(50);
		
		
		
		try {
			
			
			 File f = new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/"+"realdata"+iter+".csv");
			 
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

	
	
	
	public static void wekaClusteringWithDOExp(int nrow, int ncol, int base, int dest, int k, int radius, int dmax, int nRes, int nTargets,
			int LIMIT,int ap, HashMap<Integer,ArrayList<TargetNode>> alltargets, 
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
			//ArrayList<Integer>[] clus = allclus.get(iter);
			ArrayList<TargetNode> targets = alltargets.get(iter);//new ArrayList<TargetNode>();
			HashMap<Integer,TargetNode> targetmaps = alltargetmaps.get(iter); //new HashMap<Integer, TargetNode>();
			
			
			
			
			//printNodesWithNeighborsAndPath(targetmaps);


			Date start = new Date();
			long l1 = start.getTime();

			//ArrayList<Integer>[] clus = makeGraph(k, radius, dlim , nTargets, 2, 10, ap, targets, targetmaps);
			double res[] = wekaClusteringWithDO(base, dest, k, radius, dmax, nRes, nTargets, targets, targetmaps, iter);
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
	
	

	private static boolean isNeighbor(SuperTarget st1, SuperTarget st2) {
		
		
		// if they have nodes in neighbor
		if(st1.neighbors.values().contains(st2))
			return true;
		
		return false;
	}
	
	
	
private static boolean areBothNei(SuperTarget s1, SuperTarget s2, SuperTarget tempst) {
		
		
		if(isNeighbor(s1, tempst) && isNeighbor(s2, tempst))
			return true;
		
		return false;
	}

	
	public static HashMap<Integer, SuperTarget> clusterTargetsWeka(ArrayList<Integer> targetstocluster, 
			ArrayList<TargetNode> graph, HashMap<Integer, TargetNode> targetmaps, double dmax, int k, int radius,
			HashMap<Integer, Double> dstravel, HashMap<Integer,ArrayList<Integer>> stpaths, EM dc,
			Instances instances, HashMap<Integer,Integer> apspmap, int[][] apspmat, AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback) throws Exception
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
		
		ArrayList<Integer>[] clusters = clusterWithWeka(k, newinstance, totalcluster, targetstocluster,targetmaps, dc, apspmap, apspmat);
		
		
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
			int s1ap = -1;
			int s2ap = -1;
			double sda1a2 = -1;
			double sd1 = -1;
			double sd2 = -1;
			ArrayList<Integer> a1a2path = new ArrayList<Integer>();
			ArrayList<Integer> a2path = new ArrayList<Integer>();
			ArrayList<Integer> a1path = new ArrayList<Integer>();
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
											/*System.out.println("\n\n a1 "+ a1.getTargetid() + " , a2 "+ a2.getTargetid() +
													"\n s1 "+ s1.stid + ", s2 "+ s2.stid + 
													"\n tempst "+ tempst.stid);*/


											ArrayList<Integer> tmpa1a2spath = new ArrayList<Integer>();
											ArrayList<Integer> tmpa1path = new ArrayList<Integer>();
											ArrayList<Integer> tmpa2path = new ArrayList<Integer>();
											
											
											double[] d1 = shortestdist(a1,s1, apspmap, apspmat, tmpa1path, apsp, apspmapback); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a1,s1) "+ d1);

											double[] d2 = shortestdist(a2,s2, apspmap, apspmat, tmpa2path, apsp, apspmapback);
											//shortestdist(a2,s2); // should i use bfs for longer path other than near neighbor?
											//TODO
											//System.out.println("shortestdist(a2,s2) "+ d2);

											// next measure the intra cluster shortest traveling path using a1 and a2
											
											
											double dista1a2 = shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);
											
											
											
											
											 
											//shortestdist(a1,a2, tempst, dmax, tmpa1a2spath);

											if(dista1a2 == 0)
											{
												//throw new Exception("No path found to compute AP for st "+ tempst.stid);
												//System.out.println("Disjpint cluster need longer path "+ tempst.stid);
												
												dista1a2 = shortestdist(a1,a2, apspmap, apspmat, tmpa1a2spath, apsp, apspmapback, tempst); 
											}

											//System.out.println("shortestdist(a1,a2, tempst, dmax) "+ dista1a2);

											// if any of the dist is <0 we know that it's not possible to have a path

											if(d1[0] > 0 && d2[0] > 0 && dista1a2 > 0)
											{
												double totaldi = d1[0] + dista1a2 + d2[0];
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
													sd1 = d1[0];
													sd2= d2[0];
													s1ap = (int)d1[1];
													s2ap = (int)d2[1];
													//spath.add(tmpspath.get(0));
													a1a2path.clear();
													for(Integer in: tmpa1a2spath)
													{
														a1a2path.add(in);
													}
													
													
													a1path.clear();
													for(Integer in: tmpa1path)
													{
														a1path.add(in);
													}
													
													
													a2path.clear();
													for(Integer in: tmpa2path)
													{
														a2path.add(in);
													}

													/*System.out.println("Current mindi "+ mindi +  
															"\n a1 "+ aid1 + ", a2 "+ aid2);*/


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
			
			
			
			/*if(tempst.stid==14 && sda1a2==-1)
			{
				System.out.println("shortestdist(a1,s1) ");
				printSuperTarget(tempst);
			}*/

			if((mindi < Double.MAX_VALUE) && (mindi > 0))
			{

				//System.out.println("AP done for st "+ tempst.stid);
				
				dstravel.put(tempst.stid, sda1a2);
				stpaths.put(tempst.stid, a1a2path);
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
				
				/**
				 * add s1 and s2 as neighbor of st if they are not already
				 */
				
				addNei(tempst.stid, aid1, s1id, sd1, a1path, sts, targetmaps, s1ap);
				addNei(tempst.stid, aid2, s2id, sd2, a2path, sts, targetmaps, s2ap);
				addNei(aid1,aid2,sda1a2, a1a2path, tempst, targetmaps);
				
				
				/**
				 * add a1 a2 as neighbor if they are not already
				 */
				
				
				//printSuperTargets(sts);
				//printNodesWithNeighborsAndPath(targetmaps);
				
				
				//System.out.println("hi");
				 
			}
			
		}

		//}
		//updateNeighbors(sts);
		
		//printNodesWithNeighborsAndPath(targetmaps);
		//printSuperTargets(sts);
		
		 return sts;
	}
	

	
	private static void addNei(int aid1, int aid2, double sda1a2, ArrayList<Integer> a1a2path, SuperTarget tempst, HashMap<Integer,TargetNode> targetmaps) {
		
		
		TargetNode srcnode = tempst.nodes.get(aid1);
		TargetNode destnode = tempst.nodes.get(aid2);
		
		
		if(!srcnode.getNeighbors().contains(destnode))
		{
			//TargetNode srcnode = sts.get(stid).nodes.get(aid1);
			//TargetNode destnode = sts.get(s1id).nodes.get(a1path.get(a1path.size()-1));
			
			if(!srcnode.getNeighbors().contains(destnode))
			{
				srcnode.addNeighbor(destnode);
				destnode.addNeighbor(srcnode);
				
				// add path
				
				
				ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
				
				for(Integer n: a1a2path)
				{
					if(n != srcnode.getTargetid() && n!= destnode.getTargetid())
					{
						pathnodes.add(targetmaps.get(n));
					}
				}
				
				srcnode.setPath(destnode, pathnodes);
				srcnode.addDistance(destnode, sda1a2);
				
				
				ArrayList<TargetNode> revpathnodes = new ArrayList<TargetNode>();
				
				for(int i=0; i<a1a2path.size(); i++)
				{
					revpathnodes.add(pathnodes.get(a1a2path.size()-i-1));
				}
				
				destnode.setPath(srcnode, revpathnodes);
				destnode.addDistance(srcnode, sda1a2);
				
				
			}
		}
		
		
		
		
	}


	private static void addNei(int stid, int aid1, int s1id, double sd1, ArrayList<Integer> a1path,
			HashMap<Integer, SuperTarget> sts, HashMap<Integer,TargetNode> targetmaps, int s1ap) {
		
		
		/**
		 * a1......s1....
		 */
		
		// first check if s1 is neibor
		// first check if s1 is neibor
		
		SuperTarget src = sts.get(stid);
		SuperTarget dest = sts.get(s1id);
		
		
		
		if(!src.neighbors.keySet().contains(dest.stid))
		{
			src.neighbors.put(dest.stid,dest);
			dest.neighbors.put(src.stid, src);
			
			// add path
			
			src.path.put(dest, a1path);
			src.distances.put(dest, sd1);
			
			
			ArrayList<Integer> revpath = new ArrayList<Integer>();
			
			for(int i=0; i<a1path.size(); i++)
			{
				revpath.add(a1path.get(a1path.size()-i-1));
			}
			
			dest.path.put(src, revpath);
			dest.distances.put(src, sd1);
			
			
			
			
		}
		

		
		
		// do the same thing for TargetNode
		
		
		
		TargetNode srcnode = sts.get(stid).nodes.get(aid1);
		TargetNode destnode = targetmaps.get(s1ap);
		
		if(!srcnode.getNeighbors().contains(destnode))
		{
			srcnode.addNeighbor(destnode);
			destnode.addNeighbor(srcnode);
			
			// add path
			
			
			ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();
			
			for(Integer n: a1path)
			{
				if(n != srcnode.getTargetid() && n!= destnode.getTargetid())
				{
					pathnodes.add(targetmaps.get(n));
				}
			}
			
			srcnode.setPath(destnode, pathnodes);
			srcnode.addDistance(destnode, sd1);
			
			
			ArrayList<TargetNode> revpathnodes = new ArrayList<TargetNode>();
			
			for(int i=0; i<a1path.size(); i++)
			{
				revpathnodes.add(pathnodes.get(a1path.size()-i-1));
			}
			
			destnode.setPath(srcnode, revpathnodes);
			destnode.addDistance(srcnode, sd1);
			
			
			
			
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

	
	
	
	
	private static double shortestdist(TargetNode a1, TargetNode a2, 
			HashMap<Integer,Integer> apspmap, int[][] apspmat, ArrayList<Integer> tmpa1path, 
			AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback, SuperTarget tempst) {
		// TODO Auto-generated method stub



		double dmin = Double.MAX_VALUE;
		//boolean isnei = false;
		ArrayList<Integer> minpath = new ArrayList<Integer>();
		
		//for(TargetNode t: s1.ap.values())
		//{


			System.out.println("FInding shortest dist for target "+ a2.getTargetid());
			ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
			ArrayList<SuperTarget>  pnodes = new ArrayList<SuperTarget>();
			double distcovered = -1;
			if(a1.getNeighbors().contains(a2))
			{
				//pnodes = base.getPath(dest);
				int src = a1.getTargetid();
				int des = a2.getTargetid();

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];


			}
			else
			{




				int src = a1.getTargetid();
				int des = a2.getTargetid();


				

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];
				System.out.print("dist covered "+ distcovered+"\n");



				ArrayList<Integer>	tmppathnodes = apsp.getPathInST(src, des, apspmap, apspmapback, tempst);

				for(int k=0; k<tmppathnodes.size(); k++)
				{
					pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
				}

				

				//throw new Exception("Base to not neighbor for initial set of paths **********8");

			}
			
			
			if(dmin>distcovered)
			{

				dmin = distcovered;
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				//tmppath.add(a1.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
				//	System.out.print(dest.getTargetid()+"\n");
				//tmppath.add(a2.getTargetid());
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

				System.out.print("\n");
				
				minpath.clear();
				for(Integer p: tmppath)
				{
					minpath.add(p);
				}
			}

		//}
		
		
		
		tmpa1path.clear();
		for(Integer p: minpath)
		{
			tmpa1path.add(p);
		}
		
		
		
		return dmin;
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
			//System.out.println("Polling from queue ");
			//System.out.println("Queue size before polling "+ fringequeue.size());
			TargetNode node = fringequeue.poll();
			//System.out.println("Pulled node "+ node.getTargetid() + ", distance covered "+ node.distancecoveredyet);
			//System.out.println("Queue size after polling "+ fringequeue.size());
			if( (node.getTargetid()==a2.getTargetid()) &&
					(node.distancecoveredyet>0) && node.distancecoveredyet<=dmax)
			{
				//System.out.println("Adding node "+ node.getTargetid() +" to goals..."+ node.distancecoveredyet+ ", pcount:  "+ pathcounter);
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
				//Logger.logit("Expanding node "+ node.getTargetid()+"\n");
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


	

	
	
	private static double[] shortestdist(TargetNode a1, SuperTarget s1, 
			HashMap<Integer,Integer> apspmap, int[][] apspmat, ArrayList<Integer> tmpa1path, 
			AllPairShortestPath apsp, HashMap<Integer,Integer> apspmapback) {
		// TODO Auto-generated method stub



		double dmin = Double.MAX_VALUE;
		//boolean isnei = false;
		ArrayList<Integer> minpath = new ArrayList<Integer>();
		int tid = -1;
		
		for(TargetNode t: s1.ap.values())
		{


			//System.out.println("FInding shortest dist for target "+ s1.stid);
			ArrayList<Integer>  pathnodes = new ArrayList<Integer>();
			ArrayList<SuperTarget>  pnodes = new ArrayList<SuperTarget>();
			double distcovered = -1;
			if(a1.getNeighbors().contains(t))
			{
				//pnodes = base.getPath(dest);
				int src = a1.getTargetid();
				int des = t.getTargetid();

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];


			}
			else
			{




				int src = a1.getTargetid();
				int des = t.getTargetid();


				

				 distcovered = apspmat[apspmap.get(src)][apspmap.get(des)];
				//System.out.print("dist covered "+ distcovered+"\n");



				ArrayList<Integer>	tmppathnodes = apsp.getPath(src, des, apspmap, apspmapback);

				for(int k=0; k<tmppathnodes.size(); k++)
				{
					pathnodes.add(tmppathnodes.get(tmppathnodes.size()-k-1));
				}

				

				//throw new Exception("Base to not neighbor for initial set of paths **********8");

			}
			
			
			if(dmin>distcovered)
			{

				dmin = distcovered;
				tid = t.getTargetid();
				ArrayList<Integer> tmppath = new ArrayList<Integer>();
				//System.out.print("\n0->");
				//tmppath.add(a1.getTargetid());
				for(int k=0; k<pathnodes.size(); k++)
				{
					tmppath.add(pathnodes.get(pathnodes.size()-k-1));
					//System.out.print(pathnodes.get(pathnodes.size()-k-1)+"->");
				}
				//	System.out.print(dest.getTargetid()+"\n");
				//tmppath.add(t.getTargetid());
				//System.out.print("\n");
				/**
				 * make rev path
				 */
				for(int j=tmppath.size()-2; j>=0; j--)
				{
					tmppath.add(tmppath.get(j));
				}
				//System.out.print("complete path : \n");
				for(int k=0; k<tmppath.size(); k++)
				{

					//System.out.print(tmppath.get(k)+"->");
				}

				//System.out.print("\n");
				
				minpath.clear();
				for(Integer p: tmppath)
				{
					minpath.add(p);
				}
			}

		}
		
		
		
		tmpa1path.clear();
		for(Integer p: minpath)
		{
			tmpa1path.add(p);
		}
		
		
		
		return new double[]{dmin, tid};
	}
	
	
	private static ArrayList<Integer>[] clusterWithWeka(int k, Instances newinstance, int totalcluster,
			ArrayList<Integer> targetstocluster, HashMap<Integer, TargetNode> targetmaps, EM dc,
			HashMap<Integer,Integer> apspmap, int[][] apspmat) throws Exception {
		
		
		
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
		 
		// printClusters(clusters);
		 
		 int j = k;
		 for(Integer t: targetmaps.keySet())
		 {
			 if(!targetstocluster.contains(t))
			 {
				 clusters[j++].add(t.intValue());
			 }
		 }
		 

		// printClusters(clusters);
		
		return clusters;
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

	public static double[] wekaClusteringWithDO(int base, int dest, int ncluster, int radius, int dmax, 
			int nRes, int nTargets, ArrayList<TargetNode> targets, HashMap<Integer, TargetNode> targetmaps, int iter) throws Exception
	{
		
		
		
		
		
		
		
		
		int[][] targetssorted = SecurityGameContraction.sortTargets(targets);
		//SecurityGameContraction.printSortedTargets(targetssorted);
		
		
		//Get the list of initial targets using GCR from Tsrt, Tcur = GreedyCoverR()
		
		HashMap<Integer, Integer> apspmap = new HashMap<>();
		HashMap<Integer, Integer> apspmapback = new HashMap<>();
		int[][] apspmat = new int[nTargets+1][nTargets+1];
		AllPairShortestPath apsp = new AllPairShortestPath(nTargets);
		
		ArrayList<Integer> targetstocluster = buildGreedyCoverMultResWeka(targets, dmax, nTargets, 0, nRes,apspmap,
				apspmapback,apspmat, apsp );
		
		int currentPlace = targetstocluster.size()-1;
		
		//ArrayList<TargetNode> domindatednodes = new ArrayList<TargetNode>();

		ArrayList<TargetNode> tmpgraph = new ArrayList<TargetNode>();
		HashMap<Integer, TargetNode> tmptargetmaps = new HashMap<Integer, TargetNode>();
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
		
		
		
		
		

		 CSVLoader csvload = new CSVLoader();
		 csvload.setSource(new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/realdata"+iter+".csv"));
		 Instances data = csvload.getDataSet();
		 
		 
		 ArffSaver arf = new ArffSaver();
		 arf.setInstances(data);
		 
		 File f = new File("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/newdata"+iter+".arff");
		 
		 if(f.exists())
		 {
			 f.delete();
			 f.createNewFile();
		 }
		 
		 arf.setFile(f);
		 arf.writeBatch();
		
		
		
		
		
		
		
		FileInputStream fstream = new FileInputStream("/Users/anjonsunny/Documents/workspace/IntervalSGAbstraction/newdata"+iter+".arff");
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
		 
		/* for(int i=0; i<instances.size(); i++) //without base
		 {
			 System.out.println("instance  "+i +", cluster "+ dc.clusterInstance(instances.get(i)));
		  
		 }
		*/
		

		
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


			tmpgraph = SecurityGameContraction.getDuplicateGraph(targets);
			for(TargetNode t: tmpgraph)
			{
				tmptargetmaps.put(t.getTargetid(), t);
			}
			
			
			
			
			System.out.println();

			//SecurityGameContraction.printNodesWithNeighborsAndPath(domindatednodes, tmpgraph);


			Date start = new Date();
			long l1 = start.getTime();


			//Build an abstract graph Gt through target clustering given Tcur, Gt
			
			
			HashMap<Integer, Double> dstravel = new HashMap<Integer, Double>();
			HashMap<Integer, ArrayList<Integer>> stpaths = new HashMap<Integer, ArrayList<Integer>>();
			
			HashMap<Integer, SuperTarget> currentst = clusterTargetsWeka(targetstocluster, tmpgraph, 
					tmptargetmaps, dmax, ncluster, radius, dstravel, stpaths, dc, instances, apspmap, apspmat, apsp, apspmapback );
			
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
			assignSTValues(currentst, tmptargetmaps);
			
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
			pathseq = SecurityGameContraction.generatePathsForSuperTargetsAPSP(dmax, currentst, tmptargetmaps, targetstocluster, nRes, dstravel);
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
					probdistribution = MIPSolver4.solveForAttackerLPST(p, currentst, tmptargetmaps, nRes, attackerstrategy);

					

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

					origpmat = makeSuperOrigPMatWOMap(p, pathseq, jset, nTargets, map, mapback, 
							tmptargetmaps, currentst, stpaths);
					attackedtarget = SecurityGameContraction.findAttackTarget(origpmat, probdistribution, tmptargetmaps);
					
					//int u = getTargetNode(MIPSolver4.attackedtarget, tmpgraph).getTargetid();
					attackerv = SecurityGameContraction.expectedPayoffAtt(attackedtarget, origpmat, tmptargetmaps, probdistribution);
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

					ArrayList<ArrayList<Integer>> newpathseq = SecurityGameContraction.buildSuperGreedyCoverMultRes2(tmptargetmaps, 
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
			
			double ulimit = tmptargetmaps.get(attackedtarget).attackerreward;

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

		
		double defpayoff = SecurityGameContraction.expectedPayoffDef(attackedtarget, origpmat, tmptargetmaps, probdistribution);




		

		double[] res1 = {defpayoff, clusteringtime, solvingtime, targetsize, attackeru, slavetime, revmaptime};
		return res1;
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
						
						ArrayList<int[]> done = new ArrayList<int[]>();
						
						for(TargetNode n: currentst.get(st).nodes.values())
						{
							for(TargetNode n2: currentst.get(st).nodes.values())
							{
								if(n.getNeighbors().contains(n2)  &&  (n.getTargetid() != n2.getTargetid()) && (!isDone(done, n, n2)))
								{
									
									
									ArrayList<TargetNode> ppath = n.getPath(n2);
									
									for(TargetNode pnode: ppath)
									{
										origpmat[pnode.getTargetid()][jindex] = 1;
									}
									
									
									done.add(new int[]{n.getTargetid(), n2.getTargetid()});
									
								}
							}
						}
						
						
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


	private static boolean isDone(ArrayList<int[]> done, TargetNode n, TargetNode n2) {
		
		
		for(int[] x: done)
		{
			if((x[0] == n.getTargetid() && x[1] == n2.getTargetid()) ||(x[1] == n.getTargetid() && x[0] == n2.getTargetid()))
				return true;
		}
		
		
		return false;
	}

	
	

}
