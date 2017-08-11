package cs.Interval.Clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;








public class KmeanClustering {

	private static  boolean RAND_POINTS_FROM_OBSERVATION = false; 
	private static  boolean RAND_ACTION_INIT_TO_CLUSTERS = true;
	private static  boolean DIST_METRIC_LINE = true; //used to create dir, if this is true, then other one is false, vice versa
	private static  boolean DIST_METRIC_EUCLIDEAN =  false;
	private static  boolean MAX_DIST = false;   //make it false if euclidean
	private static  final boolean SUM_DIST = true; //make it false if euclidean
	private static  boolean MAX_DELTA = true; 
	private static  boolean AVRG_DELTA = false;


	public static void setRAND_POINTS_FROM_OBSERVATION(
			boolean rAND_POINTS_FROM_OBSERVATION) {
		RAND_POINTS_FROM_OBSERVATION = rAND_POINTS_FROM_OBSERVATION;
	}

	public static void setRAND_ACTION_INIT_TO_CLUSTERS(
			boolean rAND_ACTION_INIT_TO_CLUSTERS) {
		RAND_ACTION_INIT_TO_CLUSTERS = rAND_ACTION_INIT_TO_CLUSTERS;
	}

	public static int getBitValue(int var, int bitposition)
	{
		int x = var & (1<<bitposition);
		if(x>0)
			return 1;


		return 0; 
	}

	public static boolean isRandActionInitToClusters() {
		return RAND_ACTION_INIT_TO_CLUSTERS;
	}

	public static boolean isRandPointsFromObservation() {
		return RAND_POINTS_FROM_OBSERVATION;
	}

	public static boolean isDistMetricLine() {
		return DIST_METRIC_LINE;
	}

	public static boolean isDistMetricEuclidean() {
		return DIST_METRIC_EUCLIDEAN;
	}

	public static boolean isMaxDist() {
		return MAX_DIST;
	}

	public static boolean isSumDist() {
		return SUM_DIST;
	}

	public static boolean isMaxDelta() {
		return MAX_DELTA;
	}

	public static boolean isAvrgDelta() {
		return AVRG_DELTA;
	}




	public static List<Integer>[] clusterTargets(int numberofclusters, int[][] gamedata, boolean isequalcluster)
	{
		int numberoftargets = gamedata.length;
		int targetspercluster = gamedata.length/numberofclusters;


		//int opponent = 1^player;
		//final int   RUNNING_FOR = 20;
		//int[] numberofactions = mg.getNumActions();
		//double[] extreampayoffs = mg.getExtremePayoffs();
		int INDEX_LIMIT = 4;
		List<Integer[]>[] clusters = new List[numberofclusters]; // create an array of list. Each list will contain arrays of double
		double[][] clusterpoints = new double[numberofclusters][INDEX_LIMIT];  // cluster points for each cluster
		double[][] oldclusterpoints = new double[numberofclusters][INDEX_LIMIT];  
		double[][] distancefromcluster = new double[numberofclusters][INDEX_LIMIT];
		int runningInstance = 0;
		double[] sumofdifferences = new double[numberofclusters];  //store the sum of differences for clusters
		double[] maxdifference = new double[numberofclusters]; 
		double[] squaredrootdifference = new double[numberofclusters];
		ArrayList<Integer> alreadyassignedactions = new ArrayList<Integer>(); // for random partition method
		//boolean flagforrandompartition = false;
		for(int i=0; i< numberofclusters; i++)
		{
			clusters[i] = new ArrayList<Integer[]>(); 
		}
		if(KmeanClustering.isRandPointsFromObservation()) // implemented only for two players
		{
			ArrayList<Integer> unassignedpoints = new ArrayList<Integer>();
			int targetindex = 0;
			for(int i=0; i<numberoftargets; i++)
			{
				unassignedpoints.add(i);
			}
			for(int i=0; i< numberofclusters; i++)
			{
				targetindex = randInt(0, unassignedpoints.size()-1);
				int tmptarget = unassignedpoints.get(targetindex);
				unassignedpoints.remove(targetindex);
				int index = 0;
				while(index<INDEX_LIMIT)
				{
					clusterpoints[i][index] = gamedata[tmptarget][index];
					index++;
				}
			}
			//System.out.println("Initial cluster points...");
			for(int i =0; i< numberofclusters; i++)
			{
				//System.out.print("Cluster: "+ i + " ");
				for(int j =0; j< INDEX_LIMIT; j++)
				{
					//System.out.print(" "+clusterpoints[i][j]); 
				}
				//System.out.print("\n");

			}
		}
		if(KmeanClustering.isRandActionInitToClusters())
		{

			//	Logger.log("\n entered RandActionInitToClusters() ", false);
			/*
			 * pick an action and randomly assign it to a cluster.  
			 * then calculate the mean
			 */
			ArrayList<Integer> unassignedactions = new ArrayList<Integer>();

			for(int i=0;i<numberoftargets;i++)
			{
				unassignedactions.add(i);
			}
			// iterate through all the actions and assign to a cluster randomly
			for(int clusterid=0; clusterid<numberofclusters; clusterid++)
			{

				for(int j=0; j<(numberoftargets/numberofclusters); j++) //slot, it ensures that every cluster has equal number of actions for random partition
				{

					if(unassignedactions.size()>1)
					{
						int chosenactionindex = randInt(1, unassignedactions.size()) -1;
						int z = unassignedactions.get(chosenactionindex);
						assignToCluster(clusters, unassignedactions.get(chosenactionindex), clusterid, numberoftargets, gamedata, INDEX_LIMIT);
						//System.out.println("\n random paritiion<<<>>>>> Action is assigned to cluster "+clusterid);
						unassignedactions.remove(chosenactionindex);
					}
					else if(unassignedactions.size() ==1)
					{
						assignToCluster(clusters, unassignedactions.get(0), clusterid, numberoftargets, gamedata, INDEX_LIMIT);
						int x = unassignedactions.get(0);
						//System.out.println("\n random paritiion<<<>>>>> Action is assigned to cluster "+clusterid);
						unassignedactions.remove(0); //remove the last element
						break;
					}

				}
				/*
				 * minindex has the index for cluster
				 * make a tuple like (action, payoffs1, payoff2,...)
				 */

			}  // end of for loop

			//check if there are any actions remained unassigned
			if(unassignedactions.size() != 0)
			{
				for(Integer x: unassignedactions)
				{
					int a = numberofclusters-1;
					assignToCluster(clusters, x, numberofclusters-1, numberoftargets, gamedata, INDEX_LIMIT); // assign all the remainning actions to the last cluster.
					//System.out.println("\n random paritiion<<<>>>>> Action  is assigned to cluster "+a);

				}

			}
			// print the clusters..
			//System.out.println("\n\nInitialization to clusters \n");
			for(int i=0; i<clusters.length; i++)
			{
				//System.out.print("Cluster "+ i+ " : ");

				for(Integer[] x: clusters[i])
				{
					//System.out.print(x[0]+ " ");
				}
				//System.out.println();
			}
			
			calculateClusterMean(clusters, clusterpoints, numberofclusters, gamedata, INDEX_LIMIT);
			//System.out.println("Initial Cluster means...");
			for(int i =0; i< numberofclusters; i++)
			{

				//System.out.print("Cluster: "+ i + " ");
				//Logger.log("Cluster: "+ i + " ", false);
				for(int j =0; j< INDEX_LIMIT; j++)
				{

					//System.out.print(" "+clusterpoints[i][j]); 
					//Logger.log(" "+clusterpoints[i][j], false);

				}
				//System.out.print("\n");
				//Logger.log("\n", false);
			}
		}


		while(true)
		{
			//System.out.println("\nKmeans Iteration: "+ runningInstance);

			if(runningInstance>=100)
			{
				List<Integer>[] finalcluster = new List[numberofclusters];
				for(int i=0; i< numberofclusters; i++){

					finalcluster[i] = new ArrayList<Integer>(); 
				}
				for(int i=0; i<numberofclusters; i++)
				{
					for(Integer[] x: clusters[i])
					{
						finalcluster[i].add(x[0].intValue());
					}
				}
				return finalcluster;
			}
			//copy the cluster points to old cluster points.
			for(int i=0; i< numberofclusters; i++)
			{
				for(int j=0; j<INDEX_LIMIT; j++)
				{
					oldclusterpoints[i][j] = clusterpoints[i][j];
				}
			}
			// now clear/create e cluster object the new cluster for a new iteration. 
			for(int i=0; i< numberofclusters; i++)
			{
				clusters[i]= new ArrayList<Integer[]>(); //.clear();
			}
			/*
			 * Now iterate over all the possible action touples for player 1. 
			 * calclate the difference from cluster points
			 * assign to the cluster with the minimum difference.
			 *  
			 */
			for(int target = 0; target < numberoftargets; target++)
			{
				for(int rewardindex = 0; rewardindex < INDEX_LIMIT; rewardindex++)
				{

					double tmppayoff = gamedata[target][rewardindex]; //mg.getPayoff(outcome, player); //get the payoff for player 1 or player 2
					for(int clusterindex =0; clusterindex<numberofclusters; clusterindex++)
					{
						if(KmeanClustering.isDistMetricLine())
						{
								////System.out.println("\n entered DistMetricLine() ");

							if((tmppayoff < 0 && clusterpoints[clusterindex][rewardindex] < 0) ||  (tmppayoff >=0  && clusterpoints[clusterindex][rewardindex] >= 0))

							{
								distancefromcluster[clusterindex][rewardindex] = Math.abs((Math.abs(tmppayoff) - Math.abs(clusterpoints[clusterindex][rewardindex])));
							}
							else if((tmppayoff >= 0 && clusterpoints[clusterindex][rewardindex] < 0))
							{
								distancefromcluster[clusterindex][rewardindex] = (tmppayoff + Math.abs(clusterpoints[clusterindex][rewardindex]));
							}
							else if((tmppayoff < 0 && clusterpoints[clusterindex][rewardindex] >= 0))
							{
								distancefromcluster[clusterindex][rewardindex] = (clusterpoints[clusterindex][rewardindex]  + Math.abs(tmppayoff));
							}
						}
						/*
						 * calculate the differences of payoffs for each cluster points 
						 * calculate euclidean distance: first take the squares of difference....
						 */
						if(KmeanClustering.isDistMetricEuclidean())
						{
							////System.out.println("\n entered DistMetricEuclidean() ");
							distancefromcluster[clusterindex][rewardindex] = (clusterpoints[clusterindex][rewardindex]  - (tmppayoff));
							distancefromcluster[clusterindex][rewardindex] = distancefromcluster[clusterindex][rewardindex] * distancefromcluster[clusterindex][rewardindex];
						}
					}
				}
				//double min = Double.POSITIVE_INFINITY;
				//int minindex = 0;
				double[][] rank = new double[numberofclusters][2];
				if(KmeanClustering.isDistMetricEuclidean())
				{
					//System.out.println("Entered euclidean distance");
					/*
					 * find the squared root distances
					 */
					for(int l =0; l< numberofclusters; l++)
					{
						squaredrootdifference[l] = 0;
						for(int m =0; m< distancefromcluster[l].length; m++)
						{
							squaredrootdifference[l] += distancefromcluster[l][m];
						}
						squaredrootdifference[l] = Math.sqrt(squaredrootdifference[l]);
						//int a = target+1;
						//System.out.println("\n Target "+ target+"'s euclidean distance from cluster "+ l+ " : "+squaredrootdifference[l]);
					}
					// find the minimum squared root distance
					/*for(int n =0; n< squaredrootdifference.length; n++)
					{
						if(min > squaredrootdifference[n])
						{
							min = squaredrootdifference[n];
							minindex = n;
						}
					}*/
					//System.out.println("Ranking for target(euclid) "+target);
					rank = genRankFromMaxDistances(squaredrootdifference);

				}
				if(KmeanClustering.isMaxDist())
				{
					//System.out.println("Entered max distance");
					/*
					 * calculate the max difference instead of summing them
					 */
					for(int l =0; l< numberofclusters; l++)
					{
						double maxdiff =Double.NEGATIVE_INFINITY;
						for(int m =0; m< 4; m++)
						{
							if(maxdiff<distancefromcluster[l][m])
							{
								maxdiff = distancefromcluster[l][m];
							}
						}
						maxdifference[l] = maxdiff;
					}
					/*
					 * find the minimum difference among the maximum differences
					 */
					/*for(int n =0; n< maxdifference.length; n++)
					{
						if(min > maxdifference[n])
						{
							min = maxdifference[n];
							minindex = n;
						}
					}*/
					//System.out.println("Ranking for target (max)"+target);
					rank = genRankFromMaxDistances(maxdifference);
				}
				if(KmeanClustering.isSumDist())
				{

					//
					//System.out.println("\n entered SumDist() ");
					/*
					 * Here you have all the differences for action i 
					 * calculate the sum of the differences
					 * Then find the minimum sum
					 * 
					 */
					for(int l =0; l< numberofclusters; l++)
					{
						sumofdifferences[l] = 0;
						for(int m =0; m< distancefromcluster[l].length; m++)
						{
							sumofdifferences[l] += distancefromcluster[l][m];
						}
						//int a = i+1;
						//Logger.log("\n Action "+ a+"'s sum distance from cluster "+l+" is : "+sumofdifferences[l] , false);

					}
					/*for(int n =0; n< sumofdifferences.length; n++)
					{
						if(min > sumofdifferences[n])
						{
							min = sumofdifferences[n];
							minindex = n;
						}
					}*/
					//System.out.println("Ranking for target (sum)"+target);
					rank = genRankFromMaxDistances(sumofdifferences);

				}


				//System.out.println("\nIteration: "+ runningInstance + " \n Assigning cluster points ");
				//System.out.println("target "+target +" is assigned to cluster "+ minindex);
				//assignToCluster(clusters, target, minindex, gamedata, INDEX_LIMIT);
				assignToCluster(clusters, target, rank, gamedata, targetspercluster, INDEX_LIMIT, isequalcluster);

			}  // end of outer for loop
			/*
			 * now recalculate the cluster points
			 */

			//System.out.println("Clustered Targets : " );
			/*for(int i=0; i< clusters.length; i++)
			{
				//System.out.print("Cluster " + i + " : ");
				for(Integer[] target: clusters[i])
				{
					//System.out.print(target[0]);
					if(clusters[i].indexOf(target) < (clusters[i].size()-1) )
					{
						//System.out.print(",");
					}
				}
				//System.out.print("\n");
			}*/





			calculateClusterMean(clusters, clusterpoints, numberofclusters, gamedata, INDEX_LIMIT);
			//System.out.println("\n\nIteration: "+ runningInstance + " ");
			//System.out.println("\n\nK-mean Iteration: "+ runningInstance  +" new cluster points(mean)\n");
			/*for(int i =0; i< numberofclusters; i++)
			{

				//System.out.print("Cluster: "+ i + " ");
				//Logger.log("Cluster: "+ i + " ", false);
				for(int j =0; j< INDEX_LIMIT; j++)
				{

					//System.out.print(" "+clusterpoints[i][j]); 
					//Logger.log(" "+clusterpoints[i][j], false);

				}
				//System.out.print("\n");
				//Logger.log("\n", false);
			}*/
			//System.out.println("Checking or stop ");
			boolean checkforstop = true;
			for(int i=0; i< numberofclusters; i++)
			{

				for(int j=0; j<INDEX_LIMIT; j++)
				{
					if(clusterpoints[i][j] != oldclusterpoints[i][j])
					{
						checkforstop = false;
						break;
					}
				}
				if(checkforstop==false)
				{
					break;
				}
			}

			//System.out.println("\nIteration: "+ runningInstance + "Checking exit condition ");
			if(checkforstop == true && (!isClusterIsEmpty(clusters)))
			{
				//System.out.println("\n Exiting..." );
				break;
			}
			//Logger.log("\n\n", false);
			runningInstance++;
			//printCluster(clusters);



		}//end of outer while loop

		//printCluster(clusters);

		List<Integer>[] finalcluster = new List[numberofclusters];
		for(int i=0; i< numberofclusters; i++){

			finalcluster[i] = new ArrayList<Integer>(); 
		}
		for(int i=0; i<numberofclusters; i++)
		{
			for(Integer[] x: clusters[i])
			{
				finalcluster[i].add(x[0].intValue());
			}
		}
		return finalcluster;


	}


	private static void assignToCluster(List<Integer[]>[] clusters,
			Integer target, int clusterid, int numberoftargets, int[][] gamedata, int INDEX_LIMIT) {
		
		
		Integer[] tupleincluster = new Integer[numberoftargets+1]; // +1 for the action
		tupleincluster[0] = target; //the action in the first index
		/*

		 * now assign the payoffs

		 */
		int[] tmpoutcome = new int[2];
		for(int p = 0; p< INDEX_LIMIT; p++)
		{
			
			tupleincluster[p+1] = gamedata[target][p];

		}
		clusters[clusterid].add(tupleincluster); 
		
	}

	private static double[][] genRankFromMaxDistances(double[] maxdifference) {
		double[][] rank = new double[maxdifference.length][2]; // [clusterindex][maxdifference]
		for(int i=0; i< maxdifference.length; i++)
		{
			rank[i][0] = i;
			rank[i][1] = maxdifference[i];
		}
		sortArrayAsc(rank);
		
		for(int i=0; i<rank.length; i++)
		{
			
				//System.out.println(rank[i][0] + ",  "+rank[i][1]);
			

		}
		return rank;

	}

	private static void sortArrayAsc(double[][] ds)
	{

		double[] swap = {0.0,0.0};

		for (int i = 0; i < ds.length; i++) 
		{
			for (int d = 1; d < ds.length-i; d++) 
			{
				if (ds[d-1][1] > ds[d][1]) /* For descending order use < */
				{
					swap = ds[d];
					ds[d]  = ds[d-1];
					ds[d-1] = swap;
				}
			}
		}

	}

	/**
	 * calculates cluster mean for security games
	 * @param clusters
	 * @param clusterpoints
	 * @param numberofclusters
	 * @param gamedata
	 */
	private static void calculateClusterMean(List<Integer[]>[] clusters,
			double[][] clusterpoints, int numberofclusters, int[][] gamedata, int INDEX_LIMIT) {

		/*

		 * now recalculate the cluster mean

		 */
		//int opponent = 1^player;
		double average = 0;
		for(int clusterindex = 0; clusterindex< numberofclusters; clusterindex++)
		{
			int clustersize = clusters[clusterindex].size();
			if(clustersize==0)
			{
				//System.out.println("\n\nEmpty cluster: "+ clusterindex + " ");
				int randomtarget;
				while(true)
				{
					// if cluster is empty, assign random points from the strategy
					randomtarget = randInt(0,gamedata.length-1 );
					//check if the payoffs are same as centroid of another cluster
					//System.out.println(".....");
					break;

				}
				//Logger.log("\nAction "+ randomtarget+"'s payoffs are assigned to cluster "+ clusterindex, false);
				for(int j =0; j<INDEX_LIMIT; j++)
				{
					double reward = gamedata[randomtarget][j];  //mg.getPayoff(outcome, player);
					clusterpoints[clusterindex][j] = reward;
				}
				/**
				 * assign the target to the cluster
				 */
				
				assignToCluster(clusters, randomtarget, clusterindex, gamedata, INDEX_LIMIT);
				
				//Logger.log("\n cluster: "+ clusterindex + " is empty, points after reassignment:\n", false);


			}
			else if(clustersize>0)
			{
				for(int rewardindex = 0; rewardindex< INDEX_LIMIT; rewardindex++)
				{
					average = 0;   // corrected, average should be reset after calculating for every action
					for(Integer[] x: clusters[clusterindex])
					{
						average += x[rewardindex+1]; 
					}
					if(clustersize != 0)
					{
						clusterpoints[clusterindex][rewardindex] = average/clustersize; 
					}
				}
			}

		}

	}


	

	

	private static boolean isClusterIsEmpty(List<Integer[]>[] clusters) {
		for(int i=0; i<clusters.length; i++)
		{
			if(clusters[i].size()==0)
				return true;
		}

		return false;
	}

	public static void printCluster(List<Integer[]>[] clusters)
	{
		for(int clusterid = 0; clusterid<clusters.length; clusterid++)
		{
			System.out.println("\n*****Cluster : "+ clusterid);
			for(Integer[] x: clusters[clusterid])
			{
				for(Integer y: x)
				{
					System.out.print(y+" ");
				}
				System.out.println();

			}
		}

	}




	/**
	 * assigns targets to a cluster
	 * @param clusters
	 * @param target
	 * @param minindex
	 * @param targettoassign
	 * @param gamedata
	 */
	private static void assignToCluster(List<Integer[]>[] clusters, int target,
			int assignedcluster, int[][] gamedata, int INDEX_LIMIT) {

		//int opponent = 1^player;
		//int oppnumaction = mg.getNumActions(opponent);
		Integer[] tupleincluster = new Integer[INDEX_LIMIT+1]; // +1 for the target
		tupleincluster[0] = target; //the target in the first index
		/*
		 * now assign the rewards
		 */
		for(int p = 0; p<INDEX_LIMIT; p++)
		{
			tupleincluster[p+1] = gamedata[target][p]; //mg.getPayoff(tmpoutcome, player); 
		}
		clusters[assignedcluster].add(tupleincluster); 

	}

	/**
	 * 
	 * @param clusters
	 * @param target
	 * @param rank [clusterindex][maxdiff]
	 * @param gamedata
	 * @param targetspercluster
	 * @param INDEX_LIMIT
	 */
	private static void assignToCluster(List<Integer[]>[] clusters, int target,
			double[][] rank, int[][] gamedata, int targetspercluster, int INDEX_LIMIT, boolean isequalcluster) {


		boolean assigned = false;
		int assignedcluster = -1;
		/**
		 * find the nearest empty cluster cluster and try to assign according to rank
		 */
		for(int rankid = 0; rankid<rank.length; rankid++)
		{
			int tmpcluster = (int)rank[rankid][0]; 
			if(!isequalcluster)
			{
				assignedcluster = tmpcluster;
				assigned = true;
				//System.out.println("Target "+target+ " is assigned to cluster wo ranking "+assignedcluster);
				break;
			}
			else if(clusters[tmpcluster].size()<targetspercluster)
			{
				
				assignedcluster = tmpcluster;
				assigned = true;
				//System.out.println("Target "+target+ " is assigned to cluster "+assignedcluster);
				break;
			}

		}
		if(!assigned)
		{
			assignedcluster = (int)rank[0][0];
			//System.out.println("FOrcibly Target "+target+ " is assigned to cluster "+assignedcluster);
		}


		//int opponent = 1^player;
		//int oppnumaction = mg.getNumActions(opponent);
		Integer[] tupleincluster = new Integer[INDEX_LIMIT+1]; // +1 for the target
		tupleincluster[0] = target; //the target in the first index
		/*
		 * now assign the rewards
		 */
		for(int p = 0; p<INDEX_LIMIT; p++)
		{
			tupleincluster[p+1] = gamedata[target][p]; //mg.getPayoff(tmpoutcome, player); 
		}
		clusters[assignedcluster].add(tupleincluster); 

	}



	/**
	 * 
	 * @param clusters contains the cluster actions
	 * @return true if no cluster is empty
	 */
	private static boolean noClusterIsEmpty(List<Double[]>[] clusters) 
	{

		for(int i=0; i<clusters.length; i++)
		{
			if(clusters[i].size()==0)
				return false;
		}

		return true;
	}




	public static int randInt(int min, int max) {

		// Usually this should be a field rather than a method variable so
		// that it is not re-seeded every call.
		Random rand = new Random();

		// nextInt is normally exclusive of the top value,
		// so add 1 to make it inclusive
		int randomNum = rand.nextInt((max - min) + 1) + min;

		return randomNum;
	}

	public static List<Integer>[] getBestclusteredTargets(int numCluster,
			int[][] gamedata, boolean isequalcluster) {

		HashMap<Integer,List<Integer>[]> clusterswithdiffintervals = new HashMap<Integer, List<Integer>[]>();
		int minindex = -1;
		int minmaxinterval = Integer.MAX_VALUE;
		for(int i=0; i<10; i++)
		{
			List<Integer>[] clusteredtargets = KmeanClustering.clusterTargets(numCluster, gamedata, isequalcluster);
			int maxinterval = getMaxInterval(clusteredtargets, gamedata);
			clusterswithdiffintervals.put(i, clusteredtargets);
			//System.out.println("Clustering# "+i+", Max Interval "+maxinterval);
			if(minmaxinterval>maxinterval)
			{
				minmaxinterval = maxinterval;
				minindex = i;
			}
		}
		//System.out.println("Choosing Clustering# "+minindex+", Max Interval "+minmaxinterval);
		//printCluster(clusterswithdiffintervals.get(minindex), gamedata);
		return clusterswithdiffintervals.get(minindex);
	}
	
	
	private static void printCluster(List<Integer>[] clusters, int[][] gamedata) 
	{
		for(int clusterid = 0; clusterid<clusters.length; clusterid++)
		{
			System.out.println("\n*****Cluster : "+ clusterid);
			for(Integer x: clusters[clusterid])
			{
				System.out.print(x+" ");
				for(int i=0; i<4; i++)
				{
					System.out.print(gamedata[x][i]+" ");
				}
				System.out.println();

			}
		}
		
	}

	public static int getMaxInterval(List<Integer>[] clusteredtargets,
			int[][] gamedata) 
	{
		
		int[][][] securitygame = new int[clusteredtargets.length][4][2];
		int dminr = Integer.MAX_VALUE;
		int dmaxr = Integer.MIN_VALUE;
		
		int dminp = Integer.MAX_VALUE;
		int dmaxp = Integer.MIN_VALUE;
		
		int aminr = Integer.MAX_VALUE;
		int amaxr = Integer.MIN_VALUE;
		
		int aminp = Integer.MAX_VALUE;
		int amaxp = Integer.MIN_VALUE;
		
		int maxInterval = Integer.MIN_VALUE;
		for(int clusterindex = 0; clusterindex<clusteredtargets.length; clusterindex++)
		{
			dminr = Integer.MAX_VALUE;
			dmaxr = Integer.MIN_VALUE;

			dminp = Integer.MAX_VALUE;
			dmaxp = Integer.MIN_VALUE;

			aminr = Integer.MAX_VALUE;
			amaxr = Integer.MIN_VALUE;

			aminp = Integer.MAX_VALUE;
			amaxp = Integer.MIN_VALUE;
			for(Integer target : clusteredtargets[clusterindex])
			{
				/**
				 * min condition for defender reward
				 */
				if(gamedata[target][0]<0)
				{
					//System.out.print("OOO");
				}
				if(dminr > gamedata[target][0])
				{
					dminr = gamedata[target][0];
				}
				/**
				 * max condition for defender reward
				 */
				if(dmaxr < gamedata[target][0])
				{
					dmaxr = gamedata[target][0];
				}
				/**
				 * min condition for defender penalty
				 */
				if(dminp > gamedata[target][1])
				{
					dminp = gamedata[target][1];
				}
				/**
				 * max condition for defender penalty
				 */
				if(dmaxp < gamedata[target][1])
				{
					dmaxp = gamedata[target][1];
				}

				/**
				 * min condition for defender reward
				 */
				if(aminr > gamedata[target][2])
				{
					aminr = gamedata[target][2];
				}
				/**
				 * max condition for defender reward
				 */
				if(amaxr < gamedata[target][2])
				{
					amaxr = gamedata[target][2];
				}
				/**
				 * min condition for defender penalty
				 */
				if(aminp > gamedata[target][3])
				{
					aminp = gamedata[target][3];
				}
				/**
				 * max condition for defender penalty
				 */
				if(gamedata[target][3]>0)
				{
					//System.out.print("OOOshh");
				}
				if(amaxp < gamedata[target][3])
				{
					amaxp = gamedata[target][3];
				}
			}
			securitygame[clusterindex][0][0] = dminr;
			securitygame[clusterindex][0][1] = dmaxr;

			securitygame[clusterindex][1][0] = dminp;
			securitygame[clusterindex][1][1] = dmaxp;

			securitygame[clusterindex][2][0] = aminr;
			securitygame[clusterindex][2][1] = amaxr;

			securitygame[clusterindex][3][0] = aminp;
			securitygame[clusterindex][3][1] = amaxp;
			int clustermaxinterval = Integer.MIN_VALUE;
			for(int index=0; index<4; index++)
			{
				int tmpdist = calcLineDist(securitygame[clusterindex][index]);
				if(clustermaxinterval < tmpdist)
				{
					clustermaxinterval = tmpdist;
				}
				if(tmpdist>maxInterval)
				{
					maxInterval = tmpdist;
				}
			}
			//System.out.println("Cluster "+clusterindex + " interval "+ clustermaxinterval);

		}
		return maxInterval;
		
		
	}

	public static int calcLineDist(int[] interval) {
		// TODO Auto-generated method stub
		
		if((interval[0] < 0 && interval[1] < 0) ||  (interval[0] >=0  && interval[1] >= 0))

		{
			return Math.abs((Math.abs(interval[0]) - Math.abs(interval[1])));
		}
		else if((interval[0] >= 0 && interval[1] < 0))
		{
			return (interval[0] + Math.abs(interval[1]));
		}
		else if((interval[0] < 0 && interval[1] >= 0))
		{
			return (interval[1]  + Math.abs(interval[0]));
		}
		return -1;
		
	}
	
	public static double calcLineDist(double[] interval) {
		// TODO Auto-generated method stub
		
		if((interval[0] < 0 && interval[1] < 0) ||  (interval[0] >=0  && interval[1] >= 0))

		{
			return Math.abs((Math.abs(interval[0]) - Math.abs(interval[1])));
		}
		else if((interval[0] >= 0 && interval[1] < 0))
		{
			return (interval[0] + Math.abs(interval[1]));
		}
		else if((interval[0] < 0 && interval[1] >= 0))
		{
			return (interval[1]  + Math.abs(interval[0]));
		}
		return -1;
		
	}














}
