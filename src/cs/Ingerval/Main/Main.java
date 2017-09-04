package cs.Ingerval.Main;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;
import groupingtargets.ClusterTargets;
import groupingtargets.GroupingTargets;

public class Main {


	/**
	 * 
	 * 1. Assume that S is a 2-D array where a Sudoku grid is stored.
2. sum = S.length * (S.length + 1) / 2;
3. for r =0 to S.length - 1
4. missingNum = sum;
5. for i = 0 to S.length – 1
6. if (S[r][i] contains a known value)
7. missingNum - = S[r][i];
8. print("The missing number in row: " + r + " is " + missingNum)
	 */



	public static void NotNaive(int[][] s)
	{
		int sum = s.length * ((s.length +1)/2);

		for(int i = 0; i < s.length; i++)
		{
			int missing = sum;

			for(int j = 0; j < s.length; j++)
			{
				for(int p = 1; p <= s.length; p++)
				{
					if(s[i][j] == p);
					{
						missing -= s[i][j];
					}
				}
			}
			System.out.println("The missing number in row " + i + " is: " + missing);
		}
	}





	public static void main(String[] args) throws Exception 
	{
		
		
		
		
		//GroupingTargets.createCSVTestData();
		//GroupingTargets.testWeka();

		int lstart=0,  lend=2, hstart=8, hend=10;

		int perc = 10;//Integer.parseInt(args[3]);
		/*int nrow=2;//Integer.parseInt(args[0]), 
		int ncol=2;//Integer.parseInt(args[1]);
*/		//int ITER=1;
		//int nTargets=nrow*ncol;
		//int nRes=1;
		//double dmax = 4;//Integer.parseInt(args[2]);

	
		
		
		int nrow = 20;
		int ncol = 20;
		int nTargets = nrow*ncol;
		int dmax = 60;
		
		
		
		int base = 0;
		int dest = 0;
		int k = 8;
		int radius = 3;
		
		int nRes=2;
		int dlim = 5;
		
		int ITER = 1;
		int ap = 4; // should be </= than cluster size
		int utiliy_l=0;
		int utility_h=10;
		/*int dmaxsuper = 30;
		int dminsuper = 5;*/
		
		
		
		
		
		
		// best algo 
		/*
		SecurityGameContraction.targets.clear();
		
		double[][] density=SecurityGameContraction.generateRandomDensity( perc, ITER, lstart, lend,  hstart, hend, nTargets, false);
		//double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);



		int percentages[] ={100};//{Integer.parseInt(args[1])};
		int thresholds[] = {2};
		
		SecurityGameContraction.doubleOracleTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		*/
		
		
		
		///////////////////////
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		HashMap<Integer, ArrayList<TargetNode>> alltargets = new HashMap<Integer, ArrayList<TargetNode>>();
		HashMap<Integer, HashMap<Integer, TargetNode>> alltargetmaps = new HashMap<Integer, HashMap<Integer, TargetNode>>();
		//HashMap<Integer, ArrayList<Integer>[]> allclus = new HashMap<Integer, ArrayList<Integer>[]>();
		//HashMap<Integer, ArrayList<Integer>[]> allclus = new HashMap<Integer, ArrayList<Integer>[]>();
		//double[][] density=SecurityGameContraction.generateRandomDensity( perc, ITER, lstart, lend,  hstart, hend, nTargets, false);
		
		//double[][] density = new double[ITER][nTargets];
		int ncat = 3;
		int ranges[][] = {{0,1}, {9, 10}, {2,8}};
		int[] percforcats = {80, 10, 10};
		int[] targetsincat = getTargetsInCats(nTargets, percforcats);
		double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);
		
		
		for(int iter = 0; iter<ITER; iter++)
		{
			ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
			HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
			//int[][] gamedata = new int[nTargets][4];

			//double u [][] = new double[nrow][ncol];
			
			//for(int i=0; i<)
			
			ClusterTargets.buildcsvGraphExp(nrow,ncol,density,targets, iter );
			
			//SecurityGameContraction.assignRandomDensityZeroSum(density, gamedata, targets, iter);
			//SecurityGameContraction.buildGraph(nrow, ncol, gamedata, targets);
			//SecurityGameContraction.assignRandomDensityZeroSum(density, gamedata, targets, iter);
			alltargets.put(iter, targets);
			
			for(TargetNode t : targets)
			{
				targetmaps.put(t.getTargetid(), t);
				
			}
			
			alltargetmaps.put(iter, targetmaps);
			
			//ArrayList<Integer>[] clus = GroupingTargets.makeCluster(k , nTargets, utiliy_l, utility_h, targets, targetmaps, density, iter);
			
			ClusterTargets.buildFile(nrow,ncol,density,targets, iter );
			
			
			//allclus.put(iter, clus);
			
			
		}
		
		
		int percentages[] ={100};//{Integer.parseInt(args[1])};
		int thresholds[] = {2};
		

		

		//4 DO + GC multi + GP 3 + LP + GC multi 
		ClusterTargets.DOWithClusteringTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes, alltargets, alltargetmaps);
		SecurityGameContraction.targets.clear();
		
		
		ClusterTargets.wekaClusteringWithDOExp(nrow,ncol,base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, alltargets, alltargetmaps);
		
		
		
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.DOTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes, alltargets, alltargetmaps);
		SecurityGameContraction.targets.clear();
		
		
		//14
		//SecurityGameContraction.noContractionNoColumnGenerationTest(density, ITER, nrow, ncol, percentages, dmax, nRes, alltargets, alltargetmaps );
		SecurityGameContraction.targets.clear();
		
		
		
		
		
		
		//
		
		
		
		
		
		
		
		
		
		
		
		
		
		//GroupingTargets.testClustering();
		
		//GroupingTargets.groupingWithDOExp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps);
		
		
		//GroupingTargets.wekaClusteringWithDOExpRW(nrow,ncol,base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps);
		
		
		
		
		
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(alltargets, alltargetmaps, ITER, nTargets , dmax, nRes);
		SecurityGameContraction.targets.clear();
		
		
		
		
		
		
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		
		/*ArrayList<TargetNode> targets = new ArrayList<TargetNode>();  //createGraph();
		HashMap<Integer, TargetNode> targetmaps = new HashMap<Integer, TargetNode>();
		ArrayList<Integer>[] clus = GroupingTargets.makeGraph(k, radius, dlim , nTargets, utiliy_l, utility_h, ap, targets, targetmaps);*/
		//GroupingTargets.createGraph2(targets, targetmaps);	
		//ArrayList<Integer>[] clus = GroupingTargets.makeDummyCluster(k);
		
		
		//only one path per super targets
		// but add path incrementally
		//GroupingTargets.groupTargetBaseline3Exp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
				
		
		
		
		//only one path per super targets
		// with modified LP
		//GroupingTargets.groupTargetBaseline2ExpNewMILP(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
		
		
		
		//only one path per super targets
		//GroupingTargets.groupTargetBaseline2Exp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
		
		
		
		//GroupingTargets.groupTargetTest();
		//GroupingTargets.groupTargetBaselineExp(base, dest, k, radius, dmax, nRes, nTargets, ITER, ap, allclus,  alltargets, alltargetmaps, dmaxsuper, dminsuper);
		
		SecurityGameContraction.targets.clear();
		
		//SecurityGameContraction.baselineForGroupTest(ITER, nTargets, dmax, nRes, alltargets, alltargetmaps);
 		//SecurityGameContraction.targets.clear();
		
		
		
		
		
		


		 targetsincat = getTargetsInCats(nTargets, percforcats);
		//double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);


		/*double[][] density=SecurityGameContraction.generateRandomDensity( perc, ITER, lstart, lend,  hstart, hend, nTargets, false);
		//double[][] density=SecurityGameContraction.generateRandomDensityV2(ncat, ITER, ranges, nTargets, targetsincat);



		int percentages[] ={100};//{Integer.parseInt(args[1])};
		int thresholds[] = {2};
		
		SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();*/




		//grouping test
		//SecurityGameContraction.GroupingTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();








		/*** REAL WORLD DATA TEST
		 * 
		 * 
		 * 
		 */


		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiRealWorldDataTest(1,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		/**
		 * END
		 */



		/****
		 * 
		 * 
		 * Exact solver
		 */
		int[][] gamedata = new int[nTargets][4];

		//SecurityGameContraction.contractionWithDoubleOracle(gamedata, nTargets, nRes, density, dmax, 0, 5, 5);



		//Exact DO + GC multi + GP 3 + LP + GC multi 
			//SecurityGameContraction.doubleOracleGCMultiExactLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		//SecurityGameContraction.doubleOracleGCMultiGP3LP_OPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		//SecurityGameContraction.doubleOracleGCMultiGP3LP_lexicoOPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		//SecurityGameContraction.doubleOracleGCMultiGP3LP_TOPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		//SecurityGameContraction.doubleOracleGCMultiGP3LP_modifiedOPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		//SecurityGameContraction.targets.clear();







		/**
		 * experiments
		 */

		//1 DO + GC single + GP multi + LP + GC multi 
		//	SecurityGameContraction.doubleOracleGCSingleGPMultiLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//2 DO + GC multi + GP multi + LP + GC multi 
		//    SecurityGameContraction.doubleOracleGCMultiGPMultiLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//3 DO + GC single + GP 3 + LP + GC multi 
			//SecurityGameContraction.doubleOracleGCSingleGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		//5 DO + GC single + GP multi + LP + sample path 
		//	SecurityGameContraction.doubleOracleGCSingleGPMultiLPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//6 DO + GC multi + GP multi + LP + sample path
		//z	SecurityGameContraction.doubleOracleGCMultiGPMultiLPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//7 DO + GC single + GP 3 + LP + sample path 
		//	SecurityGameContraction.doubleOracleGCSingleGP3LPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//8 DO + GC multi + GP 3 + LP + sample path 
		//	SecurityGameContraction.doubleOracleGCMultiGP3LPSamplePathTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();


		// 9 single oracle + greedy cover + greedy path
		//	SecurityGameContraction.contractionWithSingleOracleGreedyTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//10 single oracle + 2 targets + all paths
		//	SecurityGameContraction.contractionWithSingleOracleTest1(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		// 11 basic + GC Single + ins + All paths + LP
		//	SecurityGameContraction.basicGCSingleInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//12 basic + GC Multi + ins + All paths + LP 
		//	SecurityGameContraction.basicGCMultiInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//13 basic + GC Multi + ins + GP 3 + LP 
		//	SecurityGameContraction.basicGCMultiInstantGP3LPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//14
		//SecurityGameContraction.noContractionNoColumnGenerationTest(density, ITER, nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//15
		//	SecurityGameContraction.noContractionWithColumnGenerationTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();








		/**
		 * end of experiments
		 */



		/**
		 * 
		 * GP3 vs GP3 with APSP
		 */

		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.doubleOracleGCMultiGP3APSPLPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);







		/**
		 * APSP TESTS
		 */


		//SecurityGameContraction.doubleOracleAPSPGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//4 DO + GC multi + GP 3 + LP + GC multi 
		//SecurityGameContraction.doubleOracleGCMultiGP3LPGCMultiTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);



		//12 basic + GC Multi + ins + All paths + LP 
		//SecurityGameContraction.basicGCMultiInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		//SecurityGameContraction.targets.clear();
		//SecurityGameContraction.basicAPSPGCMultiInstantAllPathsLPTest(density, ITER,nrow, ncol, percentages, dmax, nRes);

		/////

		/**
		 * APSP TESTS ENDS
		 */





		//SecurityGameContraction.testPathGenerationv2();
		//SecurityGameContraction.basicAbstraction();
		//1
		//SecurityGameContraction.noContractionNoColumnGenerationTest(density, ITER, nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//2
		//SecurityGameContraction.noContractionWithColumnGenerationTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		//3 basic + seq + single res cover +  sample path + MILP
		//SecurityGameContraction.basicSeqAbstractionAlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//4 basic + seq + single res cover +  sample path + LP
		//SecurityGameContraction.basicSeqAbstractionWithLPGreedyCover1AlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//5 basic + seq + mult res cover +  sample path + LP
		//SecurityGameContraction.basicSeqAbstractionWithLPGreedyCoverMultResAlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);



		SecurityGameContraction.targets.clear();
		//6 basic + ins + single res cover +  sample path + MILP 
		//SecurityGameContraction.basicInstantAbstractionAlgorithmTest(density, ITER,nrow, ncol, percentages, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//7 basic + ins + single res cover +  sample path + LP 




		//SecurityGameContraction.instantContractionTest(density, ITER, nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.instantContractionTestAPSP(density, ITER, nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.spanningTreeTest(density,ITER);
		//SecurityGameContraction.testPathGenerationV3(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		//SecurityGameContraction.MPI();
		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.contractionWithSingleOracleGreedyTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.contractionWithSingleOracleTest1(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();






		//9 Double oracle +  mult res cover + greedy path 2(initial paths) + MILP + path sample 
		//SecurityGameContraction.contractionWithDoubelOracleTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();
		//10 Double oracle +  mult res cover + greedy path 2(initial paths) + MILP + greedycover 2 
		//SecurityGameContraction.contractionWithDoubelOracleGreedyCoverTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		SecurityGameContraction.targets.clear();



		SecurityGameContraction.targets.clear();

		SecurityGameContraction.targets.clear();
		//13 Double oracle +  mult res cover + greedy path 3(initial paths) + LP + greedycover 2 
		//SecurityGameContraction.contractionWithDoubelOracleGreedyCover3LPTest(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);
		//14 Double oracle +  mult res cover + greedy path 3(initial paths) + LP + path sample  
		//SecurityGameContraction.contractionWithDoubelOracleGreedyPath3TestLP(density,ITER,nrow, ncol, percentages, thresholds, dmax, nRes);


		/**
		 * do tests with extreme pruning
		 */


		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.basicSeqAbstractionAlgorithmWithExtreamPrunningTest(density, ITER,nrow, ncol, percentages, dmax);
		SecurityGameContraction.targets.clear();

		//SecurityGameContraction.basicInstantAbstractionAlgorithmWithExtrmPruningTest(density, ITER,nrow, ncol, percentages, dmax);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.instantContractionWithExtreamPruningTest(density, ITER, nrow, ncol, percentages, thresholds, dmax);
		SecurityGameContraction.targets.clear();
		//SecurityGameContraction.testPathGenerationWithExtreamPruningV3(density, ITER, nrow, ncol, percentages, thresholds, dmax);




		/*int[] contractionsizes = {5 };
		ITER = 1;
		//int x = ITER-1;
		double origcoinval = 0;
		double[] vals = new double[ITER];

		int targetssize = 15;
		int row =5;
		int col = 3;           
		for(int contractionsize : contractionsizes)
		{
			System.out.println("\n contraction size "+ contractionsize + "...");
			double sumcoins = 0.0;
			long sumcontracttime = 0;
			long sumsolvetime = 0;
			double nduplicatetargets = 0;
			double valentedge = 0;
			double constrsize = 0;


			for(int iter =0; iter<(ITER); iter++)
			{

				int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
				//SecurityGameAbstraction.convertToZeroSum(gamedata);
				Date start = new Date();
				long l1 = start.getTime();

				ArrayList<TargetNode> dominatednodes =  SecurityGameContraction.testContraction(targetssize, contractionsize, gamedata, 6, row, col);
				//System.out.println("\n iter "+ iter+ ", #dominated nodes "+ dominatednodes.size());
				Date stop = new Date();
				long l2 = stop.getTime();
				System.out.println("\n iter "+ iter+ " contraction size: "+dominatednodes.size());
				constrsize+= dominatednodes.size();


				long diff = l2 - l1;
				sumcontracttime += diff;

				for(TargetNode x: dominatednodes)
				{
					System.out.print(x.getTargetid() + " ");

				}
				System.out.println();
         		SecurityGameContraction.removePathsToDominatedNodes(SecurityGameContraction.targets,dominatednodes);
				SecurityGameContraction.copyInOrigTargets(SecurityGameContraction.targets);
				SecurityGameContraction.removeDominatedTargets(SecurityGameContraction.targets,dominatednodes);

				SecurityGameContraction.printNodesWithNeighborsAndPath(dominatednodes, SecurityGameContraction.targets);

				SecurityGameContraction.addVirtualBase(0,targetssize, SecurityGameContraction.targets);
				//int q = SecurityGameContraction.calculateEdgesWithCoin();
				HashMap<Integer, Integer> nodewithcoins= SecurityGameContraction.calculateNodesWithCoin(SecurityGameContraction.targets);
				SecurityGameContraction.transformToDirectedGraph(SecurityGameContraction.targets);
				//SecurityGameContraction.printNodesWithNeighborsAndPath(dominatednodes);

				//SecurityGameContraction.printDuplicateTargets();


				//sumcoins += LpSolverIP.solve(SecurityGameContraction.targets, dominatednodes, 9);
				Date start1 = new Date();
				l1 = start1.getTime();
				double vvv= MIPSolver3.solve(SecurityGameContraction.duplicatetargets, SecurityGameContraction.targets, nodewithcoins, dominatednodes, 6);

				if(contractionsize==0)
				{
					vals[iter] = vvv;
				}
				if((contractionsize>0) && (vvv>vals[iter]))
				{
					System.out.println("\n");
				}
				sumcoins += vvv;
				Date stop1 = new Date();
				l2 = stop1.getTime();
				long diff1 = l2 - l1;
				sumsolvetime += diff1;
				//TargetNode goal = CoinCollection.coinCollection(SecurityGameContraction.targets, dominatednodes);


				long secondInMillis = 1000;
				long minuteInMillis = secondInMillis * 60;
				long hourInMillis = minuteInMillis * 60;
				long dayInMillis = hourInMillis * 24;


				long elapsedDays = diff / dayInMillis;
				diff = diff % dayInMillis;
				long elapsedHours = diff / hourInMillis;
				diff = diff % hourInMillis;
				long elapsedMinutes = diff / minuteInMillis;
				diff = diff % minuteInMillis;
				long elapsedSeconds = diff / secondInMillis;
				//sumcoins += goal.getCoinvalue();

				//System.out.println("\niter"+ iter+", coin collected "+ goal.getCoinvalue());
				//System.out.println("Running Time "+ elapsedDays+" days, "+ elapsedHours+" hrs, "+ elapsedMinutes+" mins, "+ elapsedSeconds+" sec, " + diff + "ms" );
				//double coins = CoinCollection.printPath(goal);
				//System.out.println("\niter"+ iter+", coin collected path value "+ coins);
				//SecurityGameContraction.unimportanttargets.clear();
				nduplicatetargets += SecurityGameContraction.duplicatetargets.size();
				valentedge += nodewithcoins.size();
				SecurityGameContraction.targets.clear();
				SecurityGameContraction.duplicatetargets.clear();
				SecurityGameContraction.origtargets.clear();

			}
			double avgcoin = sumcoins/ITER;
			double avgcontracttime = sumcontracttime/ITER;
			double avgsolvetime = sumsolvetime/ITER;
			double avgndup = nduplicatetargets/ITER;
			valentedge = valentedge/ITER;
			constrsize = constrsize/ITER;



			writeInFile(contractionsize ,avgcoin, avgcontracttime, avgsolvetime, avgndup, valentedge, constrsize);
		}*/

		/*//int numCluster = 10;
		double[] numResources = {20};
		//int numTarget = 200;
		int ITER = 100;
		int[] clusternumbers = {40,20,10,5};
		for(int  numCluster: clusternumbers )
		{
			//double sumsolmmr = 0;
			//double summaxintrvlmmr = 0;
			double sumsolpoly = 0;
			double summaxintrvl = 0;
			double sumsolperfectmmr = 0;
			//double summaxintrvlperfectmmr = 0;

			double sumavg = 0;
			//double summaxintrvlpoly2 = 0;
			//double sumsolperfectmmr2 = 0;
			//double summaxintrvlperfectmmr2 = 0;


			boolean isequalcluster = false;

			for(int iter=0; iter<ITER; iter++)
			{
				//int[][] gamedata = SecurityGameAbstraction.createdummyData(numTarget,1,1);
				int[][] gamedata = SecurityGameAbstraction.parseSecurityGameFile("inputr-0.700000.csv", iter);
				SecurityGameAbstraction.convertToZeroSum(gamedata);
				List<Integer>[] clusteredtargets = KmeanClustering.getBestclusteredTargets(numCluster, gamedata, isequalcluster);
				int[] clustersizes = getClusterSizes(clusteredtargets);
				int maxinterval = KmeanClustering.getMaxInterval(clusteredtargets, gamedata);
				StrategyMapping strmap = new StrategyMapping(clusteredtargets, numCluster, gamedata);

				int[][][] abstractgame = SecurityGameAbstraction.makeAbstractIntervalSecurityGame(strmap);
				double[][] avgabstractgame = SecurityGameAbstraction.makeAVGPayoffAbstractSecurityGame(strmap);
				double origexp = SecurityGameAbstraction.originalGameTesing(numResources, gamedata);
				//strmap.printSecurityGameMapping();
				//SecurityGameAbstraction.makePositive(gamedata, -125);
				//SecurityGameAbstraction.addNoise(gamedata, 5, 10, true);
		 *//**
		 * try with interval abstract game
		 *//*
				//double[] resultmmr = SecurityGameAbstraction.testingBoundedRationalMMR(numResources, gamedata,numCluster);
				//sumsolmmr += resultmmr[0];
				//summaxintrvlmmr += resultmmr[1];
				double[] resultperfectmmr = SecurityGameAbstraction.testingPerfectRationalMMR(numResources, 
						gamedata, abstractgame, strmap, numCluster, isequalcluster, origexp, clustersizes);
				sumsolperfectmmr += resultperfectmmr[0];
				summaxintrvl += maxinterval;

				double[] resultpoly = SecurityGameAbstraction.wildlifeAbstraction(numResources,
						gamedata, abstractgame, strmap, numCluster, isequalcluster, origexp, clustersizes);
				sumsolpoly += resultpoly[0];
				//summaxintrvlpoly += maxinterval;


		  *//**
		  * try with avg payoff abstract game
		  *//*

				double[] resultavgabsgm = SecurityGameAbstraction.testingAvgAbstraction(numResources, 
						gamedata, avgabstractgame, strmap, numCluster, isequalcluster, origexp, clustersizes);
				sumavg += resultavgabsgm[0];
				//summaxintrvlperfectmmr2 += maxinterval;

				double[] resultpoly2 = SecurityGameAbstraction.wildlifeAbstraction(numResources,
						gamedata, avgabstractgame, strmap, numCluster, isequalcluster, origexp);
				sumsolpoly2 += resultpoly2[0];
				//summaxintrvlpoly2 += maxinterval;
			}
			sumsolperfectmmr /= ITER;
			sumsolpoly /= ITER;
			sumavg /= ITER;
			summaxintrvl /= ITER;
			System.out.println("Sol quality mmr "+ sumsolperfectmmr+ " maxinterval "
			+ summaxintrvl + " sumsolpoly "+ sumsolpoly );
			writeInFile(numCluster, sumsolperfectmmr, sumsolpoly, sumavg, summaxintrvl);
		}
		//SecurityGameAbstraction.testingAvgAbstraction(numResources, gamedata, numCluster);
		   */
	}











	private static int[] getTargetsInCats(int nTargets, int[] percforcats) {


		int x[] = new int[percforcats.length];


		int sum = 0;

		for(int i=0; i<percforcats.length; i++)
		{
			x[i] = (int)Math.floor(nTargets*(percforcats[i]/100.00));
			sum += x[i];
		}




		if(sum<nTargets)
		{

			int max = Integer.MIN_VALUE;
			int index= 0;
			for(int i=0; i<percforcats.length; i++)
			{
				if(max>percforcats[i])
				{
					max = percforcats[i];
					index = i;
				}
			}

			x[index] += (nTargets-sum);

		}



		return x;
	}











	private static void writeInFile(int contractionsize, double avgcoin,
			double avgcontracttime, double avgsolvetime, double avgndup, double valentedge, double constrsize) 
	{

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(contractionsize+ "," + avgcoin  
					+ "," + avgcontracttime +","+avgsolvetime+
					","+avgndup+","+valentedge+","+constrsize+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}

	private static int[] getClusterSizes(List<Integer>[] clusteredtargets) {

		int[] clustersizes = new int[clusteredtargets.length];
		for(int i=0; i<clusteredtargets.length; i++)
		{
			clustersizes[i] = clusteredtargets[i].size();
		}
		return clustersizes;

	}

	/*public static void writeInFile(int numCluster, double sumsolmmr,
			double sumsopoly, double sumavg, double summaxintrvl, double valentedge, double constrsize) {
		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(new File("/Users/fake/Documents/workspace/IntervalSGAbstraction/"+"result.csv"),true));
			pw.append(numCluster+ "," + sumsolmmr  + "," + sumsopoly +"," +
			sumavg + ","+summaxintrvl+","+constrsize+"\n");
			pw.close();

		}
		catch(Exception e)
		{

		}

	}
	 */




}
