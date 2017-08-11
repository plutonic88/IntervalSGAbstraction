package cs.Interval.Abstraction;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import lpWrapper.LPSolverException;
import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import models.AttBoundStructure;
import models.PayoffStructure;
import models.SUQRAdversary;
import models.Target;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import perfectMMR.MMRZeroSumNew2;
import src.tambe.usc.DOBSS.gameRepresentation.SSGPayoffs;
import src.tambe.usc.DOBSS.gameRepresentation.SSGSolution;
import src.tambe.usc.DOBSS.gameRepresentation.StructuredSecurityGame;
import src.tambe.usc.DOBSS.solvers.FastSecurityGameSolver;
import algorithms.MMR;
import cs.Interval.Clustering.KmeanClustering;
import cs.Interval.Clustering.StrategyMapping;
import cs.Interval.Solvers.IntervalSecurityGame;
import cs.Interval.Solvers.PolynomialIntervalSolver;


//import src.tambe.usc.DOBSS.solvers.FastSecurityGameSolver.TargetData;



public class SecurityGameAbstraction {

	public static Random rand = new Random();


	
/**
 * 
 * @param numResources
 * @param gamedata 2d array containing rewards and penalties for targets
 * gamedata[t][0]= def reward
 * gamedata[t][1]= def penalty
 * gamedata[t][2]= att reward
 * gamedata[t][3]= att penalty
 * 
 * 
 * @return
 */
	public static double originalGameTesing(double[] numResources, int[][] gamedata)
	{
		//int[] numResources = {5};
		long totalTimeOrigmai = 0;
		FastSecurityGameSolver fsgs = new FastSecurityGameSolver();
		//int[][] gamedata = createdummyData();
		//int[][] gamedata = parseSecurityGameFile("inputr-0.100000.csv");
		//convertToZeroSum(gamedata);
		//makePositive(gamedata, -10);
		StructuredSecurityGame samplegame = genStructuredSecurityGame(gamedata, 1, new double[]{numResources[0]});
		//SSGPayoffs[] payoffs = samplegame.getPayoffs();
		Long startTimeOrigami = System.currentTimeMillis();
		SSGSolution sampleSolutionOrigmai = fsgs.solveGame(samplegame);
		//checkNormality(numResources[0], gamedata.length, sampleSolutionOrigmai);
		Long endTimeOrigami = System.currentTimeMillis();
		double coverageProb=0;
		sampleSolutionOrigmai.computeExpectedPayoffs(samplegame);

		//System.out.println("Result: " + sampleSolutionOrigmai);

		long time = endTimeOrigami -startTimeOrigami;
		totalTimeOrigmai += time;
		//System.out.println("Total Time Origami: " + totalTimeOrigmai);
		//System.out.println("Average time origmai:" +avgTimeOrigami + "\n" );
		return sampleSolutionOrigmai.getDefenderPayoff();


	}



	private static void checkNormality(double numResource, int numTargets,
			SSGSolution sampleSolutionOrigmai) {

		double sum =0; 
		for(int i=0; i<numTargets; i++)
		{
			double[] tmp= sampleSolutionOrigmai.getTargetProbs(i);
			for(double x: tmp)
			{
				sum += x;
			}
		}
		System.out.println("SUm prob " + sum);

	}

	public static void makePositive(int[][] gamedata, int minpayoff) {


		for(int i=0; i<gamedata.length; i++)
		{
			for(int j=0; j<4; j++)
			{
				gamedata[i][j] += Math.abs(minpayoff);
			}

		}

	}

	public static double[] testingAvgAbstraction( double[] numResources, int[][] gamedata, 
			double[][] abstractgame, StrategyMapping strmap,
			int numCluster, boolean isequalcluster, double origexp, int[] clustersizes) throws Exception
	{


		FastSecurityGameSolver fsgs = new FastSecurityGameSolver();
		StructuredSecurityGame samplegame = genStructuredSecurityGame(gamedata, 1, new double[]{numResources[0]});
		double resAbstract = getNumResourceAbstractGame(gamedata.length, numResources, numCluster);
		StructuredSecurityGame absgame = genStructuredSecurityGame(abstractgame, 1, new double[]{resAbstract} );
		Long startTimeOrigami = System.currentTimeMillis();
		SSGSolution sampleSolutionOrigmai = fsgs.solveGame(absgame, numResources[0], clustersizes);
		Long endTimeOrigami = System.currentTimeMillis();
		//checkNormality(numResources[0], clustersizes.length, sampleSolutionOrigmai);
		//printSStrategy(coverage);
		//System.out.println("Using abstractions : ");
		//willife
		//System.out.println("average ");
		SSGSolution originalstr = buildOriginalSGStrategy(sampleSolutionOrigmai, strmap, numResources[0], resAbstract, isequalcluster);
		checkNormality(numResources[0], gamedata.length, originalstr);
		originalstr.computeExpectedPayoffs(samplegame);
		//System.out.println("Result: " + originalstr);
		double absdefexppayoff = originalstr.getDefenderPayoff();
		double solquality = KmeanClustering.calcLineDist(new double[]{origexp, absdefexppayoff});
		return new double[]{solquality};


		//printSStrategy(originalstr);


		//System.out.print("hi");
	}

	private static SSGSolution buildOriginalSGStrategy(
			SSGSolution abstractsolution, StrategyMapping strmap, double numOrigResources,
			double resAbstract, boolean isequalcluster) 
	{

		int numberoforigtargets = strmap.getSecuritygamedata().length;
		double[] originalstr = new double[numberoforigtargets];
		//double restotargetratio = (double)numOrigResources/numberoforigtargets;

		List<Integer>[] clusterfortargets = strmap.getClusterfortargets();

		for(int abstarget = 0; abstarget< abstractsolution.getDefenderProbs(0).length; abstarget++)
		{
			
			for(Integer origtarget: clusterfortargets[abstarget])
			{
				if(!isequalcluster)
				{
					originalstr[origtarget] = abstractsolution.getProb(abstarget, 0);
				}
				else
				{
					originalstr[origtarget] = abstractsolution.getProb(abstarget, 0);
				}

			}
		}

		//int numResource = 5;
/*
		double sum = checkNormality(numOrigResources, originalstr);
		if(!(Math.abs(numOrigResources-sum)<.2))
		{



			*//**
			 * count the number of cov less than 1
			 *//*
			int countcovltone = 0;
			double lefttospill = 1.0;

			for(double x: originalstr)
			{
				if(x<1)
				{
					countcovltone++;
				}
			}
			double percov = 1.0/countcovltone;
			for(int index = 0; index<originalstr.length; index++)
			{
				if(originalstr[index]<1)
				{

					originalstr[index] += percov;
					lefttospill -= percov;
					countcovltone--;
					if(originalstr[index]>1)
					{
						lefttospill +=  (originalstr[index] -1);
						originalstr[index] = 1;
						percov = lefttospill / countcovltone;
					}
				}
			}
		}


		checkNormality(numOrigResources, originalstr);
*/
		SSGSolution solution = new SSGSolution(numberoforigtargets, 1);

		for (int t = 0; t < numberoforigtargets; t++) {
			solution.setProb(t, 0, originalstr[t]);
		}
		return solution;
	}

	public static double[] wildlifeAbstraction( double[] numResources, int[][] originalgame, 
			int[][][] abstractgame, StrategyMapping strmap ,
			int numCluster, boolean isequalcluster, double origdefexppayoff, int[] clustersizes) throws Exception
	{

		long totalTimeOrigmai = 0;
		FastSecurityGameSolver fsgs = new FastSecurityGameSolver();
		StructuredSecurityGame samplegame = genStructuredSecurityGame(originalgame, 1, new double[]{numResources[0]});
		double resAbstract = getNumResourceAbstractGame(originalgame.length, numResources, numCluster);
		IntervalSecurityGame issg = IntervalSecurityGame.generateAbstractIntervalSecurityGame(abstractgame,resAbstract);
		PolynomialIntervalSolver psolver = new PolynomialIntervalSolver();
		double [] coverage = psolver.solve(issg, (int)numResources[0], clustersizes);
		//System.out.println("wildlife");
		//checkNormality(resAbstract, coverage);
		//printSStrategy(coverage);
		//System.out.println("Using abstractions : ");
		SSGSolution revmappedstr = buildOriginalSGStrategy(coverage, strmap, numResources[0], isequalcluster);
		//checkNormality(numResources[0], originalgame.length, revmappedstr);
		revmappedstr.computeExpectedPayoffs(samplegame);
		//System.out.println("Result: " + revmappedstr);
		double absdefexppayoff = revmappedstr.getDefenderPayoff();
		double solquality = KmeanClustering.calcLineDist(new double[]{origdefexppayoff, absdefexppayoff});
		return new double[]{solquality};
		//printSStrategy(originalstr);


		//System.out.print("hi");
	}

	public static double[][] makeAbstractAvgSG(StrategyMapping strmap)
	{
		List<Integer>[] clusterfortargets = strmap.getClusterfortargets(); 
		double[][] abstractsecuritygame = new double[clusterfortargets.length][4];
		int[][] securitygamedata = strmap.getSecuritygamedata();
		for(int clusterindex = 0; clusterindex<clusterfortargets.length; clusterindex++)
		{
			for(int rewardindex = 0; rewardindex<4; rewardindex++)
			{
				double tmppayoffsum = 0;
				for(Integer target : clusterfortargets[clusterindex])
				{

					tmppayoffsum += securitygamedata[target][rewardindex];
				}
				tmppayoffsum = tmppayoffsum/clusterfortargets[clusterindex].size();
				abstractsecuritygame[clusterindex][rewardindex] = tmppayoffsum;
			}
		}

		return abstractsecuritygame;
	}


	public static double[][] makeAVGPayoffAbstractSecurityGame(StrategyMapping strmap)
	{
		List<Integer>[] clusterfortargets = strmap.getClusterfortargets(); 
		double[][] abssecuritygame = new double[clusterfortargets.length][4];
		int[][] securitygamedata = strmap.getSecuritygamedata();
		for(int clusterindex = 0; clusterindex<clusterfortargets.length; clusterindex++)
		{
			for(int outcome = 0; outcome< 4; outcome++)
			{
				double sum = 0;
				for(Integer target : clusterfortargets[clusterindex])
				{
					sum += securitygamedata[target][outcome];
				}
				sum /= clusterfortargets[clusterindex].size();
				abssecuritygame[clusterindex][outcome] = sum;
			}
		}
		return abssecuritygame;
	}

	public static int[][][] makeAbstractIntervalSecurityGame(StrategyMapping strmap)
	{
		List<Integer>[] clusterfortargets = strmap.getClusterfortargets(); 
		int[][][] securitygame = new int[clusterfortargets.length][4][2];
		int[][] securitygamedata = strmap.getSecuritygamedata();
		int dminr = Integer.MAX_VALUE;
		int dmaxr = Integer.MIN_VALUE;

		int dminp = Integer.MAX_VALUE;
		int dmaxp = Integer.MIN_VALUE;

		int aminr = Integer.MAX_VALUE;
		int amaxr = Integer.MIN_VALUE;

		int aminp = Integer.MAX_VALUE;
		int amaxp = Integer.MIN_VALUE;
		for(int clusterindex = 0; clusterindex<clusterfortargets.length; clusterindex++)
		{
			dminr = Integer.MAX_VALUE;
			dmaxr = Integer.MIN_VALUE;

			dminp = Integer.MAX_VALUE;
			dmaxp = Integer.MIN_VALUE;

			aminr = Integer.MAX_VALUE;
			amaxr = Integer.MIN_VALUE;

			aminp = Integer.MAX_VALUE;
			amaxp = Integer.MIN_VALUE;
			for(Integer target : clusterfortargets[clusterindex])
			{
				/**
				 * min condition for defender reward
				 */
				if(securitygamedata[target][0]<0)
				{
					//System.out.print("OOO");
				}
				if(dminr > securitygamedata[target][0])
				{
					dminr = securitygamedata[target][0];
				}
				/**
				 * max condition for defender reward
				 */
				if(dmaxr < securitygamedata[target][0])
				{
					dmaxr = securitygamedata[target][0];
				}
				/**
				 * min condition for defender penalty
				 */
				if(dminp > securitygamedata[target][1])
				{
					dminp = securitygamedata[target][1];
				}
				/**
				 * max condition for defender penalty
				 */
				if(dmaxp < securitygamedata[target][1])
				{
					dmaxp = securitygamedata[target][1];
				}

				/**
				 * min condition for defender reward
				 */
				if(aminr > securitygamedata[target][2])
				{
					aminr = securitygamedata[target][2];
				}
				/**
				 * max condition for defender reward
				 */
				if(amaxr < securitygamedata[target][2])
				{
					amaxr = securitygamedata[target][2];
				}
				/**
				 * min condition for defender penalty
				 */
				if(aminp > securitygamedata[target][3])
				{
					aminp = securitygamedata[target][3];
				}
				/**
				 * max condition for defender penalty
				 */
				if(securitygamedata[target][3]>0)
				{
					//System.out.print("OOOshh");
				}
				if(amaxp < securitygamedata[target][3])
				{
					amaxp = securitygamedata[target][3];
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
		}


		return securitygame;
	}





	public static double checkNormality(double numResource, double[] coverage)
	{
		double sum = 0;
		for(double x: coverage)
		{
			sum += x;
		}
		//System.out.println("sum  : "+ sum);
		//sum = Math.ceil(sum);
		//System.out.println("after ceil sum  : "+ sum);

		return sum;
	}


	public static double[] testingPerfectRationalMMR(double[] numResources,
			int[][] originalgame, int[][][] abstractgame, StrategyMapping strmap,
			int numCluster, boolean isequalcluster, double origexp, int[] clustersizes) throws IOException, LPSolverException, MatlabConnectionException, MatlabInvocationException 
	{
		//List<Integer> clusters[] = strmap.getClusterfortargets();

		int numTargets = originalgame.length;
		double absexp = 0;
		FastSecurityGameSolver fsgs = new FastSecurityGameSolver();
		StructuredSecurityGame samplegame = genStructuredSecurityGame(originalgame, 1, new double[]{numResources[0]});
		int numTargetsAbs = abstractgame.length;
		//strmap.printSecurityGameMapping();
		//printAbstractSecurityGame(abstractgame);
		try 
		{
			lpWrapper.Configuration.loadLibrariesCplex();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.gc();
		List<Target> targetList = new ArrayList<Target>();
		double[] advPayoffLB = new double[2 * numTargetsAbs];
		double[] advPayoffUB = new double[2 * numTargetsAbs];
		double[] weight = new double[numTargetsAbs];
		for(int i = 0; i < numTargetsAbs; i++)
		{
			weight[i] = clustersizes[i];
//			weight[i] = 1;
			//System.out.print(weight[i] + "\t");
		}
		for(int target = 0; target < numTargetsAbs; target++)
		{	
			//List<PayoffStructure> payoffList = null;
			double attRewardLB;
			double attPenaltyLB;
			double attRewardUB;
			double attPenaltyUB;

			attRewardLB = abstractgame[target][2][0]; // reward lower bound //payoff[target + 2 * numTargets];
			attRewardUB = abstractgame[target][2][1];//reward upper bound //attRewardLB + intervalRand[target];
			attPenaltyLB = abstractgame[target][3][0];//penalty lower bound //payoff[target + 3 * numTargets] - interval;
			attPenaltyUB = abstractgame[target][3][1]; //penalty upper bound //attPenaltyLB + intervalRand[target + numTargets];

			advPayoffLB[target] = attRewardLB; 
			advPayoffLB[target + numTargetsAbs] = attPenaltyLB; 					
			advPayoffUB[target] =  attRewardUB;            
			advPayoffUB[target + numTargetsAbs] =  attPenaltyUB; 
			AttBoundStructure attBoundStructure = new AttBoundStructure(advPayoffLB[target]
					, advPayoffLB[target + numTargetsAbs], advPayoffUB[target], advPayoffUB[target + numTargetsAbs]);
			Target t = new Target(target, null, attBoundStructure);
			targetList.add(t);
		}
		double numResourceAbs = getNumResourceAbstractGame(numTargets, new int[]{(int)numResources[0]}, numCluster);

		//SUQRAdversary adversary = new SUQRAdversary(0, -9.85, 0.45, 0.32, 1.0);
		//MMR mmr = new MMR(targetList, adversary, numResourceAbs, numSamples, isZeroSum);
		MMRZeroSumNew2 mmr = new MMRZeroSumNew2(numTargetsAbs, numResources[0], advPayoffLB, advPayoffUB);
		mmr.setWeight(weight);
		mmr.loadProblem();
		mmr.solve();
		//mmr.deletePayoffConstraint();
		double [] coverage = mmr.getOptDefCov();

		//double[] coverage = mmr.getOptCov();
		//checkNormality(numResourceAbs, coverage);

		//printSStrategy(coverage);
		//System.out.println("Using abstractions : ");
		SSGSolution originalstr = buildOriginalSGStrategyMMR(coverage, strmap, numResources[0], isequalcluster);
		//checkNormality(numResources[0], numTargets, originalstr);

		originalstr.computeExpectedPayoffs(samplegame);
		absexp = originalstr.getDefenderPayoff();

		//System.out.println("Result: " + originalstr);

		mmr.end();
		//	}
		//System.out.println(" original game Avg def payoff "+ origexp/NUM_ITER);
		//System.out.println(" abstraction Avg def payoff "+ absexp/NUM_ITER);
		double solquality = KmeanClustering.calcLineDist(new double[]{origexp, absexp});
		return new double[]{solquality};


	}


	private static SSGSolution buildOriginalSGStrategyMMR(double[] coverage,
			StrategyMapping strmap, double numOrigResource, boolean isequalcluster) {
		
		int numberoforigtargets = strmap.getSecuritygamedata().length;
		double[] originalstr = new double[numberoforigtargets];
		//double restotargetratio = (double)numOrigResource/numberoforigtargets;
		List<Integer>[] clusterfortargets = strmap.getClusterfortargets();
		for(int abstarget = 0; abstarget< coverage.length; abstarget++)
		{
			//double ratio = -.1;
			
			for(Integer origtarget: clusterfortargets[abstarget])
			{
				if(!isequalcluster)
				{
					originalstr[origtarget] = coverage[abstarget];
				}
				else
				{
					originalstr[origtarget] = coverage[abstarget];
				}

			}
		}
		double sum = checkNormality(numOrigResource, originalstr);
		
		SSGSolution solution = new SSGSolution(numberoforigtargets, 1);

		for (int t = 0; t < numberoforigtargets; t++) {
			solution.setProb(t, 0, originalstr[t]);
		}
		return solution;

	}



	public static double[] testingBoundedRationalMMR(double[] numResources, int[][] gamedata, int numCluster, boolean isequalcluster  ) throws Exception
	{

		int numSamples = 1;
		//double interval = 8.0;
		boolean isZeroSum = true;
		int numTargets = gamedata.length;
		//int cov = 0;
		//int payoffIndex = 0;
		double origexp = 0;
		double absexp = 0;
		int NUM_ITER = 1;
		//for(int itr = 0; itr<NUM_ITER; itr++)
		//{

		origexp = originalGameTesing(numResources, gamedata);
		//convertToZeroSum(gamedata);




		//long totalTimeOrigmai = 0;
		FastSecurityGameSolver fsgs = new FastSecurityGameSolver();
		StructuredSecurityGame samplegame = genStructuredSecurityGame(gamedata, 1, new double[]{numResources[0]});
		//List<Integer>[] clusteredtargets = KmeanClustering.clusterTargets(numCluster, gamedata);
		List<Integer>[] clusteredtargets = KmeanClustering.getBestclusteredTargets(numCluster, gamedata, true);
		int maxinterval = KmeanClustering.getMaxInterval(clusteredtargets, gamedata);

		StrategyMapping strmap = new StrategyMapping(clusteredtargets, numCluster, gamedata);
		int[][][] abstractgame = makeAbstractIntervalSecurityGame(strmap);
		//strmap.printSecurityGameMapping();
		//convertToZeroSum(abstractgame);
		//printAbstractSecurityGame(abstractgame);
		//System.out.print("hi");
		try {
			lpWrapper.Configuration.loadLibrariesCplex();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		//Random rand = new Random();
		System.gc();
		List<Target> targetList = new ArrayList<Target>();
		for(int target = 0; target < numCluster; target++)
		{	
			List<PayoffStructure> payoffList = null;
			double attRewardLB;
			double attPenaltyLB;
			double attRewardUB;
			double attPenaltyUB;
			attRewardLB = abstractgame[target][2][0];//payoff[target + 2 * numTargets];
			attRewardUB = abstractgame[target][2][1];//attRewardLB + intervalRand[target];
			attPenaltyLB = abstractgame[target][3][0];//payoff[target + 3 * numTargets] - interval;
			attPenaltyUB = abstractgame[target][3][1];//attPenaltyLB + intervalRand[target + numTargets];
			AttBoundStructure attBoundStructure = new AttBoundStructure(attRewardLB, attPenaltyLB, attRewardUB, attPenaltyUB);
			if(!isZeroSum)
			{
				payoffList = new ArrayList<PayoffStructure>();
				payoffList.add(new PayoffStructure(abstractgame[target][0][0], abstractgame[target][1][0], abstractgame[target][2][0], abstractgame[target][3][0]));
			}
			Target t = new Target(target, payoffList, attBoundStructure);
			targetList.add(t);
		}
		double numResourceAbs = getNumResourceAbstractGame(numTargets, new int[]{(int)numResources[0]}, numCluster);
		SUQRAdversary adversary = new SUQRAdversary(0, -9.85, 0.45, 0.32, 1.0);
		MMR mmr = new MMR(targetList, adversary, numResourceAbs, numSamples, isZeroSum);
		mmr.loadProblem();
		mmr.solve();
		mmr.deletePayoffConstraint();
		double [] coverage = mmr.getDefCov();

		//double[] coverage = mmr.getOptCov();
		//checkNormality(numResourceAbs, coverage);

		//printSStrategy(coverage);
		//System.out.println("Using abstractions : ");
		SSGSolution originalstr = buildOriginalSGStrategy(coverage, strmap, numResources[0], isequalcluster);
		originalstr.computeExpectedPayoffs(samplegame);
		absexp = originalstr.getDefenderPayoff();

		//System.out.println("Result: " + originalstr);

		mmr.end();
		//	}
		//System.out.println(" original game Avg def payoff "+ origexp/NUM_ITER);
		//System.out.println(" abstraction Avg def payoff "+ absexp/NUM_ITER);
		double solquality = KmeanClustering.calcLineDist(new double[]{origexp, absexp});
		return new double[]{solquality, maxinterval};

	}

	private static void convertToZeroSum(int[][][] abstractgame) {

		for(int target = 0; target< abstractgame.length; target++)
		{
			abstractgame[target][0][0] = -abstractgame[target][3][0];
			abstractgame[target][0][1] = -abstractgame[target][3][1];

			abstractgame[target][1][0] = -abstractgame[target][2][0];
			abstractgame[target][1][1] = -abstractgame[target][2][1];


		}

	}

	/**
	 * converts the attacker's payoffs to opposite of defender
	 * @param gamedata
	 */
	public static void convertToZeroSum(int[][] gamedata)
	{
		for(int target = 0; target<gamedata.length; target++)
		{
			gamedata[target][2] = -(gamedata[target][1]);
			gamedata[target][3] = -(gamedata[target][0]);

		}
	}

	public static SSGSolution buildOriginalSGStrategy(double[] coverage, 
			StrategyMapping strmap, double numOrigResource, boolean isclusterequal) 
	{

		int numberoforigtargets = strmap.getSecuritygamedata().length;
		double[] originalstr = new double[numberoforigtargets];
		//double restotargetratio = (double)numOrigResource/numberoforigtargets;
		List<Integer>[] clusterfortargets = strmap.getClusterfortargets();
		for(int abstarget = 0; abstarget< coverage.length; abstarget++)
		{
			//double ratio = -.1;
			if(!isclusterequal)
			{
				//ratio = (double)numberoforigtargets/clusterfortargets[abstarget].size();
				//ratio *= numOrigResource;
				//DecimalFormat df = new DecimalFormat("#.#"); 
				//df.setRoundingMode(RoundingMode.CEILING);
				//ratio = Double.parseDouble(df.format(ratio));
			}
			for(Integer origtarget: clusterfortargets[abstarget])
			{
				if(!isclusterequal)
				{
					originalstr[origtarget] = coverage[abstarget];
				}
				else
				{
					originalstr[origtarget] = coverage[abstarget];
				}

			}
		}
		double sum = checkNormality(numOrigResource, originalstr);
		/*if(!(Math.abs(numOrigResource-sum)<.2))
		{
			*//**
			 * count the number of cov less than 1
			 *//*
			int countcovltone = 0;
			double lefttospill = 1.0;

			for(double x: originalstr)
			{
				if(x<1)
				{
					countcovltone++;
				}
			}
			double percov = 1.0/countcovltone;
			for(int index = 0; index<originalstr.length; index++)
			{
				if(originalstr[index]<1)
				{

					originalstr[index] += percov;
					lefttospill -= percov;
					countcovltone--;
					if(originalstr[index]>1)
					{
						lefttospill +=  (originalstr[index] -1);
						originalstr[index] = 1;
						percov = lefttospill / countcovltone;
					}
				}
			}
		}*/


		checkNormality(numOrigResource, originalstr);

		SSGSolution solution = new SSGSolution(numberoforigtargets, 1);

		for (int t = 0; t < numberoforigtargets; t++) {
			solution.setProb(t, 0, originalstr[t]);
		}
		return solution;
	}


	private static void printSStrategy(double[] originalstr) {

		System.out.println();
		for(double x : originalstr)
		{
			System.out.print(x + " ");

		}
		System.out.println();

	}


	public static StructuredSecurityGame genStructuredSecurityGame(int[][] gamedata, int nDefenseTypes, double[] nDefenders)
	{
		StructuredSecurityGame ssg = new StructuredSecurityGame(gamedata.length, nDefenseTypes);
		ssg.setNDefenders(nDefenders);
		for (int target = 0; target < gamedata.length; target++) 
		{
			SSGPayoffs tmp = new SSGPayoffs(false);

			//-------------------------------------------
			//Attacker and Defender's covered payoff is set to zero here
			//-------------------------------------------
			//double r1 = Math.random() * -100;
			double r1=gamedata[target][3];
			double r3=gamedata[target][0];

			double r2 = gamedata[target][2];
			//double r3 = Math.random() * 100;
			double r4 = gamedata[target][1];

			// System.out.print("R2");
			//System.out.println( r2);

			//---------------------------------------------
			tmp.setAttackerCoveredPayoff(r1);
			tmp.setAttackerUncoveredPayoff(r2);
			tmp.setDefenderCoveredPayoff(r3);
			tmp.setDefenderUncoveredPayoff(r4);
			ssg.setPayoffs(target, tmp);
		}
		Random r = new Random();
		// first, assign at least 1 resource to each schedule (so there are no uncoverable schedules)
		for (int target = 0; target < gamedata.length; target++) 
		{
			ssg.setDefenseCapability(target, r.nextInt(nDefenseTypes), true);
		}

		// next, randomly add additional resource capabilities up to the density requested
		int cnt = gamedata.length;
		int max = (int)Math.round((double)(gamedata.length * nDefenseTypes) * 1.0d);
		while (cnt < max) 
		{
			int tmp1 = r.nextInt(gamedata.length);
			int tmp2 = r.nextInt(nDefenseTypes);
			if (ssg.getDefenseCapability(tmp1, tmp2)) continue;
			ssg.setDefenseCapability(tmp1, tmp2, true);
			cnt++;
		}
		return ssg;


	}

	public static StructuredSecurityGame genStructuredSecurityGame(double[][] gamedata, int nDefenseTypes, double[] nDefenders)
	{
		StructuredSecurityGame ssg = new StructuredSecurityGame(gamedata.length, nDefenseTypes);
		ssg.setNDefenders(nDefenders);
		for (int target = 0; target < gamedata.length; target++) 
		{
			SSGPayoffs tmp = new SSGPayoffs(false);

			//-------------------------------------------
			//Attacker and Defender's covered payoff is set to zero here
			//-------------------------------------------
			//double r1 = Math.random() * -100;
			double r1=gamedata[target][3];
			double r3=gamedata[target][0];

			double r2 = gamedata[target][2];
			//double r3 = Math.random() * 100;
			double r4 = gamedata[target][1];

			// System.out.print("R2");
			//System.out.println( r2);

			//---------------------------------------------
			tmp.setAttackerCoveredPayoff(r1);
			tmp.setAttackerUncoveredPayoff(r2);
			tmp.setDefenderCoveredPayoff(r3);
			tmp.setDefenderUncoveredPayoff(r4);
			ssg.setPayoffs(target, tmp);
		}
		Random r = new Random();
		// first, assign at least 1 resource to each schedule (so there are no uncoverable schedules)
		for (int target = 0; target < gamedata.length; target++) 
		{
			ssg.setDefenseCapability(target, r.nextInt(nDefenseTypes), true);
		}

		// next, randomly add additional resource capabilities up to the density requested
		int cnt = gamedata.length;
		int max = (int)Math.round((double)(gamedata.length * nDefenseTypes) * 1.0d);
		while (cnt < max) 
		{
			int tmp1 = r.nextInt(gamedata.length);
			int tmp2 = r.nextInt(nDefenseTypes);
			if (ssg.getDefenseCapability(tmp1, tmp2)) continue;
			ssg.setDefenseCapability(tmp1, tmp2, true);
			cnt++;
		}
		return ssg;


	}

	public static int[][] createdummyData(int numTarget, int gap, int clustersize)
	{
		int[][] data = new int[numTarget][4];
		int counter = 0;
		int val = -1;
		int increment = 0;
		int startval = randInt(1, 5);
		for(int i=0; i<numTarget; i++)
		{

			//KmeanClustering.randInt(0, 10);
			if(counter==0)
			{
				val = startval+ increment;//KmeanClustering.randInt(0, 100);
			}
			//for(int j=0; j<4; j++)
			{
				data[i][0] = val;
				data[i][1] = -(val-10);
				if(data[i][1]>0)
				{
					data[i][1] = -data[i][1];
				}
				data[i][2] = -data[i][1];
				data[i][3] = -data[i][0];
				if(data[i][2]<data[i][3])
				{
					data[i][3] -=  Math.abs(data[i][3]-data[i][2]);
				}

			}
			counter++;
			if(counter==(clustersize))
			{
				counter = 0;
				increment += gap;
			}

		}
		return data;
	}

	private static void printAbstractSecurityGame(int[][][] abstractgame) {

		for(int i=0; i<abstractgame.length; i++)
		{
			for(int j=0; j<4; j++)
			{
				System.out.print("[");
				for(int k=0; k<2; k++)
				{
					System.out.print(abstractgame[i][j][k]);
					if(k==0)
					{
						System.out.print(",");
					}
				}
				System.out.print("] ");
			}
			System.out.println();
		}

	}

	public static int[][] parseSecurityGameFile(String filename, int gamenumber)
	{
		int[][] gamedata;
		File csvtrainData = new File(filename);
		try 
		{
			CSVParser parser = CSVParser.parse(csvtrainData, StandardCharsets.US_ASCII, CSVFormat.EXCEL);
			int gamecounter = 0;
			for (CSVRecord csvRecord : parser) 
			{
				if(gamenumber==gamecounter)
				{
					int numberoftargets = csvRecord.size()/4;
					gamedata = new int[numberoftargets][4];
					Iterator<String> itr = csvRecord.iterator();
					int targetcounter = 0;
					int singlepayoffcounter = 0;
					while(itr.hasNext())
					{
						gamedata[targetcounter][singlepayoffcounter] = Integer.parseInt(itr.next());
						targetcounter++;
						//singlepayoffcounter++;
						if(targetcounter==numberoftargets)
						{
							targetcounter=0;
							singlepayoffcounter++;
							//targetcounter++;
						}
					}
					return gamedata;
				}
				gamecounter++;
			}

		}
		catch (IOException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	public static double getNumResourceAbstractGame(int numTargets,
			double[] numResources, int numCluster) {
		double originalratio = (double)numResources[0]/numTargets;
		double numResAbstract = numCluster*originalratio;
		//numResAbstract = Math.floor(numResAbstract);


		return numResAbstract;
	}

	public static double getNumResourceAbstractGame(int numTargets,
			int[] numResources, int numCluster) {
		double originalratio = (double)numResources[0]/numTargets;
		double numResAbstract = numCluster*originalratio;
		//numResAbstract = Math.floor(numResAbstract);
		return numResAbstract;
	}

	public static int randInt(int min, int max) {

		// Usually this should be a field rather than a method variable so
		// that it is not re-seeded every call.
		// rand = new Random();

		// nextInt is normally exclusive of the top value,
		// so add 1 to make it inclusive
		int randomNum = rand.nextInt((max - min) + 1) + min;

		return randomNum;
	}

	public static void addNoise(int[][] gamedata, int noiselimit, int increment, boolean rand) {


		int noiseval = 0;
		Random random = new Random();
		for(int i=0; i<gamedata.length; i++)
		{
			if(rand)
			{
				noiseval = randInt(1,noiselimit);
			}
			else
			{
				noiseval += randInt(1,increment);
			}
			for(int j=0; j<gamedata[i].length; j++)
			{
				//if(i%2==0)
				if(j==1 || j==3)
				{
					gamedata[i][j] -= noiseval;
				}
				else
				{
					gamedata[i][j] += noiseval;
				}


			}
			//noiseval = 0;
		}

	}







	public static double[] wildlifeAbstraction(double[] numResources,
			int[][] gamedata, double[][] avgabstractgame,
			StrategyMapping strmap, int numCluster, boolean isequalcluster,  double originalexppayoff) {
		// TODO Auto-generated method stub
		return null;
	}





}
