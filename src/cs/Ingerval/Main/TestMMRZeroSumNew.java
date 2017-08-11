package cs.Ingerval.Main;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.Vector;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import models.AttBoundStructure;
import models.Target;
import perfectMMR.MMRZeroSumNew2;

public class TestMMRZeroSumNew {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws MatlabInvocationException 
	 * @throws MatlabConnectionException 
	 */
	public static void main(String[] args) throws IOException, MatlabConnectionException, MatlabInvocationException {
		// TODO Auto-generated method stub
		try {
			lpWrapper.Configuration.loadLibrariesCplex();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		int nTargets = 1000;
		int nRes = 5;
		
		double interval = 4.0;
		int cov = 0;
		int nSample = 7;
//		SUQRAdversary adversary = new SUQRAdversary(0, -2.0, 0.5, 0.1, 1.0);
//		double[] intervalRand = {7.857129711
//				,5.866335365
//				,4.931078642
//				,2.031476809
//				,1.362847066
//				,2.180880732
//				,7.815141337
//				,2.757595423
//				,5.808231719
//				,2.790687848
//				,3.867027912
//				,0.961492257
//				,3.55790175
//				,6.596312162
//				,3.273376195
//				,2.470197521
//				,2.283837811
//				,0.981378746
//				,1.336572904
//				,4.691145798
//				,0.920223476
//				,6.958494397
//				,7.828782255
//				,0.596760744
//				,4.080052995
//				,2.929820923
//				,0.420656903
//				,1.742278995
//				,0.059728635
//				,2.440484569
//				,5.017987055
//				,6.933273768
//				,7.284432629
//				,6.592058505
//				,1.191099315
//				,3.07338888
//				,1.023615634
//				,4.322415086
//				,7.483382328
//				,7.88686925
//				,4.432514391
//				,5.547312476
//				,7.236880373
//				,2.624527137
//				,1.452165126
//				,7.198135591
//				,2.674455583
//				,4.525849578
//				,3.82694855
//				,2.479281361
//				,1.329583352
//				,5.724015643
//				,4.475262146
//				,1.20773423
//				,0.516407983
//				,1.444653524
//				,2.675475244
//				,5.88644558
//				,2.384248705
//				,5.763228062
//				,1.005891717
//				,0.202188243
//				,4.202306276
//				,7.850016386
//				,1.703999743
//				,5.92762561
//				,5.871887676
//				,7.545687968
//				,3.086856442
//				,6.746391437
//				,7.265020868
//				,2.274568848
//				,6.901015422
//				,7.309043534
//				,3.839906049
//				,5.716733474
//				,6.946966868
//				,4.866813348
//				,2.085197753
//				,6.324129212};
		
		for(nTargets = 25; nTargets <= 50; nTargets += 25)
		{
			System.out.println("Targets: " + nTargets);
			double[] advPayoffLB = new double[2 * nTargets];
			double[] advPayoffUB = new double[2 * nTargets];
			for(cov = 0; cov <= 10; cov += 2)
			{
				System.out.println("Cov: " + cov);
				double[][] optDefCovs = new double[nSample][];
				double[][] optAdvPayoffs = new double[nSample][];
				double[] optRegret = new double[nSample];
				double[] runtime = new double[nSample];
				double[] trueRegret = new double[nSample];
				Vector<int[]> payoffs = loadData(nTargets, cov);
		
				for(int idx = 0; idx < nSample; idx++)
				{
					System.gc();
					List<Target> targetList = new ArrayList<Target>();
					String intervalPath = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/";
					String intervalFileName = "/" + nTargets + "T" + nRes + "R/MMRInterval" + cov + "_" + idx + ".csv";
//					double[] intervalRand = loadInterval(intervalPath + intervalFileName, nTargets);
					
//					System.out.println("Bound: ");
					for(int i = 0; i < nTargets; i++)
					{
						advPayoffLB[i] = payoffs.get(idx)[i + 2 * nTargets];
						advPayoffLB[i + nTargets] = payoffs.get(idx)[i + 3 * nTargets] - interval;					
						advPayoffUB[i] = payoffs.get(idx)[i + 2 * nTargets] + interval;
						advPayoffUB[i + nTargets] = advPayoffLB[i + nTargets] + interval;
//						System.out.println(advPayoffLB[i] + "\t" + advPayoffUB[i] + "\t" + advPayoffLB[nTargets + i]
//								+ "\t" + advPayoffUB[nTargets + i]);
						AttBoundStructure attBoundStructure = new AttBoundStructure(advPayoffLB[i]
								, advPayoffLB[i + nTargets], advPayoffUB[i], advPayoffUB[i + nTargets]);
						Target t = new Target(i, null, attBoundStructure);
						targetList.add(t);
					}
					
					MMRZeroSumNew2 mmr = new MMRZeroSumNew2(nTargets, nRes, advPayoffLB, advPayoffUB);
			//		MMRZeroSumNew mmr = new MMRZeroSumNew(nTargets, nRes, advPayoffLB, advPayoffUB);
					mmr.loadProblem();
					mmr.solve();
					
					optDefCovs[idx] = mmr.getOptDefCov();
//					Map<Target, Double> defCov = new HashMap<Target, Double>();
//					for(Target t : targetList)
//					{
//						defCov.put(t, optDefCovs[idx][t.id]);
//					}
					optAdvPayoffs[idx] = mmr.getOptAdvPayoff();
					optRegret[idx] = mmr.getOptRegret();
					runtime[idx] = mmr.getRuntime();
					
//					MaxRegretZeroSum mr = new MaxRegretZeroSum(targetList, defCov, adversary, nRes, 40);
//					mr.solve();
//					System.out.println("\nMax Regret:" + mr.getMaxRegret());
//					trueRegret[idx] = mr.getMaxRegret();
//					mr.end();
					mmr.end();
//					defCov.clear();
					
				}
				// Save results
				int nRes_temp = (int)nRes;
				String runtime_name = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/PerfectMMR/NewAlg/" + nTargets + "T" + nRes_temp + 
						"R/runtime" + cov + ".csv";
				try{
					FileOutputStream runtimeFileStream = new FileOutputStream(runtime_name);
					PrintStream runtimeStream = new PrintStream(runtimeFileStream);
					for(int j = 0; j < nSample; j++)
					{
						runtimeStream.println(runtime[j]);
					}
					runtimeStream.close();
				}catch (FileNotFoundException e1)
				{
					e1.printStackTrace();
				}
				
				// Saving results
				String obj_name = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/PerfectMMR/NewAlg/" + nTargets + "T" + nRes_temp + 
						"R/obj" + cov + ".csv";
				try{
					FileOutputStream objFileStream = new FileOutputStream(obj_name);
					PrintStream objStream = new PrintStream(objFileStream);
					for(int j = 0; j < nSample; j++)
					{
						objStream.println(optRegret[j]);
					}
					objStream.close();
				}catch (FileNotFoundException e1)
				{
					e1.printStackTrace();
				}
				
//				// Saving results
//				String regret_name = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/PerfectMMR/NewAlg/" + nTargets + "T" + nRes_temp + 
//						"R/regret" + cov + ".csv";
//				try{
//					FileOutputStream regretFileStream = new FileOutputStream(regret_name);
//					PrintStream regretStream = new PrintStream(regretFileStream);
//					for(int j = 0; j < nSample; j++)
//					{
//						regretStream.println(trueRegret[j]);
//					}
//					regretStream.close();
//				}catch (FileNotFoundException e1)
//				{
//					e1.printStackTrace();
//				}
				
				// Saving results
				String def_name = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/PerfectMMR/NewAlg/" + nTargets + "T" + nRes_temp + 
						"R/strategy" + cov + ".csv";
				try{
					FileOutputStream defFileStream = new FileOutputStream(def_name);
					PrintStream defStream = new PrintStream(defFileStream);
					for(int j = 0; j < nSample; j++)
					{
						for(int k = 0; k < nTargets; k++)
						{
							defStream.print(optDefCovs[j][k] + ",");
						}
						defStream.println();
					}
					
					defStream.close();
				}catch (FileNotFoundException e1)
				{
					e1.printStackTrace();
				}
				
				String adv_name = "/Users/thanhnguyen/Documents/WORKS/UAV/JAVA/RESULTS/PerfectMMR/NewAlg/" + nTargets + "T" + nRes_temp + 
						"R/payoff" + cov + ".csv";
				try{
					FileOutputStream advFileStream = new FileOutputStream(adv_name);
					PrintStream advStream = new PrintStream(advFileStream);
					for(int j = 0; j < nSample; j++)
					{
						for(int k = 0; k < 2 * nTargets; k++)
						{
							advStream.print(optAdvPayoffs[j][k] + ",");
						}
						advStream.println();
					}
					
					advStream.close();
				}catch (FileNotFoundException e1)
				{
					e1.printStackTrace();
				}
			}
		}
		
	}
	public static double[] loadInterval(String filePath, int nTargets)
	{
		double[] interval = new double[2 * nTargets];
		FileInputStream in = null;
		FileReader fin = null;
		Scanner src = null;
		try {
			fin = new FileReader(filePath);
			src = new Scanner(fin);
			for(int i = 0; i < 2 * nTargets; i++)
				interval[i] = src.nextDouble();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.out.println("Couldn't open file for reading.");
			e.printStackTrace();
		}
		src.close();
		return interval;
	}
	public static Vector<int[]> loadData(int numTargets, int cov)
	{
		String payoffFile;
		if (cov < 10)
			payoffFile = "/Users/thanhnguyen/Documents/WORKS/UAV/GAMEGENERATION/" + numTargets + "Target/inputr-0." + cov + "00000.csv";
		else
			payoffFile = "/Users/thanhnguyen/Documents/WORKS/UAV/GAMEGENERATION/" + numTargets + "Target/inputr-1.000000.csv";
		Vector<int[]> payoffs = new Vector<int[]>();
		FileInputStream in = null;
		FileReader fin = null;
		Scanner src = null;
		try {
			fin = new FileReader(payoffFile);
			src = new Scanner(fin);
			while (src.hasNext()) {
				String line = src.nextLine();
				String[] values = line.split(",");
				int[] payoff = new int[4 * numTargets];
				for (int i = 0; i < 4 * numTargets; i++)
					payoff[i] = Integer.parseInt(values[i]);
				payoffs.add(payoff);

			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.out.println("Couldn't open file for reading.");
			e.printStackTrace();
		}
		src.close();
		return payoffs;
	}
}
