package ltarPractice;

//import java.io.FileWriter;
import java.io.FileWriter;
import java.io.IOException;

import src.tambe.usc.DOBSS.gameRepresentation.SSGPayoffs;
import src.tambe.usc.DOBSS.gameRepresentation.SSGSolution;
import src.tambe.usc.DOBSS.gameRepresentation.StructuredSecurityGame;
import src.tambe.usc.DOBSS.solvers.FastSecurityGameSolver;

public class LtarMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		test5();

	}
	
	private static void test5()
	 {

		   int nSamples = 1;
		   long totalTimeLtra = 0;
		   long totalTimeOrigmai=0;
		   int[] numTargets =   {10	 };
		   int[] numResources = {2};

		   FastSecurityGameSolver fsgs = new FastSecurityGameSolver();

		   for (int i = 0; i < numTargets.length; i++)		   
		   {
		     for (int sample = 0; sample < nSamples; sample++)
		     {
		       StructuredSecurityGame sampleGame = StructuredSecurityGame.genRandomGame(numTargets[i], 1, new double[]{numResources[i]});
		       
		       //---------------------------
		       //sampleGame.getPayoffs(i);
		       
		       SSGPayoffs[] payoffs=sampleGame.getPayoffs();
		       
		       //System.out.println("# of targets :"+numTargets[0]);
		      
		       
		       
//		       for (int j = 0; j < numTargets[0]; j++)
//		       {
//		      	 System.out.println("payoffs of target "+ j+":  " + payoffs[j].getAttackerUncoveredPayoff() );
//		      	
//		       } 
		       //---------------------------------
		       long startTimeLtra = System.currentTimeMillis(); 
		       //SSGSolution sampleSolution = fsgs.solveGame(sampleGame);
		       SSGSolution sampleSolution = fsgs.linearTimeSolveGame(sampleGame);

		       long endTimeLtra = System.currentTimeMillis();

		       
		       double coverageProb=0;
		       
		       for (int j = 0; j < numTargets[0]; j++)
		       {
		      	 //System.out.println("payoffs of target "+ j+":  " + payoffs[j].getAttackerUncoveredPayoff() );
		    	 //System.out.println ("Target "+j +  ": " + sampleSolution.getProb(3, 0));
		    	   
		    	   System.out.println ("Target "+j +  ": " + sampleSolution.getProb(j, 0));
		    	   double ci= sampleSolution.getProb(j, 0);
			       coverageProb = coverageProb + ci;
		       }
		       
		       System.out.println("Total Coverage Probability: "+ coverageProb);
		       		       
		     
		       long time = endTimeLtra-startTimeLtra;
		       totalTimeLtra += time;
		       System.out.println("Total Time LTRA : " + time);
		       
		       
		       /*---------for Origami----------------*/
		       
		       Long startTimeOrigami = System.currentTimeMillis();
		       SSGSolution sampleSolutionOrigmai = fsgs.solveGame(sampleGame);
		       Long endTimeOrigami = System.currentTimeMillis();

		       sampleSolutionOrigmai.computeExpectedPayoffs(sampleGame);
		       System.out.println("Result: " + sampleSolutionOrigmai);
		     
		       time = endTimeOrigami -startTimeOrigami;
		       totalTimeOrigmai += time;
		       
		      /*--------------------------------*/
	//		     }
	//			 System.out.println("Total time Ori: " + totalTime);
	//			 System.out.println("Ave time: " + totalTime/nSamples);
	//			 totalTime = 0;

		     }
		     
		     double avgTimeLtra=(double)totalTimeLtra/nSamples;
		     double avgTimeOrigami=(double)totalTimeOrigmai/nSamples;
		     
		     System.out.println("Total Time LTRA : " + totalTimeLtra);
		     System.out.println("Average time Ltra:" +avgTimeLtra + "\n" );
		     
		     
		     System.out.println("Total Time Origami: " + totalTimeOrigmai);
		     System.out.println("Average time origmai:" +avgTimeOrigami + "\n" );
		     
		     
		     
		    // System.out.println( Double.toString(totalTimeOrigmai));
		     
		     
		     try
		 	{
		 	    FileWriter writer = new FileWriter("C:\\Users\\csuser\\Desktop\\ltra_origami.csv");
		  
		 	    writer.append(Integer.toString(numTargets[0]));
		 	    writer.append(',');
			    writer.append(Double.toString(avgTimeLtra));
		 	    writer.append(',');
			    writer.append(Double.toString(avgTimeOrigami));
		 	    writer.append(',');
		 	    writer.append('\n');

		 	    //generate whatever data you want
		  
		 	    writer.flush();
		 	    writer.close();
		 	}
		 	catch(IOException e)
		 	{
		 	     e.printStackTrace();
		 	} 
		     }
		    
//			 System.out.println("Total time: " + totalTime);
//			 System.out.println("Ave time: " + totalTime/nSamples);
//			 totalTime = 0;
		  }

	 
	
} 


