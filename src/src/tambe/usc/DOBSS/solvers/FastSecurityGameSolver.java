package src.tambe.usc.DOBSS.solvers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import ltarPractice.Median;
import src.tambe.usc.DOBSS.gameRepresentation.SSGSolution;
import src.tambe.usc.DOBSS.gameRepresentation.StructuredSecurityGame;
//import sun.org.mozilla.javascript.internal.ObjToIntMap.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: ckiekint
 * Date: Jul 9, 2008
 * Time: 2:32:11 AM
 */
public class FastSecurityGameSolver {

	public long executionTime;


	public void printArrayList(ArrayList<Double> targetVar)
	{

		ListIterator<Double> itrTargetVar= targetVar.listIterator();

		while(itrTargetVar.hasNext())
		{		
			//System.out.println(itrTargetVar.next());		
		}

	}
	public void printHashList(HashSet<Double> targetVar)
	{

		Iterator<Double> itrTargetVar= targetVar.iterator();

		while(itrTargetVar.hasNext())
		{		
			//System.out.println(itrTargetVar.next());		
		}

	}

	public double calculateSvaluesHash(HashSet<Double> targetVar)
	{
		double sPlus=0.0;

		//ListIterator<Double> itrTargetVar= targetVar.listIterator();
		Iterator<Double> itrTargetVar=targetVar.iterator();

		while(itrTargetVar.hasNext())
		{
			//System.out.println("tplus  : " +  itrTargetVar.next());
			//System.out.println("   :  "+ itrTargetVar.next());
			//System.out.println("size  : " + targetVar.size()  );
			double temp = 1/ itrTargetVar.next();

			//System.out.println(itrTargetVar.next());
			sPlus=sPlus+temp;
		}


		return sPlus;



	}
	public double calculateSvaluesList(ArrayList<Double> targetVar)
	{
		double sPlus=0.0;
		ListIterator<Double> itrTargetVar= targetVar.listIterator();
		//Iterator<Double> itrTargetVar=targetVar.iterator();

		while(itrTargetVar.hasNext())
		{		  
			double temp = 1/ itrTargetVar.next();

			sPlus=sPlus+temp;
		}	  	 
		return sPlus;  
	}


	public SSGSolution linearTimeSolveGame(StructuredSecurityGame ssg)
	{
		long startTime = System.currentTimeMillis();
		List<TargetData> targets =  new ArrayList<TargetData>();
		//Double [] targetQarray= new Double [ssg.getNTargets()];
		ArrayList<Double> targetQ=new ArrayList<Double>();


		for (int t = 0; t < ssg.getNTargets(); t++) {
			TargetData td = new TargetData();
			td.attackerUncoveredPayoff = ssg.getPayoffs(t).getAttackerUncoveredPayoff();
			td.currentCoverage = 0;
			td.newCoverage=0;
			td.originalIndex=t;

			targets.add(td);

			//targetQarray[t]=(double)ssg.getPayoffs(t).getAttackerUncoveredPayoff();
			targetQ.add(ssg.getPayoffs(t).getAttackerUncoveredPayoff());
		}
		//<double> targetU = new List<double>();
		//ArrayList<double>= new ArrayList<>

		//List<Double> targetU=new ArrayList<Double>();

		// Double [] targetU		= new Double [ssg.getNTargets()];
		//Double [] targetC		= new Double [ssg.getNTargets()];
		//Double [] targetPlus	= new Double [ssg.getNTargets()];
		//Double [] targetMinus	= new Double [ssg.getNTargets()];

		//ArrayList<Double> targetC=new ArrayList<Double>();
		//ArrayList<Double> targetU=new ArrayList<Double>();

		HashSet<Double>targetC=new HashSet<Double>();
		HashSet<Double>targetU=new HashSet<Double>();
		//ArrayList<Double> targetQ=new ArrayList<Double>();
		//targetQ will be used in loop instead of targetQarray

		ArrayList<Double> targetPlus=new ArrayList<Double>();
		ArrayList<Double> targetMinus=new ArrayList<Double>();


		double tkPlus=Double.POSITIVE_INFINITY;
		double tkMinus=0;

		Median varMedian=new Median();
		boolean optimalFound = false;
		double q=0.0;


		/* targetQ arralist and targetQarray is the same thing.
		 * targetQarray is send as argument to find out median.
		 * targetQ arraylist is used in optimal Q finding loop for its additive facilty.
		 */


		//----------Need iteration--------------------------------------------
		while(optimalFound != true)
			//for(int j=0;j<3;j++)
		{
			System.out.println("in loop");

			//--------------SECTION 1------------------------------------

			//doube a=(Double)(targetQ.size()/2);
			double size=targetQ.size();
			size=size/2;
			int middlePoint =(int) Math.ceil(size);



			ListIterator<Double> itrTargetQ=targetQ.listIterator();
			Double [] targetQarray= new Double [targetQ.size()];
			int indexOfQarray=0;

			//inserting data to targetQarray from targetQ list
			while(itrTargetQ.hasNext())
			{
				targetQarray[indexOfQarray]=itrTargetQ.next();
				indexOfQarray++;
			}

			double median=varMedian.DeterministicSelect(targetQarray,middlePoint);
			//System.out.println("......median.....:   "+ median);

			//---------------SECTION 2-------------------------------
			//-----------Find T+,T-,tk+,tk-   -----------------------

			targetPlus.clear();
			targetMinus.clear();
			tkPlus=Double.POSITIVE_INFINITY;
			tkMinus=Double.NEGATIVE_INFINITY;


			if(targetQarray.length==2)
			{
				for(int i=0;i<targetQarray.length;i++)
				{
					if(targetQarray[i] > median) ////for T? of lenght 2 this condition doesn't has the = sign
					{
						targetPlus.add(targetQarray[i]);
						if(targetQarray[i]<tkPlus) tkPlus=targetQarray[i];
					}
					else
					{	
						targetMinus.add(targetQarray[i]);
						if(targetQarray[i]>tkMinus) tkMinus=targetQarray[i];
					}

				}
			}
			else
			{
				for(int i=0;i<targetQarray.length;i++)
				{
					if(targetQarray[i] >= median)
					{
						targetPlus.add(targetQarray[i]);
						if(targetQarray[i]<tkPlus) tkPlus=targetQarray[i];
					}
					else
					{	
						targetMinus.add(targetQarray[i]);
						if(targetQarray[i]>tkMinus) tkMinus=targetQarray[i];
					}

				}
			}


			/*System.out.println("...........T+.......................");
			    printArrayList(targetPlus);

			    System.out.println(".............T-.......................");
			    printArrayList(targetMinus);

			    System.out.println("tk+  : " + tkPlus + "  tk- : " + tkMinus);
			 */
			//-----------SECTION 3----------------------------
			//------------compute S+,s+,Sc,q-----------------

			double sBigPlus=calculateSvaluesList(targetPlus);
			double sC=calculateSvaluesHash(targetC);
			double sSmallPlus= sBigPlus + sC;
			//			     System.out.println("S+ :  "+ sBigPlus+ " Sc: "+sC+" s+: "+sSmallPlus);

			double k= targetC.size() + targetPlus.size();

			q= (k - ssg.getNumDefenders(0))/sSmallPlus;
			//System.out.println(".....q...... :  "+ q);

			////SECTION-----------4,5,6----------------------

			//--------Newly Added----------------
			//		    	if(targetQ.size()==1)
			//		    	{
			//		    		Iterator<Double> itrTargetQ2= targetQ.iterator();
			//		    		targetPlus.clear();
			//		    		while(itrTargetQ2.hasNext())
			//		    		{
			//		    			targetPlus.add(itrTargetQ2.next());
			//		    		}
			//		    		
			//		    	}
			if(q<tkMinus)
			{
				//ListIterator<Double> itrTargetC= targetC.listIterator();

				//Iterator<Double> itrTargetC=targetC.iterator();

				//we added new iterator variable to T+ & T-
				ListIterator<Double> itrTargetPlus2= targetPlus.listIterator();
				ListIterator<Double> itrTargetMinus2= targetMinus.listIterator();

				//replace Tc with (Tc U T+)
				while(itrTargetPlus2.hasNext())   targetC.add(itrTargetPlus2.next()); 
				//replace T? with T-
				targetQ.clear();
				while(itrTargetMinus2.hasNext())  targetQ.add(itrTargetMinus2.next());
				sC=sSmallPlus;

			}
			if(q>tkPlus)
			{
				//ListIterator<Double> itrTargetU= targetU.listIterator();
				ListIterator<Double> itrTargetPlus2= targetPlus.listIterator();
				ListIterator<Double> itrTargetMinus2= targetMinus.listIterator();

				//replace Tu with (Tu U T-)
				while(itrTargetMinus2.hasNext())   targetU.add(itrTargetMinus2.next());
				//replace T? with T+
				targetQ.clear();
				while(itrTargetPlus2.hasNext())  targetQ.add(itrTargetPlus2.next());
				sBigPlus=sSmallPlus;
				//Sc unchanged

			}
			if(q>=tkMinus && q<=tkPlus)
			{
				optimalFound=true;	    	

			}
			///print
			/*System.out.println("--------Tc --------------");
			     printHashList(targetC);
			     System.out.println("--------Tu --------------");
			     printHashList(targetU);
			     System.out.println("--------T? --------------");
			     printArrayList(targetQ);
			 */
		}

		//---------------------------------------------------------------------


		//	    executionTime = System.currentTimeMillis() - startTime;
		SSGSolution solution = new SSGSolution(ssg.getNTargets(), 1);

		for (int t = 0; t < ssg.getNTargets(); t++) {
			double coveringProbability=Math.max(0,1-q/(targets.get(t).attackerUncoveredPayoff));
			targets.get(t).newCoverage=coveringProbability; 

			//TargetData td = targets.get(t);
			if(coveringProbability<0)
				System.out.println(coveringProbability+" kkkk");
			solution.setProb(targets.get(t).originalIndex, 0, coveringProbability);
		}

		return solution;

	}

	public SSGSolution solveGame(StructuredSecurityGame ssg) 
	{
		long startTime = System.currentTimeMillis();

		List<TargetData> targets =  new ArrayList<TargetData>();

		// intialize target data
		for (int t = 0; t < ssg.getNTargets(); t++) 
		{
			TargetData td = new TargetData();
			td.attackerCoveredPayoff = ssg.getPayoffs(t).getAttackerCoveredPayoff();
			td.attackerUncoveredPayoff = ssg.getPayoffs(t).getAttackerUncoveredPayoff();
			td.currentCoverage = 0;
			td.newCoverage = 0;
			td.originalIndex = t;
			targets.add(td);
		}
		Collections.sort(targets, new TargetComparator());

		if (targets.size() < 2) {
			System.out.println("Trying to solve game with < 2 targets!");
			return null;
		}

		if (ssg.getNDefenseTypes() > 1) {
			System.out.println("Tried to run fast security solver on game with more than 1 resource type!");
			return null;
		}

		//    for (TargetData td : targets) {
		//      System.out.println(td);
		//    }

		boolean optimalFound = false;
		double remainingCoverage = ssg.getNumDefenders(0);
		double startingCoverage;
		int nextTarget = 1;
		while(!optimalFound) 
		{

			/*for (int t = 0; t < ssg.getNTargets(); t++) {
				TargetData td = targets.get(t);
				//if(t==5)
				{
					System.out.println("Prev coverage "+ td.currentCoverage); 
				}
			}*/
			double nextExpectedValue = targets.get(nextTarget).attackerUncoveredPayoff;
			startingCoverage = remainingCoverage;

			//      System.out.println("Starting coverage: " + startingCoverage);
			//      System.out.println("nextExpected: " + nextExpectedValue);

			for (int t = 0; t < nextTarget; t++) 
			{
				double additional = targets.get(t).getAdditionalCoverage(nextExpectedValue);
				//        System.out.println("T" + t + ": " + additional);
				remainingCoverage -= additional;
				
				//second terminal condition
				if (targets.get(t).currentCoverage + additional > 1) 
				{
					//          System.out.println("Coverage probability > 1");
					optimalFound = true;
				}


				// we have run out of coverage probability, and cannot induce the next target in the list
				if (remainingCoverage < 0) 
				{
					//          System.out.println("No coverage remaining; BREAK");
					remainingCoverage = startingCoverage;
					optimalFound = true;
					break;
				}
			}
			//System.out.println("Reamining cov" + ": " + remainingCoverage);

			if (!optimalFound) 
			{

				// update all of the coverage probabilities
				for (int t = 0; t < nextTarget; t++) 
				{
					if(targets.get(t).newCoverage<0)
					{
						//System.out.println("Negative assigning... ");
					}
					targets.get(t).currentCoverage = targets.get(t).newCoverage;
				}

				// increment next target, accounting for targets with equal uncovered payoffs
				nextTarget++;
				while (nextTarget < ssg.getNTargets() && targets.get(nextTarget).attackerUncoveredPayoff == nextExpectedValue) {
					nextTarget++;
				}

				//        System.out.println("New next target: " + nextTarget);

				// no more targets to consider...
				if (nextTarget == ssg.getNTargets()) 
					break;
			}
		}

		//    System.out.println("after loop");

		double normalizationFactor = 0d;
		for (int t = 0; t < nextTarget; t++) {
			normalizationFactor += targets.get(t).getRatio();
		}

		//    System.out.println("Normalization: " + normalizationFactor);

		for (int t = 0; t < nextTarget; t++) {
			TargetData td = targets.get(t);
			double tmp = (td.getRatio() / normalizationFactor) * remainingCoverage;
			//      System.out.println("T" + t + ": " + tmp);
			//System.out.println(tmp + "ggg");
			td.currentCoverage += tmp;
		}

		//    System.out.println("Done");

		// construct and return SSGSolution object
		SSGSolution solution = new SSGSolution(ssg.getNTargets(), 1);

		for (int t = 0; t < ssg.getNTargets(); t++) {
			TargetData td = targets.get(t);
			//if(t==5)
			{
				//System.out.println("Current coverage "+ td.currentCoverage); 
			}
			solution.setProb(td.originalIndex, 0, Math.min(1,td.currentCoverage));
		}

		executionTime = System.currentTimeMillis() - startTime;
		return solution;
	}

	private class TargetData {
		int originalIndex;
		double attackerCoveredPayoff;
		double attackerUncoveredPayoff;
		double currentCoverage;
		double newCoverage;

		public double getAdditionalCoverage(double newExpectedPayoff) 
		{
			double d1 = (newExpectedPayoff - attackerUncoveredPayoff);
			double d2 = (attackerCoveredPayoff - attackerUncoveredPayoff);
			newCoverage = d1/d2 ;
			if(newCoverage<0)
			{
				System.out.println("Negative assigning... ");
			}
			return newCoverage - currentCoverage;
		}

		public double getRatio() {
			return 1 / (attackerUncoveredPayoff - attackerCoveredPayoff);
		}

		public String toString() {
			return "[" + attackerCoveredPayoff + "," + attackerUncoveredPayoff + "," + currentCoverage + "," + newCoverage + "]";
		}
	}

	private class TargetComparator implements Comparator<TargetData> {
		public int compare(TargetData td1, TargetData td2) {
			if( td1.attackerUncoveredPayoff < td2.attackerUncoveredPayoff) return 1;
			else if( td1.attackerUncoveredPayoff > td2.attackerUncoveredPayoff ) return -1;
			else return 0;
		}
	}

	public SSGSolution solveGame(StructuredSecurityGame ssg, double numOrigRes,
			int[] clustersizes) 
	{

		long startTime = System.currentTimeMillis();

		List<TargetData> targets =  new ArrayList<TargetData>();

		// intialize target data
		for (int t = 0; t < ssg.getNTargets(); t++) 
		{
			TargetData td = new TargetData();
			td.attackerCoveredPayoff = ssg.getPayoffs(t).getAttackerCoveredPayoff();
			td.attackerUncoveredPayoff = ssg.getPayoffs(t).getAttackerUncoveredPayoff();
			td.currentCoverage = 0;
			td.newCoverage = 0;
			td.originalIndex = t;
			targets.add(td);
		}
		Collections.sort(targets, new TargetComparator());

		if (targets.size() < 2) {
			System.out.println("Trying to solve game with < 2 targets!");
			return null;
		}

		if (ssg.getNDefenseTypes() > 1) {
			System.out.println("Tried to run fast security solver on game with more than 1 resource type!");
			return null;
		}

		//    for (TargetData td : targets) {
		//      System.out.println(td);
		//    }

		boolean optimalFound = false;
		double remainingCoverage = numOrigRes;//ssg.getNumDefenders(0);
		double normremainingcov = ssg.getNumDefenders(0);
		double normstartcov = 0;
		double startingCoverage;
		int nextTarget = 1;
		//boolean isleft = true;
		//boolean isboundbroke = false;
		double covbound = Double.NEGATIVE_INFINITY;
		while(!optimalFound) 
		{

			/*for (int t = 0; t < ssg.getNTargets(); t++) {
				TargetData td = targets.get(t);
				//if(t==5)
				{
					System.out.println("Prev coverage "+ td.currentCoverage); 
				}
			}*/
			double nextExpectedValue = targets.get(nextTarget).attackerUncoveredPayoff;
			startingCoverage = remainingCoverage;
			normstartcov = normremainingcov;

			//      System.out.println("Starting coverage: " + startingCoverage);
			//      System.out.println("nextExpected: " + nextExpectedValue);

			for (int t = 0; t < nextTarget; t++) 
			{
				int origindex = targets.get(t).originalIndex;
				double additional = targets.get(t).getAdditionalCoverage(nextExpectedValue);
				//        System.out.println("T" + t + ": " + additional);
				remainingCoverage = remainingCoverage - (additional*clustersizes[origindex]);
				normremainingcov = normremainingcov - additional;
				
				//second terminal condition
				
				if ( ((targets.get(t).currentCoverage+additional)*clustersizes[origindex]) >= clustersizes[origindex]) 
				{
					//          System.out.println("Coverage probability > 1");
					optimalFound = true;
					covbound = Math.max(covbound, targets.get(t).attackerCoveredPayoff);
				}


				// we have run out of coverage probability, and cannot induce the next target in the list
				if ((remainingCoverage <= 0)) 
				{
					//          System.out.println("No coverage remaining; BREAK");
					remainingCoverage = startingCoverage;
					normremainingcov = normstartcov;
					optimalFound = true;
					
					break;
				}
			}
			//System.out.println("Reamining cov" + ": " + remainingCoverage);

			if (!optimalFound) 
			{

				// update all of the coverage probabilities
				for (int t = 0; t < nextTarget; t++) 
				{
					if(targets.get(t).newCoverage<0)
					{
						//System.out.println("Negative assigning... ");
					}
					targets.get(t).currentCoverage = targets.get(t).newCoverage;
				}

				// increment next target, accounting for targets with equal uncovered payoffs
				nextTarget++;
				while (nextTarget < ssg.getNTargets() && targets.get(nextTarget).attackerUncoveredPayoff == nextExpectedValue) {
					nextTarget++;
				}

				//        System.out.println("New next target: " + nextTarget);

				// no more targets to consider...
				if (nextTarget == ssg.getNTargets()) 
					break;
			}
		}
		
		
		/*double sum = 0; 
		for (int t = 0; t < ssg.getNTargets(); t++) 
		{
			TargetData td = targets.get(t);
			System.out.println("Traget "+t + " current cov " + td.currentCoverage);
			sum += td.currentCoverage;
			
		}

		System.out.println("sum "+ sum);*/

		double normalizationFactor = 0d;
		for (int t = 0; t < nextTarget; t++)
		{
			normalizationFactor += (targets.get(t).getRatio());
		}

		//    System.out.println("Normalization: " + normalizationFactor);

		for (int t = 0; t < nextTarget; t++) 
		{
			TargetData td = targets.get(t);
			double tmp = (td.getRatio() / normalizationFactor) * (remainingCoverage/clustersizes[td.originalIndex]);
			//      System.out.println("T" + t + ": " + tmp);
			//System.out.println(tmp + "ggg");
			td.currentCoverage += tmp;
		}
		
		/*for (int t = 0; t < nextTarget; t++) 
		{
			TargetData td = targets.get(t);
			covbound = Double.NEGATIVE_INFINITY;
			//System.out.println("Traget "+t + " current cov " + td.currentCoverage);
			if(td.currentCoverage >= 1)
			{
				covbound = Math.max(covbound, td.attackerCoveredPayoff);
			}
			if(covbound> Double.NEGATIVE_INFINITY)
			{
					//TargetData td = targets.get(t);
					double d1 = covbound - td.attackerUncoveredPayoff;
					double d2 = td.attackerCoveredPayoff-td.attackerUncoveredPayoff;
					td.currentCoverage = (d1)/(d2);
			}
		}*/
		
		
		
		
		

		//    System.out.println("Done");

		// construct and return SSGSolution object
		SSGSolution solution = new SSGSolution(ssg.getNTargets(), 1);

		for (int t = 0; t < ssg.getNTargets(); t++) {
			TargetData td = targets.get(t);
			//if(t==5)
			{
				//System.out.println("Current coverage "+ td.currentCoverage); 
			}
			//System.out.println("Traget "+t + " current cov " + td.currentCoverage);
			solution.setProb(td.originalIndex, 0, Math.min(1,td.currentCoverage));
			//solution.setProb(td.originalIndex, 0, td.currentCoverage);
		}

		executionTime = System.currentTimeMillis() - startTime;
		return solution;
		
		
		
		
	}
}


