package cs.Interval.Clustering;

import java.util.ArrayList;
import java.util.List;





/**
 * @author anjonsunny
 *
 */
public class StrategyMapping {

	private int numberofplayers;
	private int[] numberofactions;
	private int[] numberofclusters; // this can be an array to incorporate different number of cluster for different players
	private List<Integer>[][] clusterforplayerswithstrategy;// = new List[numberofplayers][numberofclusters]; //yeah, right
	private List<Integer>[] clusterfortargets; // for wildlife 
	public List<Integer>[] getClusterfortargets() {
		return clusterfortargets;
	}




	




	private int[][] securitygamedata;
	//private Game originalgame;
	private String gamename;
	//private Map<int[], MixedStrategy[]> subgameSols = new HashMap<int[], MixedStrategy[]>(); // strategy profiles for subgames given by the clustering



	



	/**
	 * constructor for wildlife security games
	 * @param clusterfortargets
	 * @param numberofcluster
	 * @param securitygamedata
	 */
	public StrategyMapping(List<Integer>[] clusterfortargets, int numberofcluster, 
			int[][] securitygamedata) 
	{
		super();
		this.clusterfortargets = new List[numberofcluster];
		for(int i=0; i< this.clusterfortargets.length; i++)
		{
			this.clusterfortargets[i] = new ArrayList<Integer>();
		}
		for(int i=0; i< numberofcluster; i++)
		{
			for(Integer x: clusterfortargets[i])
			{
				this.clusterfortargets[i].add(x);
			}
		}
		this.securitygamedata = new int[securitygamedata.length][4];
		for(int i=0; i<this.securitygamedata.length; i++)
		{
			for(int j=0; j<4; j++)
				this.securitygamedata[i][j] = securitygamedata[i][j];
		}
	}
	
	public void printSecurityGameMapping()
	{
		for(int i=0; i< this.clusterfortargets.length; i++)
		{
			System.out.print("Cluster " + i + " : ");
			for(Integer target: this.clusterfortargets[i])
			{
				System.out.print(target);
				if(this.clusterfortargets[i].indexOf(target) < (this.clusterfortargets[i].size()-1) )
				{
					System.out.print(",");
				}
			}
			System.out.print("\n");
		}
	}




	
	/**
	 * 
	 * @param abstractedaction abstracted action or pass the cluster number
	 * @param player
	 * @return a list of actions belong to the abstracted action. 
	 */
	public List<Double> getOriginalActionsFromAbstractedAction(int abstractedaction, int player)
	{
		List<Double> elm = new ArrayList<Double>();
		for(Integer x: this.clusterforplayerswithstrategy[player][abstractedaction-1])
		{
			elm.add((Math.floor(x)));
		}

		return elm;
	}



	


	/*
	 * calculates the max payoff for the actions in the original game in a cluster
	 */

	




		

	


	/*public double[][] getActionsMaxDeviation(MatrixGame absgame, int[] outcome)
	{

	}*/




	


	


	public  int getClusterSize(int action, int player)
	{

		int cluster = this.getClusternumber(action, player);
		int size = this.clusterforplayerswithstrategy[player][cluster].size();
		return size; 
	}

	public int getClusterSize1(int cluster, int player)
	{
		return this.clusterforplayerswithstrategy[player][cluster].size();

	}



	/**
	 * 
	 * @param clustermappingtoaction just one dimensional array containing the cluster number for each action
	 * @param player player number
	 */
	public void mapActions(int[] clustermappingtoaction, int player)
	{




		for(int i=0; i< this.numberofclusters[player]; i++)
		{
			for(int j=0; j < clustermappingtoaction.length; j++)
			{
				if( (i+1) == clustermappingtoaction[j]) //the clustermappingtoaction should contain the cluster number starts from 1
				{
					double[] actionwithprobability = new double[2]; // [0]<- action, [1]<- probability 
					actionwithprobability[0] = j+1;
					int action = j+1;
					this.clusterforplayerswithstrategy[player][i].add(action);
				}

			}

		}

		for(List<Integer> x: this.clusterforplayerswithstrategy[player])
		{
			for(Integer y: x)
			{
				System.out.print(y);
			}
			System.out.print("\n");
		}

	}



	/**
	 * 
	 * @param probability set the probability for an action
	 * @param player
	 */
	/*public void setStrategy(int[] probability, int player ) // there will be  n number of probability where n = number of actions in the abstracted game. 
	{
		for(int i =0; i< this.numberofclusters[player]; i++)
		{
			for(double[] x: this.clusterforplayerswithstrategy[player][i])
			{
				x[1] = probability[i];
			}
		}

	}
	 */

	


	/*
	 *//**
	 * Build the abstracted game from original game
	 *//*
	public String buildAbstractedGame()
	{





		for(int i =0; i< this.numberofplayers; i++)
		{
			System.out.print("\n\n For player: "+ i + "\n");
			for(int j =0; j< this.numberofclusters[i]; j++)
			{
				System.out.print("\nCLuster "+ j +": ");
				for(double[] x: this.clusterforplayerswithstrategy[i][j])
				{
					System.out.print(x[0] + " ");
				}
			}
			System.out.print( "\n");
		}





		int[] N = new int[this.numberofplayers];

		for(int i=0; i<N.length; i++)
		{
			N[i] = this.numberofclusters[i];
		}


		MatrixGame abstractedgame  = new MatrixGame(originalgame.getNumPlayers(), N);



		For all players, do the following steps:
	  * 
	  *1. take the abstracted game.
	  *2. take an action X
	  *    3. take the original game
	  *    4. iterate over all the actions and check whether it belongs to action X
	  *    5. SUm the payoffs and take average
	  *    6. Set the payoff for player i.
	  * 7. go to 2. 
	  *    
	  *    
	  * 




		OutcomeIterator iteratorabstractgame = abstractedgame.iterator();
		int[] outcomeabstractgame = new int[this.numberofplayers];
		int[] outcomeoriginalgame = new int[this.numberofplayers]; 

		double payoffsum =0;
		int counter = 0;

		for(int i=0; i< this.numberofplayers; i++)
		{
			iteratorabstractgame = abstractedgame.iterator();

			while(iteratorabstractgame.hasNext())
			{
				outcomeabstractgame = iteratorabstractgame.next();
				OutcomeIterator iteratororiginalgame = this.originalgame.iterator();
				payoffsum = 0;
				counter = 0;

				//now do step 6
				while(iteratororiginalgame.hasNext())
				{ 
					outcomeoriginalgame = iteratororiginalgame.next();
					//  now check
					if(checkIfMappingMatches(outcomeabstractgame, outcomeoriginalgame) == true)
					{
						payoffsum+= this.originalgame.getPayoff(outcomeoriginalgame, i);
						counter++;
					}


				}

				if(counter!=0)
				{

					abstractedgame.setPayoff(outcomeabstractgame, i, payoffsum/counter);
				}
				else if(counter==0)
				{
					abstractedgame.setPayoff(outcomeabstractgame, i, 0);
				}


			}// end of while loop
		} // end of outer for loop

		String gamename = Parameters.GAME_FILES_PATH+Referee.experimentdir+"/"+"k"+this.numberofclusters[0]+"-"+this.gamename+Parameters.GAMUT_GAME_EXTENSION;
		String gmname = "k"+this.numberofclusters[0]+"-"+this.gamename;

		try{

			PrintWriter pw = new PrintWriter(gamename,"UTF-8");
			SimpleOutput.writeGame(pw,abstractedgame);
			pw.close();
		}
		catch(Exception ex){
			System.out.println("StrategyMapping class :something went terribly wrong during clustering abstraction ");
		}


		try{
			PrintWriter pw = new PrintWriter(Parameters.GAME_FILES_PATH+Referee.experimentdir+"/"+"k"+this.numberofclusters[0]+"-"+this.gamename+".mapping","UTF-8");

			pw.write("Number of players: " + this.numberofplayers + "\n");
			pw.write("Number of actions: "+"\n");
			for(int i =0; i<this.numberofplayers; i++)
			{
				pw.write("Player "+ (i+1) +" : "+this.numberofactions[i]+ "\n");
			}


			pw.write("\n \n");

			pw.write("Number of clusters: "+"\n");
			for(int i =0; i<this.numberofplayers; i++)
			{
				pw.write("Player "+ (i+1) +" : "+this.numberofclusters[i]+ "\n");
			}


			pw.write("\n \n");
			pw.write("Mapping of actions: "+"\n\n");


			for(int i =0; i< this.numberofplayers; i++)
			{
				pw.write("Player "+ (i+1)+" : \n");
				int k =1 ;
				for(List<double[]> x: clusterforplayerswithstrategy[i])
				{
					pw.write("Cluster "+ (k++) + " Actions: ");
					for(double[] y: x)
					{
						pw.write(y[0]+", ");
					}
					pw.write("\n");
				}

				pw.write("\n \n");
			}





			//SimpleOutput.writeGame(pw,abstractedgame);
			pw.close();
		}
		catch(Exception ex){
			System.out.println("StrategyMapping class :something went terribly wrong during clustering abstraction ");
		}





		return Referee.experimentdir+"/"+gmname;




	}*/


	/**
	 * 
	 * @param abstractactions array of actions in abstracted game
	 * @param originalactions array of actions in original game
	 * @return true or false if the original actions belong to those abstracted actions
	 */
	public boolean checkIfMappingMatches(int[] abstractactions, int[] originalactions)
	{

		for(int i=0; i< this.numberofplayers; i++)
		{
			if(abstractactions[i] != (getClusternumber(originalactions[i], i) + 1)) // plus 1 because cluster number starts from 0 but action starts from 1.
			{
				return false;
			}

		}

		return true;


	}

	/**
	 * 
	 * @param action pass an action of original game
	 * @param player player number
	 * @return the cluster number the action belongs to
	 */
	public int getClusternumber(int action, int player)
	{
		int clusternumber = -1;
		for(int i=0; i<this.numberofclusters[player]; i++)
		{
			for(Integer x: this.clusterforplayerswithstrategy[player][i])
			{
				if(action == x)
				{
					clusternumber = i;
					break;
				}
			}
		}
		return clusternumber;

	}

	


	





	public void setClusterfortargets(List<Integer>[] clusterfortargets) {
		this.clusterfortargets = clusterfortargets;
	}

	public int[][] getSecuritygamedata() {
		return securitygamedata;
	}




	public void setSecuritygamedata(int[][] securitygamedata) {
		this.securitygamedata = securitygamedata;
	}




	







}
