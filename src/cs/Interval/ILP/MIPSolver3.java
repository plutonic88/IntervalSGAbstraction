package cs.Interval.ILP;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Queue;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

public class MIPSolver3 {



	/**
	 * remove dominated targets from targets array list
	 * @param duplicatetargets
	 * @param targets targets in the original graph. 
	 * @param edgewithcoins
	 * @param nodewithcoins 
	 * @param dominatedtargets
	 * @param distancelimit
	 * @param nTargets 
	 * @return
	 * @throws IloException
	 */
	public static ArrayList<Integer> solve(ArrayList<TargetNode> duplicatetargets, 
			ArrayList<TargetNode> targets, HashMap<Integer, Integer> nodewithcoins,
			ArrayList<TargetNode> dominatedtargets, double distancelimit, int nTargets) throws Exception
	{
		//edgewithcoins = 0;

		int n = duplicatetargets.size();
		int q = nodewithcoins.size()+dominatedtargets.size(); 
		double v[] = new double[q]; 
		double [][] e = new double[n][n];
		double[][] d = new double[n][n]; 
		int A[][] = new int[n][q];
		int B[][][] = new int[n][n][q];
		int startnode = duplicatetargets.size()-2;

		System.out.println("n="+duplicatetargets.size()+", q="+q);


		/*int icount =0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();

		 *//**
		 * this mapping is for v
		 *//*

		for(int i=0; i<targets.size(); i++)
		{

			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		  */


		/**
		 * fill v
		 *
		 */
		int kk=0;
		for(int i=0; i<q; i++)
		{
			/* if(i<nodewithcoins.size())
			 {
				 v[i] = SecurityGameContraction.getTargetNode(nodewithcoins.get(i), SecurityGameContraction.origtargets).getAnimaldensity();
			 }
			 else
			 {
				v[dominatedtargets.get(kk).getTargetid()] = dominatedtargets.get(kk).getAnimaldensity();
				kk++;
			 }*/
			v[i] = SecurityGameContraction.origtargets.get(i).getAnimaldensity();
		}



		/*for(Integer x : nodewithcoins.keySet())
		{
			//v[x] = targets.get(nodewithcoins.get(x)).getAnimaldensity() ;
			v[x] = SecurityGameContraction.getNodeWithCoin(nodewithcoins.get(x)).getAnimaldensity();
		}*/
		/**
		 * hashmap for mapping from edge to k.
		 * This hashmap starts from targets.size()
		 */

		/*HashMap<Integer, Integer[]> ktoedge = SecurityGameContraction.fillKToEdge(nodewithcoins.size()); // 

		for(Integer k: ktoedge.keySet())
		{
			v[k] = ktoedge.get(k)[2]; // 2nd index has the path utility
		}
		 */

		/**
		 * fill d and e
		 */
		for(int i=0; i<duplicatetargets.size(); i++)
		{
			for(int j=0; j<duplicatetargets.size(); j++)
			{
				
				
				/*if(i==2 && j==19)
				{
					System.out.println(duplicatetargets.get(i).getTargetid()+" ??? "+duplicatetargets.get(j).getTargetid());
					for(TargetNode a: duplicatetargets.get(i).getNeighbors())
					{
						System.out.print(a.getTargetid()+" ");
					}
					System.out.println();
					
				}*/

				if(duplicatetargets.get(i).getNeighbors().contains(duplicatetargets.get(j)))  // if i-> j
				{


					double pathutility = duplicatetargets.get(i).getPathUtility(duplicatetargets.get(j));
					//System.out.println(i +" --> "+ j+ " : "+ duplicatetargets.get(i).getDistance(duplicatetargets.get(j)) 
							///*+ ", pathsize "+duplicatetargets.get(i).getPath(duplicatetargets.get(j)).size()*/ );
					/*if(pathutility==0)
					{ 
					 *//**
					 * means no intermediate targets
					 *//*
						//pathutility += targets.get(i).getAnimaldensity();

					}*/

					d[i][j] = duplicatetargets.get(i).getDistance(duplicatetargets.get(j));


					e[i][j] = 1;

					//System.out.println(i +" --> "+ j+ " : "+ e[i][j] );


				}
				else if(!duplicatetargets.get(i).getNeighbors().contains(duplicatetargets.get(j)))
				{


					e[i][j] = 0;
					d[i][j] = Double.MAX_VALUE;
					//System.out.println(i +" --> "+ j+ " : "+ e[i][j] );

				}
			}
		}




		/**
		 * fil Aik and Bijk
		 */

		for(int i=0; i<duplicatetargets.size(); i++)
		{
			//for(int j=0; j<duplicatetargets.size(); j++)
			{
				//int cnt=0;
				for(int k=0; k<q; k++)
				{
					//if(k<nodewithcoins.size())
					{
						// for A...node coins.. fill B with 0
						if(duplicatetargets.get(i).getTargetid()==
								SecurityGameContraction.origtargets.get(k).getTargetid() )
						{
							A[i][k] = 0;
						}
						else if(i==startnode)
						{
							A[i][k] = 0;
						}
						//	B[i][j][k] = 0;
					}




				}
			}
		}




		for(int i=0; i<duplicatetargets.size(); i++)
		{
			for(int j=0; j<duplicatetargets.size(); j++)
			{
				if(i!=j)
				{
					for(int k=0; k<q; k++)
					{


						// for B...edge coins... fill A with 0
						TargetNode nodei= duplicatetargets.get(i);
						TargetNode nodej= duplicatetargets.get(j);

						/**
						 * if k<targets.size  k is a node
						 * else k is an edge
						 */

						//if(k<(nodewithcoins.size()))
						{
							//
							TargetNode kthnode = SecurityGameContraction.origtargets.get(k);//SecurityGameContraction.getNodeWithCoin(nodewithcoins.get(k));//targets.get(k);



							if(e[i][j]==1 && isInDominatedTarget(kthnode, dominatedtargets))
							{
								/**
								 * check if kth node is a part of path from i to j
								 */

								//System.out.println(i+"..-->"+j);
								boolean isinpath = isInPath(kthnode, nodei, nodej, targets);
								if(isinpath)
								{
									B[i][j][k] = 1;
								}
								else
								{
									B[i][j][k] = 0;
								}

							}
							else
							{
								if((nodei.getTargetid()==kthnode.getTargetid()) && 
										(nodej.getTargetid()==kthnode.getTargetid()) && 
										(e[i][j]==1) )
								{
									B[i][j][k] = 1;
									//B[j][i][k] = 1;
								}
								else
								{
									B[i][j][k] = 0;
									//B[j][i][k] = 0;
								}
							}


						}
						/*else if(k>=(nodewithcoins.size()))
						{
							Integer[] x = ktoedge.get(k);
							if( (((nodei.getTargetid()==x[0]) && (nodej.getTargetid()==x[1]))
									|| ((nodei.getTargetid()==x[1]) && (nodej.getTargetid()==x[0])))
									&& nodei.getNeighbors().contains(nodej))
							{
								B[i][j][k] = 1;
								//B[j][i][k] = 1;
							}
							else
							{
								B[i][j][k] = 0;
								//B[j][i][k] = 0;
							}

						}*/

					}
				}
			}
		}


		//printABED(A,B,e,d);


		//B[3][2][0] = 0;
		//d[3][2] = Double.MAX_VALUE;


		try
		{
			IloCplex cplex = new IloCplex();
			IloNumVar[][] x = new IloNumVar[n][];
			for(int i=0; i<n; i++)
			{
				x[i] = cplex.boolVarArray(n);
			}

			IloNumVar[] y = new IloNumVar[q];
			for(int i=0; i<q; i++)
			{
				y[i] = cplex.intVar(0, 1);

			}
			IloNumVar[] u = cplex.numVarArray(n, 1, n);

			for(int i=0; i<n ; i++)
			{
				if(i==startnode)
				{
					u[i] = cplex.numVar(1, 1);
				}
				else
				{
					u[i] = cplex.numVar(2, n);
				}
			}

			//objective
			IloLinearNumExpr obj = cplex.linearNumExpr();

			for(int k=0; k<q; k++)
			{
				obj.addTerm(v[k], y[k]);
			}
			cplex.addMaximize(obj);



			/*//constraint 17
			for(int i=0; i<n; i++)
			{
				//IloLinearNumExpr expr1 = cplex.linearNumExpr();
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						cplex.addLe(cplex.prod(1, x[i][j]), e[i][j]);


					}

				}

			}
			 */



			//constraint 19
			for(int k=0; k<q; k++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();
				for(int i=0; i<n; i++)
				{
					for(int j=0; j<n; j++)
					{
						expr1.addTerm(-(A[i][k]+B[i][j][k]), x[i][j]);
						//expr1.addTerm(-B[i][j][k], x[i][j]);
					}
				}
				expr1.addTerm(1, y[k]);
				cplex.addLe(expr1, 0);

			}

			//constraint 20
			IloLinearNumExpr exprgoal = cplex.linearNumExpr();
			for(int i=0; i<n; i++)
			{
				if(i!=startnode)
					exprgoal.addTerm(1, x[startnode][i]);
			}
			cplex.addEq(exprgoal, 1);


			IloLinearNumExpr exprgoal1 = cplex.linearNumExpr();
			for(int i=0; i<n; i++)
			{
				if(i!=startnode)
					exprgoal1.addTerm(1, x[i][startnode]);
			}
			cplex.addEq(exprgoal1, 1);

			//cplex.addEq(cplex.prod(1, y[startnode]), 1);



			//constraint 21
			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						expr1.addTerm(1.0, x[i][j]);
						expr1.addTerm(-1, x[j][i]);
					}
				}
				cplex.addLe(expr1,1);
			}




			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						expr1.addTerm(1.0, x[i][j]);

					}
				}
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						expr2.addTerm(-1, x[j][i]);
					}
				}
				cplex.addEq(cplex.sum(expr1,cplex.prod(1, expr2)), 0);

			}






			//constraint 22
			IloLinearNumExpr expr3 = cplex.linearNumExpr();


			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{

						expr3.addTerm(d[i][j], x[i][j]);

					}
				}
			}
			cplex.addLe(expr3, distancelimit);


			//constraint 25
			for(int i=0; i<n; i++)
			{

				for(int j=0; j<n; j++)
				{
					if((i!=j) &&(i!=startnode) && (j!=startnode))
					{	
						IloLinearNumExpr expr = cplex.linearNumExpr();
						expr.addTerm(1.0, u[i]);
						expr.addTerm(-1, u[j]);
						expr.addTerm(n-1, x[i][j]);
						cplex.addLe(expr, n-2);


					}

				}

			}

			cplex.solve();

			System.out.println("obj: "+ cplex.getObjValue());
			ArrayList<Integer> pathtoconsider = new ArrayList<Integer>();
			
			int cs = nTargets;
			int matadj [][] = new int[cs+1][cs+1];
			double distt [][] = new double[cs+1][cs+1];
			
			
			

			int[] f= new int[cs+1];
			
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						//System.out.print("x["+i+"]["+j+"]="+cplex.getValue(x[i][j])+ " #### ");
						if(cplex.getValue(x[i][j])==1)
						{
							int inode = duplicatetargets.get(i).getTargetid();
							int jnode = duplicatetargets.get(j).getTargetid();
							
							TargetNode in = SecurityGameContraction.getTargetNode(inode, targets);
							TargetNode jn = SecurityGameContraction.getTargetNode(jnode, targets);
							
							if(inode != jnode && jnode!=nTargets && inode!=nTargets)
							{
								//System.out.println(inode+"----"+jnode);
								f[inode]++;
								matadj[inode][jnode] = 1;
								distt[inode][jnode] = in.getDistance(jn);
							}
						}
					}
				}
			}
			
			
			  
			ArrayList<ArrayList<Integer>> pathseq = extractPath(matadj, 0, distt, distancelimit);
			
		
			// find the largest path
			
			ArrayList<Integer> lpath = findLPath(pathseq);
			
			
			
			
			///////////
			
			/*
			int current = 0; // basenode
			pathtoconsider.add(current);
			*/
			
			/*int flag = 0;
			
			while(flag!=1)
			{
				// find edge that starts with current and end with a different node y. 
				
				 * if current is not base and y is base then end the loop
				 
				
				
				//TODO build frequency table, make sure that yid maintain the frequency
				
				int yid = findPath(cplex, x, current, duplicatetargets, nTargets);
				
				if(current != 0 && yid == 0 && current !=nTargets)
				{
					flag = 1;
				}
				
				current = yid;
				if(!pathtoconsider.contains(yid) && yid != nTargets)
				{
					pathtoconsider.add(yid);
				}
				
				
				
			}
			pathtoconsider.add(0);
			
			*/
			
			//////////

			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						//System.out.print("x["+i+"]["+j+"]="+cplex.getValue(x[i][j])+ " #### ");
						if(cplex.getValue(x[i][j])==1)
						{
							TargetNode inode = duplicatetargets.get(i);
							TargetNode jnode = duplicatetargets.get(j);
							System.out.println( inode.getTargetid()+ ", in:"+inode.isIn()+",out:"+inode.isOut()+
									"-->"+jnode.getTargetid()+ ", in:"+jnode.isIn()+",out:"+jnode.isOut());


							//if( !pathtoconsider.contains(inode.getTargetid()) && inode.getTargetid()!=nTargets)
							//	pathtoconsider.add(inode.getTargetid());
							//if(!pathtoconsider.contains(jnode.getTargetid()) && jnode.getTargetid()!=nTargets)
								//pathtoconsider.add(jnode.getTargetid());
							
							
							
							//if(inode.getTargetid()==jnode.getTargetid() && !pathtoconsider.contains(inode.getTargetid()))
							{
								//pathtoconsider.add(inode.getTargetid());
								
								//  *try to find an edge where start node id is jnode.id and end node a different than jnode.id
								 
							//	int tid = findPath(cplex, x, jnode.getTargetid(), duplicatetargets);
								//pathtoconsider.add(tid);
							}
							


						}
					}
				}
				//System.out.println();
			}
			
			System.out.println();
			//System.out.println();

			/*for(int i=0; i<q; i++)
			{



				//Integer p[] = ktoedge.get(i);
				System.out.print("y["+i+"]="+cplex.getValue(y[i])+ " ");


				//System.out.println();
			}
			System.out.println();

			for(int i=0; i<v.length; i++)
			{



				System.out.print("v["+i+"]="+(v[i])+ " ");


				//System.out.println();
			}*/

			/*for(int i=0; i<n; i++)
			{

				if(i!=startnode)
				{

					System.out.print("u["+i+"]="+cplex.getValue(u[i])+ " ");
				}

				//System.out.println();
			}
			System.out.println();*/
			double objv = cplex.getObjValue(); 
			//System.out.println("coin collected  "+ objv);
			cplex.end();
			return lpath;






		}
		catch (Exception ex)
		{

		}

		return null;


	}
	
	
	private static ArrayList<Integer> findLPath(ArrayList<ArrayList<Integer>> pathseq) {
		
		int maxl= Integer.MIN_VALUE;
		ArrayList<Integer> lp = new ArrayList<Integer>();
		
		for(ArrayList<Integer> p: pathseq)
		{
			if(p.size()>maxl)
			{
				lp = p;
			}
		}
		
		return lp;
	}


	private static void printPath(TargetNode node, TargetNode goal, ArrayList<Integer> p) {

		if(node == null)
			return;
		printPath(node.parent, goal, p);
		System.out.print(  node.getTargetid() +"->");
		p.add(node.getTargetid());

	}
	
	
	public static ArrayList<ArrayList<Integer>> extractPath(int[][] mat, int src, double[][] d, double dmax)
	{
		
		
		 
		TargetNode start = new TargetNode();
		start.setTargetid(src);
		Queue<TargetNode> fringequeue = new LinkedList<TargetNode>();
		ArrayList<TargetNode> goals = new ArrayList<TargetNode>();
		//System.out.println("start node "+ start.getTargetid());
		ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
		
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
				goals.add(node);
				ArrayList<Integer> p = new ArrayList<Integer>();
				printPath(node, node, p); // dont comment this out
				paths.add(p);
				pathcounter++;
				//System.out.println();
				

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
				ArrayList<TargetNode> succs = Expand(node, dmax, mat, d);
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
		
		return paths;
		
	}
	
	
	private static ArrayList<TargetNode> Expand(  TargetNode node, double dmax, int[][] matadj, double[][] d) 
	{
		ArrayList<TargetNode> successors = new ArrayList<TargetNode>();

		/**
		 * find the index for  node
		 */
		int nodeid = node.getTargetid();
		//TargetNode tmpnode = new TargetNode(); 
		
		for(int j=0; j<matadj[nodeid].length; j++)
		{
			
			if(matadj[nodeid][j]==1)
			{
				TargetNode newnei = new TargetNode();
				newnei.setTargetid(j);
				newnei.distancecoveredyet = node.distancecoveredyet + d[nodeid][j];//tmpnode.getDistance(nei);
				newnei.parent = node;
				if(newnei.distancecoveredyet<=dmax)
				{
					successors.add(newnei);
				}
			}
			
		}
		return successors;

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
			while(tmpstart!=null)
			{
				//Pmat[map.get(tmpstart.getTargetid())][pathindex]  = 1;
				tmppathseq.add(tmpstart.getTargetid());
				tmpstart = tmpstart.parent;
				if(tmpstart.getTargetid()==tmpgoal.getTargetid())
				{
					tmppathseq.add(tmpgoal.getTargetid());
					break;
				}

			}
			pathindex++;
			pathseq.add(tmppathseq);
		}

		//return Pmat;





	}

	

	private static int findPath(IloCplex cplex, IloNumVar[][] x, int targetid, 
			ArrayList<TargetNode> duplicatetargets, int nTargets) throws UnknownObjectException, IloException {
		
		int n= duplicatetargets.size();
		int resid=-1;
		
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				if(i!=j)
				{
					TargetNode in= duplicatetargets.get(i);
					TargetNode jn= duplicatetargets.get(j);
					
					
					if(cplex.getValue(x[i][j])==1 && in.getTargetid()==targetid && jn.getTargetid()!=targetid && jn.getTargetid()!=nTargets)
					{
						
						resid = jn.getTargetid();
						return resid;
					}
				}
			}
		}
		
		
		return resid;
	}

	private static boolean isInPath(TargetNode kthnode, TargetNode nodei,
			TargetNode nodej, ArrayList<TargetNode> targets) {


		ArrayList<TargetNode> nodes = targets;

		if(nodei.getTargetid()==nodej.getTargetid())
			return false;

		/**
		 * get the path first
		 */
		ArrayList<TargetNode> pathnodes = new ArrayList<TargetNode>();

		TargetNode i = getTarget(nodes, nodei);
		TargetNode j = getTarget(nodes, nodej);
		/*System.out.println("Neighbors of "+ i.getTargetid());
		
		for(TargetNode m: i.getNeighbors())
		{
			System.out.print(m.getTargetid()+" ");
		}
		System.out.println();
		

		System.out.println("Pathnodes between "+ i.getTargetid() + ", "+ j.getTargetid() + ", dist : "+ i.getDistance(j));

		*/

		for(TargetNode x: pathnodes)
		{
			if(x.getTargetid()==kthnode.getTargetid())
			{
				return true;
			}
		}




		return false;
	}

	private static TargetNode getTarget(ArrayList<TargetNode> nodes,
			TargetNode nodei) {

		for(TargetNode x: nodes)
		{
			if(x.getTargetid()==nodei.getTargetid())
			{
				return x;
			}
		}


		return null;
	}

	private static boolean isInDominatedTarget(TargetNode kthnode,
			ArrayList<TargetNode> dominatedtargets) {


		for(TargetNode x: dominatedtargets)
		{
			if(x.getTargetid()==kthnode.getTargetid())
			{
				return true;
			}
		}

		return false;
	}

	private static void printABED(int[][] a, int[][][] b, double[][] e, double[][] d) {

		for(int i=0; i<a.length; i++)
		{
			for(int j=0; j<a[i].length; j++)
			{
				System.out.print("a["+i+"]["+j+"]="+a[i][j]+ " ");
			}
			System.out.println();

		}

		for(int i=0; i<b.length; i++)
		{
			for(int j=0; j<b[i].length; j++)
			{
				for(int k=0; k<b[i][j].length; k++)
				{

					System.out.print("b["+i+"]["+j+"]["+k+"]="+b[i][j][k]+ " ");
				}
				System.out.println();
			}
			System.out.println("\n");

		}

		for(int i=0; i<e.length; i++)
		{
			for(int j=0; j<e[i].length; j++)
			{
				System.out.print("e["+i+"]["+j+"]="+e[i][j]+ " ");
			}
			System.out.println();

		}

		DecimalFormat df = new DecimalFormat("#.00");

		for(int i=0; i<d.length; i++)
		{
			for(int j=0; j<d[i].length; j++)
			{
				if(d[i][j]==1 || d[i][j]==0 || d[i][j]==2 )
				{
					System.out.print("d["+i+"]["+j+"]="+d[i][j]+ " ");
				}
				else
				{
					System.out.print("d["+i+"]["+j+"]="+"infinity"+ " ");
				}
			}
			System.out.println();

		}





	}

}
