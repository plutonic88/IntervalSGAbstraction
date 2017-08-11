package cs.Interval.ILP;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;

public class MIPSolver2 {



	/**
	 * remove dominated targets from targets array list
	 * @param duplicatetargets
	 * @param targets targets in the original graph. 
	 * @param edgewithcoins
	 * @param nodewithcoins 
	 * @param dominatedtargets
	 * @param distancelimit
	 * @return
	 * @throws IloException
	 */
	public static double solve(ArrayList<TargetNode> duplicatetargets, ArrayList<TargetNode> targets, int edgewithcoins, HashMap<Integer, Integer> nodewithcoins, ArrayList<TargetNode> dominatedtargets, double distancelimit) throws IloException
	{
		//edgewithcoins = 0;

		int n = duplicatetargets.size();
		int q = nodewithcoins.size()+edgewithcoins; 
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
		 */

		for(Integer x : nodewithcoins.keySet())
		{
			//v[x] = targets.get(nodewithcoins.get(x)).getAnimaldensity() ;
			v[x] = SecurityGameContraction.getNodeWithCoin(nodewithcoins.get(x)).getAnimaldensity();
		}
		/**
		 * hashmap for mapping from edge to k.
		 * This hashmap starts from targets.size()
		 */

		HashMap<Integer, Integer[]> ktoedge = SecurityGameContraction.fillKToEdge(nodewithcoins.size()); // 

		for(Integer k: ktoedge.keySet())
		{
			v[k] = ktoedge.get(k)[2]; // 2nd index has the path utility
		}


		/**
		 * fill d and e
		 */
		for(int i=0; i<duplicatetargets.size(); i++)
		{
			for(int j=0; j<duplicatetargets.size(); j++)
			{
			
				if(duplicatetargets.get(i).getNeighbors().contains(duplicatetargets.get(j)))  // if i-> j
				{


					double pathutility = duplicatetargets.get(i).getPathUtility(duplicatetargets.get(j));
					//System.out.println(i +" --> "+ j+ " : "+ duplicatetargets.get(i).getDistance(duplicatetargets.get(j)));
					if(pathutility==0)
					{ 
						/**
						 * means no intermediate targets
						 */
						//pathutility += targets.get(i).getAnimaldensity();

					}

					d[i][j] = duplicatetargets.get(i).getDistance(duplicatetargets.get(j));
					e[i][j] = 1;

				}
				else if(!duplicatetargets.get(i).getNeighbors().contains(duplicatetargets.get(j)))
				{


					e[i][j] = 0;
					d[i][j] = Double.MAX_VALUE;

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
				for(int k=0; k<q; k++)
				{
					if(k<nodewithcoins.size())
					{
						// for A...node coins.. fill B with 0
						if(duplicatetargets.get(i).getTargetid()==targets.get(k).getTargetid())
						{
							A[i][k] = 0;
						}
						else
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

						if(k<(nodewithcoins.size()))
						{
							//
							TargetNode kthnode = SecurityGameContraction.getNodeWithCoin(nodewithcoins.get(k));//targets.get(k);



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
						else if(k>=(nodewithcoins.size()))
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

						}

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
							System.out.println( i+ ", in:"+inode.isIn()+",out:"+inode.isOut()+
									"-->"+j+ ", in:"+jnode.isIn()+",out:"+jnode.isOut());


						}
					}
				}
				//System.out.println();
			}
			System.out.println();
			//System.out.println();

			for(int i=0; i<q; i++)
			{


				if(i<nodewithcoins.size())
				{
					System.out.print("y["+nodewithcoins.get(i)+"]="+cplex.getValue(y[i])+ " ");
				}
				else
				{
					Integer p[] = ktoedge.get(i);
					System.out.print("y["+p[0]+"-->"+p[1]+"]="+cplex.getValue(y[i])+ " ");
				}

				//System.out.println();
			}
			System.out.println();
			
			for(int i=0; i<v.length; i++)
			{


				
					System.out.print("v["+i+"]="+(v[i])+ " ");
				

				//System.out.println();
			}

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
			cplex.end();
			return objv;






		}
		catch (Exception ex)
		{

		}

		return -1.0;


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
				if(d[i][j]==1 || d[i][j]==0)
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
