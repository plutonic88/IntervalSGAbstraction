package cs.Interval.ILP;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import cs.Interval.contraction.TargetNode;

public class LpSolverIP2 {
	
	public static double solve(ArrayList<TargetNode> targets, ArrayList<TargetNode> dominatedtargets, double distancelimit) throws IloException
	{
		int n = targets.size()-dominatedtargets.size();
		double[][] w = new double[n][n];
		double[][] d = new double[n][n]; 
		double [][] e = new double[n][n];
		double[] v = new double[n];
		int M = 10000;

		//int[] t = {0,1,2,3,4,5,6,7,8,9,10};

		Set<Integer> mySet = new HashSet<Integer>();
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		int icount =0;

		for(int i=0; i<targets.size(); i++)
		{
			if(!dominatedtargets.contains(targets.get(i)))
			{
				mySet.add( targets.get(i).getTargetid());
				map.put(targets.get(i).getTargetid(), icount);
				mapback.put(icount, targets.get(i).getTargetid());
				icount++;
			}
			
		}

		Set<Set<Integer>> sets = powerSet(mySet);






		/**
		 * fill w
		 */
		for(int i=0; i<targets.size(); i++)
		{
			for(int j=0; j<targets.size(); j++)
			{
				if(targets.get(i).getNeighbors().contains(targets.get(j)) && (!dominatedtargets.contains(targets.get(i))))  // if i-> j
				{


					double pathutility = targets.get(i).getPathUtility(targets.get(j));
					System.out.println(i +" --> "+ j+ " : "+ pathutility);
					double icoin = targets.get(i).getAnimaldensity();

					w[map.get(i)][map.get(j)] =icoin + pathutility;
					d[map.get(i)][map.get(j)] = targets.get(i).getDistance(targets.get(j));
					e[map.get(i)][map.get(j)] = 1;

				}
				else if(!targets.get(i).getNeighbors().contains(targets.get(j)) && (!dominatedtargets.contains(targets.get(i)) && !dominatedtargets.contains(targets.get(j)) ))
				{
					int b = map.get(i);
					int t = map.get(j);
					
					d[map.get(i)][map.get(j)] = Integer.MAX_VALUE;
					w[map.get(i)][map.get(j)] = Integer.MIN_VALUE;
					
					e[b][t] = 0;
				}


			}
			/*if(!dominatedtargets.contains(targets.get(i)))
			{
				double icoin = targets.get(i).getAnimaldensity();
				v[map.get(i)] = icoin;
			}*/
		}

		try
		{

			IloCplex cplex = new IloCplex();

			//variables
			IloNumVar[][] x = new IloNumVar[n][];
			for(int i=0; i<n; i++)
			{
				x[i] = cplex.intVarArray(n, 0, Integer.MAX_VALUE);
			}


			

			IloNumVar[][] z = new IloNumVar[n][];
			for(int i=0; i<n; i++)
			{
				z[i] = cplex.intVarArray(n, 0, 1);
			}


			/*IloNumVar[] y = new IloNumVar[n];
			for(int i=0; i<n; i++)
			{
				y[i] = cplex.intVar(0, 1);
			}*/


			int l = sets.size()-1;

			IloNumVar[] u = new IloNumVar[l*l];
			for(int i=0; i<(l*l); i++)
			{
				u[i] = cplex.intVar(0, 1);
			}


			//Iloi




			//objective

			IloLinearNumExpr obj = cplex.linearNumExpr();
			for(int i=0; i<(n-1); i++)
			{
				for(int j=(i+1); j<n; j++)
				{
					obj.addTerm(w[i][j], z[i][j]);
				}
			}
			
			cplex.addMaximize(obj);



			/*
			 * constraints
			 */


			
			/**
			 * constraint 3
			 */
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						cplex.addLe(cplex.sum(cplex.prod(1, x[i][j]), (-M* e[i][j])),0);

					}
				}
			}
			
			
			/**
			 * constraint 4
			 */
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						cplex.addLe(cplex.sum(cplex.prod(1, z[i][j]), cplex.sum(cplex.prod(-1, x[i][j]), cplex.prod(-1, x[j][i]))),0);

					}
				}
			}

			
			/**
			 * constraint 6
			 *//*
			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();


				for(int j=0; j<n; j++)
				{

					//expr1.addTerm(1.0, y[i]);
					if(i!=j)
						expr1.addTerm(-1.0, x[i][j]);

				}
				expr1.addTerm(1, y[i]);
				cplex.addLe(expr1, 0.0);
			}*/
			
			
			
			/**
			 * constraint 8
			 */

			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();


				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{

						expr1.addTerm(1.0, x[i][j]);
						expr1.addTerm(-1.0, x[j][i]);
					}

				}
				cplex.addEq(expr1, 0.0);
			}
			
			
			
			/**
			 * base node constraint
			 */
			/*int r=0;
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();
				IloLinearNumExpr expr2 = cplex.linearNumExpr();


				for(int j=0; j<n; j++)
				{
					if(r!=j)
					{

						expr1.addTerm(1.0, x[r][j]);
						expr2.addTerm(1.0, x[j][r]);
					}

				}
				cplex.addGe(expr1, 1);
				cplex.addGe(expr2, 1);
			}*/
			
			
			
			
			
			
			/**
			 * constraint 9
			 */

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
			
			
			
			/**
			 * constraint 10
			 */
			int ucount = 0;
			for(Set<Integer> s1 : sets)
			{

				for(Set<Integer> s2: sets)
				{

					if(!s1.isEmpty() && !s2.isEmpty())
					{

						/**
						 * s1, s1
						 */

						IloLinearNumExpr expr4 = cplex.linearNumExpr();
						//int icount = 0;
						for(Integer p: s1)
						{
							//int jcount = 0;
							for(Integer q: s1)
							{
								if(p!=q)
								{
									expr4.addTerm(1, x[map.get(p)][map.get(q)]);
								}
								//jcount++;
							}
							//icount++;
						}
						expr4.addTerm(-M,  u[ucount]);
						cplex.addLe(expr4, 0);

						/**
						 * s2 s2
						 */

						IloLinearNumExpr expr5 = cplex.linearNumExpr();
						for(Integer p: s2)
						{
							for(Integer q: s2)
							{
								if(p!=q)
								{
									expr5.addTerm(1, x[map.get(p)][map.get(q)]);
								}
							}
						}
						expr5.addTerm(-M,  u[ucount]);
						cplex.addLe(expr5, 0);




						/**
						 * s1 s2
						 */

						IloLinearNumExpr expr6 = cplex.linearNumExpr();
						for(Integer p: s1)
						{
							for(Integer q: s2)
							{
								if(p!=q)
								{
									expr6.addTerm(1, x[map.get(p)][map.get(q)]);
								}
							}
						}
						expr6.addTerm(-1,  u[ucount]);
						cplex.addGe(expr6, 0);
						ucount++;

					}

					//jcount++;
				}
				//icount++;
			}

			cplex.solve();

			System.out.println(" obj: "+ cplex.getObjValue());

			for(int i=0; i<n; i++)
			{

				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						System.out.print("x["+mapback.get(i)+"]["+mapback.get(j)+"]="+cplex.getValue(x[i][j])+ " ");
					}
				}
				System.out.println();
			}

			System.out.println();




			for(int i=0; i<n; i++)
			{

				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						System.out.print("z["+mapback.get(i)+"]["+mapback.get(j)+"]="+cplex.getValue(z[i][j])+ " ");
					}
				}
				System.out.println();

			}
			System.out.println();


			



			/*for(int i=0; i<n; i++)
			{


				System.out.print("y["+mapback.get(i)+"]="+cplex.getValue(y[i])+ " ");

				//System.out.println();
			}
			System.out.println();*/
			
			/*for(int i=0; i<(l*l); i++)
			{
				System.out.print("u["+i+"]="+cplex.getValue(u[i])+ " ");
			}*/


			return cplex.getObjValue();

		}
		catch(Exception ex)
		{
			ex.printStackTrace();




		}
		return 0;



	}

	public static Set<Set<Integer>> powerSet(Set<Integer> originalSet) {
		Set<Set<Integer>> sets = new HashSet<Set<Integer>>();
		if (originalSet.isEmpty()) {
			sets.add(new HashSet<Integer>());
			return sets;
		}
		List<Integer> list = new ArrayList<Integer>(originalSet);
		Integer head = list.get(0);
		Set<Integer> rest = new HashSet<Integer>(list.subList(1, list.size()));
		for (Set<Integer> set : powerSet(rest)) {
			Set<Integer> newSet = new HashSet<Integer>();
			newSet.add(head);
			newSet.addAll(set);
			sets.add(newSet);
			sets.add(set);
		}
		return sets;
	}

}
