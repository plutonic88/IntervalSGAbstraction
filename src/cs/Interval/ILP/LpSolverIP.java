package cs.Interval.ILP;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.IIS.Status;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;

public class LpSolverIP {

	public static double solve(ArrayList<TargetNode> targets, ArrayList<TargetNode> dominatedtargets, double distancelimit) throws IloException
	{
		int n = targets.size()-dominatedtargets.size();
		double[][] w = new double[n][n];
		double[][] d = new double[n][n]; 
		double [][] e = new double[n][n];
		double[] v = new double[n];
		int M = 1000000;
		int Q = 1000000;
		int startnode = 0;
		int basenodeid = 0;
		
		
		
		/*TargetNode basenode = targets.get(basenodeid);
		
		TargetNode virtualbase = new TargetNode(n, 0);
		
		ArrayList<TargetNode> pathtobase = new ArrayList<TargetNode>();
		
		
		
		virtualbase.setPath(basenode, pathtobase);
		virtualbase.setPathUtility(basenode, 0.0);
		virtualbase.addNeighbor(basenode);
		virtualbase.addDistance(basenode, 0.0);
		targets.add(virtualbase);
		
		basenode.addNeighbor(virtualbase);
		basenode.addDistance(virtualbase, 0.0);
		basenode.setPath(virtualbase, pathtobase);
		basenode.setPathUtility(virtualbase, 0.0);*/
		
		
		SecurityGameContraction.printNodesWithNeighborsAndPath(dominatedtargets, targets);
		
		

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
					System.out.println(i +" --> "+ j+ " : "+ targets.get(i).getDistance(targets.get(j)));
					if(pathutility==0)
					{ 
						/**
						 * means no intermediate targets
						 */
						//pathutility += targets.get(i).getAnimaldensity();
						
					}
					w[map.get(i)][map.get(j)] = pathutility;
					d[map.get(i)][map.get(j)] = targets.get(i).getDistance(targets.get(j));
					e[map.get(i)][map.get(j)] = 1;

				}
				else if(!targets.get(i).getNeighbors().contains(targets.get(j)) && (!dominatedtargets.contains(targets.get(i)) && !dominatedtargets.contains(targets.get(j)) ))
				{
					int b = map.get(i);
					int t = map.get(j);
					
					e[b][t] = 0;
					d[map.get(i)][map.get(j)] = Double.MAX_VALUE;
					w[map.get(i)][map.get(j)] = 0;
				}


			}
			if(!dominatedtargets.contains(targets.get(i)))
			{
				double icoin = targets.get(i).getAnimaldensity();
				v[map.get(i)] = icoin;
			}
		}

		try
		{

			IloCplex cplex = new IloCplex();
			
			List<IloRange> constraints = new ArrayList<IloRange>();

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
				//z[i] = cplex.boolVarArray(n);
 			}


			IloNumVar[] y = new IloNumVar[n];
			for(int i=0; i<n; i++)
			{
				y[i] = cplex.intVar(0, 1);
				//y[i] = cplex.boolVar();
			}


			int l = sets.size()-1;

			IloNumVar[] u = new IloNumVar[l*l];
			for(int i=0; i<(l*l); i++)
			{
				u[i] = cplex.intVar(0, 1);
				//u[i] = cplex.boolVar();
				
			}
			
			/*int len = 3*(n)*(n) + 2*n + 2*(sets.size()-1)*(sets.size()-1);
			
			int tcounter = 0;
			int K = 1000;*/
			
			
			/*IloNumVar[] t = new IloNumVar[len];
			int cost[] = new int[len];
			for(int i=0; i<(len); i++)
			{
				t[i] = cplex.intVar(0, 1);
				cost[i] = 30;
				//u[i] = cplex.boolVar();
				
			}*/
			
			
			


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
			for(int j=0; j<n; j++)
			{
				obj.addTerm(v[j], y[j]);
			}
			cplex.addMaximize(obj);
			
			/*for(int j=0; j<len; j++)
			{
				obj.addTerm(cost[j], t[j]);
			}
			cplex.addMinimize(obj);*/



			/*
			 * constraints
			 */
			
			IloLinearNumExpr exprgoal = cplex.linearNumExpr();
			for(int i=0; i<n; i++)
			{
				exprgoal.addTerm(1, x[startnode][i]);
			}
			constraints.add(cplex.addGe(exprgoal, 1));
			
			
			IloLinearNumExpr exprgoal1 = cplex.linearNumExpr();
			for(int i=0; i<n; i++)
			{
				exprgoal1.addTerm(1, x[i][startnode]);
			}
			constraints.add(cplex.addGe(exprgoal1, 1));



			
			/**
			 * constraint 3
			 */
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						//System.out.println("..3333... tcounter: "+ tcounter);
						//cplex.addLe(cplex.sum( cplex.prod(K, t[tcounter++]) , cplex.sum(cplex.prod(1, x[i][j]), (-M* e[i][j]))),0);
						constraints.add(cplex.addLe( cplex.sum(cplex.prod(1, x[i][j]), (-M* e[i][j])),0));
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
						//System.out.println("..4444... tcounter: "+ tcounter);
						//cplex.addLe(    cplex.sum( cplex.prod(K, t[tcounter++])  , cplex.sum(cplex.prod(1, z[i][j]), cplex.sum(cplex.prod(-1, x[i][j]), cplex.prod(-1, x[j][i])))),0);
						constraints.add(cplex.addLe(     cplex.sum(cplex.prod(1, z[i][j]), cplex.sum(cplex.prod(-1, x[i][j]), cplex.prod(-1, x[j][i]) )),0));

					}
				}
			}
			
			
			
			
			

			
			/**
			 * constraint 6
			 */
			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();


				for(int j=0; j<n; j++)
				{

					//expr1.addTerm(1.0, y[i]);
					if(i!=j)
						expr1.addTerm(-1.0, x[i][j]);

				}
				//System.out.println("...666.. tcounter: "+ tcounter + ", i :" + i);
				//expr1.addTerm(-K, t[tcounter++]);
				
				expr1.addTerm(1, y[i]);
				constraints.add(cplex.addLe(expr1, 0.0));
			}
			
			
			
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
						expr1.addTerm(1.0, x[j][i]);
						//expr1.addTerm(-1.0, x[j][i]);
					}

				}
				
				for(int k=0; k<n; k++)
				{
					if(i!=k)
					{
						expr1.addTerm(-1.0, x[i][k]);
						//expr1.addTerm(-1.0, x[j][i]);
					}
				}
				constraints.add(cplex.addEq(expr1, 0.0));
			}
			
			
			/**
			 * constraint 9
			 */

			IloLinearNumExpr expr3 = cplex.linearNumExpr();
			//IloLinearNumExpr expr31 = cplex.linearNumExpr();
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						//expr3.addTerm(-K, t[tcounter++]);
						expr3.addTerm(d[i][j], x[i][j]);
						//expr31.addTerm(d[i][j], x[i][j]);
					}
				}
			}
			
			/*int R = 1;
			IloNumVar g = cplex.numVar(0, 1);
			expr3.addTerm(-R, g);
			
			
		//	expr31.addTerm(R, g);
*/			//System.out.println(" tcounter: "+ tcounter);
			
			
			
			
			constraints.add(cplex.addLe(expr3, distancelimit));
			
		//	cplex.addLe(expr31, distancelimit+R);
			
			
			constraints.add(cplex.addEq(cplex.prod(1, y[startnode]), 1));
			
			
			
			/**
			 * constraint 10
			 */
			int ucount = 0;
			for(Set<Integer> s1 : sets)
			{

				for(Set<Integer> s2: sets)
				{

					if( !commonElement(s1,s2) && !s1.isEmpty() && !s2.isEmpty() && s1.size()>1 && s2.size()>1)
					{

						/**
						 * s1, s1
						 */
						for(int m: s1)
						{
							System.out.print(m+ ",");
							
						}
						System.out.println();
						for(int m: s2)
						{
							System.out.print(m+ ",");
							
						}
						System.out.println("\n");
						

						IloLinearNumExpr expr4 = cplex.linearNumExpr();
						
						ArrayList<Integer[]> pairs = new ArrayList<Integer[]>();
						
						//int icount = 0;
						for(Integer p: s1)
						{
							//int jcount = 0;
							for(Integer q: s1)
							{
								if((p!=q) && !donePair(p,q, pairs))
								{
									pairs.add(new Integer[]{p,q});
									expr4.addTerm(1, x[map.get(p)][map.get(q)]);
								}
								//jcount++;
							}
							//icount++;
						}
						//expr4.addTerm(-K, t[tcounter++]);
						//System.out.println("..10. ..1.. tcounter: "+ tcounter);
						expr4.addTerm(-Q,  u[ucount]);
						//constraints.add(cplex.addLe(expr4, 0));

						/**
						 * s2 s2
						 */
						pairs.clear();

						IloLinearNumExpr expr5 = cplex.linearNumExpr();
						for(Integer p: s2)
						{
							for(Integer q: s2)
							{
								if((p!=q) && !donePair(p,q, pairs))
								{
									pairs.add(new Integer[]{p,q});
									expr5.addTerm(1, x[map.get(p)][map.get(q)]);
								}
							}
						}
						//expr5.addTerm(-K, t[tcounter++]);
						//System.out.println("..10. ..2.. tcounter: "+ tcounter);
						expr5.addTerm(-Q,  u[ucount]);
					//	constraints.add(cplex.addLe(expr5, 0));




						/**
						 * s1 s2
						 */
						pairs.clear();

						IloLinearNumExpr expr6 = cplex.linearNumExpr();
						for(Integer p: s1)
						{
							for(Integer q: s2)
							{
								if((p!=q) && !donePair(p,q, pairs))
								{
									pairs.add(new Integer[]{p,q});
									expr6.addTerm(1, x[map.get(p)][map.get(q)]);
									//expr6.addTerm(1, x[map.get(q)][map.get(p)]);
								}
							}
						}
						expr6.addTerm(-1,  u[ucount]);
					//	constraints.add(cplex.addGe(expr6, 0));
						ucount++;

					}

					//jcount++;
				}
				//icount++;
			}

			cplex.solve();
			//cplex.get
			
						
			
			/*for(IloRange a: constraints)
			{
				System.out.println(cplex.getInfeasibility(a));
			}*/
			
			
			/*IloCplex.IIS q =  cplex.getIIS();
			
			Status[] s = q.getConstraintStatuses();
			
			for(int i=0; i<s.length; i++)
			{
				System.out.println(s[i]);
			}*/

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


			



			for(int i=0; i<n; i++)
			{


				System.out.print("y["+mapback.get(i)+"]="+cplex.getValue(y[i])+ " ");

				//System.out.println();
			}
			System.out.println();
			
			
			/*for(int i=0; i<len; i++)
			{


				System.out.print("t["+i+"]="+cplex.getValue(t[i])+ " ");

				//System.out.println();
			}
			System.out.println();
*/			
			
			
			
			
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

	private static boolean donePair(Integer p, Integer q,
			ArrayList<Integer[]> pairs) {
		
		
		for(Integer[] x: pairs)
		{
			if((x[0]==p && x[1]==q) || (x[1]==p && x[0]==q))
			{
				return true;
			}
		}
		return false;
	}

	private static boolean commonElement(Set<Integer> s1, Set<Integer> s2) {
		
		
		for(Integer x: s1)
		{
			for(Integer y : s2)
			{
				if(x==y)
					return true;
			}
			
		}
		
		
		return false;
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
