package cs.Interval.ILP;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import java.util.ArrayList;

import cs.Interval.contraction.TargetNode;

public class LpSolver
{
	public static double solve(ArrayList<TargetNode> targets, ArrayList<TargetNode> dominatedtargets, double distancelimit) throws IloException
	{
		int n = targets.size();
		double[][] w = new double[n][n];
		double[][] d = new double[n][n]; 
		/**
		 * fill w
		 */
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				if(targets.get(i).getNeighbors().contains(targets.get(j)))  // if i-> j
				{
					double icoin = targets.get(i).getAnimaldensity();
					double pathutility = targets.get(i).getPathUtility(targets.get(j));
					w[i][j] = icoin + pathutility;
					d[i][j] = targets.get(i).getDistance(targets.get(j));

				}
				else
				{
					w[i][j] = -1;
					d[i][j] = -1;
				}

			}
		}

		try
		{

			IloCplex cplex = new IloCplex();

			//variables
			IloNumVar[][] x = new IloNumVar[n][];
			for(int i=0; i<n; i++)
			{
				x[i] = cplex.boolVarArray(n);
			}

			IloLinearNumExpr obj = cplex.linearNumExpr();
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						obj.addTerm(w[i][j], x[i][j]);
					}
				}
			}
			cplex.addMaximize(obj);


			/**
			 * inflow equal outflow
			 * 
			 */
			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr1 = cplex.linearNumExpr();


				for(int j=0; j<n; j++)
				{

					expr1.addTerm(1.0, x[i][j]);
					expr1.addTerm(-1.0, x[j][i]);

				}
				cplex.addEq(expr1, 0.0);
			}

			for(int j=0; j<n; j++)
			{
				IloLinearNumExpr expr = cplex.linearNumExpr();
				for(int i=0; i<n; i++)
				{
					if(i!=j)
					{
						expr.addTerm(1.0, x[i][j]);
					}

				}
				cplex.addEq(expr, 1.0);
			}

			for(int i=0; i<n; i++)
			{
				IloLinearNumExpr expr = cplex.linearNumExpr();
				for(int j=0; j<n; j++)
				{
					if(i!=j)
					{
						expr.addTerm(1.0, x[i][j]);
					}

				}
				cplex.addEq(expr, 1.0);
			}

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

			cplex.solve();

			System.out.println(" obj: "+ cplex.getObjValue());
			
			
			/*for(int i=0; i<n; i++)
			{
				//IloLinearNumExpr expr = cplex.linearNumExpr();
				for(int j=0; j<n; j++)
				{
					System.out.print(x[i][j] + " ");

				}
				System.out.println();
				
			}*/

			return cplex.getObjValue();

		}
		catch(Exception ex)
		{




		}
		return 0;



	}

}
