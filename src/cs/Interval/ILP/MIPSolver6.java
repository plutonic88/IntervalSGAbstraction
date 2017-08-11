package cs.Interval.ILP;

import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloQuadNumExpr;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.HashMap;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;

public class MIPSolver6 {

	public static double[] solve(int[][] p, int[][] gamedata, ArrayList<TargetNode> targets, int nResource)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = SecurityGameContraction.targets.size();
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		
		System.out.println();
		
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{
			double dummyx[] = {0.5, 0.5, 0.0}; 

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar dd = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] a = new IloNumVar[nTargets];
			for(int i=0; i<nTargets; i++)
			{
				a[i] = cplex.boolVar();
			}


			/**
			 * objective
			 */
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1,  dd);
			
			cplex.addMaximize(obj);
			/**
			 * constraint 4
			 */
			//double sum =0;

			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);
				double prob = 0;
				double expectedpayoff=0;

				//IloLinearNumExpr expr3_1 = cplex.linearNumExpr();
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-p[target][jointpath]*(gamedata[target][0]-gamedata[target][1]), x[jointpath]);
					
					prob += (p[target][jointpath]*dummyx[jointpath]*(gamedata[target][0]-gamedata[target][1]));
					
					double tmp = (dummyx[jointpath])*gamedata[target][0]*p[target][jointpath]; //reward
					double tmp2 = (dummyx[jointpath])*gamedata[target][1]*(1-p[target][jointpath]); //penalty
					
					System.out.println("Target "+ target + ", jointpath "+jointpath+ 
							", reward "+ tmp + ", penalty "+ tmp2);
					expectedpayoff += tmp + tmp2; 


				}
				
				double exp2 = prob+ gamedata[target][1];
				
				
				System.out.println("Target "+ target + 
						", exppayoff***2 "+ exp2 +"");
				
				
				System.out.println("Target "+ target + 
						", exppayoff "+ expectedpayoff +"");
				
				/*double exp = prob*gamedata[target][0] + (1-prob)*gamedata[target][1];
				System.out.println("Target "+ target + 
						", sum exp "+ exp +"\n\n");*/
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+gamedata[target][1]);
			}

			//System.out.println("Sum constraint 4 "+sum);

			/**
			 * constraint 5
			 */

			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, kk);
				double prob =0;
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-p[target][jointpath]*(gamedata[target][3]-gamedata[target][2]), x[jointpath]);
					prob += (p[target][jointpath]*dummyx[jointpath]);
					
				}
				
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+gamedata[target][2]);
			}
			/**
			 * constraint 6
			 */
			
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				//double sum=0;
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(p[target][jointpath]*(gamedata[target][3]-gamedata[target][2]), x[jointpath]);
					
				}
				expr3.addTerm(-1, kk);
				cplex.addLe(expr3, -gamedata[target][2]);
			}

			//System.out.println("Sum constraint 6 "+sum);
			/**
			 * constraint 7
			 */
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for(int j=0; j<nJointSchedule; j++)
			{

				expr.addTerm(1.0, x[j]);
			}
			cplex.addEq(expr, 1.0);

			/**
			 * constraint for attakced target
			 */
			IloLinearNumExpr expra = cplex.linearNumExpr();
			for(int j=0; j<nTargets; j++)
			{

				expra.addTerm(1.0, a[j]);
			}
			cplex.addEq(expra, 1.0);





			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());
			double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				//if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			System.out.println();
			System.out.println();
			System.out.println();
			for(int i=0; i<nTargets; i++)
			{

				if(cplex.getValue(a[i])>0)
				{

					System.out.print("a["+mapback.get(i)+"]="+cplex.getValue(a[i])+ " \n");
				}

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}


}
