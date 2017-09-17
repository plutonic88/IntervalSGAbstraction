package cs.Interval.ILP;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import cs.Interval.contraction.TargetNode;
import groupingtargets.SuperTarget;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloQuadNumExpr;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

public class MIPSolver4 {

	public static int attackedtarget=-1;




	public static ArrayList<ArrayList<Integer>> slavePath(int nPath, int[][] gamedata,
			ArrayList<TargetNode> targets, int nResource, int nOrigtargets, double dmax)
	{

		/**
		 * make a copy node same as the base node,
		 * add the node to the neighbors of the base node
		 * After finishing the work remove the copy node
		 * 
		 * 
		 * OR add a copy node information in the matrix
		 */


		/**
		 * create the copy node
		 *//*
		TargetNode copynode = new TargetNode(nOrigtargets, targets.get(0).getAnimaldensity());


		  *//**
		  * add the node to the neighbors of the base node
		  * add distances and path
		  *//*



		for(TargetNode neighbor: basenode.getNeighbors())
		{
			// add the copy node to neighbor's neighbor
			neighbor.addNeighbor(copynode);
			// add distance
			neighbor.addDistance(copynode, basenode.getDistance(neighbor));
			// add path


		}*/

		int copynodeid = nOrigtargets;
		int nTargets = targets.size()+1;
		TargetNode basenode = targets.get(0);
		int [][] d= new int[nTargets][nTargets];
		int [] s = new int[nTargets];

		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		// add the copy node
		mapback.put(icount, copynodeid);
		map.put(copynodeid, icount);

		System.out.println("#origtargets "+ nOrigtargets);
		System.out.println("#targets(including copy) "+ nTargets);
		System.out.println("#targets Mapping");
		for(int i=0; i<nTargets; i++)
		{
			System.out.println(i+"<--"+mapback.get(i));
		}



		// distance variable d[i,j]
		for(int i=0; i<nTargets-1; i++)
		{
			for(int j=0; j<nTargets-1; j++)
			{

				TargetNode inode = targets.get(i);
				TargetNode jnode = targets.get(j);

				if(inode.getNeighbors().contains(jnode))
				{
					d[i][j] = inode.getDistance(jnode).intValue();
				}
				else
				{
					d[i][j] = Integer.MAX_VALUE;
				}
			}

		}
		/*
		 * now add the copy node info to d
		 * for all the neighbor nodes of basenode add distance to copy node
		 */

		// add distance for neighbor nodes

		for(TargetNode nei: basenode.getNeighbors())
		{
			int id = map.get(nei.getTargetid());
			d[id][d.length-1] = basenode.getDistance(nei).intValue();
			d[d.length-1][id] = d[id][d.length-1];
		}
		for(int i=0; i<nTargets; i++)
		{
			if(!(d[i][d.length-1]>0) || !(d[d.length-1][i]>0))
			{
				d[i][d.length-1] = Integer.MAX_VALUE;
				d[d.length-1][i] = Integer.MAX_VALUE;

			}
		}
		for(int i=0; i<nTargets-1; i++)
		{
			s[i] = (int)Math.floor(targets.get(i).getAnimaldensity()) ;
		}
		//copy node
		s[nTargets-1] = (int)Math.floor(basenode.getAnimaldensity());
		printArray(d);
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();



			IloNumVar[][] x_ = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				x_[i] = cplex.boolVarArray(nTargets);
			}








			/*
			 * xijp variable
			 */

			IloNumVar[][] x = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				x[i] = cplex.numVarArray(nTargets, 0, Integer.MAX_VALUE, IloNumVarType.Int);
			}



			/**
			 * yip variable
			 */
			IloNumVar[] y = cplex.boolVarArray(nTargets);

			/**
			 * uip variable
			 */
			IloNumVar[] u = cplex.numVarArray(nTargets, 0, Integer.MAX_VALUE);


			// objective

			IloLinearNumExpr obj = cplex.linearNumExpr();

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{


					if(i!=j)
						obj.addTerm(s[i], x_[i][j]);

				}
			}

			cplex.addMaximize(obj);







			// loop is not allowed
			IloLinearNumExpr expr13 = cplex.linearNumExpr();

			for(int j=0; j<nTargets; j++)
			{
				expr13.addTerm(1.0, x[j][j]);
			}
			cplex.addEq(expr13, 0);



			//constraint for start at base node 0 and end in copy node
			IloLinearNumExpr expr1 = cplex.linearNumExpr();

			for(int j=1; j<nTargets; j++)
			{
				expr1.addTerm(1.0, x[0][j]);
			}
			cplex.addEq(expr1, 1);


			IloLinearNumExpr expr9 = cplex.linearNumExpr();

			for(int j=0; j<nTargets-1; j++)
			{
				expr9.addTerm(1.0, x[j][nTargets-1]);
			}

			cplex.addEq(expr9, 1);





			IloLinearNumExpr expr12 = cplex.linearNumExpr();


			for(int j=1; j<nTargets; j++)
			{
				expr12.addTerm(1.0, x[0][j]);
			}
			for(int j=0; j<nTargets-1; j++)
			{
				expr12.addTerm(-1.0, x[j][nTargets-1]);
			}

			cplex.addEq(expr12, 0);






			/**
			 * every vertex is visited once
			 */



			/*for(int j=0; j<nTargets; j++)
			{
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				expr2.addTerm(1.0, y[j]);
				cplex.addLe(expr2, 1);
			}
			 */


			/**
			 * connectivity of path
			 */

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{

					if(i!=k)
						expr3.addTerm(1.0, x_[i][k]);
				}
				//expr3.addTerm(-1.0, y[k]);
				//	cplex.addLe(expr3, 1);


			}

			for(int k=1; k<nTargets-1; k++)
			{
				IloLinearNumExpr expr8 = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					if(j!=k)
						expr8.addTerm(1.0, x_[k][j]);
				}
				//expr8.addTerm(-1.0, y[k]);
				//	cplex.addLe(expr8, 1);
			}


			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{
					if(i!=k)
						expr3.addTerm(1.0, x_[i][k]);
				}
				//expr3.addTerm(-1.0, y[k]);

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
						expr3.addTerm(-1.0, x_[k][j]);
				}



				cplex.addEq(expr3, 0);


			}



			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{


				for(int i=0; i<nTargets-1; i++)
				{

					if(i!=k)
					{
						IloLinearNumExpr expr3 = cplex.linearNumExpr();
						expr3.addTerm(1.0, x[i][k]);
						cplex.addLe(expr3, 5);
					}
				}
				//expr3.addTerm(-1.0, y[k]);

			}

			for(int k=1; k<nTargets-1; k++)
			{

				for(int j=1; j<nTargets; j++)
				{
					if(j!=k)
					{
						IloLinearNumExpr expr8 = cplex.linearNumExpr();
						expr8.addTerm(1.0, x[k][j]);
						cplex.addLe(expr8, 5);
					}
				}
				//expr8.addTerm(-1.0, y[k]);

			}





			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{
					if(i!=k)
						expr3.addTerm(1.0, x[i][k]);
				}
				//expr3.addTerm(-1.0, y[k]);

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
						expr3.addTerm(-1.0, x[k][j]);
				}



				cplex.addEq(expr3, 0);


			}










			// x and x_

			for(int k=0; k<nTargets-1; k++)
			{

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
					{
						IloLinearNumExpr expr10 = cplex.linearNumExpr();
						expr10.addTerm(1.0, x_[k][j]);
						expr10.addTerm(-1.0, x[k][j]);
						cplex.addLe(expr10, 0);
					}
				}

				//expr10.addTerm(1.0, y[k]);



			}


			for(int k=0; k<nTargets-1; k++)
			{

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
					{
						IloLinearNumExpr expr11 = cplex.linearNumExpr();
						expr11.addTerm(1.0, x_[k][j]);
						cplex.addLe(expr11, 1);
					}
				}
				//expr11.addTerm(1.0, y[k]);


			}








			//path length constraint


			IloLinearNumExpr expr4 = cplex.linearNumExpr();
			for(int j=0; j<nTargets-1; j++)
			{

				for(int k=1; k<nTargets; k++)
				{
					if(k!=j)
						expr4.addTerm(d[j][k], x[j][k]);
				}

			}
			cplex.addEq(expr4, dmax);




			IloLinearNumExpr expr14 = cplex.linearNumExpr();
			for(int j=0; j<nTargets-1; j++)
			{

				for(int k=1; k<nTargets; k++)
				{
					if(k!=j)
						expr14.addTerm(d[j][k], x_[j][k]);
				}

			}
			//cplex.addLe(expr14, dmax);




			//subtour constraints

			for(int j=0; j<nTargets; j++)
			{

				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				expr5.addTerm(1, u[j]);
				cplex.addGe(expr5, 2);  // reason lol


			}


			for(int j=0; j<nTargets; j++)	
			{
				IloLinearNumExpr expr6 = cplex.linearNumExpr();
				expr6.addTerm(1, u[j]);
				cplex.addLe(expr6, nTargets);

			}




			for(int j=0; j<nTargets; j++)
			{

				for(int k=0; k<nTargets; k++)
				{

					if(j!=k)
					{

						IloLinearNumExpr expr7 = cplex.linearNumExpr();
						expr7.addTerm(1, u[j]);
						expr7.addTerm(-1, u[k]);
						expr7.addTerm(nTargets-1, x_[j][k]);
						cplex.addLe(expr7, nTargets-2);
					}


				}

			}

			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());
			/*double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}*/
			System.out.println();
			System.out.println();
			System.out.println();




			/*for(int i=0; i<nTargets; i++)
			{

				//if(cplex.getValue(y[0][i])>0)
				{

					System.out.print("y["+mapback.get(i)+"]="+cplex.getValue(y[i])+ " \n");
					//attackedtarget=mapback.get(i);
					//break;
				}

			}
			System.out.println();*/




			//for(int p=0; p<nPath; p++)


			//System.out.print("Path "+ p+ ": ");
			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{
					//System.out.println(i+"-->"+j);

					if(i!=j)
						if(cplex.getValue(x_[i][j])>0)
						{

							System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x_[i][j])+" \n");
						}
					//attackedtarget=mapback.get(i);
					//break;
				}

			}
			System.out.println();

			double sum =0;

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{

					if(i!=j)
						if(i!=j && cplex.getValue(x[i][j])>0)
						{


							int tmp= (int)Math.round(cplex.getValue(x[i][j]));

							System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x[i][j])/* + "  "+d[i][j] +" x " +tmp*/+" \n");

							sum += d[i][j]*tmp;
						}
					//attackedtarget=mapback.get(i);
					//break;
				}

			}







			return null;

		}
		catch(Exception ex)
		{

		}
		return null;


	}







	public static ArrayList<ArrayList<Integer>> modifiedTOP(int nPath, int[][] gamedata,
			ArrayList<TargetNode> targets, int nResource, int nOrigtargets, double dmax)
	{

		/**
		 * make a copy node same as the base node,
		 * add the node to the neighbors of the base node
		 * After finishing the work remove the copy node
		 * 
		 * 
		 * OR add a copy node information in the matrix
		 */


		/**
		 * create the copy node
		 *//*
		TargetNode copynode = new TargetNode(nOrigtargets, targets.get(0).getAnimaldensity());


		  *//**
		  * add the node to the neighbors of the base node
		  * add distances and path
		  *//*



		for(TargetNode neighbor: basenode.getNeighbors())
		{
			// add the copy node to neighbor's neighbor
			neighbor.addNeighbor(copynode);
			// add distance
			neighbor.addDistance(copynode, basenode.getDistance(neighbor));
			// add path


		}*/

		int copynodeid = nOrigtargets;
		int nTargets = targets.size()+1;
		TargetNode basenode = targets.get(0);
		int [][] d= new int[nTargets][nTargets];
		int [] s = new int[nTargets];

		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		// add the copy node
		mapback.put(icount, copynodeid);
		map.put(copynodeid, icount);

		System.out.println("#origtargets "+ nOrigtargets);
		System.out.println("#targets(including copy) "+ nTargets);
		System.out.println("#targets Mapping");
		for(int i=0; i<nTargets; i++)
		{
			System.out.println(i+"<--"+mapback.get(i));
		}



		// distance variable d[i,j]
		for(int i=0; i<nTargets-1; i++)
		{
			for(int j=0; j<nTargets-1; j++)
			{

				TargetNode inode = targets.get(i);
				TargetNode jnode = targets.get(j);

				if(inode.getNeighbors().contains(jnode))
				{
					d[i][j] = inode.getDistance(jnode).intValue();
				}
				else
				{
					d[i][j] = Integer.MAX_VALUE;
				}
			}

		}
		/*
		 * now add the copy node info to d
		 * for all the neighbor nodes of basenode add distance to copy node
		 */

		// add distance for neighbor nodes

		for(TargetNode nei: basenode.getNeighbors())
		{
			int id = map.get(nei.getTargetid());
			d[id][d.length-1] = basenode.getDistance(nei).intValue();
			d[d.length-1][id] = d[id][d.length-1];
		}
		for(int i=0; i<nTargets; i++)
		{
			if(!(d[i][d.length-1]>0) || !(d[d.length-1][i]>0))
			{
				d[i][d.length-1] = Integer.MAX_VALUE;
				d[d.length-1][i] = Integer.MAX_VALUE;

			}
		}
		for(int i=0; i<nTargets-1; i++)
		{
			s[i] = (int)Math.floor(targets.get(i).getAnimaldensity()) ;
		}
		//copy node
		s[nTargets-1] = (int)Math.floor(basenode.getAnimaldensity());
		printArray(d);
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();



			IloNumVar[][] x_ = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				x_[i] = cplex.boolVarArray(nTargets);
			}

			IloNumVar[][] x = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				x[i] = cplex.boolVarArray(nTargets);
			}

			IloNumVar[][] z = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				z[i] = cplex.boolVarArray(nTargets);
			}




			/**
			 * uip variable
			 */
			IloNumVar[] u = cplex.numVarArray(nTargets, 0, Integer.MAX_VALUE);

			IloNumVar[] u_ = cplex.numVarArray(nTargets, 0, Integer.MAX_VALUE);


			// objective

			IloLinearNumExpr obj = cplex.linearNumExpr();

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{


					if(i!=j)
					{
						obj.addTerm(s[i], z[i][j]);
						//obj.addTerm(s[i], x[i][j]);
					}


				}
			}

			cplex.addMaximize(obj);
			
			
			for(int i=0; i<nTargets-1; i++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					if(i!=j)
					{
						expr3.addTerm(1.0, z[i][j]);
					}
				}
				cplex.addLe(expr3, 1);
			}
			
			for(int i=0; i<nTargets-1; i++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					if(i!=j)
					{
						expr3.addTerm(1.0, z[i][j]);
						expr3.addTerm(-1.0, x[i][j]);
						expr3.addTerm(-1.0, x_[i][j]);
						
					}
				}
				cplex.addLe(expr3, 0);
			}
			
			


			//constraint for start at base node 0 and end in copy node
			IloLinearNumExpr expr1 = cplex.linearNumExpr();

			for(int j=1; j<nTargets; j++)
			{
				expr1.addTerm(1.0, x_[0][j]);
			}
			cplex.addEq(expr1, 1);

			IloLinearNumExpr ex1 = cplex.linearNumExpr();

			for(int j=1; j<nTargets; j++)
			{
				ex1.addTerm(1.0, x[0][j]);
			}
			cplex.addEq(ex1, 1);




			IloLinearNumExpr expr9 = cplex.linearNumExpr();

			for(int j=0; j<nTargets-1; j++)
			{
				expr9.addTerm(1.0, x_[j][nTargets-1]);
			}

			cplex.addEq(expr9, 1);

			IloLinearNumExpr ex9 = cplex.linearNumExpr();

			for(int j=0; j<nTargets-1; j++)
			{
				ex9.addTerm(1.0, x[j][nTargets-1]);
			}

			cplex.addEq(ex9, 1);





			IloLinearNumExpr expr12 = cplex.linearNumExpr();


			for(int j=1; j<nTargets; j++)
			{
				expr12.addTerm(1.0, x_[0][j]);
			}
			for(int j=0; j<nTargets-1; j++)
			{
				expr12.addTerm(-1.0, x_[j][nTargets-1]);
			}

			cplex.addEq(expr12, 0);


			IloLinearNumExpr ex12 = cplex.linearNumExpr();


			for(int j=1; j<nTargets; j++)
			{
				ex12.addTerm(1.0, x[0][j]);
			}
			for(int j=0; j<nTargets-1; j++)
			{
				ex12.addTerm(-1.0, x[j][nTargets-1]);
			}

			cplex.addEq(ex12, 0);






			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{

					if(i!=k)
					{

						expr3.addTerm(1.0, x_[i][k]);

					}
				}
				cplex.addLe(expr3, 1);
				//expr3.addTerm(-1.0, y[k]);

			}

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{

					if(i!=k)
					{

						expr3.addTerm(1.0, x[i][k]);

					}
				}
				cplex.addLe(expr3, 1);
				//expr3.addTerm(-1.0, y[k]);

			}




			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr8 = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					if(j!=k)
					{

						expr8.addTerm(1.0, x_[k][j]);

					}
				}
				cplex.addLe(expr8, 1);
				//expr8.addTerm(-1.0, y[k]);

			}

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr8 = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					if(j!=k)
					{

						expr8.addTerm(1.0, x[k][j]);

					}
				}
				cplex.addLe(expr8, 1);
				//expr8.addTerm(-1.0, y[k]);

			}






			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{
					if(i!=k)
						expr3.addTerm(1.0, x_[i][k]);
				}
				//expr3.addTerm(-1.0, y[k]);

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
						expr3.addTerm(-1.0, x_[k][j]);
				}
				cplex.addEq(expr3, 0);

			}

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{
					if(i!=k)
						expr3.addTerm(1.0, x[i][k]);
				}
				//expr3.addTerm(-1.0, y[k]);

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
						expr3.addTerm(-1.0, x[k][j]);
				}
				cplex.addEq(expr3, 0);

			}





			IloLinearNumExpr expr14 = cplex.linearNumExpr();
			for(int j=0; j<nTargets-1; j++)
			{

				for(int k=1; k<nTargets; k++)
				{
					if(k!=j)
						expr14.addTerm(d[j][k], x_[j][k]);
				}

			}
			cplex.addLe(expr14, dmax);


			IloLinearNumExpr ex14 = cplex.linearNumExpr();
			for(int j=0; j<nTargets-1; j++)
			{

				for(int k=1; k<nTargets; k++)
				{
					if(k!=j)
						ex14.addTerm(d[j][k], x[j][k]);
				}

			}
			cplex.addLe(ex14, dmax);




			//subtour constraints

			for(int j=0; j<nTargets; j++)
			{

				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				expr5.addTerm(1, u[j]);
				cplex.addGe(expr5, 2);  // reason lol


			}


			for(int j=0; j<nTargets; j++)
			{

				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				expr5.addTerm(1, u_[j]);
				cplex.addGe(expr5, 2);  // reason lol


			}





			for(int j=0; j<nTargets; j++)	
			{
				IloLinearNumExpr expr6 = cplex.linearNumExpr();
				expr6.addTerm(1, u[j]);
				cplex.addLe(expr6, nTargets);

			}

			for(int j=0; j<nTargets; j++)	
			{
				IloLinearNumExpr expr6 = cplex.linearNumExpr();
				expr6.addTerm(1, u_[j]);
				cplex.addLe(expr6, nTargets);

			}





			for(int j=0; j<nTargets; j++)
			{

				for(int k=0; k<nTargets; k++)
				{

					if(j!=k)
					{

						IloLinearNumExpr expr7 = cplex.linearNumExpr();
						expr7.addTerm(1, u[j]);
						expr7.addTerm(-1, u[k]);
						expr7.addTerm(nTargets-1, x_[j][k]);
						cplex.addLe(expr7, nTargets-2);
					}


				}

			}

			for(int j=0; j<nTargets; j++)
			{

				for(int k=0; k<nTargets; k++)
				{

					if(j!=k)
					{

						IloLinearNumExpr expr7 = cplex.linearNumExpr();
						expr7.addTerm(1, u_[j]);
						expr7.addTerm(-1, u_[k]);
						expr7.addTerm(nTargets-1, x[j][k]);
						cplex.addLe(expr7, nTargets-2);
					}


				}

			}


			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());

			System.out.println();
			System.out.println();
			System.out.println();







			double sum =0;

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{

					if(i!=j)
						if(i!=j && cplex.getValue(x_[i][j])>0)
						{


							int tmp= (int)Math.round(cplex.getValue(x_[i][j]));

							System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x_[i][j])/* + "  "+d[i][j] +" x " +tmp*/+" \n");

							sum += d[i][j]*tmp;
						}
					//attackedtarget=mapback.get(i);
					//break;
				}

			}


			System.out.println("Another path ");

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{

					if(i!=j)
						if(i!=j && cplex.getValue(x[i][j])>0)
						{


							int tmp= (int)Math.round(cplex.getValue(x[i][j]));

							System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x[i][j])/* + "  "+d[i][j] +" x " +tmp*/+" \n");

							sum += d[i][j]*tmp;
						}
					//attackedtarget=mapback.get(i);
					//break;
				}

			}


			ArrayList<Integer> pathtoconsider = new ArrayList<Integer>();
			int current = 0; // basenode
			pathtoconsider.add(current);
			int flag = 0;

			while(flag!=1)
			{
				// find edge that starts with current and end with a different node y. 
				/*
				 * if current is not base and y is base then end the loop
				 */

				int yid = findPath(cplex, x_, current, nTargets, copynodeid, mapback);


				if(current != 0 && yid == 0)
				{
					flag = 1;
				}

				current = yid;
				if(!pathtoconsider.contains(yid))
				{
					pathtoconsider.add(yid);
				}



			}
			pathtoconsider.add(0);
			
			
			
			
			ArrayList<Integer> pathtoconsider2 = new ArrayList<Integer>();
			current = 0; // basenode
			pathtoconsider2.add(current);
			flag = 0;

			while(flag!=1)
			{
				// find edge that starts with current and end with a different node y. 
				/*
				 * if current is not base and y is base then end the loop
				 */

				int yid = findPath(cplex, x, current, nTargets, copynodeid, mapback);


				if(current != 0 && yid == 0)
				{
					flag = 1;
				}

				current = yid;
				if(!pathtoconsider2.contains(yid))
				{
					pathtoconsider2.add(yid);
				}



			}
			pathtoconsider2.add(0);




			//////////



			ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();

			res.add(pathtoconsider);
			res.add(pathtoconsider2);	
			return res;

		}
		catch(Exception ex)
		{

		}
		return null;


	}





	public static ArrayList<ArrayList<Integer>> lexicographicOP(int nPath, int[][] gamedata,
			ArrayList<TargetNode> targets, int nResource, int nOrigtargets, double dmax)
	{

		/**
		 * make a copy node same as the base node,
		 * add the node to the neighbors of the base node
		 * After finishing the work remove the copy node
		 * 
		 * 
		 * OR add a copy node information in the matrix
		 */


		/**
		 * create the copy node
		 *//*
		TargetNode copynode = new TargetNode(nOrigtargets, targets.get(0).getAnimaldensity());


		  *//**
		  * add the node to the neighbors of the base node
		  * add distances and path
		  *//*



		for(TargetNode neighbor: basenode.getNeighbors())
		{
			// add the copy node to neighbor's neighbor
			neighbor.addNeighbor(copynode);
			// add distance
			neighbor.addDistance(copynode, basenode.getDistance(neighbor));
			// add path


		}*/

		int copynodeid = nOrigtargets;
		int nTargets = targets.size()+1;
		TargetNode basenode = targets.get(0);
		int [][] d= new int[nTargets][nTargets];
		int [] s = new int[nTargets];

		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		// add the copy node
		mapback.put(icount, copynodeid);
		map.put(copynodeid, icount);

		System.out.println("#origtargets "+ nOrigtargets);
		System.out.println("#targets(including copy) "+ nTargets);
		System.out.println("#targets Mapping");
		for(int i=0; i<nTargets; i++)
		{
			System.out.println(i+"<--"+mapback.get(i));
		}



		// distance variable d[i,j]
		for(int i=0; i<nTargets-1; i++)
		{
			for(int j=0; j<nTargets-1; j++)
			{

				TargetNode inode = targets.get(i);
				TargetNode jnode = targets.get(j);

				if(inode.getNeighbors().contains(jnode))
				{
					d[i][j] = inode.getDistance(jnode).intValue();
				}
				else
				{
					d[i][j] = Integer.MAX_VALUE;
				}
			}

		}
		/*
		 * now add the copy node info to d
		 * for all the neighbor nodes of basenode add distance to copy node
		 */

		// add distance for neighbor nodes

		for(TargetNode nei: basenode.getNeighbors())
		{
			int id = map.get(nei.getTargetid());
			d[id][d.length-1] = basenode.getDistance(nei).intValue();
			d[d.length-1][id] = d[id][d.length-1];
		}
		for(int i=0; i<nTargets; i++)
		{
			if(!(d[i][d.length-1]>0) || !(d[d.length-1][i]>0))
			{
				d[i][d.length-1] = Integer.MAX_VALUE;
				d[d.length-1][i] = Integer.MAX_VALUE;

			}
		}
		for(int i=0; i<nTargets-1; i++)
		{
			s[i] = (int)Math.floor(targets.get(i).getAnimaldensity()) ;
		}
		//copy node
		s[nTargets-1] = (int)Math.floor(basenode.getAnimaldensity());
		printArray(d);
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();


			// start a loop 

			for(int r=0; r<nResource; r++)
			{

				IloCplex cplex = new IloCplex();



				IloNumVar[][] x_ = new IloNumVar[nTargets][];
				for(int i=0; i<nTargets; i++)
				{
					x_[i] = cplex.boolVarArray(nTargets);
				}


				/**
				 * uip variable
				 */
				IloNumVar[] u = cplex.numVarArray(nTargets, 0, Integer.MAX_VALUE);


				// objective

				IloLinearNumExpr obj = cplex.linearNumExpr();

				for(int i=0; i<nTargets-1; i++)
				{

					for(int j=1; j<nTargets; j++)
					{


						if(i!=j)
							obj.addTerm(s[i], x_[i][j]);

					}
				}

				cplex.addMaximize(obj);


				//constraint for start at base node 0 and end in copy node
				IloLinearNumExpr expr1 = cplex.linearNumExpr();

				for(int j=1; j<nTargets; j++)
				{
					expr1.addTerm(1.0, x_[0][j]);
				}
				cplex.addEq(expr1, 1);


				IloLinearNumExpr expr9 = cplex.linearNumExpr();

				for(int j=0; j<nTargets-1; j++)
				{
					expr9.addTerm(1.0, x_[j][nTargets-1]);
				}

				cplex.addEq(expr9, 1);





				IloLinearNumExpr expr12 = cplex.linearNumExpr();


				for(int j=1; j<nTargets; j++)
				{
					expr12.addTerm(1.0, x_[0][j]);
				}
				for(int j=0; j<nTargets-1; j++)
				{
					expr12.addTerm(-1.0, x_[j][nTargets-1]);
				}

				cplex.addEq(expr12, 0);




				//connectivity

				for(int k=1; k<nTargets-1; k++)
				{

					IloLinearNumExpr expr3 = cplex.linearNumExpr();
					for(int i=0; i<nTargets-1; i++)
					{

						if(i!=k)
						{

							expr3.addTerm(1.0, x_[i][k]);

						}
					}
					cplex.addLe(expr3, 1);
					//expr3.addTerm(-1.0, y[k]);

				}

				for(int k=1; k<nTargets-1; k++)
				{

					IloLinearNumExpr expr8 = cplex.linearNumExpr();
					for(int j=1; j<nTargets; j++)
					{
						if(j!=k)
						{

							expr8.addTerm(1.0, x_[k][j]);

						}
					}
					cplex.addLe(expr8, 1);
					//expr8.addTerm(-1.0, y[k]);

				}





				//connectivity

				for(int k=1; k<nTargets-1; k++)
				{

					IloLinearNumExpr expr3 = cplex.linearNumExpr();
					for(int i=0; i<nTargets-1; i++)
					{
						if(i!=k)
							expr3.addTerm(1.0, x_[i][k]);
					}
					//expr3.addTerm(-1.0, y[k]);

					for(int j=1; j<nTargets; j++)
					{
						if(k!=j)
							expr3.addTerm(-1.0, x_[k][j]);
					}



					cplex.addEq(expr3, 0);


				}





				IloLinearNumExpr expr14 = cplex.linearNumExpr();
				for(int j=0; j<nTargets-1; j++)
				{

					for(int k=1; k<nTargets; k++)
					{
						if(k!=j)
							expr14.addTerm(d[j][k], x_[j][k]);
					}

				}
				cplex.addLe(expr14, dmax);




				//subtour constraints

				for(int j=0; j<nTargets; j++)
				{

					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					expr5.addTerm(1, u[j]);
					cplex.addGe(expr5, 2);  // reason lol


				}


				for(int j=0; j<nTargets; j++)	
				{
					IloLinearNumExpr expr6 = cplex.linearNumExpr();
					expr6.addTerm(1, u[j]);
					cplex.addLe(expr6, nTargets);

				}




				for(int j=0; j<nTargets; j++)
				{

					for(int k=0; k<nTargets; k++)
					{

						if(j!=k)
						{

							IloLinearNumExpr expr7 = cplex.linearNumExpr();
							expr7.addTerm(1, u[j]);
							expr7.addTerm(-1, u[k]);
							expr7.addTerm(nTargets-1, x_[j][k]);
							cplex.addLe(expr7, nTargets-2);
						}


					}

				}

				cplex.solve();
				System.out.println("obj: "+ cplex.getObjValue());

				System.out.println();
				System.out.println();
				System.out.println();







				double sum =0;

				for(int i=0; i<nTargets-1; i++)
				{

					for(int j=1; j<nTargets; j++)
					{

						if(i!=j)
							if(i!=j && cplex.getValue(x_[i][j])>0)
							{


								int tmp= (int)Math.round(cplex.getValue(x_[i][j]));

								System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x_[i][j])/* + "  "+d[i][j] +" x " +tmp*/+" \n");

								sum += d[i][j]*tmp;
							}
						//attackedtarget=mapback.get(i);
						//break;
					}

				}


				ArrayList<Integer> pathtoconsider = new ArrayList<Integer>();
				int current = 0; // basenode
				pathtoconsider.add(current);
				int flag = 0;

				while(flag!=1)
				{
					// find edge that starts with current and end with a different node y. 
					/*
					 * if current is not base and y is base then end the loop
					 */

					int yid = findPath(cplex, x_, current, nTargets, copynodeid, mapback);

					if(current != 0 && yid == 0)
					{
						flag = 1;
					}
					if(yid==0 && current==0)
					{
						break;
					}

					current = yid;
					if(!pathtoconsider.contains(yid))
					{
						pathtoconsider.add(yid);
					}



				}
				if(flag==1)
				{
					pathtoconsider.add(0);
					res.add(pathtoconsider);
				}

				// modifiy the utility


				for(int i=0; i<nTargets-1; i++)
				{

					for(int j=1; j<nTargets; j++)
					{

						if(i!=j)
							if(i!=j && cplex.getValue(x_[i][j])>0)
							{


								s[i] = 0;
							}
						//attackedtarget=mapback.get(i);
						//break;
					}

				}

			}

			//end loop
			return res;

		}
		catch(Exception ex)
		{

		}
		return null;


	}



	public static ArrayList<ArrayList<Integer>> originalOP(int nPath, int[][] gamedata,
			ArrayList<TargetNode> targets, int nResource, int nOrigtargets, double dmax)
	{

		/**
		 * make a copy node same as the base node,
		 * add the node to the neighbors of the base node
		 * After finishing the work remove the copy node
		 * 
		 * 
		 * OR add a copy node information in the matrix
		 */


		/**
		 * create the copy node
		 *//*
		TargetNode copynode = new TargetNode(nOrigtargets, targets.get(0).getAnimaldensity());


		  *//**
		  * add the node to the neighbors of the base node
		  * add distances and path
		  *//*



		for(TargetNode neighbor: basenode.getNeighbors())
		{
			// add the copy node to neighbor's neighbor
			neighbor.addNeighbor(copynode);
			// add distance
			neighbor.addDistance(copynode, basenode.getDistance(neighbor));
			// add path


		}*/

		int copynodeid = nOrigtargets;
		int nTargets = targets.size()+1;
		TargetNode basenode = targets.get(0);
		int [][] d= new int[nTargets][nTargets];
		int [] s = new int[nTargets];

		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		// add the copy node
		mapback.put(icount, copynodeid);
		map.put(copynodeid, icount);

		System.out.println("#origtargets "+ nOrigtargets);
		System.out.println("#targets(including copy) "+ nTargets);
		System.out.println("#targets Mapping");
		for(int i=0; i<nTargets; i++)
		{
			System.out.println(i+"<--"+mapback.get(i));
		}



		// distance variable d[i,j]
		for(int i=0; i<nTargets-1; i++)
		{
			for(int j=0; j<nTargets-1; j++)
			{

				TargetNode inode = targets.get(i);
				TargetNode jnode = targets.get(j);

				if(inode.getNeighbors().contains(jnode))
				{
					d[i][j] = inode.getDistance(jnode).intValue();
				}
				else
				{
					d[i][j] = Integer.MAX_VALUE;
				}
			}

		}
		/*
		 * now add the copy node info to d
		 * for all the neighbor nodes of basenode add distance to copy node
		 */

		// add distance for neighbor nodes

		for(TargetNode nei: basenode.getNeighbors())
		{
			int id = map.get(nei.getTargetid());
			d[id][d.length-1] = basenode.getDistance(nei).intValue();
			d[d.length-1][id] = d[id][d.length-1];
		}
		for(int i=0; i<nTargets; i++)
		{
			if(!(d[i][d.length-1]>0) || !(d[d.length-1][i]>0))
			{
				d[i][d.length-1] = Integer.MAX_VALUE;
				d[d.length-1][i] = Integer.MAX_VALUE;

			}
		}
		for(int i=0; i<nTargets-1; i++)
		{
			s[i] = (int)Math.floor(targets.get(i).getAnimaldensity()) ;
		}
		//copy node
		s[nTargets-1] = (int)Math.floor(basenode.getAnimaldensity());
		printArray(d);
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();



			IloNumVar[][] x_ = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				x_[i] = cplex.boolVarArray(nTargets);
			}


			/**
			 * uip variable
			 */
			IloNumVar[] u = cplex.numVarArray(nTargets, 0, Integer.MAX_VALUE);


			// objective

			IloLinearNumExpr obj = cplex.linearNumExpr();

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{


					if(i!=j)
						obj.addTerm(s[i], x_[i][j]);

				}
			}

			cplex.addMaximize(obj);


			//constraint for start at base node 0 and end in copy node
			IloLinearNumExpr expr1 = cplex.linearNumExpr();

			for(int j=1; j<nTargets; j++)
			{
				expr1.addTerm(1.0, x_[0][j]);
			}
			cplex.addEq(expr1, 1);


			IloLinearNumExpr expr9 = cplex.linearNumExpr();

			for(int j=0; j<nTargets-1; j++)
			{
				expr9.addTerm(1.0, x_[j][nTargets-1]);
			}

			cplex.addEq(expr9, 1);





			IloLinearNumExpr expr12 = cplex.linearNumExpr();


			for(int j=1; j<nTargets; j++)
			{
				expr12.addTerm(1.0, x_[0][j]);
			}
			for(int j=0; j<nTargets-1; j++)
			{
				expr12.addTerm(-1.0, x_[j][nTargets-1]);
			}

			cplex.addEq(expr12, 0);




			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{

					if(i!=k)
					{

						expr3.addTerm(1.0, x_[i][k]);

					}
				}
				cplex.addLe(expr3, 1);
				//expr3.addTerm(-1.0, y[k]);

			}

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr8 = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					if(j!=k)
					{

						expr8.addTerm(1.0, x_[k][j]);

					}
				}
				cplex.addLe(expr8, 1);
				//expr8.addTerm(-1.0, y[k]);

			}





			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{

				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				for(int i=0; i<nTargets-1; i++)
				{
					if(i!=k)
						expr3.addTerm(1.0, x_[i][k]);
				}
				//expr3.addTerm(-1.0, y[k]);

				for(int j=1; j<nTargets; j++)
				{
					if(k!=j)
						expr3.addTerm(-1.0, x_[k][j]);
				}



				cplex.addEq(expr3, 0);


			}





			IloLinearNumExpr expr14 = cplex.linearNumExpr();
			for(int j=0; j<nTargets-1; j++)
			{

				for(int k=1; k<nTargets; k++)
				{
					if(k!=j)
						expr14.addTerm(d[j][k], x_[j][k]);
				}

			}
			cplex.addLe(expr14, dmax);




			//subtour constraints

			for(int j=0; j<nTargets; j++)
			{

				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				expr5.addTerm(1, u[j]);
				cplex.addGe(expr5, 2);  // reason lol


			}


			for(int j=0; j<nTargets; j++)	
			{
				IloLinearNumExpr expr6 = cplex.linearNumExpr();
				expr6.addTerm(1, u[j]);
				cplex.addLe(expr6, nTargets);

			}




			for(int j=0; j<nTargets; j++)
			{

				for(int k=0; k<nTargets; k++)
				{

					if(j!=k)
					{

						IloLinearNumExpr expr7 = cplex.linearNumExpr();
						expr7.addTerm(1, u[j]);
						expr7.addTerm(-1, u[k]);
						expr7.addTerm(nTargets-1, x_[j][k]);
						cplex.addLe(expr7, nTargets-2);
					}


				}

			}

			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());

			System.out.println();
			System.out.println();
			System.out.println();







			double sum =0;

			for(int i=0; i<nTargets-1; i++)
			{

				for(int j=1; j<nTargets; j++)
				{

					if(i!=j)
						if(i!=j && cplex.getValue(x_[i][j])>0)
						{


							int tmp= (int)Math.round(cplex.getValue(x_[i][j]));

							System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x_[i][j])/* + "  "+d[i][j] +" x " +tmp*/+" \n");

							sum += d[i][j]*tmp;
						}
					//attackedtarget=mapback.get(i);
					//break;
				}

			}


			ArrayList<Integer> pathtoconsider = new ArrayList<Integer>();
			ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();
			int current = 0; // basenode
			pathtoconsider.add(current);
			int flag = 0;

			while(flag!=1)
			{
				// find edge that starts with current and end with a different node y. 
				/*
				 * if current is not base and y is base then end the loop
				 */

				int yid = findPath(cplex, x_, current, nTargets, copynodeid, mapback);


				if(current != 0 && yid == 0)
				{
					flag = 1;
				}
				else if(current==0 && yid==0)
				{
					break;
				}

				current = yid;
				if(!pathtoconsider.contains(yid))
				{
					pathtoconsider.add(yid);
				}



			}
			if(flag==1)
			{
				pathtoconsider.add(0);
				res.add(pathtoconsider);
			}



			//////////






			return res;

		}
		catch(Exception ex)
		{

		}
		return null;


	}




	public static ArrayList<ArrayList<Integer>> originalTOP(int nPath, int[][] gamedata,
			ArrayList<TargetNode> targets, int nResource, int nOrigtargets, double dmax)
	{

		/**
		 * make a copy node same as the base node,
		 * add the node to the neighbors of the base node
		 * After finishing the work remove the copy node
		 * 
		 * 
		 * OR add a copy node information in the matrix
		 */


		/**
		 * create the copy node
		 *//*
		TargetNode copynode = new TargetNode(nOrigtargets, targets.get(0).getAnimaldensity());


		  *//**
		  * add the node to the neighbors of the base node
		  * add distances and path
		  *//*



		for(TargetNode neighbor: basenode.getNeighbors())
		{
			// add the copy node to neighbor's neighbor
			neighbor.addNeighbor(copynode);
			// add distance
			neighbor.addDistance(copynode, basenode.getDistance(neighbor));
			// add path


		}*/

		int copynodeid = nOrigtargets;
		int nTargets = targets.size()+1;
		TargetNode basenode = targets.get(0);
		int [][] d= new int[nTargets][nTargets];
		int [] s = new int[nTargets];
		//int nPath = nResource;

		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
		for(int i=0; i<targets.size(); i++)
		{
			map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		// add the copy node
		mapback.put(icount, copynodeid);
		map.put(copynodeid, icount);

		System.out.println("#origtargets "+ nOrigtargets);
		System.out.println("#targets(including copy) "+ nTargets);
		System.out.println("#targets Mapping");
		for(int i=0; i<nTargets; i++)
		{
			System.out.println(i+"<--"+mapback.get(i));
		}



		// distance variable d[i,j]
		for(int i=0; i<nTargets-1; i++)
		{
			for(int j=0; j<nTargets-1; j++)
			{

				TargetNode inode = targets.get(i);
				TargetNode jnode = targets.get(j);

				if(inode.getNeighbors().contains(jnode))
				{
					d[i][j] = inode.getDistance(jnode).intValue();
				}
				else
				{
					d[i][j] = Integer.MAX_VALUE;
				}
			}

		}
		/*
		 * now add the copy node info to d
		 * for all the neighbor nodes of basenode add distance to copy node
		 */

		// add distance for neighbor nodes

		for(TargetNode nei: basenode.getNeighbors())
		{
			int id = map.get(nei.getTargetid());
			d[id][d.length-1] = basenode.getDistance(nei).intValue();
			d[d.length-1][id] = d[id][d.length-1];
		}
		for(int i=0; i<nTargets; i++)
		{
			if(!(d[i][d.length-1]>0) || !(d[d.length-1][i]>0))
			{
				d[i][d.length-1] = Integer.MAX_VALUE;
				d[d.length-1][i] = Integer.MAX_VALUE;

			}
		}
		for(int i=0; i<nTargets-1; i++)
		{
			s[i] = (int)Math.floor(targets.get(i).getAnimaldensity()) ;
		}
		//copy node
		s[nTargets-1] = (int)Math.floor(basenode.getAnimaldensity());
		printArray(d);
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();



			IloNumVar[][][] x = new IloNumVar[nTargets][nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				x[i] = new IloNumVar[nTargets][nPath];
			}

			for(int i=0; i<nTargets; i++)
			{
				for(int j=0; j<nTargets; j++)
					x[i][j] = cplex.boolVarArray(nPath); 
			}

			IloNumVar[][] y = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				y[i] = cplex.boolVarArray(nPath);
			}


			/**
			 * uip variable
			 */
			IloNumVar[][] u = new IloNumVar[nTargets][];
			for(int i=0; i<nTargets; i++)
			{
				u[i] = cplex.numVarArray(nPath, 0, Integer.MAX_VALUE);
			}




			// objective

			IloLinearNumExpr obj = cplex.linearNumExpr();

			for(int p=0; p<nPath; p++)
			{

				for(int j=1; j<nTargets-1; j++)
				{


					if(p!=j)
						obj.addTerm(s[j], y[j][p]);

				}
			}

			cplex.addMaximize(obj);


			//constraint for start at base node 0 and end in copy node
			IloLinearNumExpr expr1 = cplex.linearNumExpr();

			for(int p=0; p<nPath; p++)
			{
				for(int j=1; j<nTargets; j++)
				{
					expr1.addTerm(1.0, x[0][j][p]);
				}
			}
			cplex.addEq(expr1, nPath);


			IloLinearNumExpr expr9 = cplex.linearNumExpr();

			for(int p=0; p<nPath; p++)
			{
				for(int j=0; j<nTargets-1; j++)
				{
					expr9.addTerm(1.0, x[j][nTargets-1][p]);
				}
			}

			cplex.addEq(expr9, nPath);





			/////



			for(int p=0; p<nPath; p++)
			{
				IloLinearNumExpr exp = cplex.linearNumExpr();
				for(int j=1; j<nTargets; j++)
				{
					exp.addTerm(1.0, x[0][j][p]);
				}
				cplex.addEq(exp, 1);
			}





			for(int p=0; p<nPath; p++)
			{
				IloLinearNumExpr exp9 = cplex.linearNumExpr();
				for(int j=0; j<nTargets-1; j++)
				{
					exp9.addTerm(1.0, x[j][nTargets-1][p]);
				}
				cplex.addEq(exp9, 1);
			}




			////







			IloLinearNumExpr expr12 = cplex.linearNumExpr();

			for(int p=0; p<nPath; p++)
			{
				for(int j=1; j<nTargets; j++)
				{
					expr12.addTerm(1.0, x[0][j][p]);
				}
			}

			for(int p=0; p<nPath; p++)
			{
				for(int j=0; j<nTargets-1; j++)
				{
					expr12.addTerm(-1.0, x[j][nTargets-1][p]);
				}
			}

			cplex.addEq(expr12, 0);





			for(int i=1;i<nTargets-1; i++)
			{
				IloLinearNumExpr exprr = cplex.linearNumExpr();
				for(int p=0; p<nPath; p++)
				{
					exprr.addTerm(1, y[i][p]);
				}
				cplex.addLe(exprr, 1);

			}




			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{
				for(int p=0; p<nPath; p++)
				{
					IloLinearNumExpr expr3 = cplex.linearNumExpr();
					for(int i=0; i<nTargets-1; i++)
					{

						if(i!=k)
						{
							expr3.addTerm(1.0, x[i][k][p]);

						}
					}
					expr3.addTerm(-1, y[k][p]);
					cplex.addEq(expr3, 0);
				}


			}



			for(int k=1; k<nTargets-1; k++)
			{


				for(int p=0; p<nPath; p++)
				{
					IloLinearNumExpr expr8 = cplex.linearNumExpr();
					for(int j=1; j<nTargets; j++)
					{
						if(j!=k)
						{

							expr8.addTerm(1.0, x[k][j][p]);

						}
					}
					expr8.addTerm(-1, y[k][p]);
					cplex.addEq(expr8, 0);
				}
				//expr8.addTerm(-1.0, y[k]);

			}





			//connectivity

			for(int k=1; k<nTargets-1; k++)
			{

				for(int p=0; p<nPath; p++)
				{
					IloLinearNumExpr expr3 = cplex.linearNumExpr();
					for(int i=0; i<nTargets-1; i++)
					{
						if(i!=k)
							expr3.addTerm(1.0, x[i][k][p]);
					}
					//expr3.addTerm(-1.0, y[k]);

					for(int j=1; j<nTargets; j++)
					{
						if(k!=j)
							expr3.addTerm(-1.0, x[k][j][p]);
					}
					cplex.addEq(expr3, 0);
				}


			}







			for(int p=0; p<nPath; p++)
			{
				IloLinearNumExpr expr14 = cplex.linearNumExpr();
				for(int j=0; j<nTargets-1; j++)
				{

					for(int k=1; k<nTargets; k++)
					{
						if(k!=j)
							expr14.addTerm(d[j][k], x[j][k][p]);
					}

				}
				cplex.addLe(expr14, dmax);
			}




			//subtour constraints

			for(int j=0; j<nTargets; j++)
			{
				for(int p=0; p<nPath; p++)
				{

					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					expr5.addTerm(1, u[j][p]);
					cplex.addGe(expr5, 1);  // reason lol
				}


			}


			for(int j=1; j<nTargets; j++)
			{
				for(int p=1; p<nPath; p++)	
				{
					IloLinearNumExpr expr6 = cplex.linearNumExpr();
					expr6.addTerm(1, u[j][p]);
					cplex.addLe(expr6, nTargets);

				}
			}




			for(int j=1; j<nTargets; j++)
			{

				for(int k=1; k<nTargets; k++)
				{


					for(int p=0; p<nPath; p++)
					{

						if(j!=k)
						{

							IloLinearNumExpr expr7 = cplex.linearNumExpr();
							expr7.addTerm(1, u[j][p]);
							expr7.addTerm(-1, u[k][p]);
							expr7.addTerm(nTargets-1, x[j][k][p]);
							cplex.addLe(expr7, nTargets-2);
						}
					}


				}

			}

			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());

			System.out.println();
			System.out.println();
			System.out.println();


			/*int[][] pp = new int[nTargets][nPath];

			for(int p=0; p<nPath; p++)
			{


				System.out.println("path "+ p);

				for(int i=1; i<nTargets; i++)
				{
					System.out.print(cplex.getValue(u[i][p])+"->");
					pp[(int)cplex.getValue(u[i][p])][p]=mapback.get(i);
				}
			}

			System.out.println();
			 */



			for(int p=0; p<nPath; p++)
			{

				double sum =0;
				System.out.println("path "+ p);

				for(int i=0; i<nTargets-1; i++)
				{

					for(int j=1; j<nTargets; j++)
					{

						if(i!=j)
							if(i!=j && cplex.getValue(x[i][j][p])>0)
							{
								int tmp= (int)Math.round(cplex.getValue(x[i][j][p]));

								System.out.print(mapback.get(i)+"-->"+mapback.get(j) +":"+cplex.getValue(x[i][j][p]) + "  "+d[i][j] +" x " +tmp+" \n");

								sum += d[i][j]*tmp;
							}

					}

				}
				System.out.println();
			}


			ArrayList<ArrayList<Integer>> path = new ArrayList<ArrayList<Integer>>();
			for(int p=0; p<nPath; p++)
			{
				ArrayList<Integer> pathtoconsider = new ArrayList<Integer>();
				int current = 0; // basenode
				pathtoconsider.add(current);
				int flag = 0;

				while(flag!=1)
				{
					// find edge that starts with current and end with a different node y. 
					/*
					 * if current is not base and y is base then end the loop
					 */

					System.out.println("Searching...current "+ current);
					int yid = findPath(cplex, x, current, nTargets, copynodeid,  p,mapback);

					System.out.println("Searching..done. yid "+ yid);


					if(current != 0 && yid == 0)
					{
						flag = 1;
					}
					else if(current==0 && yid==0)
					{
						break;
					}

					current = yid;
					if(!pathtoconsider.contains(yid))
					{
						pathtoconsider.add(yid);
					}



				}
				if(flag==1)
				{
					pathtoconsider.add(0);
					path.add(pathtoconsider);
				}

			}

			//////////



			//ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();

			//res.add(pathtoconsider);		
			return path;

		}
		catch(Exception ex)
		{

		}
		return null;


	}





	private static int findPath(IloCplex cplex, IloNumVar[][][] x, int targetid, 
			int nTargets, int copynodeid, int path, HashMap<Integer,Integer> mapback) throws UnknownObjectException, IloException {

		int n= nTargets;

		for(int i=0; i<n-1; i++)
		{
			for(int j=1; j<n; j++)
			{
				if(i!=j)
				{
					int in= mapback.get(i);
					int jn= mapback.get(j);
					//System.out.println(in+"...."+jn);

					if(cplex.getValue(x[i][j][path])==1 && in==targetid && jn!=targetid && jn!=copynodeid)
					{
						System.out.println(in+"..FOUND.."+jn);
						return jn;
					}
				}
			}
		}


		return 0;
	}

	private static int findPath(IloCplex cplex, IloNumVar[][] x, int targetid, 
			int nTargets, int copynodeid, HashMap<Integer,Integer> mapback) throws UnknownObjectException, IloException {

		int n= nTargets;

		for(int i=0; i<n-1; i++)
		{
			for(int j=1; j<n; j++)
			{
				if(i!=j)
				{
					int in= mapback.get(i);
					int jn= mapback.get(j);
					//System.out.println(in+"...."+jn);

					if(cplex.getValue(x[i][j])==1 && in==targetid && jn!=targetid && jn!=copynodeid)
					{
						//System.out.println(in+"..FOUND.."+jn);
						return jn;
					}
				}
			}
		}


		return 0;
	}




	/*private static int findPath(IloCplex cplex, IloNumVar[][] x, int targetid, 
			ArrayList<TargetNode> duplicatetargets, int nTargets) throws UnknownObjectException, IloException {

		int n= duplicatetargets.size();

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

						return jn.getTargetid();
					}
				}
			}
		}


		return 0;
	}
	 */



	private static void printArray(int[][] d) {


		for(int i=0; i<d.length; i++)
		{
			for(int j=0; j<d[0].length; j++)
			{
				System.out.print(d[i][j]+" ");
			}
			System.out.println();
		}

	}






	public static double[] solveForAttacker(int[][] p, int[][] gamedata, ArrayList<TargetNode> targets, int nResource)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = targets.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		for(int target =0; target<nTargets; target++)
		{
			int targetid= mapback.get(target);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[target][target] = gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[target][target] = gamedata[targetid][3] - gamedata[targetid][2];
		}
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(int target=0; target<nTargets; target++)
		{
			int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[target] = gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[target] = gamedata[targetid][2];

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
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
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}
			/**
			 * constraint 5
			 */

			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)
			{

				for(int target=0; target<nTargets; target++)
				{
					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
					{
						expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
					}
					expr5.addTerm(-1, kk);
					cplex.addLe(expr5, -Ua[target]);
				}
			}
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
				if(cplex.getValue(x[i])>0)
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
					attackedtarget=mapback.get(i);
					break;
				}

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}
	




	public static double[] solveForSuperAttackerNewMILP(int[][] p, HashMap<Integer, SuperTarget> sts, HashMap<String, double[]>  supertable, int nResource, 
			HashMap<Integer,Double> attackerstrategy, HashMap<String,Integer> subgamemap,
			HashMap<Integer,String> subgamemapback, HashMap<String,Integer> subgamecontainer, HashMap<Integer,Double> superuncovattckr)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = supertable.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		/*int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
*/

		/*for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}*/
		
		/*for(String starget : supertable.keySet())
		{
			String [] x = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}
			
			int targetindex= subgamemap.get(starget);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[targetindex][targetindex] = (int)supertable.get(starget)[0]0 - (-(int)supertable.get(starget)[1]);  //gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[targetindex][targetindex] = -(int)supertable.get(starget)[0]0 -((int)supertable.get(starget)[1]); //gamedata[targetid][2];
		}*/
		
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		//double[] Ud = new double[nTargets];
		//double[] Ua = new double[nTargets];

		/*for(String starget : supertable.keySet())
		{
			String [] x = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}
			
			int targetindex= subgamemap.get(starget);
			
			//int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[targetindex] = -(int)supertable.get(starget)[1];//gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[targetindex] = (int)supertable.get(starget)[1]; //gamedata[targetid][2];

		}*/
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			IloNumVar[] a = new IloNumVar[sts.size()];
			IloNumVar[] c = new IloNumVar[subgamecontainer.size()];
			for(int i=0; i<sts.size(); i++)
			{
				a[i] = cplex.boolVar();
			}
			
			for(int i=0; i<subgamecontainer.size(); i++)
			{
				c[i] = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			}
			
			/**
			 * objective
			 */
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1,  dd);
			cplex.addMaximize(obj);
			/**
			 * constraint 1
			 */
			
			// for every supertarget
			for(int target=0; target<sts.size(); target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);
				// for every subgame of a supertarget

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					//expr3.addTerm(-1.0*gamedata[target][0]*p[target][jointpath], x[jointpath]); reward
					//expr3.addTerm(-1.0*gamedata[target][1]*(1-p[target][jointpath]), x[jointpath]); penalty
					//expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);

					int maxsubgameid = -1;
					String maxsubgame = "none";
					double maxu = Double.NEGATIVE_INFINITY;
					for(String subgame: subgamecontainer.keySet())
					{

						int subgameid = subgamemap.get(subgame);
						if(subgamecontainer.get(subgame)==target && p[subgameid][jointpath]==1) // if that's the super target where subgame belongs
						{
							if((-(int)supertable.get(subgame)[1]) > maxu)
							{
								maxsubgameid = subgameid;
								maxsubgame = subgame;
								maxu = (-(int)supertable.get(subgame)[1]);

							}

						}
					}


					expr3.addTerm(-1.0* (0 /*-(int)supertable.get(subgame)[1]*/)  * p[maxsubgameid][jointpath], x[jointpath]); //reward
					expr3.addTerm(-1.0* (-(int)supertable.get(maxsubgame)[1])  * (p[maxsubgameid][jointpath]), x[jointpath]); //penalty

				}

				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M);
			}

			
			
			for(int target=0; target<sts.size(); target++)
			{
				
				//expr3.addTerm(1, dd);
				// for every subgame of a supertarget
				for(String subgame: subgamecontainer.keySet())
				{
					
					if(subgamecontainer.get(subgame)==target) // if that's the super target where subgame belongs
					{
						IloLinearNumExpr expr3 = cplex.linearNumExpr();
						int subgameid = subgamemap.get(subgame);
						expr3.addTerm(1, c[subgameid]);
						for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
						{
							
							//expr3.addTerm(-1.0* (0 /*-(int)supertable.get(subgame)[1]*/)  * p[subgameid][jointpath], x[jointpath]); //reward
							expr3.addTerm(-1.0* /*((int)supertable.get(subgame)[1])*/1*(p[subgameid][jointpath]), x[jointpath]); //penalty

						}
						cplex.addEq(expr3, 0);
					}
				}
				//expr3.addTerm(M, a[target]);
				
			}



			
			
			
			/**
			 * constraint 2
			 */

			// for every supertarget
			for(int target=0; target<sts.size(); target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, kk);
				// for every subgame of a supertarget
				
						for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
						{
							//expr3.addTerm(-1.0*gamedata[target][0]*p[target][jointpath], x[jointpath]); reward
							//expr3.addTerm(-1.0*gamedata[target][1]*(1-p[target][jointpath]), x[jointpath]); penalty
							//expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);
							
							int maxsubgameid = -1;
							String maxsubgame = "none";
							double maxu = Double.NEGATIVE_INFINITY;
							for(String subgame: subgamecontainer.keySet())
							{

								int subgameid = subgamemap.get(subgame);
								if(subgamecontainer.get(subgame)==target && p[subgameid][jointpath]==1) // if that's the super target where subgame belongs
								{
									if(((int)supertable.get(subgame)[1]) > maxu)
									{
										maxsubgameid = subgameid;
										maxsubgame = subgame;
										maxu = ((int)supertable.get(subgame)[1]);

									}

								}
							}

						
							expr3.addTerm(-1.0* (0 /*-(int)supertable.get(subgame)[1]*/)  * p[maxsubgameid][jointpath], x[jointpath]); //penalty
							expr3.addTerm(-1.0* ((int)supertable.get(maxsubgame)[1])  * (p[maxsubgameid][jointpath]), x[jointpath]); //reward

						}
					
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M);
			}


			
			
			
			

			/**
			 * constraint 3
			 */
			/*//for(int res=0; res<nResource; res++)
			{

				for(int target=0; target<nTargets; target++)
				{
					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
					{
						expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
					}
					expr5.addTerm(-1, kk);
					cplex.addLe(expr5, -Ua[target]);
				}
			}
			*/
			
			
			/**
			 * constraint 3
			 */
			// for every supertarget
			for(int target=0; target<sts.size(); target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				//expr3.addTerm(1, kk);
				// for every subgame of a supertarget
				for(String subgame: subgamecontainer.keySet())
				{

					if(subgamecontainer.get(subgame)==target) // if that's the super target where subgame belongs
					{
						for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
						{
							//expr3.addTerm(-1.0*gamedata[target][0]*p[target][jointpath], x[jointpath]); reward
							//expr3.addTerm(-1.0*gamedata[target][1]*(1-p[target][jointpath]), x[jointpath]); penalty
							//expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);

							int subgameid = subgamemap.get(subgame);
							expr3.addTerm(1.0* (0 /*-(int)supertable.get(subgame)[1]*/)  * p[subgameid][jointpath], x[jointpath]); //penalty
							expr3.addTerm(1.0* ((int)supertable.get(subgame)[1])  * (p[subgameid][jointpath]), x[jointpath]); //reward

						}
					}
				}
				expr3.addTerm(-1, kk);
				cplex.addLe(expr3,/*-superuncovattckr.get(target)*/0);
			}

			
			
			
			
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
			for(int j=0; j<sts.size(); j++)
			{

				expra.addTerm(1.0, a[j]);
			}

			cplex.addEq(expra, 1.0);
			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());
			double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			
			for(int i: subgamecontainer.values())
			{

				if(cplex.getValue(a[i])>0)
				{

					System.out.print("a["+i+"]="+cplex.getValue(a[i])+ " \n");
					//attackedtarget=subgamemapback.get(i);
					break;
				}

			}
			
			
			for(int target=0; target<sts.size(); target++)
			{
				System.out.println("\nST "+target);
				for(String subgame: subgamecontainer.keySet())
				{
					
					if(subgamecontainer.get(subgame)==target) // if that's the super target where subgame belongs
					{
						int subgameid = subgamemap.get(subgame);
						System.out.println("subgame  "+subgame + ", c "+cplex.getValue(c[subgameid]));
							
						
					}
				}
				
				
			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}
	
	
	public static double[] solveForSuperAttacker(int[][] p, HashMap<Integer, SuperTarget> sts, HashMap<String, double[]>  supertable, int nResource, 
			HashMap<Integer,Double> attackerstrategy, HashMap<String,Integer> subgamemap, HashMap<Integer,String> subgamemapback)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = supertable.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		/*int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
*/

		/*for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}*/
		
		for(String starget : supertable.keySet())
		{
			
			
			
			
			/*String [] x = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}*/
			
			int targetindex= subgamemap.get(starget);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[targetindex][targetindex] = /*(int)supertable.get(starget)[0]*/0 - (-(int)supertable.get(starget)[1]);  //gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[targetindex][targetindex] = /*-(int)supertable.get(starget)[0]*/0 -((int)supertable.get(starget)[1]); //gamedata[targetid][2];
		}
		
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(String starget : supertable.keySet())
		{
			/*String [] x = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}*/
			
			int targetindex= subgamemap.get(starget);
			
			//int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[targetindex] = -(int)supertable.get(starget)[1];//gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[targetindex] = (int)supertable.get(starget)[1]; //gamedata[targetid][2];

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
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
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}
			/**
			 * constraint 5
			 */

			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)
			{

				for(int target=0; target<nTargets; target++)
				{
					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
					{
						expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
					}
					expr5.addTerm(-1, kk);
					cplex.addLe(expr5, -Ua[target]);
				}
			}
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
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			System.out.println();
			System.out.println();
			System.out.println();
			for(int i: subgamemapback.keySet())
			{

				if(cplex.getValue(a[i])>0)
				{

					//System.out.print("a["+subgamemapback.get(i)+"]="+cplex.getValue(a[i])+ " \n");
					//attackedtarget=subgamemapback.get(i);
					break;
				}

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}
	
	
	
	
	
	public static double[] solveForSuperAttackerForOneSubgame(int[][] p, HashMap<Integer, SuperTarget> sts,
			HashMap<String, double[]>  supertable, int nResource, 
			HashMap<Integer,Double> attackerstrategy, HashMap<String,Integer> subgamemap,
			HashMap<Integer,String> subgamemapback, HashMap<String,Integer> subgamecontainer)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = supertable.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		/*int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();
*/

		/*for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}*/
		
		for(String starget : supertable.keySet())
		{
			
			
			
			
			/*String [] x = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}*/
			
			int targetindex= subgamemap.get(starget);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[targetindex][targetindex] = /*(int)supertable.get(starget)[0]*/0 - (-(int)supertable.get(starget)[1]);  //gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[targetindex][targetindex] = /*-(int)supertable.get(starget)[0]*/0 -((int)supertable.get(starget)[1]); //gamedata[targetid][2];
		}
		
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(String starget : supertable.keySet())
		{
			/*String [] x = starget.split(",");
			StringBuilder stringid = new StringBuilder();
			for(String y: x)
			{
				stringid.append(y);
			}*/
			
			int targetindex= subgamemap.get(starget);
			
			//int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[targetindex] = -(int)supertable.get(starget)[1];//gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[targetindex] = (int)supertable.get(starget)[1]; //gamedata[targetid][2];

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			IloNumVar dd1 = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			IloNumVar[] a = new IloNumVar[nTargets];
			for(int i=0; i<nTargets; i++)
			{
				a[i] = cplex.boolVar();
			}
			
			IloNumVar[] w = new IloNumVar[nTargets];
			for(int i=0; i<nTargets; i++)
			{
				w[i] = cplex.boolVar();
			}
			
			/**
			 * for every super target 
			 *  get the subgames
			 *  then for every subgames apply the constraints
			 */
			
			
			
			
			
			/**
			 * objective
			 */
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1,  dd);
			cplex.addMaximize(obj);
			
			
			
			//IloQuadNumExpr obj = cplex.quadNumExpr();
			//obj.addTerm(1, dd, dd1);
			
			

			/*IloLinearNumExpr expr71 = cplex.linearNumExpr();
			//for(int target=0; target<nTargets; target++)
			{
				expr71.addTerm(1, dd1);
			}
			cplex.addEq(expr71, 1);
			*/
			
			
			
			for(Integer stid: sts.keySet())
			{
				// get the subgames
				ArrayList<String> sbgms = getSubgames(stid, subgamecontainer);
				// for every subgame 
				IloQuadNumExpr exp = cplex.quadNumExpr();
				for(String s: sbgms)
				{
					int subid = subgamemap.get(s);
					//IloLinearNumExpr expr = cplex.linearNumExpr();
					for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
					{
						exp.addTerm(p[subid][jointpath], x[jointpath], w[subid]);
						//exp.addTerm(p[subid][jointpath], x[jointpath]);


					}
					
				}
				
				cplex.addLe(exp, 1);
			}
			
			
			
			for(Integer stid: sts.keySet())
			{
				// get the subgames
				ArrayList<String> sbgms = getSubgames(stid, subgamecontainer);
				// for every subgame 
				IloLinearNumExpr expr = cplex.linearNumExpr();
				for(String s: sbgms)
				{
					int subid = subgamemap.get(s);
					expr.addTerm(1, w[subid]);
				}
				cplex.addEq(expr, 1);
			}
			
			
			
			
			IloLinearNumExpr expr7 = cplex.linearNumExpr();
			for(int target=0; target<nTargets; target++)
			{
				expr7.addTerm(1, w[target]);
			}
			cplex.addEq(expr7, sts.size());
			
			
			
			
			
			
			
			/**
			 * constraint 4
			 */
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}
			/**
			 * constraint 5
			 */

			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)
			{

				for(int target=0; target<nTargets; target++)
				{
					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
					{
						expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
					}
					expr5.addTerm(-1, kk);
					cplex.addLe(expr5, -Ua[target]);
				}
			}
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
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			System.out.println();
			System.out.println();
			System.out.println();
			for(int i: subgamemapback.keySet())
			{

				if(cplex.getValue(w[i])>0)
				{

					System.out.print("w["+subgamemapback.get(i)+"]="+cplex.getValue(w[i])+ " \n");
					//attackedtarget=subgamemapback.get(i);
					//break;
				}

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}
	
	
	

	private static ArrayList<String> getSubgames(Integer stid, HashMap<String, Integer> subgamecontainer) {
		ArrayList<String> sbgms = new ArrayList<String>();
		for(String s: subgamecontainer.keySet())
		{
			Integer tmpstid = subgamecontainer.get(s);
			if(tmpstid==stid)
			{
				sbgms.add(s);
			}
		}
		return sbgms;
	}







	public static double[] solveForAttackerLP(int[][] p, int[][] gamedata, ArrayList<TargetNode> targets, double nResource, 
			HashMap<Integer,Double> attackerstrategy)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = targets.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction mmmmm
		 */
		int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		for(int target =0; target<nTargets; target++)
		{
			int targetid= mapback.get(target);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[target][target] = gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[target][target] = gamedata[targetid][3] - gamedata[targetid][2];
		}
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(int target=0; target<nTargets; target++)
		{
			int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[target] = gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[target] = gamedata[targetid][2];

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			//IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] a = new IloNumVar[nTargets];
			/*for(int i=0; i<nTargets; i++)
			{
				a[i] = cplex.boolVar();
			}*/
			/**
			 * objective
			 */
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1,  kk);
			cplex.addMinimize(obj);
			/**
			 * constraint 4
			 */
			/*for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}*/
			/**
			 * constraint 5
			 */

			/*for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}*/

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)

			List<IloRange> constraints = new ArrayList<IloRange>();
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
				}
				expr5.addTerm(-1, kk);
				constraints.add(cplex.addLe(expr5, -Ua[target]));
			}

			/**
			 * constraint 7
			 */
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for(int j=0; j<nJointSchedule; j++)
			{

				expr.addTerm(1.0, x[j]);
			}
			cplex.addEq(expr, 1);

			/**
			 * constraint for attakced target
			 */
			/*IloLinearNumExpr expra = cplex.linearNumExpr();
			for(int j=0; j<nTargets; j++)
			{

				expra.addTerm(1.0, a[j]);
			}

			cplex.addEq(expra, 1.0);*/



			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());
			double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			System.out.println();
			System.out.println();
			System.out.println();
			/*for(int i=0; i<nTargets; i++)
			{

				if(cplex.getValue(a[i])>0)
				{

					System.out.print("a["+mapback.get(i)+"]="+cplex.getValue(a[i])+ " \n");
					 attackedtarget=mapback.get(i);
					 break;
				}


			}*/

			for(int i=0; i<nTargets; i++)
			{
				attackerstrategy.put(mapback.get(i), -(cplex.getDual(constraints.get(i))));
				System.out.println(" constraint "+ mapback.get(i) + " dual "+ -(cplex.getDual(constraints.get(i))));

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}
	
	


	public static double[] solveForAttackerLPST(int[][] p, HashMap<Integer, SuperTarget> sts, 
			HashMap<Integer, TargetNode> targetmaps, double nResource, 
			HashMap<Integer,Double> attackerstrategy)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = sts.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(Integer stid: sts.keySet())
		{

			map.put(stid, icount);
			mapback.put(icount, stid);
			icount++;

		}
		for(Integer stid: sts.keySet())
		{
			// find supertarget's maximum value target for attacker
			
			
			
			
			int targetid= map.get(stid);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			//D[targetid][targetid] = gamedata[targetid][0] - gamedata[targetid][1];
			D[targetid][targetid] = (int)sts.get(stid).defenderreward - (int)sts.get(stid).defenderpenalty;
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			//A[targetid][targetid] = gamedata[targetid][3] - gamedata[targetid][2];
			A[targetid][targetid] = (int)sts.get(stid).attackerpenalty - (int)sts.get(stid).attackerreward;
		}
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(Integer stid: sts.keySet())
		{
			int targetid = map.get(stid);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[targetid] = (int)sts.get(stid).defenderpenalty;
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[targetid] = (int)sts.get(stid).attackerreward;

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			//IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] a = new IloNumVar[nTargets];
			/*for(int i=0; i<nTargets; i++)
			{
				a[i] = cplex.boolVar();
			}*/
			/**
			 * objective
			 */
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1,  kk);
			cplex.addMinimize(obj);
			/**
			 * constraint 4
			 */
			/*for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}*/
			/**
			 * constraint 5
			 */

			/*for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}*/

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)

			List<IloRange> constraints = new ArrayList<IloRange>();
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
				}
				expr5.addTerm(-1, kk);
				constraints.add(cplex.addLe(expr5, -Ua[target]));
			}

			/**
			 * constraint 7
			 */
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for(int j=0; j<nJointSchedule; j++)
			{

				expr.addTerm(1.0, x[j]);
			}
			cplex.addEq(expr, 1);

			/**
			 * constraint for attakced target
			 */
			/*IloLinearNumExpr expra = cplex.linearNumExpr();
			for(int j=0; j<nTargets; j++)
			{

				expra.addTerm(1.0, a[j]);
			}

			cplex.addEq(expra, 1.0);*/



			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());
			double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			System.out.println();
			System.out.println();
			System.out.println();
			/*for(int i=0; i<nTargets; i++)
			{

				if(cplex.getValue(a[i])>0)
				{

					System.out.print("a["+mapback.get(i)+"]="+cplex.getValue(a[i])+ " \n");
					 attackedtarget=mapback.get(i);
					 break;
				}


			}*/

			for(int i=0; i<nTargets; i++)
			{
				attackerstrategy.put(mapback.get(i), -(cplex.getDual(constraints.get(i))));
				System.out.println(" constraint "+ mapback.get(i) + " dual "+ -(cplex.getDual(constraints.get(i))));

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}






	public static double[] solveForAttackerSuperLP(int[][] p, HashMap<Integer, SuperTarget> sts, HashMap<String, double[]>  supertable, int nResource, 
			HashMap<Integer,Double> attackerstrategy, HashMap<Integer,Integer> subgamemap, HashMap<Integer,Integer> subgamemapback)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = supertable.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		/*int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}*/
		for(String starget : supertable.keySet())
		{
			int target = Integer.parseInt(starget);
			
			int targetindex= subgamemap.get(target);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[targetindex][targetindex] = 0 - (-(int)supertable.get(starget)[1]);  //gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[targetindex][targetindex] = 0 -((int)supertable.get(starget)[1]); //gamedata[targetid][2];
		}
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(String starget : supertable.keySet())
		{
			int target = Integer.parseInt(starget);
			
			int targetindex= subgamemap.get(target);
			
			//int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[targetindex] = -(int)supertable.get(starget)[1];//gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[targetindex] = (int)supertable.get(starget)[1]; //gamedata[targetid][2];

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			//IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] a = new IloNumVar[nTargets];
			/*for(int i=0; i<nTargets; i++)
			{
				a[i] = cplex.boolVar();
			}*/
			/**
			 * objective
			 */
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1,  kk);
			cplex.addMinimize(obj);
			/**
			 * constraint 4
			 */
			/*for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}*/
			/**
			 * constraint 5
			 */

			/*for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}*/

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)

			List<IloRange> constraints = new ArrayList<IloRange>();
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr5 = cplex.linearNumExpr();
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
				}
				expr5.addTerm(-1, kk);
				constraints.add(cplex.addLe(expr5, -Ua[target]));
			}

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
			/*IloLinearNumExpr expra = cplex.linearNumExpr();
			for(int j=0; j<nTargets; j++)
			{

				expra.addTerm(1.0, a[j]);
			}

			cplex.addEq(expra, 1.0);*/



			cplex.solve();
			System.out.println("obj: "+ cplex.getObjValue());
			double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				if(cplex.getValue(x[i])>0)
				{

					System.out.print("x["+i+"]="+cplex.getValue(x[i])+ " ");
					result[i] = cplex.getValue(x[i]);
				}
			}
			System.out.println();
			System.out.println();
			System.out.println();
			/*for(int i=0; i<nTargets; i++)
			{

				if(cplex.getValue(a[i])>0)
				{

					System.out.print("a["+mapback.get(i)+"]="+cplex.getValue(a[i])+ " \n");
					 attackedtarget=mapback.get(i);
					 break;
				}


			}*/

			for(String starget : supertable.keySet())
			{
				int target = Integer.parseInt(starget);
				
				int targetindex= subgamemap.get(target);
				
				
				attackerstrategy.put(target, -(cplex.getDual(constraints.get(targetindex))));
				System.out.println(" constraint "+ target + " dual "+ -(cplex.getDual(constraints.get(targetindex))));

			}
			return result;

		}
		catch(Exception ex)
		{

		}
		return null;


	}
	
	
	
	
	


	public static double[] solve(int[][] p, int[][] gamedata, ArrayList<TargetNode> targets, int nResource)
	{
		/**
		 * make the D matrix and A matrix
		 */
		int nTargets = targets.size();
		int[][] D = new int[nTargets][nTargets];
		int[][] A = new int[nTargets][nTargets];
		int nJointSchedule= p[0].length;
		/**
		 * include targetid when using contraction
		 */
		int icount = 0;
		//HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		HashMap<Integer, Integer> mapback = new HashMap<Integer, Integer>();


		for(int i=0; i<targets.size(); i++)
		{

			//map.put(targets.get(i).getTargetid(), icount);
			mapback.put(icount, targets.get(i).getTargetid());
			icount++;

		}
		for(int target =0; target<nTargets; target++)
		{
			int targetid= mapback.get(target);
			//System.out.println("D["+target+"]["+target+"]=("+ gamedata[targetid][0] + ")-("+ gamedata[targetid][1]+")") ;
			D[target][target] = gamedata[targetid][0] - gamedata[targetid][1];
			//System.out.println("A["+target+"]["+target+"]=("+ gamedata[targetid][3] + ")-("+ gamedata[targetid][2]+")") ;
			A[target][target] = gamedata[targetid][3] - gamedata[targetid][2];
		}
		/**
		 * the U vectors
		 * Ud
		 * Ua
		 */
		System.out.println();
		double[] Ud = new double[nTargets];
		double[] Ua = new double[nTargets];

		for(int target=0; target<nTargets; target++)
		{
			int targetid = mapback.get(target);
			//System.out.println("Ud["+target+"] = "+gamedata[targetid][1]);
			Ud[target] = gamedata[targetid][1];
			//System.out.println("Ua["+target+"] = "+gamedata[targetid][2]);
			Ua[target] = gamedata[targetid][2];

		}
		int M =100000;

		/**
		 * cplex variables
		 */
		try
		{

			IloCplex cplex = new IloCplex();
			IloNumVar[] x = new IloNumVar[nJointSchedule];//cplex.numVarArray(nPath, 0, 1);
			for(int i=0; i<nJointSchedule ; i++)
			{
				x[i] = cplex.numVar(0, 1);
			}
			//IloNumVar[] d = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar dd = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
			//IloNumVar[] k = cplex.numVarArray(nTargets, 0, Double.MAX_VALUE);
			IloNumVar kk = cplex.numVar(Double.NEGATIVE_INFINITY, Double.MAX_VALUE);
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
			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1, dd);

				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr3.addTerm(-1.0*D[target][target]*p[target][jointpath], x[jointpath]);


				}
				expr3.addTerm(M, a[target]);
				cplex.addLe(expr3, M+Ud[target]);
			}
			/**
			 * constraint 5
			 */

			for(int target=0; target<nTargets; target++)
			{
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				expr4.addTerm(1, kk);
				for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
				{
					expr4.addTerm(-1.0*A[target][target]*p[target][jointpath], x[jointpath]);
				}
				//expr4.addTerm(-1.0, Ud[target]);
				expr4.addTerm(M, a[target]);
				cplex.addLe(expr4, M+Ua[target]);
			}

			/**
			 * constraint 6
			 */
			//for(int res=0; res<nResource; res++)
			{

				for(int target=0; target<nTargets; target++)
				{
					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					for(int jointpath=0; jointpath<nJointSchedule; jointpath++)
					{
						expr5.addTerm(A[target][target]*p[target][jointpath], x[jointpath]);
					}
					expr5.addTerm(-1, kk);
					cplex.addLe(expr5, -Ua[target]);
				}
			}
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
			System.out.println();
			System.out.println("Solution status="+ cplex.getStatus());

			System.out.println("obj: "+ cplex.getObjValue());
			double[] result = new double[nJointSchedule];
			for(int i=0; i<nJointSchedule; i++)
			{
				if(cplex.getValue(x[i])>0)
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
					attackedtarget=mapback.get(i);
				}

			}
			return result;

		}
		catch(Exception ex)
		{
			System.out.println(ex.getMessage().toString());

		}
		return null;


	}


}
