package cs.Interval.Coincollection;

import java.util.ArrayList;

import cs.Interval.contraction.SecurityGameContraction;
import cs.Interval.contraction.TargetNode;

public class CoinCollection 
{
	public static TargetNode coinCollection(ArrayList<TargetNode> targets, ArrayList<TargetNode> dominatednodes)
	{
		ArrayList<Integer> maxnodeids = new ArrayList<Integer>();
		//ArrayList<Integer> deadendnodeids = new ArrayList<Integer>();
		for(int nodeid = 0; nodeid<SecurityGameContraction.targets.size(); nodeid++)
		{
			TargetNode node = SecurityGameContraction.targets.get(nodeid); 

			if(!dominatednodes.contains(node))
			{
				//System.out.println("In node "+ nodeid);
				//if(node.getRow() >0 && node.getCol()>0)
				{
					double maxcoinval = Double.NEGATIVE_INFINITY;
					Integer maxnodeid = -1;
					if(node.getTargetid()==9)
					{
						//System.out.println();
					}

					for(TargetNode neighbor : SecurityGameContraction.targets.get(nodeid).getNeighbors())
					{
						if( node.getTargetid()>0 && (node.getTargetid()>neighbor.getTargetid())/*((node.getRow()==0 && neighbor.getRow()==0) && 
								((node.getCol()-1)==neighbor.getCol()))   
								||  ((node.getCol()==0 && neighbor.getCol()==0) && 
										((node.getRow()-1)==neighbor.getRow())) 
										|| (node.getRow() >0 && node.getCol()>0)*/)
						{
							double neighborcoincolected = 0;

							if(!maxnodeids.contains(neighbor.getTargetid()))
							{
								neighborcoincolected += neighbor.getCoinvalue();
								//System.out.println("neighbor coin value "+ neighbor.getCoinvalue());
								/**
								 * sum the utility along the path with neighborcoincolected
								 */
								neighborcoincolected += node.getPathUtility(neighbor);
								//System.out.println("path coin value "+ node.getPathUtility(neighbor)+"\n");
							}

							if((neighborcoincolected) > maxcoinval)
							{
								maxcoinval = neighborcoincolected;
								maxnodeid = neighbor.getTargetid();
							}
						}

					}
					if((maxcoinval>=0) && (maxnodeid != -1))
					{
						//System.out.println("Max parent id "+ maxnodeid + " , coinval " + maxcoinval);

						node.setCoinvalue(maxcoinval+node.getAnimaldensity());
						//System.out.println("Coin collected yet "+ node.getCoinvalue());
						node.setMaxCoinParent(SecurityGameContraction.targets.get(maxnodeid));
						maxnodeids.add(maxnodeid);
					}
					else
					{
						node.setCoinvalue(node.getAnimaldensity());
						node.setMaxCoinParent(null);
						//deadendnodeids.add(node.getTargetid());
						//System.out.println("No max parent");
					}


					if(node.isGoal())
					{

						return node;
					}

				}
			}
		}
		return null;
	}


	public static double printPath(TargetNode node)
	{

		/*if(node.getRow() ==0 && node.getCol()==0)
		{
			System.out.print(node.getTargetid() + " ");
			return;
		}
		printPath(node.getMaxCoinParent());
		System.out.print(node.getTargetid() + " ");
		if(node.getMaxCoinParent().getPath(node).size()>0)
		{
			for(TargetNode n: node.getMaxCoinParent().getPath(node))
			{
				System.out.print(n.getTargetid() + " ");
			}
		}*/

		TargetNode n = new TargetNode();
		n = node;
		double pathutility = 0;
		ArrayList<Integer> visited = new ArrayList<Integer>();

		while(n.getTargetid()>=0)
		{
			System.out.print(n.getTargetid() + " ");

			if(n.getRow() ==0 && n.getCol()==0)
			{
				if(!visited.contains(n.getTargetid()))
				{
					pathutility += n.getAnimaldensity();
				}

				//System.out.print(node.getTargetid() + " ");
				break;
			}
			//System.out.print(n.getTargetid() + "... ");
			//System.out.print(n.getMaxCoinParent().getTargetid() + "++ ");

			/*for(TargetNode p: n.getPath().keySet())
			{
				System.out.println(p.getTargetid() + "** ");
			}*/

			if(!(n.getMaxCoinParent()==null))
			{
				ArrayList<TargetNode> x = n.getPath(n.getMaxCoinParent());
				if(!visited.contains(n.getTargetid()))
				{
					pathutility += n.getAnimaldensity();
				}
				if(!x.equals(null))
				{
					for(TargetNode y: x)
					{
						if(!visited.contains(y.getTargetid()))
						{
							pathutility += y.getAnimaldensity();
						}
						visited.add(y.getTargetid());
						System.out.print(y.getTargetid() + " ");
					}
				}


				visited.add(n.getTargetid());
				n = n.getMaxCoinParent();
			}
			else
			{
				return pathutility;
			}
		}
		return pathutility;

	}

}
