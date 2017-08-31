package groupingtargets;

import java.awt.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import cs.Interval.contraction.TargetNode;

public class SuperTarget {


	public int stid;
	public String combid;
	public double animaldensity;
	public double defenderreward;
	public double defenderpenalty;
	public double attackerreward;
	public double attackerpenalty;
	public SuperTarget parent;
	public double distcoveredyet = 0;
	public double allocateddistance = 0;
	public double currentprob = 1;
	public int currentindex = 0;

	public HashMap<Integer, TargetNode> nodes = new HashMap<Integer, TargetNode>();
	public ArrayList<Integer> pathnodesid = new ArrayList<Integer>(); 
	public HashMap<Integer, TargetNode> ap = new HashMap<Integer, TargetNode>();
	public HashMap<Integer, TargetNode> entrypoints = new HashMap<Integer, TargetNode>();
	public HashMap<Integer, TargetNode> exitpoints = new HashMap<Integer, TargetNode>();
	public HashMap<Integer, SuperTarget> neighbors = new  HashMap<Integer, SuperTarget>();
	public HashMap<SuperTarget, Double> distances = new HashMap<SuperTarget, Double>();
	public HashMap<SuperTarget, ArrayList<Integer>> path = new HashMap<SuperTarget, ArrayList<Integer>>();


	public SuperTarget(SuperTarget superTarget) {

		this.animaldensity= superTarget.animaldensity;
		this.defenderreward = superTarget.defenderreward;
		this.defenderpenalty = superTarget.defenderpenalty;
		this.attackerreward = superTarget.attackerreward;
		this.attackerpenalty = superTarget.attackerpenalty;
		for(TargetNode t: superTarget.nodes.values())
		{
			this.nodes.put(t.getTargetid(), t);
		}


	}


	public SuperTarget(SuperTarget superTarget, TargetNode entrypoint, TargetNode exitpoint) {

		this.animaldensity= superTarget.animaldensity;
		this.defenderreward = superTarget.defenderreward;
		this.defenderpenalty = superTarget.defenderpenalty;
		this.attackerreward = superTarget.attackerreward;
		this.attackerpenalty = superTarget.attackerpenalty;
		for(TargetNode t: superTarget.nodes.values())
		{
			this.nodes.put(t.getTargetid(), t);
		}

		this.entrypoints.put(entrypoint.getTargetid(), entrypoint);
		this.exitpoints.put(exitpoint.getTargetid(), exitpoint);


	}


	public SuperTarget() {
		// TODO Auto-generated constructor stub
	}

	public void addAP(TargetNode ap)
	{
		this.ap.put(ap.getTargetid(), ap);
	}

	public void addEntryPoint(TargetNode ep)
	{
		this.entrypoints.put(ep.getTargetid(), ep);
	}

	public TargetNode getEntryPoint(Integer tid)
	{
		return this.entrypoints.get(tid);
	}

	public TargetNode getExitPoint(Integer tid)
	{
		return this.exitpoints.get(tid);
	}


	public void addNeighbor(SuperTarget st)
	{
		this.neighbors.put(st.stid, st);
	}





	public void addExitPoint(TargetNode ep)
	{
		this.exitpoints.put(ep.getTargetid(), ep);
	}
	
	
	public static HashMap<Integer, SuperTarget> buildSuperTargets(ArrayList<Integer>[] clusters,
			HashMap<Integer, TargetNode> targetmaps) {

		HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();

		for(int i=0; i<clusters.length; i++)
		{
			
			if(clusters[i].size()>0)
			{
				// create a st object
				SuperTarget st = new SuperTarget();
				// set id
				st.stid = i;
				// list all the nodes in a hashmap
				for(Integer t: clusters[i])
				{
					TargetNode tmp = targetmaps.get(t);
					st.nodes.put(t, tmp);
					//System.out.println("Inserting target "+ t + " in super target "+ i);
				}

				// list all the accesspoints : entry and exit points
				for(TargetNode node: st.nodes.values()) // key is the target id
				{
					// get the node, check if it has any neighbor outside of this super target
					boolean hasneioutside = hasNeiOutSide(node, st.nodes);
					if(hasneioutside)
					{
						//System.out.println("Inserting target "+ node.getTargetid() + " as accesspoint of supertarget "+ i);
						st.ap.put(node.getTargetid(), node);
					}
					else
					{
						//System.out.println("target "+ node.getTargetid() + " has no outside neighbor ");
					}
				}
				sts.put(st.stid, st);
			}
		}


		// now connect the neighbors of super targets

		int[] stdids = new int[sts.keySet().size()];


		int index=0;
		for(int i: sts.keySet())
		{
			stdids[index++] = i;
		}

		for(int i=0; i<stdids.length-1; i++)
		{
			for(int j=i+1; j<stdids.length; j++)
			{
				SuperTarget sti = sts.get(stdids[i]);
				SuperTarget stj = sts.get(stdids[j]);
				// now check if they are neighbor
				boolean isneighbor = isNeighbor(sti, stj);
				if(isneighbor)
				{
					//System.out.println("Super target "+ stdids[i] + " is neighbor of super target "+ stdids[j]);
					// now add the connections
					sti.neighbors.put(stj.stid, stj);
					stj.neighbors.put(sti.stid, sti);
					//break;
				}

			}
		}


		return sts;

	}


	public static HashMap<Integer, SuperTarget> buildSuperTargets2(ArrayList<Integer>[] clusters,
			HashMap<Integer, TargetNode> targetmaps) {

		HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();

		for(int i=0; i<clusters.length; i++)
		{
			// create a st object
			SuperTarget st = new SuperTarget();
			// set id
			st.stid = i;
			// list all the nodes in a hashmap
			for(Integer t: clusters[i])
			{
				TargetNode tmp = targetmaps.get(t);
				st.nodes.put(t, tmp);
				//System.out.println("Inserting target "+ t + " in super target "+ i);
			}

			// list all the accesspoints : entry and exit points
			int apcounter = 0;
			
			if(st.nodes.size()==1)
			{
				st.ap.put(st.nodes.get(0).getTargetid(), st.nodes.get(0));
			}
			else
			{

				HashMap<Integer, Integer> countoutnei = new HashMap<Integer, Integer>();
				for(TargetNode node: st.nodes.values()) // key is the target id
				{
					// get the node, check if it has any neighbor outside of this super target
					int countneioutside = countNeiOutSide(node, st.nodes);
					countoutnei.put(node.getTargetid(), countneioutside);
					// pick 2 ap who has most connectins outside

					/*if(countneioutside)
				{
					//System.out.println("Inserting target "+ node.getTargetid() + " as accesspoint of supertarget "+ i);
					st.ap.put(node.getTargetid(), node);
					apcounter++;
					if(apcounter==2)
					{
						break;
					}
				}
				else
				{
					//System.out.println("target "+ node.getTargetid() + " has no outside neighbor ");
				}*/
				}
				// find two most connected nodes as ap

				countoutnei = sortHashMapByValues(countoutnei);

				int count = 0;

				for(int j: countoutnei.keySet())
				{
					
					if(count==(countoutnei.size()-1))
					{
						st.ap.put(j, targetmaps.get(j));
					}
					else if(count==(countoutnei.size()-2))
					{
						st.ap.put(j, targetmaps.get(j));

					}
					count++;


				}
			}



			
			
			sts.put(st.stid, st);
		}


		// now connect the neighbors of super targets

		int[] stdids = new int[sts.keySet().size()];


		int index=0;
		for(int i: sts.keySet())
		{
			stdids[index++] = i;
		}

		for(int i=0; i<stdids.length-1; i++)
		{
			for(int j=i+1; j<stdids.length; j++)
			{
				SuperTarget sti = sts.get(stdids[i]);
				SuperTarget stj = sts.get(stdids[j]);
				// now check if they are neighbor
				boolean isneighbor = isNeighbor(sti, stj);
				if(isneighbor)
				{
					//System.out.println("Super target "+ stdids[i] + " is neighbor of super target "+ stdids[j]);
					// now add the connections
					sti.neighbors.put(stj.stid, stj);
					stj.neighbors.put(sti.stid, sti);
					//break;
				}

			}
		}


		return sts;

	}
	
	
	public static HashMap<Integer, SuperTarget> buildSuperTargets3(ArrayList<Integer>[] clusters,
			HashMap<Integer, TargetNode> targetmaps, HashMap<Integer,HashMap<Integer,Integer>> stcountoutnei) {

		HashMap<Integer, SuperTarget> sts = new HashMap<Integer, SuperTarget>();

		for(int i=0; i<clusters.length; i++)
		{
			// create a st object
			SuperTarget st = new SuperTarget();
			// set id
			st.stid = i;
			// list all the nodes in a hashmap
			for(Integer t: clusters[i])
			{
				TargetNode tmp = targetmaps.get(t);
				st.nodes.put(t, tmp);
				//System.out.println("Inserting target "+ t + " in super target "+ i);
			}

			// list all the accesspoints : entry and exit points
			int apcounter = 0;
			
			if(st.nodes.size()==1)
			{
				st.ap.put(st.nodes.get(0).getTargetid(), st.nodes.get(0));
			}
			else
			{

				HashMap<Integer, Integer> countoutnei = new HashMap<Integer, Integer>();
				for(TargetNode node: st.nodes.values()) // key is the target id
				{
					// get the node, check if it has any neighbor outside of this super target
					int countneioutside = countNeiOutSide(node, st.nodes);
					countoutnei.put(node.getTargetid(), countneioutside);
					// pick 2 ap who has most connectins outside

					/*if(countneioutside)
				{
					//System.out.println("Inserting target "+ node.getTargetid() + " as accesspoint of supertarget "+ i);
					st.ap.put(node.getTargetid(), node);
					apcounter++;
					if(apcounter==2)
					{
						break;
					}
				}
				else
				{
					//System.out.println("target "+ node.getTargetid() + " has no outside neighbor ");
				}*/
				}
				// find two most connected nodes as ap

				countoutnei = sortHashMapByValues(countoutnei);
				stcountoutnei.put(i, countoutnei);

				int count = 0;

				for(int j: countoutnei.keySet())
				{
					
					if(count==(countoutnei.size()-1))
					{
						st.ap.put(j, targetmaps.get(j));
					}
					else if(count==(countoutnei.size()-2))
					{
						st.ap.put(j, targetmaps.get(j));

					}
					count++;


				}
			}



			
			
			sts.put(st.stid, st);
		}


		// now connect the neighbors of super targets

		int[] stdids = new int[sts.keySet().size()];


		int index=0;
		for(int i: sts.keySet())
		{
			stdids[index++] = i;
		}

		for(int i=0; i<stdids.length-1; i++)
		{
			for(int j=i+1; j<stdids.length; j++)
			{
				SuperTarget sti = sts.get(stdids[i]);
				SuperTarget stj = sts.get(stdids[j]);
				// now check if they are neighbor
				boolean isneighbor = isNeighbor(sti, stj);
				if(isneighbor)
				{
					//System.out.println("Super target "+ stdids[i] + " is neighbor of super target "+ stdids[j]);
					// now add the connections
					sti.neighbors.put(stj.stid, stj);
					stj.neighbors.put(sti.stid, sti);
					//break;
				}

			}
		}


		return sts;

	}
	
	
	
	public static LinkedHashMap<Integer, Integer> sortHashMapByValues(
	        HashMap<Integer,Integer> countoutnei) {
	    ArrayList<Integer> mapKeys = new ArrayList<>(countoutnei.keySet());
	    ArrayList<Integer> mapValues = new ArrayList<>(countoutnei.values());
	    Collections.sort(mapValues);
	    Collections.sort(mapKeys);

	    LinkedHashMap<Integer, Integer> sortedMap =
	        new LinkedHashMap<>();

	    Iterator<Integer> valueIt = mapValues.iterator();
	    while (valueIt.hasNext()) {
	        Integer val = valueIt.next();
	        Iterator<Integer> keyIt = mapKeys.iterator();

	        while (keyIt.hasNext()) {
	            Integer key = keyIt.next();
	            Integer comp1 = countoutnei.get(key);
	            Integer comp2 = val;

	            if (comp1==comp2) {
	                keyIt.remove();
	                sortedMap.put(key, val);
	                break;
	            }
	        }
	    }
	    return sortedMap;
	}


	private static boolean isNeighbor(SuperTarget sti, SuperTarget stj) {

		for(Integer id: sti.nodes.keySet())
		{
			TargetNode n = sti.nodes.get(id);
			// for every node check if its neighbor is in stj
			for(TargetNode nei: n.getNeighbors())
			{
				if(stj.nodes.keySet().contains(nei.getTargetid()))
				{
					return true;
				}
			}
		}
		return false;
	}


	private static boolean hasNeiOutSide(TargetNode targetid, HashMap<Integer, TargetNode> nodes) {


		TargetNode t = nodes.get(targetid.getTargetid());
		for(TargetNode nei: t.getNeighbors())
		{
			// check if neibor is not in supertarget
			if(!nodes.containsKey(nei.getTargetid()))
			{
				// if neibor is not in supertarget that means it has neighbor outside
				return true;
			}
		}
		return false;
	}
	

	private static int countNeiOutSide(TargetNode targetid, HashMap<Integer, TargetNode> nodes) {


		TargetNode t = nodes.get(targetid.getTargetid());
		int count = 0;
		for(TargetNode nei: t.getNeighbors())
		{
			// check if neibor is not in supertarget
			if(!nodes.containsKey(nei.getTargetid()))
			{
				// if neibor is not in supertarget that means it has neighbor outside
				count++;
			}
		}
		return count;
	}
	

	public static void printPath(SuperTarget node) 
	{

		if(node == null)
			return;
		printPath(node.parent);

		for(TargetNode entry: node.entrypoints.values())
			for(TargetNode exit: node.exitpoints.values())
			{
				/*System.out.print( "("+ node.stid+","+entry.getTargetid()+","+exit.getTargetid()
				+","+node.allocateddistance+""+")->");*/

				System.out.print("("+node.stid+","+node.combid+")->");
			}

	}

	public static void printFinalPathWithProb(SuperTarget node, ArrayList<Double> fpath) 
	{

		if(node == null)
			return;
		printFinalPathWithProb(node.parent, fpath);
		for(Integer tnode: node.pathnodesid)
		{
			/*System.out.print( "("+ node.stid+","+entry.getTargetid()+","+exit.getTargetid()
				+","+node.allocateddistance+""+")->");*/

			//System.out.print("("+tnode+","+node.currentprob+")->");
			fpath.add(Math.ceil(tnode));
		}

	}


	public static SuperTarget mergeSuperTargets(SuperTarget st1, SuperTarget st2) {
		
		
		SuperTarget st = new SuperTarget();
		
		//new combid
		/*if(st1.nodes.size()==1)
			st.combid += st1.stid;
		else
		{
			st.combid += ","+st1.combid;
		}
		
		if(st2.nodes.size()==1)
			st.combid += ","+st1.stid;
		else
		{
			st.combid += ","+st2.combid;
		}*/
		st.stid = 10 + st1.stid + st2.stid;
		
		
		
		//new nodes
		for(TargetNode t: st1.nodes.values())
		{
			st.nodes.put(t.getTargetid(), t);
		}
		
		for(TargetNode t: st2.nodes.values())
		{
			st.nodes.put(t.getTargetid(), t);
		}
		
		
		//new potential APS
		for(TargetNode ap: st1.ap.values())
		{
			st.ap.put(ap.getTargetid(), ap);
		}
		for(TargetNode ap: st2.ap.values())
		{
			st.ap.put(ap.getTargetid(), ap);
		}
		
		
		
		//new neighbors
		for(SuperTarget n: st1.neighbors.values())
		{
			st.neighbors.put(n.stid , n);
		}

		for(SuperTarget n: st2.neighbors.values())
		{
			st.neighbors.put(n.stid , n);

		}
		
		
		
		
		return st;
	}


	public static SuperTarget mergeSuperTargets(SuperTarget st1, SuperTarget st2, int aid1, int aid2, HashMap<Integer,TargetNode> targetmaps) {

		SuperTarget st = new SuperTarget();
		
		
		//st.stid = 10 + st1.stid + st2.stid;
		
		 st.stid = 200 + st1.stid + st2.stid;
		
		//new nodes
		for(TargetNode t: st1.nodes.values())
		{
			st.nodes.put(t.getTargetid(), t);
		}
		
		for(TargetNode t: st2.nodes.values())
		{
			st.nodes.put(t.getTargetid(), t);
		}
		
		
		//new potential APS
		st.ap.put(aid1, targetmaps.get(aid1));
		st.ap.put(aid2, targetmaps.get(aid2));
		
		
		
		
//		//new neighbors
//		for(SuperTarget n: st1.neighbors.values())
//		{
//			st.neighbors.put(n.stid , n);
//		}
//
//		for(SuperTarget n: st2.neighbors.values())
//		{
//			st.neighbors.put(n.stid , n);
//
//		}
//		
		
		
		
		return st;
	}




}
