package cs.com.allpair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import groupingtargets.SuperTarget;




public class AllPairShortestPath
{
	public int distancematrix[][];
	public int next[][];
	private int numberofvertices;
	public static final int INFINITY = 9999999;
	public static final int MY_NULL = 999999;

	public AllPairShortestPath(int numberofvertices)
	{
		distancematrix = new int[numberofvertices + 1][numberofvertices + 1];
		this.numberofvertices = numberofvertices;
		next = new int[numberofvertices + 1][numberofvertices + 1];

		for (int source = 1; source <= numberofvertices; source++)
		{
			for (int destination = 1; destination <= numberofvertices; destination++)
			{
				//distancematrix[source][destination] = adjacencymatrix[source][destination];
				next[source][destination] = MY_NULL;
			}
		}




	}

	public static ArrayList<Integer> getPath(int i, int j, int[][] next)
	{
		ArrayList<Integer> path = new ArrayList<Integer>();

		if(next[i][j]==MY_NULL)
		{
			return path;
		}

		path.add(i);

		while(i!=j)
		{
			i = next[i][j];
			path.add(i);
		}
		return path;

	}


	public  ArrayList<Integer> getPath(int src, int dest, HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback)
	{
		ArrayList<Integer> path = new ArrayList<Integer>();
		ArrayList<Integer> origpath = new ArrayList<Integer>();

		int i = map.get(src);
		int j = map.get(dest);


		if(next[i][j]==MY_NULL)
		{
			return origpath;
		}

		path.add(i);

		while(i!=j)
		{
			i = next[i][j];
			path.add(i);
		}

		for(Integer x: path)
		{
			origpath.add(mapback.get(x));
		}
		if(origpath.size()>2)
		{
			origpath.remove(0);
			origpath.remove(origpath.size()-1);
		}
		return origpath;

	}
	
	
	
	public  ArrayList<Integer> getPathInST(int src, int dest, HashMap<Integer, Integer> map, HashMap<Integer, Integer> mapback, SuperTarget tempst)
	{
		ArrayList<Integer> path = new ArrayList<Integer>();
		ArrayList<Integer> origpath = new ArrayList<Integer>();

		int i = map.get(src);
		int j = map.get(dest);


		if(next[i][j]==MY_NULL)
		{
			return origpath;
		}

		path.add(i);

		while(i!=j)
		{
			i = next[i][j];
			path.add(i);
		}

		for(Integer x: path)
		{
			origpath.add(mapback.get(x));
		}
		if(origpath.size()>2)
		{
			origpath.remove(0);
			origpath.remove(origpath.size()-1);
		}
		return origpath;

	}



	public int[][] allPairShortestPathWithPathConstruction(int adjacencymatrix[][], ArrayList<Integer> domint)
	{
		for (int source = 1; source <= numberofvertices; source++)
		{
			for (int destination = 1; destination <= numberofvertices; destination++)
			{
				distancematrix[source][destination] = adjacencymatrix[source][destination];
				next[source][destination] = destination;
			}
		}

		for (int intermediate = 1; intermediate <= numberofvertices; intermediate++)
		{
			if(domint.contains(intermediate))
			{

				for (int source = 1; source <= numberofvertices; source++)
				{
					for (int destination = 1; destination <= numberofvertices; destination++)
					{
						if (distancematrix[source][intermediate] + distancematrix[intermediate][destination]
								< distancematrix[source][destination])
						{
							distancematrix[source][destination] = distancematrix[source][intermediate] 
									+ distancematrix[intermediate][destination];
							next[source][destination] = next[source][intermediate];
						}
					}
				}
			}
		}

		/*for (int source = 1; source <= numberofvertices; source++)
			System.out.print("\t" + source);

		System.out.println();
		for (int source = 1; source <= numberofvertices; source++)
		{
			System.out.print(source + "\t");
			for (int destination = 1; destination <= numberofvertices; destination++)
			{
				System.out.print(distancematrix[source][destination] + "\t");
			}
			System.out.println();
		}*/
		return distancematrix;
	}


	public int[][] allPairShortestPath(int adjacencymatrix[][])
	{
		for (int source = 1; source <= numberofvertices; source++)
		{
			for (int destination = 1; destination <= numberofvertices; destination++)
			{
				distancematrix[source][destination] = adjacencymatrix[source][destination];
				next[source][destination] = destination;
			}
		}

		for (int intermediate = 1; intermediate <= numberofvertices; intermediate++)
		{
			for (int source = 1; source <= numberofvertices; source++)
			{
				for (int destination = 1; destination <= numberofvertices; destination++)
				{
					if (distancematrix[source][intermediate] + distancematrix[intermediate][destination]
							< distancematrix[source][destination])
					{
						distancematrix[source][destination] = distancematrix[source][intermediate] 
								+ distancematrix[intermediate][destination];
						next[source][destination] = next[source][intermediate];
					}
				}
			}
		}

		/*for (int source = 1; source <= numberofvertices; source++)
			//System.out.print("\t" + source);

		//System.out.println();
		for (int source = 1; source <= numberofvertices; source++)
		{
			//System.out.print(source + "\t");
			for (int destination = 1; destination <= numberofvertices; destination++)
			{
				//System.out.print(distancematrix[source][destination] + "\t");
			}
			//System.out.println();
		}*/
		return distancematrix;
	}


	/*
	public static void main(String... arg)
	{
		int adjacency_matrix[][];
		int numberofvertices;

		Scanner scan = new Scanner(System.in);
		System.out.println("Enter the number of vertices");
		numberofvertices = scan.nextInt();

		adjacency_matrix = new int[numberofvertices + 1][numberofvertices + 1];
		System.out.println("Enter the Weighted Matrix for the graph");
		for (int source = 1; source <= numberofvertices; source++)
		{
			for (int destination = 1; destination <= numberofvertices; destination++)
			{
				adjacency_matrix[source][destination] = scan.nextInt();
				if (source == destination)
				{
					adjacency_matrix[source][destination] = 0;
					continue;
				}
				if (adjacency_matrix[source][destination] == 0)
				{
					adjacency_matrix[source][destination] = INFINITY;
				}
			}
		}

		System.out.println("The Transitive Closure of the Graph");
		AllPairShortestPath allPairShortestPath= new AllPairShortestPath(numberofvertices);
		allPairShortestPath.allPairShortestPathWithPathConstruction(adjacency_matrix);
		ArrayList<Integer> path = getPath(1, 2, allPairShortestPath.next);

		System.out.println("Path from 1 to 2: " );
		for(Integer x: path)
		{
			System.out.print(x+" ");
		}


		scan.close();
	}*/
}
