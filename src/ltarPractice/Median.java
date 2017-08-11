package ltarPractice;

import java.util.HashSet;

public class Median {
Double median(Double[] array)
	{
		if (array.length == 1)
		{
			return array[0];
		}
		int smallerCount = 0;
		for (int i = 0 ; i < array.length ; i++)
		{
			for (int j = 0 ; j < array.length ; j++)
			{
				if (array[i] < array[j])
				{
					smallerCount++;
				}
			}
			if (smallerCount == (array.length - 1)/2)
			{
				return array[i];
			}
			smallerCount = 0;
		}
		return 1.0; //should never happen
	}
	// finds pivot element of a given array recursively using DeterministicSelect
	Double Pivot(Double[] array)
	{
		if (array.length == 1)
		{
			return array[0];
		}
		//Divide A into n/5 groups of 5 elements each
		int numGroups = array.length / 5;
		if (array.length % 5 > 0)
		{
			numGroups++;
		}
		Double[] setOfMedians = new Double[numGroups];
		for (int i = 0 ; i < numGroups ; i++)
		{
			Double[] subset;
			if(array.length % 5 > 0)
			{
				if (i == numGroups - 1)
				{
					subset = new Double[array.length % 5];
				}
				else
				{
					subset = new Double[5];
				}
			}
			else
			{
				subset = new Double[5];
			}
			for (int j = 0; j < subset.length ; j++)
			{
				subset[j] = array[5*i+j];
			}
			//Find the median of each group
			setOfMedians[i] = median(subset);
		}
		//Use DeterministicSelect to find the median, p, of the n/5 medians
		Double goodPivot = DeterministicSelect(setOfMedians, setOfMedians.length/2 );
		return goodPivot;
	}

//The algorithm in words:
//1. Divide n elements into groups of 5
//2. Find median of each group (How? How long?)
//3. Use Select() recursively to find median x of the n/5? medians
//4. Partition the n elements around x. Let k = rank(x)
//5. if (i == k) then return x else (i > k) use Select() recursively to find (i-k)th
//smallest element in last partition
//source
//Lecture PDF mentioned in the blog post
//and MIT Lecture 6 order statistics.

public	Double DeterministicSelect(Double[] array, int k)
	{
		if (array.length == 1)
		{
			return array[0];
		}
		Double pivot = Pivot(array);
		//set is used to ignore duplicate values
		HashSet<Double> A1 = new HashSet<Double>();
		HashSet<Double> A2 = new HashSet<Double>();
		for (int i = 0; i < array.length ; i++)
		{
			if (array[i] < pivot)
			{
				A1.add(array[i]);
			}
			else if (array[i] > pivot)
			{
				A2.add(array[i]);
			}
		}
		if (A1.size() >= k)
		{
			return DeterministicSelect(A1.toArray(new Double[0]) ,k);
		}
		else if (A1.size() == k-1)
		{
			return pivot;
		}
		else
		{
			return DeterministicSelect(A2.toArray(new Double[0]) , k - A1.size() - 1);
		}
	}

}
