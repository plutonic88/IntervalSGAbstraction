package perfectMMR;

import java.util.Vector;

class Result{
	Vector<Integer> left_idx;
	Vector<Integer> right_idx;
	int pos;
	int start_left;
	int end_left;
	int start_right;
	int end_right;
	public Result()
	{
		left_idx = new Vector<Integer>();
		right_idx = new Vector<Integer>();
		pos = -1;
		start_left = -1;
		end_left = -1;
		start_right = -1;
		end_right = -1;
	}
}
public class QuickSort {
	int n;
	double[] values;
	int[] id;
	//Outcome
	int[] sort_idx;
	public int[] sort_id;
	public QuickSort(double[] values, int[] id, int n)
	{
		this.values = values;
		this.n = n;
		this.id = id;
		sort_idx = new int[n];
		sort_id = new int[n];
	}
	
	public Result partition(Vector<Integer> sub_idx, int start, int end)
	{
		Result my_result = new Result();
		if(start == end)
		{
			my_result.pos = start;
			return my_result;
		}
		int pos = end;
		int num = end - start + 1;
		double r = values[sub_idx.get(num - 1)];
		for(int i = 0; i < num - 1; i++)
		{
			int temp_idx = sub_idx.get(i);
			if(values[temp_idx] > r)
			{
				my_result.right_idx.add(temp_idx);
				pos = pos - 1;
			}
			else my_result.left_idx.add(temp_idx);
		}
		my_result.start_left = start;
		my_result.end_left = pos - 1;
		my_result.start_right = pos + 1;
		my_result.end_right = end;
		my_result.pos = pos;
		return my_result;
	}
	
	public void quick_sort(Vector<Integer> sub_idx, int start, int end)
	{
		Result my_result = partition(sub_idx, start, end);
		sort_idx[my_result.pos] = sub_idx.get(sub_idx.size() - 1);
		if(!my_result.left_idx.isEmpty())
			quick_sort(my_result.left_idx, my_result.start_left, my_result.end_left);
		if(!my_result.right_idx.isEmpty())
			quick_sort(my_result.right_idx, my_result.start_right, my_result.end_right);
	}
	public void get_sort_id()
	{
		Vector<Integer> sub_idx = new Vector<Integer>();
		for(int i = 0; i < n; i++)
			sub_idx.add(i);
		quick_sort(sub_idx, 0, n - 1);
		for(int i = 0; i < n; i++)
			sort_id[i] = id[sort_idx[i]];
	}
}
