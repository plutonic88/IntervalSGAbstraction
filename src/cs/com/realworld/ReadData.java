package cs.com.realworld;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class ReadData {

	/*public static void main(String[] args) throws Exception 
	{
		int row = 560;
		int col = 560;
		double [][] u = new double[row][col];
		double [][] d = new double[row][col];
		ReadData.readData(row, col, u, d);


	}*/
	
	
	public static void getChunk(int row, int col, double[][] utility, double[][] elevation,
			int rstart, int cstart, int rend, int cend, double[][] u, double[][] e)
			{
				for(int i= rstart; i<rend; i++)
				{
					for(int j= cstart; j<cend; j++)
					{
						u[i][j] = (int)Math.floor(utility[i][j]);
						e[i][j] = elevation[i][j];
					}
				}
		
			}





	public static void readData(int row, int col, double[][] utility, double[][] elevation)
	{
		try{
			FileInputStream fstream = new FileInputStream("tn_utility_density_RasterPiece.asc");
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;

			int r = 0;

			while ((strLine = br.readLine()) != null)  
			{
				String[] tokens = strLine.split(" ");
				//System.out.println("hi");

				int j=0;
				for(String x: tokens)
				{
					utility[r][j]= Double.parseDouble(x);

					j++;
				}




				r++;

			}
			in.close();
		}catch (Exception e)
		{
			System.err.println("Error: " + e.getMessage());
		}
		
		try{
			FileInputStream fstream = new FileInputStream("tn_elevation.asc");
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;

			int r = 0;

			while ((strLine = br.readLine()) != null)  
			{
				String[] tokens = strLine.split(" ");
				//System.out.println("hi");

				int j=0;
				for(String x: tokens)
				{
					elevation[r][j]= Double.parseDouble(x);

					j++;
				}




				r++;

			}
			in.close();
		}catch (Exception e)
		{
			System.err.println("Error: " + e.getMessage());
		}
		
		
		
		//System.out.println("done");
	}

}
