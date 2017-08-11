package perfectMMR;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import lpWrapper.Configuration;
import lpWrapper.MIProblem;

public class MRZeroSumLP extends MIProblem{
	double x;
	double xstar;
	double c;
	double cstar;
	double rLB;
	double rUB;
	double pLB;
	double pUB;
	int t = -1;
	HashMap<String, Integer> varMap; // for MILP, column variables
	HashMap<String, Integer> rowMap;	// for MILP, row constraints
	
	public MRZeroSumLP()
	{
		varMap = new HashMap<String, Integer>();
		rowMap = new HashMap<String, Integer>();
	}
	public void setCoeff(double x, double xstar, double c, double cstar
			, double rLB, double rUB, double pLB, double pUB, int t)
	{
		this.x = x;
		this.xstar = xstar;
		this.c = c;
		this.cstar = cstar;
		this.rLB = rLB;
		this.rUB = rUB;
		this.pLB = pLB;
		this.pUB = pUB;
		if(this.t != t)
		{
			this.updateObj();
			this.updateCoeffStar();
			this.updateCoeff();
			this.updateBound();
		}
		else
			this.updateCoeff();
		this.t = t;
		
	}
	@Override
	protected void setProblemType() {
		// TODO Auto-generated method stub
		this.setProblemName("MRZeroSumLP");
		this.setProblemType(PROBLEM_TYPE.LP, OBJECTIVE_TYPE.MAX);
	}
	

	@Override
	protected void setColBounds() {
		int index = 1;
		// TODO Auto-generated method stub
		addAndSetColumn("v", BOUNDS_TYPE.FREE, -Configuration.MM, Configuration.MM, VARIABLE_TYPE.CONTINUOUS, 1.0);
		varMap.put("v", index++);
		addAndSetColumn("r", BOUNDS_TYPE.FREE, -Configuration.MM, Configuration.MM, VARIABLE_TYPE.CONTINUOUS, 0.0);
		varMap.put("r", index++);
		addAndSetColumn("p", BOUNDS_TYPE.FREE, -Configuration.MM, Configuration.MM, VARIABLE_TYPE.CONTINUOUS, 0.0);
		varMap.put("p", index++);
	}
	public void updateObj()
	{
		int index = varMap.get("r");
		this.setObjectiveCoef(index, 1 - xstar);
		index = varMap.get("p");
		this.setObjectiveCoef(index, xstar);
	}
	public void updateCoeffStar()
	{
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		ja.add(varMap.get("r"));
		ar.add(1 - xstar);
		ja.add(varMap.get("p"));
		ar.add(xstar);
		if(rowMap.containsKey("CStar"))
		{
			setMatRow(rowMap.get("CStar"), ja, ar);
			resetRowBound(rowMap.get("CStar"), BOUNDS_TYPE.LOWER, cstar, Configuration.MM);
		}
		else
		{
			addAndSetRow("CStar", BOUNDS_TYPE.LOWER, cstar, Configuration.MM);
			rowMap.put("CStar", this.getNumRows());
			setMatRow(this.getNumRows(), ja, ar);
		}
		ja.clear();
		ar.clear();
	}
	public void updateCoeff()
	{
		List<Integer> ja = new ArrayList<Integer>();
		List<Double> ar = new ArrayList<Double>();
		ja.add(varMap.get("v"));
		ar.add(1.0);
		ja.add(varMap.get("r"));
		ar.add(1 - x);
		ja.add(varMap.get("p"));
		ar.add(x);
		if(rowMap.containsKey("C1"))
		{
			setMatRow(rowMap.get("C1"), ja, ar);
			resetRowBound(rowMap.get("C1"), BOUNDS_TYPE.UPPER, -Configuration.MM, 0.0);
		}
		else
		{
			addAndSetRow("C1", BOUNDS_TYPE.UPPER, -Configuration.MM, 0.0);
			rowMap.put("C1", this.getNumRows());
			setMatRow(this.getNumRows(), ja, ar);
		}
		ar.clear();
		ja.clear();
		ja.add(varMap.get("v"));
		ar.add(1.0);
		if(rowMap.containsKey("C2"))
		{
			setMatRow(rowMap.get("C2"), ja, ar);
			resetRowBound(rowMap.get("C2"), BOUNDS_TYPE.UPPER, -Configuration.MM, c);
		}
		else
		{
			addAndSetRow("C2", BOUNDS_TYPE.UPPER, -Configuration.MM, c);
			rowMap.put("C2", this.getNumRows());
			setMatRow(this.getNumRows(), ja, ar);
		}
		
		ja.clear();
		ar.clear();
	}
	
	public void updateBound()
	{
		this.resetColumnBound(varMap.get("r"), BOUNDS_TYPE.DOUBLE, rLB, rUB);
		this.resetColumnBound(varMap.get("p"), BOUNDS_TYPE.DOUBLE, pLB, pUB);
	}
	public double getObjValue()
	{
		return this.getLPObjective();
	}
	public double getReward()
	{
		return this.getColumnPrimal(varMap.get("r"));
	}
	public double getPenalty()
	{
		return this.getColumnPrimal(varMap.get("p"));
	}
	public double getDefEU()
	{
		return this.getColumnPrimal(varMap.get("v"));
	}
	@Override
	protected void setRowBounds() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void generateData() {
		// TODO Auto-generated method stub
		
	}

}
