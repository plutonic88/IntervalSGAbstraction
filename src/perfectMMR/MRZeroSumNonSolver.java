package perfectMMR;

import java.io.File;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.MatlabProxyFactory;
import matlabcontrol.MatlabProxyFactoryOptions;
import matlabcontrol.extensions.MatlabNumericArray;
import matlabcontrol.extensions.MatlabTypeConverter;

public class MRZeroSumNonSolver {
	int nTargets;
	int nRes;
	int nIter = 1;
	
	double[][] xStar;
	double[][] advPayoffLB;
	double[][] advPayoffUB;
	
	public double maxRegret;
	public double[] optDefCov;
	public double[] optAdvPayoff;
	
	MatlabProxy proxy;
	public MRZeroSumNonSolver(int nTargets, int nRes, double[] xStar, double[] advPayoffLB, double[] advPayoffUB)
	{
		this.nTargets = nTargets;
		this.nRes = nRes;
		
		optDefCov = new double[nTargets];
		optAdvPayoff = new double[2 * nTargets];
		
		this.xStar = new double[nTargets][];
    	this.advPayoffLB = new double[2 * nTargets][];
    	this.advPayoffUB = new double[2 * nTargets][];
    	for(int i = 0; i < nTargets; i++)
    	{
    		if(xStar != null)
    		{
	    		double[] tempDefCov = {xStar[i]};
	    		this.xStar[i] = tempDefCov;
    		}
    		double[] tempAttRewardLB = {advPayoffLB[i]};
    		double[] tempAttRewardUB = {advPayoffUB[i]};
    		double[] tempAttPenaltyLB = {advPayoffLB[nTargets + i]};
    		double[] tempAttPenaltyUB = {advPayoffUB[nTargets + i]};
    		this.advPayoffLB[i] = tempAttRewardLB;
    		this.advPayoffLB[i + nTargets] = tempAttPenaltyLB;
    		this.advPayoffUB[i] = tempAttRewardUB;
    		this.advPayoffUB[i + nTargets] = tempAttPenaltyUB;
    	}
		
		MatlabProxyFactoryOptions options = new MatlabProxyFactoryOptions.Builder().
	    		setMatlabStartingDirectory(new File("/Users/thanhnguyen/Documents/WORKS/UAV/AAMAS15CODES/THANH/MATLAB")).
	    		setUsePreviouslyControlledSession(true).build();
	    MatlabProxyFactory factory = new MatlabProxyFactory(options);
	    try {
			proxy = factory.getProxy();
		} catch (MatlabConnectionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public void setXStar(double[] xStar)
	{
		for(int i = 0; i < nTargets; i++)
    	{
    		double[] tempDefCov = {xStar[i]};
    		this.xStar[i] = tempDefCov;
    	}
	}
	public void setIter(int nIter)
	{
		this.nIter = nIter;
	}
	public void solve() throws MatlabInvocationException
	{
		MatlabTypeConverter processor = new MatlabTypeConverter(proxy);
	    processor.setNumericArray("xStar", new MatlabNumericArray(xStar, null));
	    processor.setNumericArray("attLB", new MatlabNumericArray(advPayoffLB, null));
	    processor.setNumericArray("attUB", new MatlabNumericArray(advPayoffUB, null));
	    proxy.eval("nRes = " + nRes);
	    proxy.eval("nIter = " + nIter);
	    proxy.eval("[x, maxRegret] = computeMRZeroSumPerfectAtt(nRes, xStar, attLB, attUB, nIter)");
//	    proxy.eval("display(defCov)");
//	    proxy.eval("display(attLB)");
//	    proxy.eval("display(attUB)");
	    double[] x = (double[]) proxy.getVariable("x");
	    double[] maxRegret = (double[])proxy.getVariable("maxRegret");
	    for(int t = 0; t < nTargets; t++)
	    {
	    	optDefCov[t] = x[t];
	    	optAdvPayoff[t] = x[nTargets + t];
	    	optAdvPayoff[t + nTargets] = x[2 * nTargets + t];
	    }
	    this.maxRegret = maxRegret[0];
	}
	public void end()
	{
		proxy.disconnect();
	}
}
