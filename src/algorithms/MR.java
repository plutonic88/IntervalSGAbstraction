package algorithms;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.MatlabProxyFactory;
import matlabcontrol.MatlabProxyFactoryOptions;
import matlabcontrol.extensions.MatlabNumericArray;
import matlabcontrol.extensions.MatlabTypeConverter;
import models.SUQRAdversary;
import models.Target;

public class MR {
    int nIter = 10;
    double[][] adversary;
    double numRes;
    int numTargets;
    double[][] defCov;
    double[][] attLB;
    double[][] attUB;
    List<Target> targetList;
    
    // Output
    double[] attReward;
    double[] attPenalty;
    double[] optDefCov;
    double maxRegret;
    
    MatlabProxy proxy;
    
    int nConstraint = 0;
    Vector<double[]> coeffList;
    Vector<Double> constantList;
    public MR(List<Target> targetList, Map<Target, Double> defCov, SUQRAdversary adversary, double numRes, int nIter) 
    {
    	this.targetList = targetList;
    	this.numTargets = targetList.size();
    	this.defCov = new double[numTargets][];
    	this.attLB = new double[2 * numTargets][];
    	this.attUB = new double[2 * numTargets][];
    	this.adversary = new double[3][];
    	setDefCov(targetList, defCov);
    	setAttBound(targetList);
    	double[] tempDefCov = {adversary.w[adversary.DEF_COV]};
    	double[] tempAttReward = {adversary.w[adversary.ATT_REWARD]};
    	double[] tempAttPenalty = {adversary.w[adversary.ATT_PENALTY]};
    	this.adversary[0] = tempDefCov;
    	this.adversary[1] = tempAttReward;
    	this.adversary[2] = tempAttPenalty;
    	this.numRes = numRes;
    	this.nIter = nIter;
    	
    	attReward = new double[numTargets];
    	attPenalty = new double[numTargets];
    	optDefCov = new double[numTargets];
    	
    	coeffList = new Vector<double[]>();
    	constantList = new Vector<Double>();
    	
    	MatlabProxyFactoryOptions options = new MatlabProxyFactoryOptions.Builder().
	    		setMatlabStartingDirectory(new File("/Users/fake/research/col/MMR/MATLAB/PAWS")).
	    		setUsePreviouslyControlledSession(true).build();
	    MatlabProxyFactory factory = new MatlabProxyFactory(options);
	    try {
			proxy = factory.getProxy();
		} catch (MatlabConnectionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    public void addConstraint(double[] coeff, double constant)
    {
    	nConstraint++;
    	coeffList.add(coeff);
    	constantList.add(constant);
    }
    public void setDefCov(List<Target> targetList, Map<Target, Double> defCov)
    {
//    	int numTargets = targetList.size();
    	if(defCov != null)
    	for(Target t : targetList)
    	{
    		double[] tempDefCov = {defCov.get(t)};
    		this.defCov[t.id] = tempDefCov;
//    		double[] tempAttRewardLB = {t.getAttBoundStructure().getAttRewardLB()};
//    		double[] tempAttRewardUB = {t.getAttBoundStructure().getAttRewardUB()};
//    		double[] tempAttPenaltyLB = {t.getAttBoundStructure().getAttPenaltyLB()};
//    		double[] tempAttPenaltyUB = {t.getAttBoundStructure().getAttPenaltyUB()};
//    		this.attLB[t.id] = tempAttRewardLB;
//    		this.attLB[t.id + numTargets] = tempAttPenaltyLB;
//    		this.attUB[t.id] = tempAttRewardUB;
//    		this.attUB[t.id + numTargets] = tempAttPenaltyUB;
    	}
    }
    public void setAttBound(List<Target> targetList)
    {
    	int numTargets = targetList.size();
    	for(Target t : targetList)
    	{
    		double[] tempAttRewardLB = {t.getAttBoundStructure().getAttRewardLB()};
    		double[] tempAttRewardUB = {t.getAttBoundStructure().getAttRewardUB()};
    		double[] tempAttPenaltyLB = {t.getAttBoundStructure().getAttPenaltyLB()};
    		double[] tempAttPenaltyUB = {t.getAttBoundStructure().getAttPenaltyUB()};
    		this.attLB[t.id] = tempAttRewardLB;
    		this.attLB[t.id + numTargets] = tempAttPenaltyLB;
    		this.attUB[t.id] = tempAttRewardUB;
    		this.attUB[t.id + numTargets] = tempAttPenaltyUB;
    	}
    }
    public double[][] initSamples()
    {
    	double[][] samples = new double[nIter][3 * numTargets];
    	
    	return samples;
    }
    public void solve() throws MatlabConnectionException, MatlabInvocationException
    {
//    	int numTargets = defCov.length;
    	//Create a proxy, which we will use to control MATLAB
	    
//	    System.out.println(defCov.toString());
	    MatlabTypeConverter processor = new MatlabTypeConverter(proxy);
	    processor.setNumericArray("w", new MatlabNumericArray(adversary, null));
	    processor.setNumericArray("defCov", new MatlabNumericArray(defCov, null));
	    processor.setNumericArray("attLB", new MatlabNumericArray(attLB, null));
	    processor.setNumericArray("attUB", new MatlabNumericArray(attUB, null));
	    proxy.eval("nIter = " + nIter);
	    proxy.eval("nRes = " + numRes);
	    if(nConstraint == 0)
	    	proxy.eval("[x maxRegret] = computeMaxRegretZeroSumNoConstraint(nRes, w, defCov, attLB, attUB, nIter)");
	    else
	    {
	    	double[][] coeffArray = new double[coeffList.size()][];
	        double[][] constantArray = new double[coeffList.size()][];
	        int idx = 0;
	    	for(double[] coeff : coeffList)
	    	{
	    		coeffArray[idx] = coeff;
	    		double[] tempConsant = {constantList.get(idx)};
	    		constantArray[idx] = tempConsant;
	    		idx++;
	    	}
	    	
	    	processor.setNumericArray("coeffArray", new MatlabNumericArray(coeffArray, null));
	    	processor.setNumericArray("constantArray", new MatlabNumericArray(constantArray, null));
	    	proxy.eval("[x maxRegret] = computeMaxRegretZeroSum(nRes, w, defCov, attLB, attUB, nIter, coeffArray, constantArray)");
	    }
//	    proxy.eval("display(defCov)");
//	    proxy.eval("display(attLB)");
//	    proxy.eval("display(attUB)");
	    double[] x = (double[]) proxy.getVariable("x");
	    double[] maxRegret = (double[])proxy.getVariable("maxRegret");
	    for(int t = 0; t < numTargets; t++)
	    {
	    	optDefCov[t] = x[t];
	    	attReward[t] = x[numTargets + t];
	    	attPenalty[t] = x[2 * numTargets + t];
	    }
	    this.maxRegret = maxRegret[0];
    }
 
    public double getAttReward(Target t)
    {
    	return attReward[t.id];
    }
    public double getAttPenalty(Target t)
    {
    	return attPenalty[t.id];
    }
    public double getOptDefCov(Target t)
    {
    	return optDefCov[t.id];
    }
    public double getAttReward(int targetID)
    {
    	return attReward[targetID];
    }
    public double getAttPenalty(int targetID)
    {
    	return attPenalty[targetID];
    }
    public double getOptDefCov(int targetID)
    {
    	return optDefCov[targetID];
    }
    public double getMaxRegret()
    {
    	return maxRegret;
    }
    public void end()
    {
    	proxy.disconnect();
    	coeffList.clear();
    	constantList.clear();
    }
    
    public void setNumIter(int nIter)
    {
    	this.nIter = nIter;
    }
 
}
