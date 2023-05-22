import numpy as np

def SimulateGenotypes(maf, N):
    """
    Simulate genotypes for N samples with minor 
    allele frequency maf (assuming HWE)
    
    Parameters
    ----------
    maf : float
        Minor allele frequency of the SNP
    N : int
        Sample size (number of people)
        
    Returns
    -------
    gts : np.array of floats
        Genotypes (scaled to have mean 0, variance 1)
        of each person
    """
    gts = []
    for i in range(N):
        gt = sum(random.random() < maf)+sum(random.random() < maf)
        gts.append(gt)
    # Technical note: scale to have mean 0 var 1
    gts = np.array(gts)
    gts = (gts-np.mean(gts))/np.sqrt(np.var(gts))
    return gts

def SimulatePhenotype(gts, Beta):
    """
    Simulate phenotypes under a linear model Y=beta*X+error
    
    Parameters
    ----------
    gts : np.array of floats
        Genotypes of each person
    Beta : float
        Effect size
        
    Returns
    -------
    pts : np.array och personf floats
        Simulated phenotype value of each person
    """
    if Beta<-1 or Beta>1:
        print("Error: Beta should be between -1 and 1")
        return [None]*len(pts)
    pts = Beta*gts + np.random.normal(0, np.sqrt(1-Beta**2), size=len(gts))
    return pts
============================================================================================================================================

%pylab inline
import seaborn as sns

####### You can play around with N, Beta, maf in this cell to see how the plots change
# Simulate an example association
N = 1000 # sample size (number of people)
Beta = 0.8 # Effect size
maf = 0.2 # Minor allele frequency

gts = SimulateGenotypes(maf, N)
pts = SimulatePhenotype(gts, Beta)
