import numpy as np

def DataGenerator(N_Samples=100, N_PeriodicGenes=20, N_NonPeriodicGenes=0, noise=True, Means=False,NoiseSigma=1):
   
    phi = np.random.rand((N_Samples))*2*np.pi
    coeff=np.random.normal(0, 1, (N_PeriodicGenes,2) )
    zeta=np.matrix([np.cos(phi),np.sin(phi)])
    E=np.matmul(coeff,zeta)
    if N_NonPeriodicGenes>0:
        NP=np.zeros((N_NonPeriodicGenes, N_Samples))
        E=np.vstack([E,NP])
  
    if Means:
        means = np.random.uniform( size=(N_PeriodicGenes+N_NonPeriodicGenes,1) )
        E=E+means

    if noise:
        noise_matrix=np.random.normal(0,NoiseSigma, size=(N_PeriodicGenes+N_NonPeriodicGenes,N_Samples) )
        E=E+noise_matrix

    return E,phi,coeff

def NonUniformDataGenerator(N_Samples=100, N_PeriodicGenes=20, N_NonPeriodicGenes=0, noise=True, Means=False):
   
    phi = np.random.rand((N_Samples))*np.pi
    coeff=np.random.normal(0, 1, (N_PeriodicGenes,2) )
    zeta=np.matrix([np.cos(phi),np.sin(phi)])
    E=np.matmul(coeff,zeta)
    if N_NonPeriodicGenes>0:
        NP=np.zeros((N_NonPeriodicGenes, N_Samples))
        E=np.vstack([E,NP])
  
    if noise:
        noise_matrix=np.random.normal( size=(N_PeriodicGenes+N_NonPeriodicGenes,N_Samples) )
        E=E+noise_matrix

    return E,phi,coeff


