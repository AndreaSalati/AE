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


def optimal_shift(phi, phi0, N=200 ):

    mad=12 #median abs deviation, super high, its a starting value, could be +inf
    #there are two symetries, rotational and also multiplication by -1 (other side of the circle)
    for i in range(N):
        
        offset=i/N*2*np.pi #creates many offsetts around the cycle
        theta=(phi-offset)%(2*np.pi) #creates the shifted phi vector
        delta=np.abs(theta-phi0)%(2*np.pi) #difference between the shifted

        for j in range(len(phi)):
            delta[j]=min(delta[j], 2*np.pi-delta[j]) #this checks if we should move in clockwise or counterclockwise direction
        
        if(np.median(delta)<mad):
            mad=np.median(delta)
            bestphi=theta
            j_temp=j
            sdel=delta
        #checking the same things but on the inverted situation
        phi=-phi%(2*np.pi)
        theta=(phi-offset)%(2*np.pi) #creates the shifted phi vector
        delta=np.abs(theta-phi0)%(2*np.pi) #difference between the shifted

        for j in range(len(phi)):
                delta[j]=min(delta[j], 2*np.pi-delta[j]) #this checks if we should move in clockwise or counterclockwise direction
            
        if(np.median(delta)<mad):
            mad=np.median(delta)
            bestphi=theta
            j_temp=j
            sdel=delta
    
    return bestphi,mad