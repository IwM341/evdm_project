from .utils import *
import pickle

import numpy as np

try:
    from ._cpp_lib import pyevdm as evdm
except ImportError:
    pass
try:
    import scipy.sparse as sprs
    import scipy.sparse.linalg as sprslalg
except ImportError:
    pass

class FormatFloat:
    def __mul__(self,other):
        return f'{other:.2e}'
    def __rmul__(self,other):
        return f'{other:.2e}'
    def __add__(self,other):
        return f'{other:.2e}'
    def __radd__(self,other):
        return f'{other:.2e}'
ffmt = FormatFloat()

class SparseRop:
    def __init__(self,AP,P,m_chain,tau,
                 Nsample = 10000,MaxScat = 10000):
        self.m_chain = m_chain
        self.Nsample = Nsample
        self.MaxScat = MaxScat
        self.tau = tau
        self.AP = AP
        self.P = P
    def __matmul__(self,X):
        x1 = self.m_chain.evolute(X,self.tau,self.Nsample,self.MaxScat)
        x1 *= X.sum()
        Sl_0 = sprslalg.lgmres(self.AP,X,x0=(1/self.P.diagonal())*x1,
                #M=spsp.diags(1/Rmat.diagonal()), 
                inner_m = 30,maxiter=200)
        result = np.abs(self.P@Sl_0[0])
        if(result.sum() > 0):
            result*=(X.sum()/result.sum())
        return result
    
def PreCondition(m_mat : sprs.csc_matrix):
    if(isinstance(m_mat, sprs.csc_matrix)):
        Precond = sprs.diags(1/m_mat.diagonal())
        return (m_mat@Precond,Precond)
    else:
        Precond = m_mat.diagonal()
        return (m_mat*Precond,Precond)
def PreCondInv(Rmat : sprs.csc_matrix,MarkovChain,tau):
    if(isinstance(Rmat, sprs.csc_matrix)):
        (Rp,Pr) = PreCondition(Rmat)
        return SparseRop(Rp,Pr,MarkovChain,tau)
    else:
        return np.linalg.inv(Rmat)

def RInvMat(scat_mat,tau):
    scat_mat*=(-tau)
    if(isinstance(scat_mat,sprs.csc_matrix)):
        scat_mat += sprs.identity(scat_mat.shape[0])
    else:
        scat_mat += np.identity(scat_mat.shape[0])
    return scat_mat

def RInvPrecond(scat_mat,tau):
    scat_mat*=(-tau)
    if(isinstance(scat_mat,sprs.csc_matrix)):
        scat_mat += sprs.identity(scat_mat.shape[0])
    else:
        scat_mat += np.identity(scat_mat.shape[0])
    return PreCondition(scat_mat)

def ROperator(scat_mat,tau):
    if(isinstance(scat_mat,sprs.csc_matrix)):
        m_markov = evdm.MarkovChain(scat_mat)
        scat_mat*=(-tau)
        scat_mat += sprs.identity(scat_mat.shape[0])
    else:
        m_markov = None
        scat_mat*=(-tau)
        scat_mat += np.identity(scat_mat.shape[0])
    return PreCondInv(scat_mat,m_markov,tau)

def make_state(smat,capt,ann = None,evap = None,elastic_factor = 1):
    grid = capt.grid
    if(smat.grid.size != capt.grid.size):
        raise RuntimeError("sizes doesn't matches")
    mmat = to_diag64(smat,True)#.astype('float64'),
    N = capt.grid.size
    mmat[0:N,0:N] *= elastic_factor
    
    return {
        'grid':grid,
        'mat':mmat,
        'evap':smat.evap_histo.to_numpy().astype('float64'),
        'capt':capt.to_numpy().astype('float64'),
        'ann':ann,
        }

def load_state(filenameMat,filenameCapt,filenameAnn = None,elastic_factor = 1):
    smat = pickle.load(open(filenameMat,'rb'))
    capt = pickle.load(open(filenameCapt,'rb'))
    if(filenameAnn is None):
        ann = None
    else:
        ann = pickle.load(open(filenameAnn,'rb'))
    return make_state(smat,capt,ann,elastic_factor )

def CalcR(smatrix_np,tau,erase):
    if(erase):
        m_mat = smatrix_np
    else:
        m_mat = smatrix_np.copy()
    
    if(False):
        m_mat *= -0.5*tau
        m_diag = 1 + m_mat.diagonal()
        np.fill_diagonal(m_mat,m_diag)
        #print( np.abs(m_mat-np.identity(m_mat.shape[0])).sum())
        R2 = np.linalg.inv(m_mat)
        m_mat *= 2
        m_diag = m_mat.diagonal() - 1
        np.fill_diagonal(m_mat,m_diag)
        R1 = np.linalg.inv(m_mat)
        R2@=R2
        R2*=2
        R2-=R1 
        return R2
    else:
        return ROperator(m_mat,tau)
def CalcRFD(smatrix_np,tau,erase):
    if(erase):
        m_mat = smatrix_np
    else:
        m_mat = smatrix_np.copy()
    
    Nx = m_mat.shape[0]
    Nx_half = int(Nx/2)

    m_mat = m_mat[Nx_half:,0:Nx_half]
    m_diag = m_mat.diagonal()
    if(isinstance(m_mat,np.array)):
        np.fill_diagonal(m_mat,m_diag*0)
        np.fill_diagonal(m_mat,-np.sum(m_mat,axis = 0))
    else:
        m_mat.setdiag(0)
        m_mat.setdiag(-np.array(m_mat.sum(axis=0)).flatten())
    if(False):
        m_mat *= -0.5*tau
        m_diag = 1 + m_mat.diagonal()
        np.fill_diagonal(m_mat,m_diag)
        #print( np.abs(m_mat-np.identity(m_mat.shape[0])).sum())
        R2 = np.linalg.inv(m_mat)
        m_mat *= 2
        m_diag = m_mat.diagonal() - 1
        np.fill_diagonal(m_mat,m_diag)
        R1 = np.linalg.inv(m_mat)
        R2@=R2
        R2*=2
        R2-=R1
        return R2
    return ROperator(m_mat,tau)

def EvToTask(EvolveInfo,T_final,N,erase = False):
    tau = T_final/N
    R2 = CalcR(EvolveInfo['mat'],tau,erase)
    X = EvolveInfo['capt'].copy()
    Ann = EvolveInfo['ann']
    if(erase):
        EvolveInfo['mat'] = None
        EvolveInfo['ann'] = None
    
    
    return {'grid': EvolveInfo['grid'],'R':R2,'X':X,'A':Ann,'N':N,'tau':tau}



def make_task(capt,rmat,ann,N,tau):
    return {'grid': capt.grid,
            'R':rmat,'X':capt.to_numpy().astype('float64'),
            'A':ann,'N':N,'tau':tau}
def make_taskFD(capt,rmat,ann,N,tau):
    return {'grid': capt.grid,'R':rmat,
            'X':capt.to_numpy(1).astype('float64'),
            'A':ann,'N':N,'tau':tau}


def load_task(capt_fname,rmat_fname,ann_fname,N,tau):
    def unpickle(fname):
        return pickle.load(open(fname,'rb'))
    if(ann_fname is None):
        ann = None
    else:
        ann = unpickle(rmat_fname)
    return make_task(unpickle(capt_fname),ann,unpickle(ann_fname),N,tau)

def load_taskFD(capt_fname,rmat_fname,ann_fname,N,tau):
    def unpickle(fname):
        return pickle.load(open(fname,'rb'))
    if(ann_fname is None):
        ann = None
    else:
        ann = unpickle(rmat_fname)
    return make_taskFD(unpickle(capt_fname),ann,unpickle(ann_fname),N,tau)

def EvToTaskFastDecay(EvolveInfo,T_final,N,erase = False):
    tau = T_final/N
    Nx = EvolveInfo['mat'].shape[0]
    Nx_half = int(Nx/2)
    R2 = CalcRFD(EvolveInfo['mat'],tau,erase)
    #print('Rmat: ',R2.sum(axis=0).min()," ", R2.sum(axis=0).max())
    X = EvolveInfo['capt'][Nx_half:].copy()
    #print(np.sum(X))
    Ann = EvolveInfo['ann']
    if(erase):
        EvolveInfo['mat'] = None
        EvolveInfo['ann'] = None
    
    
    return {'grid': EvolveInfo['grid'],'R':R2,'X':X,'A':Ann,'N':N,'tau':tau}


def GetEvolveVector(DictTask,verbose = False):
    X = DictTask['X']
    R = DictTask['R']
    X1 = 0.0*X
    X0 = 1.0*X
    N = DictTask['N']
    Sum = X*0 # X*0.5

    #print( np.abs(R-np.identity(R.shape[0])).sum())
    for i in range(N):
        X1 = R@X0
        if(i==0 and i < N-1):
          Sum += 0.5*X1
        elif(i==N-1 and i > 0):
          Sum += 0.5*X1
        else:
          Sum += X1
        Y = X1
        X1 = X0
        X0 = Y
    Sum -= X0*0.5
    Sum *= DictTask['tau']
    return (X0,Sum)

def GetEvolveVectorAnn(DictTask,a_gamma,verbose = False):
    X = DictTask['X']
    R = DictTask['R']
    X1 = 0.0*X
    X0 = 1.0*X
    N = DictTask['N']
    Ann = DictTask['A'].A0
    tau = DictTask['tau']

    def AnnStep(AnnMat,N0_half,C0_half,tau,Tmp0,Tmp1,Tmp2,Tmp3,Tmp4):
        N0_half = np.fmax(0,N0_half,out=Tmp4)
        #print(Tmp0.shape,' ',AnnMat.shape,' ',N0_half.shape)
        AN_diag = np.dot(AnnMat,N0_half,out = Tmp0)
        AN_diag *= (-tau)
        ExpAnn_1 = np.exp(AN_diag,out=Tmp1) # = exp(-D tau)@C
        Tmp2[:] = C0_half
        Tmp2*=ExpAnn_1 # = exp(-D tau)@C

        AdotC = np.dot(AnnMat,Tmp2,out = Tmp3) # A*tau@exp(-D tau)@C 
        AdotC *=-tau
        AdotC*=N0_half
        AdotC += Tmp2 # (1 + A*tau)@exp(-D tau)@C 
        
        Tmp2 *= AN_diag # D*tau@C
        #Tmp2 *= 2

        AdotC -= Tmp2 # (1 + A*tau - D*tau)@exp(-D tau)@C
        AdotC*=ExpAnn_1 # exp(-D tau)@(1 + A*tau - D*tau)@exp(-D tau)@C
        return AdotC
    
    
    
    #print( np.abs(R-np.identity(R.shape[0])).sum())
    #Tmp0 = X0.copy()

    Nhlf = Ann.shape[0]
    TmpSum = X0.copy()
    Tmp0 = X0[0:Nhlf].copy()
    Tmp1 = Tmp0.copy()
    Tmp2 = Tmp0.copy()
    Tmp3 = Tmp0.copy()
    Tmp4 = Tmp0.copy()

    Sum = 0*X0 
    for i in range(N):
        R.dot(X0,out = X1)

        TmpSum[:] =  X1
        TmpSum += X0
        TmpSum *= (tau/2)
        Sum += TmpSum

        X1[0:Nhlf] = AnnStep(Ann,Sum[0:Nhlf],X1[0:Nhlf],tau*a_gamma,Tmp0,Tmp1,Tmp2,Tmp3,Tmp4)
        
        if(verbose):
            print(i,": ", 
                  "N = ", ffmt+Sum.sum(), 
                  ", C = ", ffmt+X1.sum(),
                  ' A = ', ffmt+a_gamma*np.dot(Sum[0:Nhlf],np.dot(Ann,Sum[0:Nhlf],out = Tmp0))
                  )
        Y = X1
        X1 = X0
        X0 = Y
    return (X0,Sum)

def GetEvolveVectorAnnNaive(DictTask,a_gamma,verbose = False,a_gamma_debug = None):
    if(a_gamma_debug is None):
        a_gamma_debug =a_gamma 
    C = DictTask['X']
    R = DictTask['R']
    

    N = DictTask['N']
    Ann = DictTask['A'].A0
    tau = DictTask['tau']

    def AnnStep(AnnMat,N0_half,tau,Tmp0,Tmp1,Tmp2):
        N0_half = np.fmax(0,N0_half,out=Tmp0)
        #print(Tmp0.shape,' ',AnnMat.shape,' ',N0_half.shape)
        AN_diag = np.dot(AnnMat,N0_half,out = Tmp1)
        AN_diag *= (-tau)
        ExpAnn_1 = np.exp(AN_diag,out=Tmp2) # = exp(-D tau)@C
        N0_half *= ExpAnn_1
        return N0_half
    
    
    
    #print( np.abs(R-np.identity(R.shape[0])).sum())
    #Tmp0 = X0.copy()

    Nhlf = Ann.shape[0]

    Tmp0 = C[0:Nhlf].copy()
     

    dNC = C*tau

    Tmp1 = Tmp0.copy()
    Tmp2 = Tmp0.copy()


    Sum = 0*C 
    TmpSum= Sum.copy()

    for i in range(N):
        Sum += dNC
        R.dot(Sum,out = TmpSum)
        TmpSum[0:Nhlf] = AnnStep(Ann,TmpSum[0:Nhlf],tau*a_gamma,Tmp0,Tmp1,Tmp2)
        
        if(verbose):
            print(i,": ", 
                  "N = ", ffmt+Sum.sum(), 
                  ' A = ', ffmt+a_gamma_debug*np.dot(TmpSum[0:Nhlf],np.dot(Ann,TmpSum[0:Nhlf],out = Tmp0))
                  )
        Y = Sum
        Sum = TmpSum
        TmpSum = Y

    return (C,Sum)


def GetResults(DictTask,CaptFinalD):
    D = {}
    D['Ne'] = DictTask['grid'].size_e

    Ann = DictTask['A']
    
    (CaptOut,Dout) = CaptFinalD
    Nsize = Dout.shape[0]
    
    try:
        DoutSlice = Dout[0:Ann.A0.shape[0]]
        ann_s = DoutSlice.dot(Ann.A0.dot(DoutSlice))
    except:
        ann_s = None
    m_distrib = evdm.Distrib(DictTask['grid'],Dout)

    D['Dout'] = m_distrib
    D['Cout'] = CaptOut
    D['Ann'] = ann_s 
    D['E'] = -DictTask['grid'].Epoints[0]+m_distrib.avarage(lambda e,l:e)
    return D

