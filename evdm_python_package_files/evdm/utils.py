from . import ff
from ._interfaces import *
import numpy as np

try:
    from ._cpp_lib import pyevdm as evdm
except ImportError:
    pass

try:
    import scipy as sp
except ImportError:
    pass

def pkl_save(obj,fname):
    import pickle
    pickle.dump(obj,open(fname,'wb'))
def pkl_load(fname):
    import pickle
    return pickle.load(open(fname,'rb')) 

class Cout_t:
    def __lshift__(self,value):
        print(value)
        return self
    def __repr__(self):
        return ""
cout = Cout_t()

def MakeLogGrid(Body : evdm.Body,h0,Ne,Nl = None,r = 1.1,dtype= 'float',lpoints = None):
    

    Emax = Body.phi[1][0]
    if(h0 >= Emax/Ne):
        raise ValueError(f"h0 ({h0}) >= H0 ({ Emax/Ne})")
    def logksi(x,ksi):
        return int(np.log(x)/np.log(ksi))
    def Phi(H0n,L,N,ksi,h0):
        return  L/(1/(ksi-1)*(1-h0/H0n) + N - logksi(H0n/h0,ksi))
    
    assert( Nl is not None or lpoints is not None)
    H0 = abs(Emax)/Ne
    H1 = Phi(H0,Emax,Ne,r,h0)
    sch = 0
    eq = lambda x,y: abs(x-y)/abs(x+y) < 1e-6
    while(not eq(H0,H1)):
        print(f"H0 = {H0}, h0 = {h0}, log = ",logksi(H0/h0,r))
        H0 = H1
        H1 = Phi(H0,Emax,Ne,r,h0)
        sch = sch + 1
        if(sch == 30):
            raise ValueError("bad recursion")
    N0 = logksi(H0/h0,r)
    Nr = Ne - N0
    #print(h0,', ',H0)
    #print(f"Ns = {N0}, {Nr}")
    powers = [r**k for k in range(N0)]
    msums = h0*np.array([sum(powers[0:i]) for i in range(len(powers))])
    #print(f"E: 0,{msums[1]},...,{msums[-1]} ,{Emax}")
    #print(msums)
    H0eff  = (Emax-msums[-1])/(Nr+1)
    #print(f"H0eff = {H0eff} vs {H0}")
    Epoints = np.array(list(msums) + [msums[-1] + H0eff*i for i in range(1,Nr+2)])
    Epoints[-1] = Emax
    
    Epoints = np.array(Epoints) - Emax
    #print('Epoints (normed)= ', (Epoints-Epoints[0])/(Epoints[-1]-Epoints[0]))
    
    if(lpoints is not None):
        Lpoints = lpoints
    else:
        Lpoints = np.linspace(0,1,Nl+1)
    #print(Lpoints)
    return evdm.MakeGrid(Body,Epoints,Lpoints,ptypes = 2,value_type = dtype)

def Grid_sqrt_L(N = 100,alpha = 0.04,delta=0.04):

    xi = np.linspace(0,1,N+1,endpoint=True)
    psy = ( ((1+alpha)**2-(1+alpha-xi)**2)/(1+2*alpha))
    m_l = np.sqrt(((1+2*delta)*psy+delta**2))-delta
    return m_l
    
def GridPreparaion(DoutPre : evdm.Distrib,CaptPre,Ne_pre,Ne_full,Nl,
                   alpha=0.04,delta = 0.04,width = 2,nmax_in_bin=0.1):

    import scipy as sp
    with  np.errstate(all='raise'):
        ED = DoutPre.Edistrib()
        

        Egrid = ED[0]
        Evals = ED[1]/ED[1].sum()
        E0 = Egrid[0]

        HD = sp.stats.rv_histogram((Evals,Egrid),density=True)
        
        Epoints = np.linspace(0,1,2001)

        def rho_des(x):
            index = min(np.searchsorted(Egrid,(1-x)*E0, side='right') - 1,Evals.size-1)
            h0 = (Egrid[index+1]-Egrid[index])/np.abs(E0)
            return max( (Evals[index])/(h0*nmax_in_bin),Ne_pre)
        
        Rho0 = np.vectorize(rho_des)(Epoints)

        dE = HD.ppf(0.94) - HD.ppf(0.06)
        Emin = HD.ppf(0.02)
        #print(Emin,", ",dE)

        def smooh_parab(X0,Rho0,H,Xmin):
            i_max = Rho0.argmax()
            xmax = X0[i_max]
            Rmax= Rho0[i_max]
            return np.where(X0 >= Xmin, np.maximum( Rmax*(1-((X0-xmax)/H)**2),Rho0  ), Rho0)
        Rho1 = smooh_parab(Epoints,Rho0,width*dE/np.abs(E0),1 - (Emin/E0))


        if(CaptPre):
            EC = CaptPre.Edistrib()
            EgridC = EC[0]
            EvalsC = EC[1]/EC[1].sum()

            def rho_des1(x):
                index = min(np.searchsorted(EgridC,(1-x)*E0, side='right') - 1,EvalsC.size-1)
                h0 = (EgridC[index+1]-EgridC[index])/np.abs(E0)
                return max( (EvalsC[index])/(h0*0.1),Ne_pre)
            Rho2 = np.vectorize(rho_des1)(Epoints)

            Rho_final = np.maximum(Rho1,Rho2)
        else:
            Rho_final = Rho1
        
        from scipy.integrate import cumulative_trapezoid as cumtrapz
        from scipy.interpolate import interp1d

        CDF_Nx = cumtrapz(Rho_final,Epoints,initial=0)
        import math
        Nfull = Ne_full or math.ceil(CDF_Nx[-1])
        #print(CDF_Nx)
        Probs = np.linspace(0,1,Nfull+1)
        CDF_Nx *= (1/CDF_Nx[-1])
        CDF_Nx[-1] = 1
        Func_1 = interp1d(CDF_Nx,Epoints,'linear',fill_value = (0,1))
        Epoints1 = Func_1(Probs)

        Epoints_dim = Epoints1*np.abs(E0)-np.abs(E0)
        Epoints_dim[0] = -np.abs(E0)
        Epoints_dim[-1] = 0
    
    def gen_l_grid(alpha,delta,N):
        xi = np.linspace(0,1,N+1,endpoint=True)
        psy = ( ((1+alpha)**2-(1+alpha-xi)**2)/(1+2*alpha))
        m_l = np.sqrt(((1+2*delta)*psy+delta**2))-delta
        return m_l

    grid0 = DoutPre.grid
    return evdm.MakeGrid(grid0.body,Epoints_dim,gen_l_grid(alpha,delta,Nl),grid0.ptypes)


def ConditionRefineGrid(grid,Ne,Nl,ECondition,LCondition = None):

    grid_dict = grid.__getstate__().copy()

    if(LCondition  == None):
        LCondition = lambda e,l:ECondition(e)

    el_g = grid_dict['grid']['InnerGrids']['value']
    m_grid_list = el_g['InnerGrids']
    m_e_g = el_g['Grid']
    e_final_list = []
    l_grid_final_list = []

    #print('el_g', m_e_g)
    #raise KeyboardInterrupt

    for i in range(len(m_e_g)-1):
        Ne_eff = 1
        if(ECondition(m_e_g[i+1])):
            Ne_eff = Ne
        he = (m_e_g[i+1]-m_e_g[i])/Ne_eff
        for j in range(Ne_eff):
            m_e = m_e_g[i] + he*j
            e_final_list.append(m_e_g[i] + he*j)
            m_l_g = m_grid_list[i]
            l_refine = []
            for k in range(len(m_l_g) - 1):
                Nl_eff = 1
                if(LCondition(m_e,m_l_g[k+1])):
                    Nl_eff = Nl
                hl = (m_l_g[k+1]-m_l_g[k])/Nl_eff
                for m in range(Nl_eff):
                    l_refine.append(m_l_g[k] + m*hl)
            l_refine.append(m_l_g[-1])
            l_grid_final_list.append(np.array(l_refine))

    e_final_list.append(m_e_g[-1])
    el_g['Grid'] = np.array(e_final_list)
    el_g['InnerGrids'] = l_grid_final_list
    grid_dict.pop('LE')
    grid_dict.pop('Trajs')
    return evdm.GridEL(grid_dict)


def CoarsenGrid(grid,Ne,Nl):
    import numpy as np
    grid_dict = grid.__getstate__().copy()

    el_g = grid_dict['grid']['InnerGrids']['value']
    m_grid_list = el_g['InnerGrids']
    m_e_g = el_g['Grid']
    

    l_grid_final_list = []
    #print('m_e_g', m_e_g)
    #raise KeyboardInterrupt
    if(isinstance(m_e_g,dict)):
        e_final = m_e_g
        sz = e_final['size']
        e_final['size'] = (sz-1)//Ne + 1
        for i in range( e_final['size']-1):
            m_l_g = m_grid_list[Ne*i]
            sz = m_l_g['size']
            m_l_g['size'] = (sz-1)//Nl + 1
            l_grid_final_list.append(m_l_g)
    else:
        if((len(m_e_g)-1)%Ne == 0):
            e_final = np.array(m_e_g[0::Ne])
        else:
            e_final = np.concat([m_e_g[0::Ne],m_e_g[-1]])
        
        for i in range(len(e_final)-1):
            m_l_g = m_grid_list[i]
            if((len(m_l_g)-1)%Nl == 0):
                l_final = np.array(m_l_g[0::Nl])
            else:
                l_final = np.concat([m_l_g[0::Nl],m_l_g[-1]])
            l_grid_final_list.append(l_final)

    el_g['Grid'] = e_final
    el_g['InnerGrids'] = l_grid_final_list
    grid_dict.pop('LE')
    grid_dict.pop('Trajs')
    return evdm.GridEL(grid_dict)


def get_solar_abundance(Name,ElementA = None):
    # Фотосферные данные (солнечная распространенность)
    print(f"Warning: using photospheric abundance for {Name} (inner solar abundance not available)")
    photoshperic_dens = {
        "H": 12.00, "He": 10.93, "Li": 1.05, "Be": 1.38, "B": 2.70,
        "C": 8.43, "N": 7.83, "O": 8.69, "F": 4.56, "Ne": 7.93,
        "Na": 6.24, "Mg": 7.60, "Al": 6.45, "Si": 7.51, "P": 5.41,
        "S": 7.12, "Cl": 5.50, "Ar": 6.40, "K": 5.03, "Ca": 6.34,
        "Sc": 3.15, "Ti": 4.95, "V": 3.93, "Cr": 5.64, "Mn": 5.43,
        "Fe": 7.50, "Co": 4.99, "Ni": 6.22, "Cu": 4.19, "Zn": 4.56,
        "Ga": 3.04, "Ge": 3.65, 
        "Kr": 3.25, "Rb": 2.52, "Sr": 2.87, "Y": 2.21, "Zr": 2.58,
        "Nb": 1.46, "Mo": 1.88, "Ru": 1.75, "Rh": 0.91, "Pd": 1.57,
        "Ag": 0.94,             "In": 0.80, "Sn": 2.04, 
                                            "Xe": 2.24,  "Ba": 2.18,
        "La": 1.10, "Ce": 1.58, "Pr": 0.72, "Nd": 1.42, "Sm": 0.96,
        "Eu": 0.52, "Gd": 1.07, "Tb": 0.30, "Dy": 1.10, "Ho": 0.48,
        "Er": 0.92, "Tm": 0.10, "Yb": 0.84, "Lu": 0.10, "Hf": 0.85,
                    "W" : 0.85,             "Os": 1.40, "Ir": 1.38, 
                    "Au": 0.92,             "Tl": 0.90, "Pb": 1.75, 
                    "Th": 0.02, 
    }
    
    # Метеоритные данные (распространенность в метеоритах)
    meteor_dens = {
        "H": 8.22, "He": 1.29, "Li": 3.26, "Be": 1.30, "B": 2.79,
        "C": 7.39, "N": 6.26, "O": 8.40, "F": 4.42, "Ne": -1.12,
        "Na": 6.27, "Mg": 7.53, "Al": 6.43, "Si": 7.51, "P": 5.43,
        "S": 7.15, "Cl": 5.23, "Ar": -0.50, "K": 5.08, "Ca": 6.29,
        "Sc": 3.05, "Ti": 4.91, "V": 3.96, "Cr": 5.64, "Mn": 5.48,
        "Fe": 7.45, "Co": 4.87, "Ni": 6.20, "Cu": 4.25, "Zn": 4.63,
        "Ga": 3.08, "Ge": 3.58, "As": 2.30, "Se": 3.34, "Br": 2.54,
        "Kr": -2.27, "Rb": 2.36, "Sr": 2.88, "Y": 2.17, "Zr": 2.53,
        "Nb": 1.41, "Mo": 1.94, "Ru": 1.76, "Rh": 1.06, "Pd": 1.65,
        "Ag": 1.20, "Cd": 1.71, "In": 0.76, "Sn": 2.07, "Sb": 1.01,
        "Te": 2.18, "I": 1.55, "Xe": -1.95, "Cs": 1.08, "Ba": 2.18,
        "La": 1.17, "Ce": 1.58, "Pr": 0.76, "Nd": 1.45, "Sm": 0.94,
        "Eu": 0.51, "Gd": 1.05, "Tb": 0.32, "Dy": 1.13, "Ho": 0.47,
        "Er": 0.92, "Tm": 0.12, "Yb": 0.92, "Lu": 0.09, "Hf": 0.71,
        "Ta": -0.12, "W": 0.65, "Re": 0.26, "Os": 1.35, "Ir": 1.32,
        "Pt": 1.62, "Au": 0.80, "Hg": 1.17, "Tl": 0.77, "Pb": 2.04,
        "Bi": 0.65, "Th": 0.06, "U": -0.54
    }
    
    # Сначала проверяем фотосферные данные
    if Name in photoshperic_dens:
        return 10**(photoshperic_dens[Name] - 12)
    
    # Если нет в фотосферных, проверяем метеоритные данные
    if Name in meteor_dens:
        print(f"Warning: Using meteoritic abundance for {Name} (solar abundance not available)")
        return 10**(meteor_dens[Name] - 12)
    
    # Если элемента нет ни в одном списке
    raise KeyError(f"Element {Name} not found in solar or meteoritic abundance tables")

def GetElementDense(m_body_model, element:ff.Nucleus):
    el_name_num = (element.name+str(element.A))
    if(el_name_num in m_body_model.columns):
        rho_e = m_body_model[el_name_num]
    elif(element.name in m_body_model.columns):
        rho_e = m_body_model[element.name]
    else:
        calibration = np.array(m_body_model["Fe"])/m_body_model["Fe"][1]
        try:
            abond = get_solar_abundance(element.name,element.A)*element.abondonce
            rho_e = calibration*abond*element.A
        except Exception as e:
            raise e
    
    return np.array(rho_e)*np.array(m_body_model['Rho']/element.A)



def CaptureNuc(
        ptype_in,ptype_out,
        m_grid,m_wimp_model : evdm.WimpModel,
        m_elements : list[ff.Nucleus],
        m_body_table,
        m_operator ,m_norm_operator,Nmk,seed,constrain = True,rpow=1,form_factor=None): 
    from concurrent.futures import ThreadPoolExecutor
    def heavy_computation(m_element):
        n_e = GetElementDense(m_body_table,m_element)
        m_wimp_params = m_wimp_model(ptype_in,ptype_out)
        Capt = evdm.Capture(m_grid)
        scat_mod = ff.ScatterModel(m_wimp_params,m_element,m_operator,m_norm_operator,2.06e-3,form_factor)
        CalcC = CaptureCalc(Capt,scat_mod,n_e,0.73e-3,0.423e-3,1.78e-3,
                         Nmk,r_pow= rpow, seed = seed,constrain = constrain )
        return (Capt,CalcC)
    
    Capt = evdm.Capture(m_grid)

    with ThreadPoolExecutor() as executor:
        results = list(executor.map(heavy_computation, m_elements))
    for (dCapt,CalcC) in results:
        Capt += dCapt
    return Capt

def ScatterNuc(
        ptype_in,ptype_out,
        m_grid,m_wimp_model : WimpModel,
        m_elements : list[ff.Nucleus],
        m_body_table,
        m_operator ,m_norm_operator,Nmk,seed,method = "naive",
        algol = 'naive',measure = (1,2),Nmk_traj = 10,ScatterMatrix = None,zero = 0,**kwargs)->evdm.Matrix: 
    
    print('start:')
    if(ScatterMatrix is None):
        ScatterMatrix = evdm.Matrix(m_grid)
    
    for m_element in m_elements:
        print(f'scatter for {m_element}')
        n_e = GetElementDense(m_body_table,m_element)
        m_wimp_params = m_wimp_model(ptype_in,ptype_out)
        form_factor = kwargs.get("form_factor",None)
        scat_mod = ff.ScatterModel(m_wimp_params,m_element,m_operator,m_norm_operator,2.06e-3,**kwargs)
        print()
        ScatterCalc(ScatterMatrix,scat_mod,n_e,Nmk,Nmk_traj = Nmk_traj,
                         seed=seed,method=method,measure = measure,algol=algol,**kwargs)
        if(np.isnan(ScatterMatrix.to_numpy().sum())):
          raise RuntimeError(f"ScatterMatrix have nan!, at ")
    return ScatterMatrix


def plot_distrib(distrib,ptypes = [0],mbounds = np.linspace(0, 1, 11,endpoint = True),norm = None):

    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import matplotlib.colors as pltcolors
    from matplotlib import cm
    from matplotlib.collections import LineCollection
    import matplotlib.image as mpimg

    if(len(ptypes) == 1):
        fig, _axes = plt.subplots(1,len(ptypes))
        axes = [_axes]
    else:
        fig, axes = plt.subplots(1,len(ptypes))
    
    for (ptype,ax) in zip(ptypes,axes):
        (X01,Y01,trs01,vals01) = distrib.plot(ptype,'dEdL')
        triang01 = mtri.Triangulation(X01, Y01, trs01)
        m_norm = vals01.max()
        
        bounds_01 = mbounds
        ax.tricontourf(triang01,(vals01/m_norm),levels = bounds_01,cmap='viridis')
def plot_distribs(
        distribs,ptype,mbounds = np.linspace(0, 1, 11,endpoint = True),
        norm = None,xlim = None,ylim = None):
    
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import matplotlib.colors as pltcolors
    from matplotlib import cm
    from matplotlib.collections import LineCollection
    import matplotlib.image as mpimg

    fig, _axes = plt.subplots(1,len(distribs))
    if(len(distribs) == 1):
        axes = [_axes]
    else:
        axes = _axes
    norms = None
    for (distr,ax) in zip(distribs,axes):
        (X01,Y01,trs01,vals01) = distr.plot(ptype,'dEdL')
        triang01 = mtri.Triangulation(X01, Y01, trs01)
        bounds_01 = mbounds

        if norm == None:
            m_norm = vals01.max()
        else:
            if(norms == None):
                norms = vals01.max()
            m_norm = norms

        ax.tricontourf(triang01,(vals01/m_norm),levels = bounds_01,cmap='viridis')
    if(xlim):
        plt.xlim(xlim)
    if(ylim):
        plt.ylim(ylim)



def get_deltas(distribs):
    dEdL = evdm.Metric(p_deg=1)
    return np.array([ [ dEdL(c,c1)for c1 in distribs] for c in distribs])
    
def ToPlotly(X_s,Y_s,Triangles,Values)->tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    shift = X_s.shape[0]
    X_new = np.concatenate([X_s,X_s])
    Y_new = np.concatenate([Y_s,Y_s])
    Values_new = np.concatenate([Values,Values*0])

    N_i = Triangles[:,0]
    M_i = Triangles[:,1]
    K_i = Triangles[:,2]

    Nd_i = Triangles[:,0] + shift
    Md_i = Triangles[:,1] + shift
    Kd_i = Triangles[:,2] + shift

    Triangles_new = np.concatenate(
        [
            Triangles,
            np.array([N_i, Nd_i,K_i]).transpose(),
            np.array([K_i,Nd_i, Kd_i]).transpose(),

            np.array([M_i,Md_i,K_i]).transpose(),
            np.array([K_i, Md_i,Kd_i]).transpose(),

            np.array([ M_i,Md_i,N_i]).transpose(),
            np.array([N_i,Md_i, Nd_i]).transpose(),
        ]
    )
    return (X_new,Y_new,Triangles_new,Values_new)
def MakeMesh(X_s,Y_s,Triangles,Values):
    import plotly.graph_objects as go
    return go.Mesh3d(
        x=X_s,
        y=Y_s,
        z=Values,
        i=Triangles[:,0],  # Индексы первого вершинного треугольника
        j=Triangles[:, 1],  # Индексы второго вершинного треугольника
        k=Triangles[:, 2],  # Индексы третьего вершинного треугольника
        intensity=Values,  # Цветовая карта на основе нормализованных значений
        colorscale='Viridis',   # Цветовая карта
        showscale=True          # Отображение цветовой шкалы
    )
def DistribPlot3D(distrib,mes = 'dEdL',space = "",ptype=0):
    (X_s,Y_s,Triangles,Values) = distrib.plot(ptype,mes,space)
    Values = Values/Values.max()
    return MakeMesh(*ToPlotly(X_s,Y_s,Triangles,Values))

def DistribPlot3D(distrib,mes = 'dEdL',space = "",ptype=0,max_value = None):
    (X_s,Y_s,Triangles,Values) = distrib.plot(ptype,mes,space)
    if(max_value is None):
        max_value = Values.max()
    Values = Values/max_value
    return MakeMesh(*ToPlotly(X_s,Y_s,Triangles,Values))

def Step1(S : np.ndarray,X : np.ndarray,tau):
    return X + tau*(S @ X)

def Step2(S : np.ndarray,X : np.ndarray,tau):
    X1 = (S @ X)
    return X + tau*X1 + (0.5*tau**2) * (S @ X1)

def S1_mat(S : np.ndarray,tau):
    return {'inv_1_plus_tau_S' : np.identity(S.shape[0]) - S*(0.5*tau)}

def Step2s(S1 : dict,S : np.ndarray,X : np.ndarray,tau):
    _tau = tau/2
    return np.linalg.solve(S1['inv_1_plus_tau_S'], X + _tau*(S @ X) )


def to_diag(scat_matr):
    matr = scat_matr.copy()
    matr.calc_diag()
    return matr.to_numpy()
def to_diag64(scat_matr,countEvap = False):
    matr = scat_matr.as_type('double')
    matr.calc_diag(countEvap)
    return matr.to_numpy()

def S1norm(S1 : np.ndarray):
    return np.abs(S1.diagonal()).max()

def Evolve1(Smat,X0,t,tau):
    if(tau * S1norm(Smat) > 1):
        raise Exception("Too big tau")
    N = int(t/tau + 1)
    tau_1 = t/N
    X1 = X0
    pbar = evdm.pbar.progress_bar_jupiter('steps')
    pbar.update(0,N-1)
    for i in range(N-1):
        X1 = Step1(Smat,X0,tau_1)
        pbar.update(i,N)
    return X1

def Evolve2(Smat,X0,t,tau):
    if(tau * S1norm(Smat) > 1):
        raise Exception("Too big tau")
    N = int(t/tau + 1)
    tau_1 = t/N
    X1 = X0
    pbar = evdm.pbar.progress_bar_jupiter('steps')
    pbar.update(0,N-1)
    for i in range(N):
        X1 = Step2(Smat,X1,tau_1)
        pbar.update(i,N-1)
    return X1

def Evolve2s(Smat,X0,t,tau,pbar = None):
    N = int(t/tau + 1)
    tau_1 = t/N
    S1_m = S1_mat(Smat,tau_1) 
    X1 = X0
    if(pbar == 'new'):
        pbar = evdm.pbar.progress_bar_jupiter('steps')
    if(pbar):
        pbar.update(0,N-1)
    for i in range(N):
        X1 = Step2s(S1_m,Smat,X1,tau_1)
        if(pbar):
            pbar.update(i,N-1)
    return X1

def AvarageEnergy(distrib : evdm.Distrib,ptype = -1):
    e0 = distrib.grid.body.phi[1][0]
    return distrib.avarage(lambda e,l:(e0+e),ptype)

def get_nv(m_svd,index):
    m_vector = np.abs(np.transpose(m_svd[2])[:,index])
    m_vector /= np.sum(m_vector)
    return m_vector

def get_best_svd(SVDs,grid,target_e):
    for i in range(1,5):
        m_distrib =  evdm.Distrib(grid,get_nv(SVDs,-i))
        me = AvarageEnergy(m_distrib)
        if ( (me-target_e)/(abs(me)+abs(target_e)) < 0.5):
            return m_distrib

def get_first_svd(SVDs,grid):
    m_distrib =  evdm.Distrib(grid,get_nv(SVDs,-1))
    me = AvarageEnergy(m_distrib)
    return m_distrib
        
def Richardson(X1,X2,X4):
    pow2p = (X1-X2)/(X2-X4)
    p = np.log(pow2p)/np.log(2)
    print("p = ", p)
    return (X1*X4-X2**2)/(X1-2*X2+X4)

def Richardson_p(X1,X2,X4):
    pow2p = (X1-X2)/(X2-X4)
    p = np.log(pow2p)/np.log(2)
    print("p = ", p)
    return ((X1*X4-X2**2)/(X1-2*X2+X4),p)

def plt_grid():
    from matplotlib import pyplot as plt
    plt.grid(which='major',linestyle='-')
    plt.grid(which='minor',linestyle='--')
    plt.minorticks_on()

def DelayedAssign(varname,getterfunc):
    from threading import Thread
    if(varname in globals()):
        globals().pop(varname)
    def target_func():
        globals()[varname] = getterfunc()
    T = Thread(target = target_func)
    T.start()


def set_vars(var_dict:dict):
    for key,item in var_dict.items():
        globals()[key] = item

