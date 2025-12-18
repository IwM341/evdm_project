from . import ff
from . import pbar
from ._cpp_lib.pyevdm import *
from .dm_model import *
from typing import List
K_to_GeV = 8.61732814974056e-14
def MakeGrid(mBody,Evalues,Lvalues,ptypes = 1,value_type = 'float'):
    import numpy as np
    """
        Makes Grid from Evalues and Lvalues
        
        Parameters:
        -----------
        mBody: evdm.Body
            Body model
        Evalues : list
            Energy array from Emin to Emax
        Evalues : list | list[list]
            l array from 0 to 1 or arrays for each E bin
        ptypes : int
        value_type : str
            'float' or 'double'
    """

    if(value_type == 'float'):
        dtype = np.float32
    elif(value_type == 'double'):
        dtype = np.float64
    else:
        raise TypeError(f'Unknown dtype "{value_type}"')

    Evalues = np.array(Evalues,dtype=dtype)
    not_list = (isinstance(Lvalues,np.ndarray))

    def chek_array(X,name):
        if(not np.all(np.diff(X)> 0)):
            raise ValueError(f'{name} array is not monotonic')
    def check_l(Lv):
        chek_array(Lv,'Lvalues')
        if( Lv[0] > 1e-6 or Lv[0] < 0):
            raise ValueError('Lvalues[0] should be 0')
        if( Lv[-1] < 1-1e-5 or Lv[-1] > 1):
            raise ValueError('Lvalues[-1] should be 1')
    
    chek_array(Evalues,'Evalues')
    if(not_list):
        check_l(Lvalues)
        Lvalues = [Lvalues for i in range(len(Evalues)-1)]
    else:
        for i in range(len(Lvalues)):
            check_l(Lvalues[i])
        
        Lvalues = [Lvalues[i] for i in range(len(Evalues)-1)]
    


    
    
    
    mgrid_dict = {
        'type' : 'evdm.GridEL',
        'body_t' : 'float',
        'body':mBody,
        'grid_t' : value_type,
        'gtype' : 'CVV',
        'grid': {
            'Grid':{'size':ptypes},
            'InnerGrids' : {
                'size': ptypes,
                'value' : {
                    'Grid' : Evalues,
                    'InnerGrids' : Lvalues
                }
            }
        }
    }
    return GridEL(mBody,mgrid_dict)

def CaptureCalc(capt_vector,scat_mod : ff.ScatterModel,
                n_dense,Vbody,Vdisp,Vmax,Nmk,r_pow = 1,weight = 1,
                seed = 0,constrain = False):
    """
    Calculates capture, add event to capture vector,
	returns tuple (capture,mk sigma)

	Parameters:
    -----------
    capt_vector : Capture
        Capture histogramm.
	scatter_model : ScatterModel
        instance of ScatterModel class, containing all scatter info.
	n_dense : array
        relative concentration of nucleus (n_i(r)/<n_p>).
    Vbody : float
        speed of body relative to halo.
    Vdisp : float
        dispersion of DM speed in halo.
    Nmk : int
        number of monte-carlo steps.
    r_pow : float
        impact on r distribution: r = (xi)^(r_pow), where xi uniforemly distributed.
    weight : float
        additional scale factor, default is 1.
    seed : int
        for generator
    """
    if(not getattr(scat_mod,'Zero',False)):
        wimp = scat_mod.wimp
        nuc = scat_mod.nucleus
        sc_event = ScatterEvent(n_dense,scat_mod.factor(),scat_mod.str_char())
        return CalcCaptureImpl(capt_vector,
            wimp.In,wimp.Out,wimp.mass,wimp.delta,nuc.mass,
            sc_event,Vbody,Vdisp,Vmax,Nmk,r_pow,weight,seed,constrain)
    else:
        return (0,0)
def ScatterCalc(sc_matrix,scatter_model : ff.ScatterModel,n_dense,
                Nmk,**kwargs):
    """
    Calculates scatter matrix part, add event to matrix class
    
    Parameters:
    -----------
    sc_matrix : matrix
         matrix containing scatter and evaporation info
    scatter_model : ScatterModel
        containing all scatter info.
    n_dense : array
        relative concentration of nucleus (n_i(r)/<n_p>).
    Nmk : int
        number of monte-carle steps.
    method : string
        'notherm', 'naive','soft','soft_tresh','full'
    algol : string
        'naive', 'shift', 'diffuse'
    measure: tuple
        how E,L distributed in bin.
        should be a tuple (p,q). 
        E and l would be uniformly distributed by dE^pdl^q.
        default is (1,2)
    zero: float
        replace elements of scatter matrix sij < zero to zero
    Nmk_traj : int
        number of monte-carle steps on each trajectory.
        default is 1
    weight : float
        additional scale factor.
        default is 1.
    seed : int
        for generator
    bar : any
        optional progress bar update function.
    """
    if(not getattr(scatter_model,'Zero',False)):
        wimp = scatter_model.wimp
        nuc = scatter_model.nucleus
        sc_event = ScatterEvent(n_dense,scatter_model.factor(),scatter_model.str_char())
        return CalcScatterImpl(sc_matrix,wimp.In,wimp.Out,wimp.mass,wimp.delta,nuc.mass,sc_event,
                            Nmk,**kwargs)

class Group:
    def __init__(self,input_object : dict | object):
        '''
        construct gruop of evdm objects:
            if input_object is serialized dict of evdm.Group 
                restore dict to Group
            else 
                make object Group, which could be saved/loaded
        (it is nesessary, because grid and body shouldn't be copied by serialization)
        '''

        from ._serialize import _evdm_type_
        if(_evdm_type_(input_object)):
            self.obj = input_object
        elif('type' in input_object):
            if(input_object['type'] == "evdm.Group"):
                self._from_object(input_object)
            elif(type(input_object['type']) == str 
                and input_object['type'].startswith('evdm.')):
                raise TypeError("couldn't construct evdm.Group from serialized"+
                    "to dict evdm object except evdm.Group")
        elif(type(input_object) == dict):
            self.__dict__ = input_object
        else:
            if(hasattr(input_object,'__dict__')):
                self.__dict__ = input_object.__dict__
            else:
                self.obj = input_object

    def _from_object(self,m_dict ):
        from ._serialize import deep_move
        m_res = deep_move(m_dict)
        self.__dict__ = m_res

    def to_object(self):
        from ._serialize import serialize_Object
        m_grids = []
        m_bodies = []
        ret = serialize_Object(self,m_grids,m_bodies)
        ret['bodies'] = [bd.to_object() for bd in m_bodies]
        ret['grids'] = [
            {'value':gr['value'].to_object(),
             'bodyref':gr['bodyref'],
             'object_paths':gr['object_paths']} 
             for gr in m_grids]
        ret['type'] = 'evdm.Group'
        return ret


def _plot_grid(self,is_internal,ax = None):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    mlc = LineCollection(self.plot(is_internal))
    if(not ax):
        fig, ax = plt.subplots()
    ax.add_collection(mlc)
    ax.autoscale()

GridEL.pyplot = _plot_grid



def SolarAbundance(Name):
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


class SolarDensityGetter:
    def __init__(self,table_of_element):
        self.tab = table_of_element
    def __call__(self, element : ff.Nucleus):
        import numpy as np
        el_name_num = (element.name+str(element.A))
        if(el_name_num in self.tab.columns):
            rho_e = self.tab[el_name_num]
        elif(element.name in self.tab.columns):
            rho_e = self.tab[element.name]
        else:
            calibration = np.array(self.tab["Fe"])/self.tab["Fe"][1]
            try:
                abond = SolarAbundance(element.name)
                rho_e = calibration*abond*element.A
            except Exception as e:
                raise e
        return np.array(rho_e)*np.array(self.tab['Rho']/element.A)

class FF_Provider:
    def __init__(
            self,m_wimp : WimpModel,elements : List[ff.Nucleus],
            density_getter,Operator,NormOperator = None, NormVelocity = 1e-3):
        self.wimp = m_wimp
        self.elements = elements
        self.Operator = Operator
        self.NormOperator = Operator if NormOperator is None else NormOperator
        self.NormVelocity = NormVelocity
        self.density_getter = density_getter
    def construct(self,pin,pout):
        delta = self.wimp.delta(pin,pout)
        data = []
        for nuc in self.elements:
            scat_mod = ff.ScatterModel( 
                self.wimp(pin,pout),nuc,
                self.Operator,self.NormOperator,
                self.NormVelocity)
            sc_event = ScatterEvent(
                self.density_getter(nuc),
                scat_mod.factor(),
                scat_mod.str_char())
            
            data.append({
                'A':nuc.A,
                'Z':nuc.Z,
                'mN':nuc.mass,
                'mX' : self.wimp.mass + delta,
                'delta':delta,
                'event':sc_event
            })
        return data
