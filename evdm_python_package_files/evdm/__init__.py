from . import ff
from . import pbar
from ._cpp_lib.pyevdm import *
from .dm_model import *

K_to_GeV = 8.61732814974056e-14

 
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
    if(not scat_mod.Zero):
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
        'notherm', 'naive','soft','soft_tresh'
    Nmk_traj : int
        number of monte-carle steps on each trajectory.
    weight : float
        additional scale factor, default - 1.
    seed : int
        for generator
    bar : any
        optional progress bar update function.
    """
    if(not scatter_model.Zero):
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