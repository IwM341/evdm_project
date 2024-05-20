from . import ff
from . import pbar
from ._cpp_lib.pyevdm import *
from .dm_model import *

K_to_GeV = 8.61732814974056e-14

 
def CaptureCalc(capt_vector,scat_mod : ff.ScatterModel,n_dense,Vbody,Vdisp,Nmk,r_pow = 1,weight = 1):
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
    """
    if(not scat_mod.Zero):
        wimp = scat_mod.wimp
        nuc = scat_mod.nucleus
        sc_event = ScatterEvent(n_dense,scat_mod.factor(),scat_mod.str_char())
        return CalcCaptureImpl(capt_vector,
            wimp.In,wimp.Out,wimp.mass,wimp.delta,nuc.mass,
            sc_event,Vbody,Vdisp,Nmk,r_pow,weight)
    else:
        return (0,0)
def ScatterCalc(sc_matrix,scatter_model : ff.ScatterModel,n_dense,
                Nmk,Nmk_traj = 1,count_evap = False,weight = 1,bar = None):
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
    Nmk_traj : int
        number of monte-carle steps on each trajectory.
    weight : float
        additional scale factor, default - 1.
    bar : any
        optional progress bar update function.
    """
    if(not scatter_model.Zero):
        wimp = scatter_model.wimp
        nuc = scatter_model.nucleus
        sc_event = ScatterEvent(n_dense,scatter_model.factor(),scatter_model.str_char())
        return CalcScatterImpl(sc_matrix,wimp.In,wimp.Out,wimp.mass,wimp.delta,nuc.mass,sc_event,
                            Nmk,Nmk_traj,count_evap,weight,bar)
