import numpy
from typing import Callable, overload

class Body:
    def __init__(self, Rho: object, size: int | None = ..., velocity: float = ..., Temp: object = ..., dtype: str = ...) -> None:
        """__init__(self: pyevdm.Body, Rho: object, size: Optional[int] = None, velocity: float = 0.0007, Temp: object = None, dtype: str = '') -> None

        constructs Body from Rho
        Parameters:
        ___________
        Rho : array or F
        \tvalues of mass density or density function [0,1]->real
        size : int
        \tif Rho is function then size is the number of points
        velocity : float
        \tvelocity of body relative to halo
        dtype : string
        \tfloat or double
        Temp : array
        \t optional temperature array of body

        """
    def setTemp(self, Temp: object = ...) -> None:
        """setTemp(self: pyevdm.Body, Temp: object = None) -> None"""
    @property
    def M(self) -> tuple: ...
    @property
    def Q(self) -> tuple: ...
    @property
    def dtype(self) -> str: ...
    @property
    def phi(self) -> tuple: ...
    @property
    def rho(self) -> tuple: ...
    @property
    def size(self) -> int: ...

class Capture(Distrib):
    @overload
    def __init__(self, ELGrid: GridEL, dtype: str = ..., init: object = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Capture(Distrib) class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        dtype : string
        \t'float' or 'double'
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution

        2. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, values: numpy.ndarray) -> None

        constructor of Capture(Distrib) class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        values : array
        \tnumpy array of values

        """
    @overload
    def __init__(self, ELGrid: GridEL, values: numpy.ndarray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Capture(Distrib) class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        dtype : string
        \t'float' or 'double'
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution

        2. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, values: numpy.ndarray) -> None

        constructor of Capture(Distrib) class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        values : array
        \tnumpy array of values

        """
    def append(self, second_capture: Capture) -> Capture:
        """append(self: pyevdm.Capture, second_capture: pyevdm.Capture) -> pyevdm.Capture

        adds to capture extra events capture
        """
    def as_type(self, dtype: str) -> Distrib:
        """as_type(self: pyevdm.Capture, dtype: str) -> pyevdm.Distrib

        creating Capture with another dtype
        """
    def copy(self) -> Capture:
        """copy(self: pyevdm.Capture) -> pyevdm.Capture"""
    def __add__(self, arg0: Capture) -> Capture:
        """__add__(self: pyevdm.Capture, arg0: pyevdm.Capture) -> pyevdm.Capture"""
    @property
    def events(self) -> list[scatter_event_info]: ...

class Distrib:
    @overload
    def __init__(self, ELGrid: GridEL, dtype: str = ..., init: object = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Distrib class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        dtype : string
        \t'float' or 'double'
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution

        2. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, values: numpy.ndarray) -> None

        constructor of Distrib class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        values : array
        \tnumpy array of values

        """
    @overload
    def __init__(self, ELGrid: GridEL, values: numpy.ndarray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Distrib class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        dtype : string
        \t'float' or 'double'
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution

        2. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, values: numpy.ndarray) -> None

        constructor of Distrib class
        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created
        values : array
        \tnumpy array of values

        """
    def as_type(self, dtype: str) -> Distrib:
        """as_type(self: pyevdm.Distrib, dtype: str) -> pyevdm.Distrib

        creating Distrib with another dtype
        """
    def copy(self) -> Distrib:
        """copy(self: pyevdm.Distrib) -> pyevdm.Distrib"""
    def count(self, ptype: int = ...) -> float:
        """count(self: pyevdm.Distrib, ptype: int = -1) -> float

        calculates number of particles
        """
    def plot(self, ptype: int) -> tuple:
        """plot(self: pyevdm.Distrib, ptype: int) -> tuple

        returns tuple: (vertexes,triangles,values)
        where vertexes - (N,2) shape array of vert coords,
        triangles - triangles (N,3) shape array of indexes
        values - values corresponding to vertexesParameters:
        ___________
        ptype : int
        \twimp type
        """
    def rdens(self, ptypes: object = ..., rmin: float = ..., rmax: float = ..., Nr: int = ..., Nb: int = ..., rden: object = ..., bar: object = ...) -> tuple:
        """rdens(self: pyevdm.Distrib, ptypes: object = -1, rmin: float = 0, rmax: float = 1, Nr: int = 100, Nb: int = 10000, rden: object = None, bar: object = None) -> tuple

        return density function, i.e. d^N/d^3r
        Parameters:
        ___________
        ptypes : array or int
        \t considered ptypes
        rmin : float
        \tmin radius in r distribution
        rmax : float
        \tmax radius in r distribution
        Nr : int
        \tnumber of points in r grid from rmin to rmax
        Nb : int
        \tnumber of MK iterations per bin for integration, default - 10000
        rden : function
        \toptional density function of r points location
        bar : progressbar
        \toptional progress bar update function
        """
    def to_numpy(self, ptype: int = ..., is_raw: bool = ...) -> numpy.ndarray:
        """to_numpy(self: object, ptype: int = -1, is_raw: bool = False) -> numpy.ndarray

        gives numpy array view to distribution
        Parameters:
        ___________
        ptype : int
        \twimp type, default -1, meaning all
        is_raw : bool
        \tgives distribution array with padding
        """
    @property
    def grid(self) -> GridEL: ...

class GridEL:
    def __init__(self, body: Body, ptypes: int, Ne: int, Nl_func: object, RhoE: object = ..., RhoL: object = ..., dtype: str = ..., **kwargs) -> None:
        """__init__(self: pyevdm.GridEL, body: pyevdm.Body, ptypes: int, Ne: int, Nl_func: object, RhoE: object = None, RhoL: object = None, dtype: str = 'float', **kwargs) -> None

        constructs E-L Grid
        Parameters:
        ___________
        body : Body
        \tinstance of Body class
        ptypes : int
        \tnumber of particle types
        Ne : int
        \t number of bins of E axis
        Nl_func : int or function
        \tif number then number of bins of L axis,else a function from [0.0,1.0]->int,indicating number of bins in L grid depending on e,where 0.0 correspond to Emin = phi(0), and 1.0 -- to Emax = 0
        RhoE : function
        \t [optional] bin density of e axis RhoE : e in [0.0,1.0]->float
        RhoL : function
        \t [optional] bin density of l axis RhoL : (e,l) in [0.0,1.0]x[0.0,1.0]->float
        dtype : string
        \tfloat or double
        """
    def LE(self) -> tuple:
        """LE(self: pyevdm.GridEL) -> tuple

        return tuple (E array,L(E) array)
        """
    def le_functor(self) -> Callable[[float], float]:
        """le_functor(self: pyevdm.GridEL) -> Callable[[float], float]

        returns functor of lmax(e) function
        """
    def plot(self, is_internal: bool = ...) -> numpy.ndarray:
        """plot(self: pyevdm.GridEL, is_internal: bool = False) -> numpy.ndarray

        return array of shape (N,2,2):
        [ [[x_start_i,y_start_i],[x_end_i,y_end_i]],...]
        arrays could be plottted with matplotlib:
        lc = matplotlib.collections.LineCollection(result of this function)
        fig, ax = plt.subplots()
        ax.add_collection(lc)
        if is_internal true, then l from 0 to 1, else from 0, l(e)Parameters:
        ___________
        is_internal : bool
        if true, plot in internal representation
        """
    def plot_trajp(self, pname: object = ..., ld: bool = ...) -> list:
        """plot_trajp(self: pyevdm.GridEL, pname: object = 'T', ld: bool = True) -> list

        returns list L of dicts able to be plot in plotly of someparametrs of trajectory (Tin,Tout,Tfull,rmin,rmax,theta)
        pname name of plotted parameter: 'Tin', 'Tout' or 'T'
        (if ld = False, L would be plotted in hidden space)
        to plot list in plotly use go.Figure(data = [go.Surface(P) for P in RetList])
        """
    def print_debug(self) -> str:
        """print_debug(self: pyevdm.GridEL) -> str"""
    def rmp(self, e: float, l: float) -> tuple:
        """rmp(self: pyevdm.GridEL, e: float, l: float) -> tuple

        returns (rmin,rmax)(e,l_undim)
        """
    def traj_functor(self, pname: object = ...) -> Callable[[float, float], float]:
        """traj_functor(self: pyevdm.GridEL, pname: object = 'T') -> Callable[[float, float], float]

        return functor of required param:
        T0,T1,Tin0,Tin1,Tout,rmin,rmax,theta,Tinth
        """
    @property
    def body(self) -> Body: ...
    @property
    def dtype(self) -> str: ...
    @property
    def ptypes(self) -> int: ...
    @property
    def size(self) -> int: ...

class Matrix:
    def __init__(self, ELGrid: GridEL, dtype: str = ...) -> None:
        """__init__(self: pyevdm.Matrix, ELGrid: pyevdm.GridEL, dtype: str = 'float') -> None

        constructor of Matrix class
        Parameters:
        ___________
        ELGrid : GridEL
        \tCreated EL Grid
        dtype : string
        \tfloat or double

        """
    def as_type(self, dtype: str) -> Matrix:
        """as_type(self: pyevdm.Matrix, dtype: str) -> pyevdm.Matrix

        creating Matrix with another dtype
        """
    def to_numpy(self, ptype_in: int = ..., ptype_out: int = ..., is_raw: bool = ...) -> numpy.ndarray:
        """to_numpy(self: object, ptype_in: int = -1, ptype_out: int = -1, is_raw: bool = False) -> numpy.ndarray

        gives numpy array view to scatter matrix
        Parameters:
        ___________
        ptype_in : int
        \tin wimp type, default -1, meaning all
        ptype_out : int
        \tout wimp type, default -1, meaning all
        is_raw : bool
        \tgives raw array with padding
        """
    @property
    def events(self) -> list[scatter_event_info]: ...
    @property
    def grid(self) -> GridEL: ...

class Metric:
    def __init__(self, Mes: str = ..., p_deg: float = ..., ptype: int = ...) -> None:
        """__init__(self: pyevdm.Metric, Mes: str = 'dEdL', p_deg: float = 1, ptype: int = -1) -> None

        Creation Lp Metric of distributions
        use: M(D1,D2), constructor arguments:
        Parameters:
        ___________
        Mes : Measure
        \t'dEdl' - measure differential is dE and dl, where l = L/Lmax(E)
        \t'dEdL' - measure differential is dE and dL
        \t'dEdL2' - measure differential is dE and dL^2 = 2LdL
        p_deg : float
        \tvalue of p parameter of Lp norm
        ptype : int
        \tindex of particles to compare, if -1, than summ of all particles
        """
    def __call__(self, arg0: Distrib, arg1: Distrib) -> float:
        """__call__(self: pyevdm.Metric, arg0: pyevdm.Distrib, arg1: pyevdm.Distrib) -> float

        compare two distribs
        """

class ScatterEvent:
    def __init__(self, n_e: object, sf: ScatterFactor, name: str = ..., unique: bool = ...) -> None:
        """__init__(self: pyevdm.ScatterEvent, n_e: object, sf: pyevdm.ScatterFactor, name: str = '__unnamed__', unique: bool = False) -> None

        creating ScatterEvent
        Parameters:
        ___________
        n_e : array 
        \t relative to n_p concentration of targets, where n_p - avarage concentration of nuclons
        sf : ScatterFactor
        \tscatter factor class instancename : string
        \toptional name of event
        unique : bool
        \tif true, sum of capture with same event names will be blocked
        """

class ScatterFactor:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""
    def __call__(self, y: float, v2T: float = ...) -> float:
        """__call__(self: pyevdm.ScatterFactor, y: float, v2T: float = 0) -> float

        evaluates form factor for demonstration
        """

class scatter_event_info:
    def __init__(self, *args, **kwargs) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

def CalcCaptureImpl(capt_vector: Capture, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: ScatterEvent, Vbody: float, Vdisp: float, Nmk: int, r_pow: float = ..., weight: float = ...) -> tuple:
    """CalcCaptureImpl(capt_vector: pyevdm.Capture, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: pyevdm.ScatterEvent, Vbody: float, Vdisp: float, Nmk: int, r_pow: float = 2.0, weight: float = 1.0) -> tuple

    Calculates capture, add event to capture vector,returns tuple (capture, sigma)
    Parameters:
    ___________
    capt_vector : Capture
    \tCapture histogramm
    ptype_in : int
    \tindex of input particle type
    ptype_out : int
    \tindex of output particle type
    m_wimp : float
    \tdm particle mass, GeV
    delta : float
    \tdelta mass GeV: output mass - input mass
    m_nuc : float
    \tnuclear mass, GeV
    sc_event : ScatterEvent
    Vbody : float
    \tspeed of body relative to halo
    Vdisp : float
    \tdispersion of DM speed in halo
    Nmk : int
    \tnumber of monte-carle steps
    r_pow :float
    \timpact on r distribution: r = (xi)^(r_pow), where xi uniforemly distributed
    weight :float
    \t[optional] scale factor, default - 1
    """
def CalcScatterImpl(matrix: Matrix, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: ScatterEvent, Nmk: int, Nmk_traj: int = ..., count_evap: bool = ..., weight: float = ..., bar: object = ...) -> None:
    """CalcScatterImpl(matrix: pyevdm.Matrix, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: pyevdm.ScatterEvent, Nmk: int, Nmk_traj: int = 10, count_evap: bool = False, weight: float = 1, bar: object = None) -> None

    Calculates scatter matrix part, add event to matrix class
    Parameters:
    ___________
    matrix : Matrix
    \tScatter matrix histo
    ptype_in : int
    \tindex of input particle type
    ptype_out : int
    \tindex of output particle type
    m_wimp : float
    \tdm particle mass, GeV
    delta : float
    \tdelta mass GeV: output mass - input mass
    m_nuc : float
    \tnuclear mass, GeV
    sc_event : ScatterEvent
    Nmk : int
    \tnumber of monte-carle steps
    Nmk_traj: int
    \tnumber of monte-carle steps on each trajectory
    weight : float
    \t[optional] scale factor, default - 1
    bar : object
    \t[optional] progress bar update function
    """
def func_factr(func: object) -> ScatterFactor:
    """func_factr(func: object) -> pyevdm.ScatterFactor

    create elastic form fractor from function
    Input should contain pointer to funcion with signature float ScatterFunc(float q_2,float v2T)
    where q_2 = q^2 - transferred momentum in GeV,v2T - squared norm of inelastic transfer velocity
    Parameters:
    ___________
    func : function
    \tfloat ScatterFunc(float q_2,float v2T)
    """
@overload
def qexp_factor(b: float, y_inv: bool, P_0: numpy.ndarray[numpy.float32], P_V: numpy.ndarray[numpy.float32]) -> ScatterFactor:
    """qexp_factor(*args, **kwargs)
    Overloaded function.

    1. qexp_factor(b: float, y_inv: bool, P_0: numpy.ndarray[numpy.float32], P_V: numpy.ndarray[numpy.float32]) -> pyevdm.ScatterFactor

    create exponential form factor exp(-2y)(p0_i y^i + pv_i y^i*v_{perp}^2  )/y^t
    where y = b^2*q^2/4
    Parameters:
    ___________
    b : float
    \t size of nuclei in GeVn
    y_inv : bool
    \tif false than t = 0 and y^t = 1, else t = 1
    P_0 : array
    \tcoefficients of y^i
    P_V : array
    \toptional coefficints of y^i v_perp^2

    2. qexp_factor(b: float, y_inv: bool, P_0: numpy.ndarray[numpy.float32]) -> pyevdm.ScatterFactor

    create exponential form factor exp(-2y)(p0_i y^i)/y^t
    where y = b^2*q^2/4
    Parameters:
    ___________
    b : float
    \t size of nuclei in GeVn
    y_inv : bool
    \tif false than t = 0 and y^t = 1, else t = 1
    P_0 : array
    \tcoefficients of y^i
    """
@overload
def qexp_factor(b: float, y_inv: bool, P_0: numpy.ndarray[numpy.float32]) -> ScatterFactor:
    """qexp_factor(*args, **kwargs)
    Overloaded function.

    1. qexp_factor(b: float, y_inv: bool, P_0: numpy.ndarray[numpy.float32], P_V: numpy.ndarray[numpy.float32]) -> pyevdm.ScatterFactor

    create exponential form factor exp(-2y)(p0_i y^i + pv_i y^i*v_{perp}^2  )/y^t
    where y = b^2*q^2/4
    Parameters:
    ___________
    b : float
    \t size of nuclei in GeVn
    y_inv : bool
    \tif false than t = 0 and y^t = 1, else t = 1
    P_0 : array
    \tcoefficients of y^i
    P_V : array
    \toptional coefficints of y^i v_perp^2

    2. qexp_factor(b: float, y_inv: bool, P_0: numpy.ndarray[numpy.float32]) -> pyevdm.ScatterFactor

    create exponential form factor exp(-2y)(p0_i y^i)/y^t
    where y = b^2*q^2/4
    Parameters:
    ___________
    b : float
    \t size of nuclei in GeVn
    y_inv : bool
    \tif false than t = 0 and y^t = 1, else t = 1
    P_0 : array
    \tcoefficients of y^i
    """
