import numpy
from typing import Any, Callable, overload

class AnnPre:
    def __init__(self, grid: GridEL, Nmk: int, dtype: str = ..., bar: object = ...) -> None:
        """__init__(self: pyevdm.AnnPre, grid: pyevdm.GridEL, Nmk: int, dtype: str = 'float', bar: object = None) -> None

        constructor.

        Parameters:
        \tgrid : el grid.
        \tNmk : number of MK integratons per bin.
        \tdtype : type of values.
        \tbar : update progress bar function.

        """
    def add_to_matrix(self, matrix: Matrix, ptype0: int, ptype1: int, a0: float, av: float) -> None:
        """add_to_matrix(self: pyevdm.AnnPre, matrix: pyevdm.Matrix, ptype0: int, ptype1: int, a0: float, av: float) -> None

        build part of annihilation matrix of ptype0, ptype1, the matrix would by Aij += a0*a_0_ij+av*a_v_ijwhere annihilation determines by dN_i/dt = (n_{\\infty} \\sigma_{ann} v_esc) (A_{ij}  N_j) N_iParameters:
        \tmatrix : input matrix to build (modify)
        \tptype0, ptype1: scattering types
        \ta0 : coefficient of \\sigma_{ann} v_esc = const annihilation
        \tav : coefficient of \\sigma_{ann} v_esc = const * v^2 annihilation

        """

class Body:
    @overload
    def __init__(self, arg0: dict) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Body, arg0: dict) -> None

        constructs Body from saved dict



        2. __init__(self: pyevdm.Body, Rho: object, velocity: float, Temp: object = None, dtype: str = '') -> None

        constructs Body from Rho

        Parameters:
        ___________
        Rho : array of pho(x_i)
        \tor pair (F,size), where F - density function [0,1]->real., size - number of points
        size : int
        \tif Rho is function then size is the number of points.
        velocity : float
        \tvelocity of body relative to halo.
        dtype : string
        \tfloat or double.
        Temp : array
        \t optional temperature array of body.

        """
    @overload
    def __init__(self, Rho: object, velocity: float, Temp: object = ..., dtype: str = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Body, arg0: dict) -> None

        constructs Body from saved dict



        2. __init__(self: pyevdm.Body, Rho: object, velocity: float, Temp: object = None, dtype: str = '') -> None

        constructs Body from Rho

        Parameters:
        ___________
        Rho : array of pho(x_i)
        \tor pair (F,size), where F - density function [0,1]->real., size - number of points
        size : int
        \tif Rho is function then size is the number of points.
        velocity : float
        \tvelocity of body relative to halo.
        dtype : string
        \tfloat or double.
        Temp : array
        \t optional temperature array of body.

        """
    def setTemp(self, Temp: object = ...) -> None:
        """setTemp(self: pyevdm.Body, Temp: object = None) -> None"""
    def to_object(self) -> dict:
        """to_object(self: object) -> dict"""
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
        \tgrid, where distribution is created.
        dtype : string
        \t'float' or 'double'.
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution.

        2. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, values: dict) -> None

        constructor of Capture(Distrib) class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        values : array
        \tnumpy array of values.


        3. __init__(self: pyevdm.Capture, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None
        """
    @overload
    def __init__(self, ELGrid: GridEL, values: dict) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Capture(Distrib) class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        dtype : string
        \t'float' or 'double'.
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution.

        2. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, values: dict) -> None

        constructor of Capture(Distrib) class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        values : array
        \tnumpy array of values.


        3. __init__(self: pyevdm.Capture, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None
        """
    @overload
    def __init__(self, arg0: GridEL, arg1: numpy.ndarray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Capture(Distrib) class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        dtype : string
        \t'float' or 'double'.
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution.

        2. __init__(self: pyevdm.Capture, ELGrid: pyevdm.GridEL, values: dict) -> None

        constructor of Capture(Distrib) class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        values : array
        \tnumpy array of values.


        3. __init__(self: pyevdm.Capture, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None
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
    def to_object(self) -> dict:
        """to_object(self: object) -> dict"""
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
        \tgrid, where distribution is created.
        dtype : string
        \t'float' or 'double'.
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution.

        2. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, object: dict) -> None

        constructor of Distrib class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        object : dict representation of distrib

        3. __init__(self: pyevdm.Distrib, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None
        """
    @overload
    def __init__(self, ELGrid: GridEL, object: dict) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Distrib class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        dtype : string
        \t'float' or 'double'.
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution.

        2. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, object: dict) -> None

        constructor of Distrib class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        object : dict representation of distrib

        3. __init__(self: pyevdm.Distrib, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None
        """
    @overload
    def __init__(self, arg0: GridEL, arg1: numpy.ndarray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, dtype: str = 'float', init: object = None) -> None

        constructor of Distrib class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        dtype : string
        \t'float' or 'double'.
        init : lambda
        \tinitialiser fuinction (density)of 3 args: (ptype,e,l), default - None, meaning zero dirtribution.

        2. __init__(self: pyevdm.Distrib, ELGrid: pyevdm.GridEL, object: dict) -> None

        constructor of Distrib class

        Parameters:
        ___________
        ELGrid : GridEL
        \tgrid, where distribution is created.
        object : dict representation of distrib

        3. __init__(self: pyevdm.Distrib, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None
        """
    def Edistrib(self, ptype: int = ...) -> tuple:
        """Edistrib(self: pyevdm.Distrib, ptype: int = -1) -> tuple

        reduce distribution to energy distribution
        """
    def as_type(self, dtype: str) -> Distrib:
        """as_type(self: pyevdm.Distrib, dtype: str) -> pyevdm.Distrib

        creating Distrib with another dtype
        """
    def avarage(self, F: object, ptype: int = ...) -> float:
        """avarage(self: pyevdm.Distrib, F: object, ptype: int = -1) -> float

        calc avarage over distribution of F(e,l), where l \\in [0,1]
        """
    def copy(self) -> Distrib:
        """copy(self: pyevdm.Distrib) -> pyevdm.Distrib"""
    def count(self, ptype: int = ...) -> float:
        """count(self: pyevdm.Distrib, ptype: int = -1) -> float

        calculates number of particles
        """
    def plot(self, ptype: int, mes: str = ..., space: str = ...) -> tuple:
        """plot(self: pyevdm.Distrib, ptype: int, mes: str = 'dEdL', space: str = '') -> tuple

        returns tuple: (X,Y,triangles,values)
        where X,Y - arrays of vertises x and y coords,
        triangles - triangles (N,3) shape array of indexes
        values - values corresponding to vertexes.

        Parameters:
        ___________
        ptype : int
        \twimp type.
        mes : str
        \t measure of bins (standart is dEdL, no measure - '1' or '')
        space: str
        \tif space = 'latent', L would be from 0 to 1, otherwise - from 0 to Lmax(E)

        to plot using matplotlib type:
        \timport matplotlib.tri as mtri
        \tplt.tricontourf(triang,Z)
        \tplt.triplot(triang,color = 'black')
        """
    def plot2(self, ptype: int) -> tuple:
        '''plot2(self: pyevdm.Distrib, ptype: int) -> tuple

        returns tuple: (vertexes,triangles,values)
        where vertexes - (N,2) shape array of vert coords,
        triangles - triangles (N,3) shape array of indexes
        values - values corresponding to vertexes.
        \tParameters:
        ___________
        ptype : int
        \twimp type.

        to plot using plotly type:

        fig = figf.create_trisurf(vert[:,0],vert[:,1], np.log(np.abs(-vals)),
        \tcolormap="Portland",
        \tsimplices=trs,
        \ttitle="title")
        '''
    def rdens(self, ptypes: object = ..., rmin: float = ..., rmax: float = ..., Nr: int = ..., Nb: object = ..., rden: object = ..., bar: object = ...) -> tuple:
        """rdens(self: pyevdm.Distrib, ptypes: object = -1, rmin: float = 0, rmax: float = 1, Nr: int = 100, Nb: object = 10000, rden: object = None, bar: object = None) -> tuple

        git density function, i.e. d^N/d^3r

        Parameters:
        ___________
        ptypes : array or int
        \t considered ptypes.
        rmin : float
        \tmin radius in r distribution.
        rmax : float
        \tmax radius in r distribution.
        Nr : int
        \tnumber of points in r grid from rmin to rmax.
        Nb : int
        \tnumber of MK iterations per bin for integration, default - 10000.
        rden : function
        \toptional density function of r points location.
        bar : progressbar
        \toptional progress bar update function.
        """
    def to_numpy(self, ptype: int = ..., is_raw: bool = ...) -> numpy.ndarray:
        """to_numpy(self: object, ptype: int = -1, is_raw: bool = False) -> numpy.ndarray

        gives numpy array view to distribution

        Parameters:
        ___________
        ptype : int
        \twimp type, default -1, meaning all.
        is_raw : bool
        \tgives distribution array with padding.
        """
    @overload
    def to_object(self) -> dict:
        """to_object(self: object) -> dict

        return numpy array, so Distrib/Capture could be restoredfrom grid and this object:
        \tm_array = m_distrib.to_object()
        \tm_distrib1 = Distrib(m_distrib.grid(),m_array)

        """
    @overload
    def to_object(self) -> Any:
        """to_object(self: object) -> dict

        return numpy array, so Distrib/Capture could be restoredfrom grid and this object:
        \tm_array = m_distrib.to_object()
        \tm_distrib1 = Distrib(m_distrib.grid(),m_array)

        """
    @property
    def grid(self) -> GridEL: ...

class GridEL:
    @overload
    def __init__(self, body: Body, ptypes: int, Ne: int, Nl_func: object, RhoE: object = ..., RhoL: object = ..., dtype: str = ..., **kwargs) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.GridEL, body: pyevdm.Body, ptypes: int, Ne: int, Nl_func: object, RhoE: object = None, RhoL: object = None, dtype: str = 'float', **kwargs) -> None

        constructs E-L Grid

        Parameters:
        ___________
        body : Body
        \tinstance of Body class
        ptypes : int
        \tnumber of particle types.
        Ne : int
        \t number of bins of E axis.
        Nl_func : int | function
        \tif number then number of bins of L axis,else a function from [0.0,1.0]->int,indicating number of bins in L grid depending on e,where 0.0 correspond to Emin = phi(0), and 1.0 -- to Emax = 0.
        RhoE : function
        \t[optional] bin density of e axis RhoE : e in [0.0,1.0]->float.
        RhoL : function
        \t[optional] bin density of l axis RhoL : (e,l) in [0.0,1.0]x[0.0,1.0]->float.
        dtype : string
        \tfloat or double.

        2. __init__(self: pyevdm.GridEL, arg0: pyevdm.Body, arg1: dict) -> None
        """
    @overload
    def __init__(self, arg0: Body, arg1: dict) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.GridEL, body: pyevdm.Body, ptypes: int, Ne: int, Nl_func: object, RhoE: object = None, RhoL: object = None, dtype: str = 'float', **kwargs) -> None

        constructs E-L Grid

        Parameters:
        ___________
        body : Body
        \tinstance of Body class
        ptypes : int
        \tnumber of particle types.
        Ne : int
        \t number of bins of E axis.
        Nl_func : int | function
        \tif number then number of bins of L axis,else a function from [0.0,1.0]->int,indicating number of bins in L grid depending on e,where 0.0 correspond to Emin = phi(0), and 1.0 -- to Emax = 0.
        RhoE : function
        \t[optional] bin density of e axis RhoE : e in [0.0,1.0]->float.
        RhoL : function
        \t[optional] bin density of l axis RhoL : (e,l) in [0.0,1.0]x[0.0,1.0]->float.
        dtype : string
        \tfloat or double.

        2. __init__(self: pyevdm.GridEL, arg0: pyevdm.Body, arg1: dict) -> None
        """
    def Epoints(self) -> numpy.ndarray:
        """Epoints(self: object) -> numpy.ndarray"""
    def LE(self) -> tuple:
        """LE(self: pyevdm.GridEL) -> tuple

        return tuple (E array,L(E) array)
        """
    def Lpoints(self, index: int) -> numpy.ndarray:
        """Lpoints(self: object, index: int) -> numpy.ndarray"""
    def le_functor(self) -> Callable[[float], float]:
        """le_functor(self: pyevdm.GridEL) -> Callable[[float], float]

        returns functor of lmax(e) function
        """
    def lpoints(self, index: int) -> numpy.ndarray:
        """lpoints(self: object, index: int) -> numpy.ndarray"""
    def plot(self, is_internal: bool = ...) -> numpy.ndarray:
        """plot(self: pyevdm.GridEL, is_internal: bool = False) -> numpy.ndarray

        gives array of shape (N,2,2):
        [ [[x_start_i,y_start_i],[x_end_i,y_end_i]],...]
        arrays could be plottted with matplotlib:
        lc = matplotlib.collections.LineCollection(result of this function)
        fig, ax = plt.subplots()
        ax.add_collection(lc)
        if is_internal true, then l from 0 to 1, else from 0, l(e).

        Parameters:
        ___________
        is_internal : bool
        \tif true, plot in internal representation.
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
    def refine(self, Er: int, Lr: int) -> GridEL:
        """refine(self: pyevdm.GridEL, Er: int, Lr: int) -> pyevdm.GridEL

        make refeined grid by dividing E step by Erand L step by Lr
        """
    def rmp(self, e: float, l: float) -> tuple:
        """rmp(self: pyevdm.GridEL, e: float, l: float) -> tuple

        returns (rmin,rmax)(e,l_undim)
        """
    def to_object(self) -> dict:
        """to_object(self: object) -> dict"""
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
    def full_size(self) -> int: ...
    @property
    def ptypes(self) -> int: ...
    @property
    def size(self) -> int: ...

class Matrix:
    @overload
    def __init__(self, ELGrid: GridEL, dtype: str = ...) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Matrix, ELGrid: pyevdm.GridEL, dtype: str = 'float') -> None

        constructor of Matrix class

        Parameters:
        ___________
        ELGrid : GridEL
        \tCreated EL Grid.
        dtype : string
        \tfloat or double.


        2. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray, arg2: numpy.ndarray) -> None

        3. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None

        4. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: dict) -> None
        """
    @overload
    def __init__(self, arg0: GridEL, arg1: numpy.ndarray, arg2: numpy.ndarray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Matrix, ELGrid: pyevdm.GridEL, dtype: str = 'float') -> None

        constructor of Matrix class

        Parameters:
        ___________
        ELGrid : GridEL
        \tCreated EL Grid.
        dtype : string
        \tfloat or double.


        2. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray, arg2: numpy.ndarray) -> None

        3. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None

        4. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: dict) -> None
        """
    @overload
    def __init__(self, arg0: GridEL, arg1: numpy.ndarray) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Matrix, ELGrid: pyevdm.GridEL, dtype: str = 'float') -> None

        constructor of Matrix class

        Parameters:
        ___________
        ELGrid : GridEL
        \tCreated EL Grid.
        dtype : string
        \tfloat or double.


        2. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray, arg2: numpy.ndarray) -> None

        3. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None

        4. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: dict) -> None
        """
    @overload
    def __init__(self, arg0: GridEL, arg1: dict) -> None:
        """__init__(*args, **kwargs)
        Overloaded function.

        1. __init__(self: pyevdm.Matrix, ELGrid: pyevdm.GridEL, dtype: str = 'float') -> None

        constructor of Matrix class

        Parameters:
        ___________
        ELGrid : GridEL
        \tCreated EL Grid.
        dtype : string
        \tfloat or double.


        2. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray, arg2: numpy.ndarray) -> None

        3. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: numpy.ndarray) -> None

        4. __init__(self: pyevdm.Matrix, arg0: pyevdm.GridEL, arg1: dict) -> None
        """
    def (self, *args, **kwargs):
        """(self: pyevdm.Matrix) -> pyevdm.Distrib

        make new evolution matrix (S_{ii} = -\\sum_j{S_{ji}})
        """
    def append(self, arg0: Matrix) -> Matrix:
        """append(self: pyevdm.Matrix, arg0: pyevdm.Matrix) -> pyevdm.Matrix

        adds to capture extra events capture
        """
    def as_type(self, dtype: str) -> Matrix:
        """as_type(self: pyevdm.Matrix, dtype: str) -> pyevdm.Matrix

        creating Matrix with another dtype
        """
    def calc_diag(self) -> None:
        """calc_diag(self: pyevdm.Matrix) -> None

        make diag values -summ of scatter probabilities
        """
    def copy(self) -> Matrix:
        """copy(self: pyevdm.Matrix) -> pyevdm.Matrix"""
    def diag_distrib(self) -> Distrib:
        """diag_distrib(self: pyevdm.Matrix) -> pyevdm.Distrib

        givesdistribution of diag scatter matrix.

        Parameters:
        ___________
        ptype : int
        \twimp type, default -1, meaning all.
        """
    def to_numpy(self, ptype_in: int = ..., ptype_out: int = ..., is_raw: bool = ...) -> numpy.ndarray:
        """to_numpy(self: object, ptype_in: int = -1, ptype_out: int = -1, is_raw: bool = False) -> numpy.ndarray

        gives numpy array view to scatter matrix.

        Parameters:
        ___________
        ptype_in : int
        \tin wimp type, default -1, meaning all.
        ptype_out : int
        \tout wimp type, default -1, meaning all.
        is_raw : bool
        \tgives raw array with padding.
        """
    def to_numpy_diag(self, ptype: int = ..., is_raw: bool = ...) -> numpy.ndarray:
        """to_numpy_diag(self: object, ptype: int = -1, is_raw: bool = False) -> numpy.ndarray

        gives numpy array view to diag scatter matrix.

        Parameters:
        ___________
        ptype : int
        \twimp type, default -1, meaning all.
        is_raw : bool
        \tgives raw array with padding.
        """
    def to_object(self) -> dict:
        """to_object(self: object) -> dict

        serialization into object dict
        """
    def __add__(self, arg0: Matrix) -> Matrix:
        """__add__(self: pyevdm.Matrix, arg0: pyevdm.Matrix) -> pyevdm.Matrix"""
    @property
    def evap_histo(self) -> Distrib: ...
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
        \t'dEdL2' - measure differential is dE and dL^2 = 2LdL.
        p_deg : float
        \tvalue of p parameter of Lp norm.
        ptype : int
        \tindex of particles to compare, if -1, than summ of all particles.
        """
    def __call__(self, arg0: Distrib, arg1: Distrib) -> float:
        """__call__(self: pyevdm.Metric, arg0: pyevdm.Distrib, arg1: pyevdm.Distrib) -> float

        compare two distribs
        """

class ScatterEvent:
    def __init__(self, n_e: object, sf: ScatterFactor, name: str = ..., unique: bool = ...) -> None:
        """__init__(self: pyevdm.ScatterEvent, n_e: object, sf: pyevdm.ScatterFactor, name: str = '__unnamed__', unique: bool = False) -> None

        creating ScatterEvent.

        Parameters:
        ___________
        n_e : array 
        \t relative to n_p concentration of targets, where n_p - avarage concentration of nuclons.
        sf : ScatterFactor
        \tscatter factor class instance.
        name : string
        \toptional name of event.
        unique : bool
        \tif true, sum of capture with same event names will be blocked.
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

def CalcCaptureImpl(capt_vector: Capture, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: ScatterEvent, Vbody: float, Vdisp: float, Nmk: int, r_pow: float = ..., weight: float = ..., seed: int = ...) -> tuple:
    """CalcCaptureImpl(capt_vector: pyevdm.Capture, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: pyevdm.ScatterEvent, Vbody: float, Vdisp: float, Nmk: int, r_pow: float = 2.0, weight: float = 1.0, seed: int = 2305843009213693945) -> tuple

    Calculates capture, add event to capture vector,returns tuple (capture, sigma)

    Parameters:
    ___________
    capt_vector : Capture
    \tCapture histogramm.
    ptype_in : int
    \tindex of input particle type.
    ptype_out : int
    \tindex of output particle type.
    m_wimp : float
    \tdm particle mass, GeV.
    delta : float
    \tdelta mass GeV: output mass - input mass.
    m_nuc : float
    \tnuclear mass, GeV.
    sc_event : ScatterEvent.
    Vbody : float
    \tspeed of body relative to halo.
    Vdisp : float
    \tdispersion of DM speed in halo.
    Nmk : int
    \tnumber of monte-carle steps.
    r_pow :float
    \timpact on r distribution: r = (xi)^(r_pow), where xi uniforemly distributed.
    weight :float
    \t[optional] scale factor, default - 1.
    \tseed : int
    \t[optional] for random generator
    """
def CalcScatterImpl(matrix: Matrix, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: ScatterEvent, Nmk: object, **kwargs) -> None:
    """CalcScatterImpl(matrix: pyevdm.Matrix, ptype_in: int, ptype_out: int, m_wimp: float, delta: float, m_nuc: float, sc_event: pyevdm.ScatterEvent, Nmk: object, **kwargs) -> None

    Calculates scatter matrix part, add event to matrix class

    Parameters:
    ___________
    matrix : Matrix
    \tScatter matrix histo.
    ptype_in : int
    \tindex of input particle type.
    ptype_out : int
    \tindex of output particle type.
    m_wimp : float
    \tdm particle mass, GeV.
    delta : float
    \tdelta mass GeV: output mass - input mass.
    m_nuc : float
    \tnuclear mass, GeV.
    sc_event : ScatterEvent.
    Nmk : int | function | vector
    \tnumber of monte-carle steps.
    \tMay depends on (e,l) or on (e0,e1,l0,l1)
    method : str
    \t method of generating therm velocity of nuclei.can be: 
    \t 'notherm', 'naive','soft' (more probability of high velocities), 'soft_tresh' (same as soft, but considering inelastic treshold)
    Nmk_traj: int
    \tnumber of monte-carle steps on each trajectory.
    weight : float
    \t[optional] scale factor, default - 1.
    bar : object
    \t[optional] progress bar update function.
    """
def func_factor(func: object) -> ScatterFactor:
    """func_factor(func: object) -> pyevdm.ScatterFactor

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
