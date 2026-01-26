import math
import sympy
from sympy.parsing.sympy_parser import parse_expr
from ._form_factors import _form_factors as _ff_bind
from ._form_factors import _latex_names as _ltx_nms
from ._form_factors import _symbols as _symv 
from ._form_factors import _nuc_info as _nuc_info
from . import dm_model as dm_m
from ._cpp_lib import pyevdm as _evdm



class O:
    _Expression_Latex_O_Array = _ltx_nms._make_latex_array()
    """
    Non relativistic operators
    """
    def __init__(self,n : int,t : int = 0):
        """
        symbol of non relativistic operator
        :parametr n: number of operator
        :parametr t: isospin number 0 or 1
        """
        self.symbol = _symv._O_impl(n,t)
        self.n = n
        self.t = t
    def __str__(self):
        return self.symbol.__str__()
    def latex(self):
        return 'O_{'+str(self.n)+'}^'+str(self.t) +" = "+ O._Expression_Latex_O_Array[self.n-1]
    def __add__(self,other):
        if(isinstance(other,O)):
            return self.symbol + other.symbol
        else:
            return self.symbol + other
    def __radd__(self,other):
        return self.__add__(other)
    def __mul__(self,other):
        if(isinstance(other,O)):
            return self.symbol * other.symbol
        else:
            return self.symbol * other
    def __rmul__(self,other):
        return self.__mul__(other)

class Nucleus:
    """
    class representing Nucleus
    """
    def __init__(self,A,Z = None,name = None):
        """
        A - number of nuclons
        Z - charge
        name - could be instead of charge
        """
        if(Z != None):
            self.name = _nuc_info._mindeleev_get_name(Z)
            self.Z = Z
        else:
            self.Z = _nuc_info._mindeleev_get_Z(name)
            self.name = name
        self.A = A
        self.spin = _nuc_info._get_spin(Z,A)
        self.mass = A*0.938
        try:
            self.factors = _ff_bind.FormFactors[self.name][self.A]
        except:
            print(f"Warning: element {(Z,name,A)} hasn't form factors")
            self.factors = None

        if (A == 1):
            self.b = 1e-4
        else:
            self.b = float(sympy.sqrt(41.467/(45*A**(-1.0/3) - 25*A**(-2.0/3) ))*5.067730716548338)
    def __repr__(self):
        return f"Nucleus(A = {self.A}, Z = {self.Z})"
    def __str__(self):
        return self.__repr__()
    Hydrogen = None
Nucleus.Hydrogen = Nucleus(1,1,'H')

Q2 = _symv.Q2
"""
symbol of q^2 operator (transition momentum)
"""
Delta = _symv.Delta
"""
symbol of mass difference
""" 
M_N = _symv.M_N
"""
symbol of nucleus mass
"""
M_X = _symv.M_X
"""
symbol of WIMP mass
""" 

MU_NX = M_X*M_N/(M_X+M_N)
"""
symbol of WIMP nucleus reduced mass
"""
class FormFactor_Helm:
    def __init__(self,
            wimp_pars: dm_m.WimpScatterParams,
            nucleus: Nucleus,
            R = None,
            S2 = None
        ):
        self.Zero = False 
        self.wimp = wimp_pars
        self.nucleus = nucleus
        fermi_GeV= 5
        s2 : float =  S2 or (fermi_GeV*0.9)**2
        b = (1.23*nucleus.A**(1.0/3)-0.6)*fermi_GeV
        a = 0.52*fermi_GeV
        R : float = R or math.sqrt(b*b+7*math.pi**2*a*a/3-5*s2)
        mp = Nucleus.Hydrogen.mass
        cns_fac : float = nucleus.A**4*( (wimp_pars.mass+mp)/(wimp_pars.mass+nucleus.A*mp) )**2

        if(cns_fac == 0):
            self.Zero = True 
        
        self.s2 = s2
        self.R = R  
        self.cfac = cns_fac
        def myBessel(x : float)->float:
            if(x<0.01):
                return 1.0/3-x*x*(1-x*x/28)/10
            else:
                return (math.sin(x)-x*math.cos(x))/(x*x*x)
        def ScatterFactor(q2:float,v2T:float)->float:
            bf = 3*myBessel(math.sqrt(q2)*R)*math.exp(-q2*s2/2)
            return bf*bf*cns_fac
        self.py_func = ScatterFactor
    def as_func(self):
        K = 2/(3/self.R**2+self.s2/2)
        return lambda y: self.py_func(y*K,0)
    def factor(self):
        return _evdm.helm_factor(self.R,self.s2,self.cfac)
    def str_char(self):
        return f'W_plus_{self.nucleus.name}_{self.nucleus.A}_{self.wimp.In}{self.wimp.Out}'
    

class FormFactor_Standard:
    """
    class with info of scattering model with wimp+nucleus,
    contains Wimp in - out parametrs
    
    """
    def __init__(self,
            wimp_pars: dm_m.WimpScatterParams,
            nucleus: Nucleus,
            operator,
            operator_norm,
            norm_dv,
        ):
        """
        wimp_pars: instance of class WimpScatterParams\n
        nucleus: instance of class Nucleus, contain nucleus information\n
        operator: a linear composition of O_i with coeffs.\n
        operator_norm: same as operator, but used to normilize\n 
        cross section to Hydrogen (if None, then same as operator)\n
        norm_dv: delta velocity in scatter process with hydrogen to normalize.\n
        """
        

        if(operator_norm == None):
            operator_norm = operator
        self.wimp = wimp_pars
        self.nucleus = nucleus
        self.operator = operator
        self.norm_op = operator_norm

        _H : Nucleus = Nucleus.Hydrogen
        
        sympyficate = lambda x: x if(isinstance(x,(sympy.Expr))) else x.symbol

        m_mat_el = _symv.GetMatrixElement(sympyficate(operator),nucleus.factors)
        m_mat_el_h = _symv.GetMatrixElement(sympyficate(operator_norm),_H.factors)
        self.matel = m_mat_el
        self.norm_matel = m_mat_el_h

        

        normd_arrays = _symv.FormFactorArrays(m_mat_el,
                m_mat_el_h,_H.b,nucleus.b,
                wimp_pars.mass,_H.mass,nucleus.mass,
                0,wimp_pars.delta,wimp_pars.spin,
                nucleus.spin,norm_dv
            )
        
        #self.Zero indicates that all coeffs are zero
        if(len(normd_arrays)== 3 and len(normd_arrays[2]) == 1 and normd_arrays[2][0] ==0):
            print(f"Warning: form factor of {nucleus} is zero")
            self.Zero = True
        else:
            self.Zero = False
        self.coeffs = normd_arrays
    def as_func(self):
        def eval_poly(is_inv,Coeffs,y):
            factor = 1
            if(is_inv):
                factor = 1/y
            _sum = 0
            for c in Coeffs:
                _sum = _sum + factor*c
                factor *= y
            return math.exp(-2*y)*_sum
        return lambda y: eval_poly(self.coeffs[1],self.coeffs[2],y)
    
    def __repr__(self):
        return f'Scatter({self.wimp} + {self.nucleus}, op = {self.operator})'
    def __str__(self):
        return self.__repr__()
    def str_char(self):
        return f'W_plus_{self.nucleus.name}_{self.nucleus.A}_{self.wimp.In}{self.wimp.Out}'
    def factor(self):
        return _evdm.qexp_factor(*self.coeffs)

class ScatterModel:
    """
    class with info of scattering model with wimp+nucleus,
    contains Wimp in - out parametrs
    
    """
    def __init__(self,
            wimp_pars: dm_m.WimpScatterParams,
            nucleus: Nucleus,
            operator,
            operator_norm,
            norm_dv,
            form_factor = None
        ):
        """
        wimp_pars: instance of class WimpScatterParams\n
        nucleus: instance of class Nucleus, contain nucleus information\n
        operator: a linear composition of O_i with coeffs.\n
        operator_norm: same as operator, but used to normilize\n 
        cross section to Hydrogen (if None, then same as operator)\n
        norm_dv: delta velocity in scatter process with hydrogen to normalize.\n
        """
        self.wimp = wimp_pars
        self.nucleus = nucleus

        if(form_factor == "helm"):
            self.ff = FormFactor_Helm(wimp_pars,nucleus)
        elif(form_factor == "exp"):
            self.ff = FormFactor_Helm(wimp_pars,nucleus,0)
        else:
            try:
                self.ff = FormFactor_Standard ( wimp_pars, nucleus, operator, operator_norm, norm_dv)
            except Exception as e:
                print("can't create form factors", e)
                print("fallback to helm form factors")
                self.ff = FormFactor_Helm(wimp_pars,nucleus)

    def as_func(self):
        return self.ff.as_func()
    
    def __repr__(self):
        return f'Scatter({self.wimp} + {self.nucleus}, op = {self.operator})'
    def __str__(self):
        return self.__repr__()
    def str_char(self):
        return f'W_plus_{self.nucleus.name}_{self.nucleus.A}_{self.wimp.In}{self.wimp.Out}'
    def factor(self):
        return self.ff.factor()


    
class ScatterModel_TestFF:
    def __init__(self,
            wimp_pars: dm_m.WimpScatterParams,
            nucleus: Nucleus
        ):
        import math
        import llvmlite.ir
        import llvmlite.binding as llvm
        self.wimp = wimp_pars
        self.nucleus = nucleus
        self.Zero = False
        #llvm.initialize()
        llvm.initialize_native_target()
        llvm.initialize_native_asmprinter()
        llvm_ir_str = ScatterModel_TestFF._src_llvm_ir
        delta_str = "fdiv float _delta_value_, %0"
        no_delta_str = "fsub float 0.0, 0.0"
        if(self.wimp.delta == 0):
            _delta_q = no_delta_str
        else:
            _delta_q = delta_str

        

        llvm_ir_str = llvm_ir_str.replace("_define_delta_q_",_delta_q )
        llvm_ir_str = llvm_ir_str.replace("_delta_value_",self._float_to_llvm(self.wimp.delta) )
        llvm_ir_str = llvm_ir_str.replace("_1_mu_",self._float_to_llvm((self.wimp.mass+nucleus.mass)/(2*self.wimp.mass*nucleus.mass)))

        self._engine = self._create_execution_engine()
        self._mod = self._compile_ir(self._engine, llvm_ir_str)
        self._func_ptr = self._engine.get_function_address("scatfunc")


    def factor(self):
        return _evdm.func_factor(type('_scat_factor_fptr',(),{'address':self._func_ptr}))
    def str_char(self):
        return f'W_plus_{self.nucleus.name}_{self.nucleus.A}_{self.wimp.In}{self.wimp.Out}'
    
    @staticmethod
    def _float_to_llvm(x:float):
        import llvmlite.ir
        return str(llvmlite.ir.Constant(llvmlite.ir.FloatType(),x)).replace('float',"")
    
    _src_llvm_ir = """
        define float @scatfunc(float %0, float %1){
            %3 = _define_delta_q_
            %4 = fadd float %3, _1_mu_
            %5 = fmul float %4, %4
            %6 = tail call float @llvm.fmuladd.f32(float %5, float %0, float %1)
            %7 = tail call float @sqrtf(float %6)
            %8 = fdiv float 1.000000e+00, %7
            ret float %8
        }
        declare dso_local float @sqrtf(float) local_unnamed_addr #1
        declare float @llvm.fmuladd.f32(float, float, float)
        """
    @staticmethod
    def _create_execution_engine():
        """
        Create an ExecutionEngine suitable for JIT code generation on
        the host CPU.  The engine is reusable for an arbitrary number of
        modules.
        """
        import llvmlite.binding as llvm
        # Create a target machine representing the host
        target = llvm.Target.from_default_triple()
        target_machine = target.create_target_machine()
        # And an execution engine with an empty backing module
        backing_mod = llvm.parse_assembly("")
        engine = llvm.create_mcjit_compiler(backing_mod, target_machine)
        return engine
    @staticmethod
    def _compile_ir(engine, llvm_ir):
        """
        Compile the LLVM IR string with the given engine.
        The compiled module object is returned.
        """
        import llvmlite.binding as llvm
        # Create a LLVM module object from the IR
        mod = llvm.parse_assembly(llvm_ir)
        mod.verify()
        # Now add the module and make sure it is ready for execution
        engine.add_module(mod)
        engine.finalize_object()
        engine.run_static_constructors()
        return mod
    
class ScatterModel_SimpleFF:
    def __init__(self,
            wimp_pars: dm_m.WimpScatterParams,
            nucleus: Nucleus
        ):
        import math
        import llvmlite.ir
        import llvmlite.binding as llvm
        self.Zero = False 

        self.wimp = wimp_pars
        self.nucleus = nucleus
        fermi_GeV= 5

        s2 : float =  (fermi_GeV*0.9)**2
        b = (1.23*nucleus.A**(1.0/3)-0.6)*fermi_GeV
        a = 0.52*fermi_GeV
        R : float = math.sqrt(b*b+7*math.pi**2*a*a/3-5*s2)
        mp = Nucleus.Hydrogen.mass
        cns_fac : float = nucleus.A**4*( (wimp_pars.mass+mp)/(wimp_pars.mass+nucleus.A*mp) )**2

        self.s2 = s2
        self.R = R  
        self.cfac = cns_fac
        def myBessel(x : float)->float:
            if(x<0.01):
                return 1.0/3-x*x*(1-x*x/28)/10
            else:
                return (math.sin(x)-x*math.cos(x))/(x*x*x)
        def ScatterFactor(q2:float,v2T:float)->float:
            bf = 3*myBessel(math.sqrt(q2)*R)*math.exp(-q2*s2/2)
            return bf*bf*cns_fac
        self.py_func = ScatterFactor

        #llvm.initialize()
        llvm.initialize_native_target()
        llvm.initialize_native_asmprinter()
        llvm_ir_str = ScatterModel_SimpleFF._src_llvm_ir
        llvm_ir_str = llvm_ir_str.replace("_R_value_",self._float_to_llvm(R) )
        llvm_ir_str = llvm_ir_str.replace("_s_squad_vaule_",self._float_to_llvm(-s2))
        llvm_ir_str = llvm_ir_str.replace("_const_factor_vaule_", self._float_to_llvm(cns_fac))

        self._engine = self._create_execution_engine()
        self._mod = self._compile_ir(self._engine, llvm_ir_str)
        self._func_ptr = self._engine.get_function_address("scatfunc")


    def factor(self):
        return _evdm.func_factor(type('_scat_factor_fptr',(),{'address':self._func_ptr}))
    def str_char(self):
        return f'W_plus_{self.nucleus.name}_{self.nucleus.A}_{self.wimp.In}{self.wimp.Out}'
    
    @staticmethod
    def _float_to_llvm(x:float):
        import llvmlite.ir
        return str(llvmlite.ir.Constant(llvmlite.ir.FloatType(),x)).replace('float',"")
    
    _src_llvm_ir = """
        define float @scatfunc(float %0, float %1){
            %3 = tail call float @sqrtf(float %0) #2
            %4 = fmul float %3,  _R_value_ 
            %5 = fcmp olt float %4, 1.0
            br i1 %5, label %6, label %13

        6:                                                ; preds = %2
            %7 = fmul float %4, %4
            %8 = fdiv float %7, 2.800000e+01
            %9 = fsub float 1.000000e+00, %8
            %10 = fmul float %7, %9
            %11 = fdiv float %10, 1.000000e+01
            %12 = fsub float 0x3FD5555560000000, %11
            br label %21

        13:                                               ; preds = %2
            %14 = tail call float @sinf(float %4) #2
            %15 = tail call float @cosf(float %4) #2
            %16 = fmul float %4, %15
            %17 = fsub float %14, %16
            %18 = fmul float %4, %4
            %19 = fmul float %4, %18
            %20 = fdiv float %17, %19
            br label %21

        21:                                               ; preds = %6, %13
            %22 = phi float [ %12, %6 ], [ %20, %13 ]
            %23 = fmul float %22, 3.000000e+00
            %24 = fmul float %0, _s_squad_vaule_
            %25 = fmul float %24, 5.000000e-01
            %26 = tail call float @expf(float %25) #2
            %27 = fmul float %23, %26
            %28 = fmul float %27, _const_factor_vaule_
            %29 = fmul float %27, %28
            ret float %29
        }
        declare dso_local float @sqrtf(float) local_unnamed_addr #1

        ; Function Attrs: mustprogress nofree nounwind willreturn
        declare dso_local float @expf(float) local_unnamed_addr #1

        ; Function Attrs: mustprogress nofree nounwind willreturn
        declare dso_local float @sinf(float) local_unnamed_addr #1

        ; Function Attrs: mustprogress nofree nounwind willreturn
        declare dso_local float @cosf(float) local_unnamed_addr #1
        """
    @staticmethod
    def _create_execution_engine():
        """
        Create an ExecutionEngine suitable for JIT code generation on
        the host CPU.  The engine is reusable for an arbitrary number of
        modules.
        """
        import llvmlite.binding as llvm
        # Create a target machine representing the host
        target = llvm.Target.from_default_triple()
        target_machine = target.create_target_machine()
        # And an execution engine with an empty backing module
        backing_mod = llvm.parse_assembly("")
        engine = llvm.create_mcjit_compiler(backing_mod, target_machine)
        return engine
    @staticmethod
    def _compile_ir(engine, llvm_ir):
        """
        Compile the LLVM IR string with the given engine.
        The compiled module object is returned.
        """
        import llvmlite.binding as llvm
        # Create a LLVM module object from the IR
        mod = llvm.parse_assembly(llvm_ir)
        mod.verify()
        # Now add the module and make sure it is ready for execution
        engine.add_module(mod)
        engine.finalize_object()
        engine.run_static_constructors()
        return mod
            
        
        
        


