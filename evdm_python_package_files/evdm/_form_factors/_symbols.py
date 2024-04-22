
from sympy.parsing.sympy_parser import parse_expr
import sympy

"""
MatrixElementPart = parse_expr(
    "4*pi/(2*J+1)*"
    "(R_M*W_M + R_Sigma_pp*W_Sigma_pp+R_Sigma_p*W_Sigma_p+"
    "(q2/m_N**2)*(R_Phi_pp*W_Phi_pp + R_MPhi_pp*W_MPhi_pp+"
    "R_Phitild_pp*W_Phitild_pp+R_Delta*W_Delta + R_DeltaSigma_p*W_DeltaSigma_p))")
"""

MatrixElementFact = parse_expr("4*pi/(2*J+1)")

Q2 = sympy.symbols('q2',real = True)
Delta = sympy.symbols('delta',real = True)
M_N = sympy.symbols('m_N',real = True)
M_X = sympy.symbols('m_X',real = True)
MU_NX = M_X*M_N/(M_X+M_N)
Q2MN = Q2/M_N**2
J_X =  sympy.symbols('jx',real = True)
J_N =  sympy.symbols('J',real = True)
V2_T = sympy.symbols('v2_T',real = True)

XiResponses = {"R_M":parse_expr("c1*ctc1 + jx*(jx+1)/3*( (q2/m_N**2*c5*ctc5+c8*ctc8)*VT2 + q2/m_N**2 * c11*ctc11)"),
               "R_Phi_pp":Q2MN*parse_expr("1/4*q2/m_N**2 * c3*ctc3 + jx*(jx+1)/12*(c12-q2/m_N**2*c15)*(ctc12-q2/m_N**2*ctc15)"),
               "R_MPhi_pp" : Q2MN*parse_expr("re(  c3*ctc1+ jx*(jx+1)/3*(c12-q2/m_N**2*c15)*ctc11)"),
               "R_Phi_t_p" : Q2MN*parse_expr("jx*(jx+1)/12*(c12*ctc12 + q2/m_N**2*c13*ctc13)"),
               "R_Sigma_pp" : parse_expr("q2/(4*m_N**2)*c10*ctc10 + jx*(jx+1)/12*"
                    "(c4*ctc4+q2/m_N**2 *(c4*ctc6+c6*ctc4) ) + q2**2/m_N**4*c6*ctc6 + VT2*(c12*ctc12 + q2/m_N**2*c13*ctc13) "),
               "R_Sigma_p":parse_expr("1/8*(q2/m_N**2*c3*ctc3+c7*ctc7)*VT2 + jx*(jx+1)/12*"
                    "(c4*ctc4+q2/m_N**2*c8*ctc9+VT2/2*(  (c12-q2/m_N**2*c15)*(ctc12-q2/m_N**2*ctc15)+ q2/m_N**2*c14*ctc14) )"),
               "R_Delta":Q2MN*parse_expr("jx*(jx+1)/3*(q2/m_N**2*c5*ctc5 + c8*ctc8)"),
               "R_DeltaSigma_p" : Q2MN*parse_expr("jx*(jx+1)/3*re(q2/m_N**2*c5*ctc4 - c8*ctc9)")}
pre_substritutions = [("jx",sympy.symbols("jx",real = True)),
                      ("m_N",sympy.symbols("m_N",real = True)),
                      ("VT2",sympy.symbols("VT2",real = True)),
                      ("q2",sympy.symbols("q2",real = True)),
                      ("delta",sympy.symbols("delta",real = True)),
                      ("v2_T",sympy.symbols("v2_T",real = True))]
VT2 = parse_expr("v2_T - ((m_N+m_X)/(2*m_N*m_X) + delta/q2)**2*q2").subs(pre_substritutions)

def rsym(symb):
    return sympy.Symbol(symb,real =True)

def _O_impl(n : int,t : int = 0):
    """
    symbol of non relativistic potential
    :parametr n: number of potential
    :parametr t: isospin number 0 or 1
    """
    if(type(n) != int or n <=0 or n > 15):
        raise ValueError(f"unexpected potential index {n}")
    if (t == 0):
        return sympy.symbols(f'O_{n}')
    elif (t == 1):
        return sympy.symbols(f'O_t{n}')
    else:
        raise ValueError(f"unexpected isospin number {n}")

'''
def W_names():
    return sum([ [f"W_Delta{i}{j}",f"W_DeltaSigma_p{i}{j}",
                f"W_M{i}{j}", f"W_MPhi_pp{i}{j}", f"W_Phi_pp{i}{j}",
                f"W_Phitild_pp{i}{j}", f"W_Sigma_p{i}{j}",
                f"W_Sigma_pp{i}{j}"] for i in range(1) for j in range(1)],[])
def W_name():
    return ["W_Delta","W_DeltaSigma_p",
                 "W_M", "W_MPhi_pp", "W_Phi_pp",
                "W_Phitild_pp", "W_Sigma_p",
                "W_Sigma_pp"]
def R_names():
    return sum([ [f"R_Delta{i}{j}",f"R_DeltaSigma_p{i}{j}",
                f"R_M{i}{j}", f"R_MPhi_pp{i}{j}", f"R_Phi_pp{i}{j}",
                f"R_Phitild_pp{i}{j}", f"R_Sigma_p{i}{j}",
                f"R_Sigma_pp{i}{j}"] for i in range(1) for j in range(1)],[])
'''

def GetOCoeffs(OExpression):
    return ( [OExpression.coeff(_O_impl(n,0),1).simplify() for n in range(1,16)],
            [OExpression.coeff(_O_impl(n,1),1).simplify() for n in range(1,16)])

def C_Subs(t1,t2,C_t0,C_t1):
    C1 = C_t0 if(t1 == 0) else C_t1
    C2 = C_t0 if(t2 == 0) else C_t1
    Arr1= [ (f"c{i+1}",_c) for (i,_c) in enumerate(C1)]
    Arr2 = [ (f"ctc{i+1}",_c.conjugate()) for (i,_c) in enumerate(C2)]
    return  Arr1+Arr2

def ParseWFactor(W_fact : str):
    m_body = W_fact[1:-2]
    m_t1 = int(W_fact[-2:-1])
    m_t2 = int(W_fact[-1:])
    return (m_body,m_t1,m_t2)

def GetMatrixElement(OExpression,W_dict_poly_part):
    (C_t0,C_t1) = GetOCoeffs(OExpression.expand())
    PList = []
    for (W_name,_W_val) in W_dict_poly_part.items():
        (InfoName,t1,t2) = ParseWFactor(W_name)
        W_val = parse_expr(_W_val) if(type(_W_val) == str) else _W_val

        W_pol_val = W_val if(W_val.limit("y","+oo") != 0) else W_val*parse_expr("exp(2*y)")
        R_name = "R"+InfoName
        #print("Rname = ",R_name)
        R_val = XiResponses[R_name]
        #print("R_val = ",R_val)
        PList.append(
            R_val.subs(C_Subs(t1,t2,C_t0,C_t1))
            .subs(pre_substritutions)*W_pol_val*MatrixElementFact
            )
    return sum(PList)

def scatter_model_subs(m_X,m_N,delta,j_X,j_N,b):
    """
    :parametr m_X: WIMP mass, GeV
    :parametr m_N: Nuclear Mass, GeV
    :parametr delta: WIMP delta mass, GeV
    :parametr j_X: WIMP angular momentum
    :parametr j_N: Nuclear angular momentum
    :parametr b: b oscillator parametr
    """
    return {"m_N":m_N,"m_X":m_X,"delta": delta,"jx":j_X,"J":j_N,"b":b,
            M_N : m_N,M_X : m_X,Delta : delta,J_X : j_X,J_N:j_N,sympy.Symbol('b',real = True):b}

def NormalizeMatEl(MatrixElementExpression,Vdiff,b,mx,mT,mN,delta,j_X):
    """
    normalizator object:
    makes the resulting expression normalised to Hydrogen Avarage Matrix Element
    """
    VT2_subed = VT2.subs(V2_T,Vdiff**2)
    
    mu = mx*mT/(mx+mT)
    V_0 = Vdiff
    V_1 = float(sympy.sqrt(Vdiff**2-2*delta/mu))
    p_0 = mu*V_0
    p_1 = mu*V_1
    hydrogen_subs = scatter_model_subs(mx,mT,delta,j_X,0.5,b)
    FullExpr = MatrixElementExpression.subs(
        rsym('VT2'),VT2_subed).subs(rsym("q2"),"4*y/b**2").subs(hydrogen_subs)

    #print("V",V_0,V_1)
    #print("delta",delta)
    #print("FullExpr",FullExpr)
    #print("d",2*delta/mu)
    IngegralWeight = 1/(p_0*p_1*b**2)
    #print("IngegralWeight",IngegralWeight)
    Y0 = b**2*(p_0-p_1)**2/4
    Y1 = b**2*(p_0+p_1)**2/4
    #print("Y0",Y0)
    #print("Y1",Y1)
    av_mat_el = sympy.re(
        FullExpr.integrate(("y",Y0,Y1)).evalf()
        )*IngegralWeight*(V_1/V_0)
    return (mT+mx)**2/(mN+mx)**2*(mN/mT)**2/av_mat_el


def FormFactorArrays(MatElNuc : sympy.Expr ,MatEl_norm : sympy.Expr ,
                     b_H,b_N,mx,mH,mN,delta_H,delta,j_X,j_N,
                     Vdiff_Norm):
    
    norm_coeff = NormalizeMatEl(MatEl_norm,Vdiff_Norm,b_H,mx,mH,mN,delta_H,j_X)
    pot_scatter_model = scatter_model_subs(mx,mN,delta,j_X,j_N,b_N)
    
    pre_ff= (MatElNuc.subs({Q2:"4*y/b**2"}).subs(pot_scatter_model)*
            norm_coeff).expand()
    _y = sympy.symbols("y")
    is_y_1 = (pre_ff.coeff(_y,-1) != 0)
    y_1_di = 1 if is_y_1 else 0

    m_deg = sympy.degree(pre_ff*_y,_y)-1
    V0_poly = pre_ff.coeff(rsym('VT2'),0)
    V1_poly = pre_ff.coeff(rsym('VT2'),1)
    has_vT = (V1_poly != 0)
    try:
        V0_coeffs = [float(V0_poly.coeff(_y,i - is_y_1)) for i in range(m_deg+1)]
        V1_coeffs = [float(V1_poly.coeff(_y,i - is_y_1)) for i in range(m_deg+1)]
    except Exception as exeption:
        raise Exception(f"couldn't get coeffs of V0_poly = {V0_poly} and V1_poly = {V1_poly}") from exeption
    
    b = pot_scatter_model['b']
    return (b,is_y_1,V0_coeffs,V1_coeffs) if has_vT else (b,is_y_1,V0_coeffs)