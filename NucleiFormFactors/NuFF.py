import re
import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.sympy_parser import standard_transformations,implicit_multiplication_application
#import pkg_resources
import os
class ff_expr:
    def __init__(self, coeffs, is_exp=False):
        """
        Initialize the ff_expr object.

        :param coefficients: An np.array of float coefficients for the polynomial.
        :param is_exp: A boolean indicating whether the expression is exponential.
        """
        self.coeffs = np.array(coeffs, dtype=float)
        self.is_exp = is_exp

    def __call__(self, x):
        """
        Evaluate the expression at a given value of x when the instance is called.

        :param x: The value at which to evaluate the expression.
        :return: The result of evaluating the expression.
        """
        poly_result = np.polyval(self.coeffs, x)
        exp_result = np.exp(-2 * x) if(self.is_exp  ) else 1 
        return exp_result * poly_result

    def __str__(self):
        """
        Return a string representation of the expression.

        :return: A string representation of the expression.
        """
        if self.is_exp:
            return f"exp(-2y)*poly({self.coeffs},y)"  # Display as exp(-2x) * polynomial
        else:
            return f"poly({self.coeffs},y)"  # Display as a polynomial

__form_factor_dict = None
def __init__form_factors__():
    module_directory = os.path.dirname(__file__)
    ff_filename = os.path.join(module_directory, 'ntEFT_pre.txt')
    replace_map = {
        
        "\\Sigma^{\\prime\\prime}" : "Sigma_pp",
        "\\Sigma^{\\prime}" : "Sigma_pp",
        "\\Sigma^\\prime" : "Sigma_p",
        "\\Phi^{\\prime}":"Phi_p",
        "\\Phi^\\prime":"Phi_p",
        "\\Phi^{\\prime\\prime}":"Phi_pp",
        "\\Delta" : "Delta",
        "\\tilde{\\Phi}^\\prime": "Phi_t_p"
    }

    with open(ff_filename) as mfile:
        ff_content = mfile.read()

    for old_pattern, new_value in replace_map.items():
        ff_content = ff_content.replace(old_pattern, new_value)
    #print(ff_content)
    content_array = re.split(r'\\subsection\*',ff_content)

    element_masses = []
    element_full = []
    element_names = []
    element_ff_text = []

    #formfactor = (formfactor ij (00,01,10,11), formfactor_name (some latex text), formula)

    extract_pattern = r'{([^()]+)\s*\(\$\^{?(\d+)}?\$([A-Za-z]+)\)}\s*\\begin\{flalign\}([\s\S]+?)\\end\{flalign\}'
    for content in content_array:
        #print(f'doing for content: {content}')
        match = re.match(extract_pattern, content)
        full_name = "NoName"
        mass = 0
        short_name = "NN"
        ff_text = ""

        if match:
            full_name = match.group(1)
            mass = match.group(2)
            short_name = match.group(3)
            ff_text = match.group(4)

        element_masses.append(mass)
        element_full.append(full_name)
        element_names.append(short_name)
        element_ff_text.append(ff_text)
    content = ""

    element_form_factors = []
    element_form_factors_num = []

    ff_expression_pattern = r'W\^{(\d{2})}_{([\w]*)}\s*\(\s*y\s*\)\s*&=\s*([ye\+\-\*/\.\d\s\^{}\(\)]*)\s*&'
    #pattern = r'W\^{(\d{2})}_{(.*?)}\(y\)&=\s*e\^{-2y}\s*\(([ye\+\-\*/\.\d\s\^{}\(\)]*)\)\s*&\\nonumber\\\\'

    #ff_expression_pattern = r'W\^{(\d{2})}_{(\s\S)} [\s\S]* &\\nonumber\\\\'
    for ff_text in element_ff_text:
        #print(f'doing for content: {ff_text}')
        matches = re.finditer(ff_expression_pattern, ff_text, re.DOTALL)
        form_factors = []
        mtch_sz = 0
        for match in matches:
            Wnum = match.group(1)
            Wname = match.group(2)
            Wexpression = match.group(3)
            form_factors.append((Wnum,Wname,Wexpression.replace('e^{-2y}',"C").replace('\n',"").replace('^',"**")))
            mtch_sz +=1
        element_form_factors.append(form_factors)
        element_form_factors_num.append(mtch_sz)

    #print(element_masses)
    #print(element_full)
    #print(element_names) 
    #print(element_form_factors_num)
    #print(element_form_factors)
    global __form_factor_dict
    __form_factor_dict =dict( zip([(element_names[i],element_masses[i]) for i in range(len(element_names))],
                    [{"name" : element_names[i], "fill_name" : element_full[i], "mass" : int(element_masses[i]), 
                                "ffnum" : element_form_factors_num[i],"FF":element_form_factors[i]} for i in range(len(element_names))]))

__init__form_factors__()

def get_elements():
    return __form_factor_dict.keys()

def get_iter_elements():
    return __form_factor_dict.items()

def get_form_factors(element_name,element_mass = None):
    for key,vals in __form_factor_dict.items():
        if(element_name == key[0] and (element_mass == None or element_mass == key[1])):
            return vals['FF']
    return None

def ff_make_expr(form_factor_expr_string):
    transformations = (standard_transformations + (implicit_multiplication_application,))
    expr = parse_expr(form_factor_expr_string, transformations=transformations)
    return sp.simplify(expr.subs(sp.symbols('C'),parse_expr('exp(-2*y)')))

def ff_get_poly(form_factor_expr_string):
    transformations = (standard_transformations + (implicit_multiplication_application,))
    expr = parse_expr(form_factor_expr_string, transformations=transformations)
    return ff_expr(
            sp.Poly(sp.simplify(expr.subs(sp.symbols('C'),1)), sp.symbols('y')).all_coeffs(),
            sp.symbols('C') in expr.free_symbols
        )