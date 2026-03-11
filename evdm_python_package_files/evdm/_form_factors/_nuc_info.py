from ._isotopes import isotope_table

def _get_spin(Z,A):

    try:
        return isotope_table.query(f"A == {A} and Z == {Z}")['spin'].values[0]
    except:
        raise ValueError(f"not supported element with Z = '{Z}' and A = '{A}'")
    
    specials = {
        (92,235):7/2,
        (13,27):5/2,
        (11,23):3/2,
        (79,197):3/2
    }
    if (Z,A) in specials:
        return specials[(Z,A)]

    if(A%4 == 0):
        return 0
    elif(A%2 == 0):
        if(A == 14):
            return 1
        else:
            return 0
    else:
        return 1/2
    #raise Exception(f'unsupported element with A = {A}')

def _mindeleev_get_name(Z):
    try:
        return isotope_table.query(f"Z == {Z}")['name'].values[0]
    except:
        raise ValueError(f"not supported element with Z = '{Z}'") 
    
def _mindeleev_get_Z(name):
    try:
        return isotope_table.query(f"name == '{name}'")['Z'].values[0]
    except:
        raise ValueError(f"not supported element with name '{name}'") 

def _get_abondonce(A,Z):
    try:
        return isotope_table.query(f"A == {A} and Z == {Z}")['abundance'].values[0]
    except:
        raise ValueError(f"not supported element with Z = '{Z}' and A = '{A}'")


if False:
    Mindeleev_table = {
        0:'NN',
        1:'H',
        2:'He',
        6:'C',
        7:'N',
        8:'O',
        10:'Ne',
        11:'Na',
        12:'Mg',
        13:'Al',
        14:'Si',
        16:'S',
        18:'Ar',
        20: 'Ca', 
        26 : 'Fe',
        28 : 'Ni',
        90 : 'Th',
        82 : 'Pb',
        92 : 'U',
        81 : 'Ti',
        197 : 'Au'
    }