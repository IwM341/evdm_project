def _get_spin(Z,A):

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
def _mindeleev_get_name(Z):
    try:
        return Mindeleev_table[Z] 
    except:
        raise ValueError(f"not supported element with Z = '{Z}'") 
    
def _mindeleev_get_Z(name):
    for  (Z,_name) in Mindeleev_table.items():
        if(_name == name):
            return Z
    raise ValueError(f"not supported element with name '{name}'") 


        