def _get_spin(Z,A):
    if(A%4 == 0):
        return 0
    elif(A%2 == 0):
        if(A == 14):
            return 1
        else:
            return 0
    else:
        if(A == 1):
            return 1/2
        if(A == 3):
            return 1/2
        if(A == 23):
            return 3/2
        if(A == 27):
            return 5/2
    raise Exception(f'unsupported element with A = {A}')
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
        28 : 'Ni' 
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


        