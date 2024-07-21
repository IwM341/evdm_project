def _make_latex_array():
    V = "\\vec{v}^{\\perp}"
    br = lambda _str: f"\\left({_str}\\right)"
    S_N = "\\vec{S}_N"
    S_W = "\\vec{S}_{\\chi}"
    Q = "\\frac{\\vec{q}}{m_N}"
    mlt = " \\cdot "
    tms = " \\times "
    return [
        "1",#1
        f'({V})^2',#2
        f'i {S_N} {mlt} {br(f"{Q} {tms} {V}")}',#3
        f'{S_W} {mlt} {S_N}',#4 
        f'i {S_W} {mlt} {br(f"{Q} {tms} {V}")}',#  5
        f'{br(f"{S_W}{mlt}{Q}")}{br(f"{S_N}{mlt}{Q}")}',#6
        f'{S_N}{mlt}{V}',#7
        f'{S_W}{mlt}{V}',#8
        f'{S_W}{mlt}{br(f"{S_N}{tms}{Q}")}',#9
        f'i{S_N} {mlt} {Q}',#10
        f'i{S_W} {mlt} {Q}',#11
        f'{S_W} {mlt} {br(f"{S_N} {tms} {V}")}',#12
        f'i{br(f"{S_W}{mlt}{V}")}{br(f"{S_N} {mlt} {Q}")}',#13
        f'i{br(f"{S_W}{mlt}{Q}")}{br(f"{S_N} {mlt} {V}")}',#14
        f'-{br(f"{S_W}{mlt}{Q}")}{br(f"({S_N} {tms} {V}) {mlt} {Q}")}'#15
    ]
