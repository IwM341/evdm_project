class WimpScatterParams:
    def __init__(self,ptype_in,ptype_out,mass_in,mass_delta =0 ,wimp_spin = 0.5):
        self.mass = mass_in
        self.delta = mass_delta
        self.spin = wimp_spin
        self.In  = ptype_in
        self.Out  = ptype_out
    def __repr__(self):
        return f'WimpScatter(m = {self.mass}, delta = {self.delta}, j = {self.spin})'
    def __str__(self):
        return self.__repr__()
    
class WimpModel:
    def __init__(self,wimp_mass,wimp_spin,wimp_deltas = [0]):
        wimp_deltas_new = wimp_deltas.copy()
        wimp_deltas_new.sort()
        min_delta = wimp_deltas_new[0]
        for i in range(len(wimp_deltas_new)):
            wimp_deltas_new[i] -=min_delta
        self.mass = wimp_mass
        self.deltas = wimp_deltas_new
        self.spin = wimp_spin
    def delta(self,ptype_in,ptype_out):
        return self.deltas[ptype_out]-self.deltas[ptype_in]
    def __call__(self,ptype_in = 0,ptype_out = 0) -> WimpScatterParams:
        return WimpScatterParams(
            ptype_in,ptype_out,
            self.mass,self.delta(ptype_in,ptype_out),
            self.spin)

    def __repr__(self):
        return f'Wimps(mass = {self.mass}, deltas = {self.deltas},spin = {self.spin})'
    def __str__(self):
        return self.__repr__()