
class Droplet:
    def __init__(self, T, dia, fuel_names, Y = None, vel = None, Reyn = None):
        self.T = T
        self.dia = dia
        self.fuel_names = fuel_names
        if (Y is None):
            self.Y = [1., 0.]
        else:
            self.Y = Y
        self.fixed = True
        if (abs(sum(self.Y) - 1.) > 1.E-8):
            errorstatement = "Liquid mass fractions must sum to 1"
            raise ValueError(errorstatement)
        if (vel is None):
            self.vel = 0.
        else:
            self.vel = vel
        if (Reyn is None):
            self.Reyn = 0.
        else:
            self.Reyn = Reyn

class GasPhase:
    def __init__(self, T, P, vel = None):
        self.T = T
        self.P = P
        if (vel is None):
            self.vel = 0.
        else:
            self.vel = vel
class CaseInfo:
    def __init__(self, name, droplet: Droplet, gas: GasPhase, xyunits, domain = None, cell_num = None, reftype = None):
        self.name = name
        self.droplet = droplet
        self.gas = gas
        # If reference is experimental or computational results
        if (reftype is None):
            self.reftype = "exp"
        else:
            self.reftype = "comp"
        if (xyunits[0] == "s/mm2"):
            self.xconv = 1. / (self.droplet.dia * 1.E3)**2
            self.xlabel = "$t/d_0^2$ [s/mm$^2$]"
        elif (xyunits[0] == "ms"):
            self.xconv = 1.E3
            self.xlabel = "$t$ [ms]"
        elif (xyunits[0] == "runge"):
            # This is nu_gas / r_0**2 * 1E-2
            self.xconv = 1.3465E-7 * (self.droplet.dia / 2.)**-2
            self.xlabel = "$t \nu / r_0^2 10^{-2}$"
        else:
            self.xconv = 1.
            self.xlabel = "$t$ [s]"
        if ("2" in xyunits[1]):
            self.yexp = 2.
        else:
            self.yexp = 1.
        if ("dd0" in xyunits[1]):
            self.yconv = 1. / self.droplet.dia
            if (self.yexp == 2):
                self.ylabel = "$(d/d_0)^2$"
            else:
                self.ylabel = "$d/d_0$"
        elif(xyunits[1] == "r2_mm"):
            self.yconv = 1.E3 * 0.5
            self.ylabel = "$r^2$ [mm$^2$]"

        if (domain is None):
            self.domain = [2., 2.]
        else:
            self.domain = domain
        if (cell_num is None):
            self.cell_num = [32, 32]
        else:
            self.cell_num = cell_num
        if (abs(self.cell_num[0] / self.domain[0] - self.cell_num[1] / self.domain[1]) > 0.):
            error = "Uniform grid spacing required"
            raise ValueError(error)
    def __str__():
        return name
    def set_end_time(self, time):
        self.time = time * 1.02
        self.plot_per = self.time / 100.

def Tonini_4_33():
    drop = Droplet(300., 2.E-3, ["NC8H18", "NC10H22"], [0.5, 0.5], Reyn = 110.)
    gas = GasPhase(1001., 101325.)
    case = CaseInfo("Tonini_4_33", drop, gas, xyunits = ["s/mm2", "dd02"], reftype = "comp")
    return case
def Daif():
    drop = Droplet(294., 1.334E-3, ["NC7H16", "NC10H22"], [0.667, 0.333])
    gas = GasPhase(348., 101325., vel = 3.10)
    case = CaseInfo("Daif", drop, gas, xyunits = ["s", "r2_mm"])
    return case
def Abramzon():
    drop = Droplet(300., 1.E-4, ["NC10H22", "NC8H18"], vel = 15.)
    gas = GasPhase(1500., 1013250.)
    case = CaseInfo("Abramzon", drop, gas, domain = [0.08, 0.005], cell_num = [256, 16], xyunits = ["ms", "dd0"], reftype="comp")
    return case
def RungeMix():
    drop = Droplet(272., 550.E-6, ["NC7H16", "NC10H22"], [0.5, 0.5])
    gas = GasPhase(273., 101325., vel = 2.5)
    case = CaseInfo("RungeMix", drop, gas, domain = [40., 40.], xyunits = ["runge", "dd02"])
    return case
def RungeDec():
    drop = Droplet(272., 500.E-6, ["NC7H16", "NC10H22"], [0., 1.])
    gas = GasPhase(273., 101325., vel = 2.5)
    case = CaseInfo("RungeDec", drop, gas, domain = [40., 40.], xyunits = ["runge", "dd02"])
    return case
def RungeHep():
    drop = Droplet(272., 570.E-6, ["NC7H16", "NC10H22"], [1., 0.])
    gas = GasPhase(273., 101325., vel = 2.5)
    case = CaseInfo("RungeHep", drop, gas, domain = [10., 10.], xyunits = ["runge", "dd02"])
    return case

def CreateInputParams(case):
    dom_lo = [-case.domain[0] / 2, -case.domain[1] / 2]
    dom_hi = [case.domain[0] / 2, case.domain[1] / 2]
    prob_params = "geometry.prob_lo = {:g} {:g} ".format(dom_lo[0], dom_lo[1])
    prob_params += "geometry.prob_hi = {:g} {:g} ".format(dom_hi[0], dom_hi[1])
    prob_params += "amr.plot_file = {}/plt ".format(case.name)
    prob_params += "amr.n_cell = {:d} {:d} ".format(case.cell_num[0], case.cell_num[1])
    prob_params += "amr.plot_per = {:.1g} ".format(case.plot_per)
    prob_params += "amr.stop_time = {:.1g} ".format(case.time)
    prob_params += "prob.P_mean = {:g} ".format(case.gas.P)
    prob_params += "prob.T_gas = {:g} ".format(case.gas.T)
    fixed_part = True
    if (case.droplet.Reyn > 0.):
        prob_params += "prob.re_d = {:g} ".format(case.droplet.Reyn)
    elif (case.gas.vel > 0.):
        prob_params += "prob.vel_gas = {:g} ".format(case.gas.vel)
    elif (case.droplet.vel > 0.):
        prob_params += "prob.vel_drop = {:g} 0. 0. ".format(case.droplet.vel)
        prob_params += "prob.loc_drop = {:g} 0. 0. ".format(-case.domain[0] * 0.4)
        fixed_part = False
    prob_params += "prob.T_drop = {:g} ".format(case.droplet.T)
    prob_params += "prob.dia_drop = {:g} ".format(case.droplet.dia)
    prob_params += "prob.Y_drop = {:g} {:g} ".format(case.droplet.Y[0], case.droplet.Y[1])
    prob_params += "particles.fuel_species = " + case.droplet.fuel_names[0] + " " + case.droplet.fuel_names[1] + " "
    if (fixed_part):
        prob_params += "particles.fixed_parts = 1 "
        prob_params += "peleLM.lo_bc = Inflow Interior "
        prob_params += "peleLM.hi_bc = Outflow Interior "
        prob_params += "geometry.is_periodic = 0 1 1 "
    else:
        prob_params += "particles.fixed_parts = 0 "
        prob_params += "peleLM.lo_bc = Interior Interior "
        prob_params += "peleLM.hi_bc = Interior Interior "
        prob_params += "geometry.is_periodic = 1 1 1 "
    return prob_params
