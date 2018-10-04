"""A basic colletion of input data sets parameters and variables for the aacopf."""

import numpy as np

from pyomo.environ import NonNegativeIntegers
from pyomo.environ import NonNegativeReals
from pyomo.environ import Param
from pyomo.environ import Reals
from pyomo.environ import Set
from pyomo.environ import Var


def set_and_parameter_from_dat(model):
    """Set and Param definitions.

    The Sets and Params needed for the algebraic OPF are created and initialized from the input
    data "file".dat
    """
    mod = model.model()
    mod.Bus = Set(within=NonNegativeIntegers)
    mod.PV = Set(within=NonNegativeIntegers)
    mod.PQ = Set(within=NonNegativeIntegers)
    mod.Branch = Set(within=mod.Bus * mod.Bus)
    mod.Gen = Set(within=mod.Bus)
    mod.Time = Set()

    mod.RATE_A = Param(mod.Branch, domain=Reals)
    mod.G = Param(mod.Branch, domain=Reals)
    mod.B = Param(mod.Branch, domain=Reals)
    mod.R = Param(mod.Branch, domain=Reals)
    mod.Y = Param(mod.Branch, domain=Reals)
    mod.SHIFT = Param(mod.Branch, domain=Reals)
    mod.B_ch = Param(mod.Branch, domain=Reals)
    mod.G_sh = Param(mod.Bus, domain=Reals)
    mod.B_sh = Param(mod.Bus, domain=Reals)
    mod.TAP = Param(mod.Branch, domain=Reals)
    mod.P_dem = Param(mod.Bus, mod.Time, domain=Reals)
    mod.Q_dem = Param(mod.Bus, mod.Time, domain=Reals)

    mod.V_set = Param(mod.Gen, domain=Reals)
    mod.P_gen_max = Param(mod.Gen, domain=Reals)
    mod.P_gen_min = Param(mod.Gen, domain=Reals)
    mod.Q_gen_max = Param(mod.Gen, domain=Reals)
    mod.Q_gen_min = Param(mod.Gen, domain=Reals)
    mod.V_max = Param(mod.Bus, domain=Reals)
    mod.V_min = Param(mod.Bus, domain=Reals)
    mod.cost_p = Param(
        mod.Gen, mod.Time, domain=Reals,
        doc="cost for active power (€/kW per time step)")
    mod.cost_q = Param(
        mod.Gen, mod.Time, domain=Reals,
        doc="cost for active power (€/kW per time step)")

    return(model)


def set_and_parameter_from_ppci(model, ppci):
    """Set and Param definitions.

    The Sets and Params needed for the algebraic OPF are created and initialized.
    Input data for initialization is either read from ppci via functions.
    """
    def cnvrt(index, values, pu):
        """Convert index and values to a input data dict in pu."""
        return(dict(zip(index, values / pu)))

    BUS = ppci["bus"][:, 0].astype(int)
    BUS_TYPE = ppci["bus"][:, 1].astype(int)
    PQ_BUS = [x for x in BUS if BUS_TYPE[x] in [1]]
    PV_BUS = [x for x in BUS if BUS_TYPE[x] in [2]]
    SLACK_BUS = [x for x in BUS if BUS_TYPE[x] in [3]]
    BRANCH = list(zip(
        ppci["branch"][:, 0].real.astype(int),
        ppci["branch"][:, 1].real.astype(int)))
    BRANCH_INVERSE = list(zip(
        ppci["branch"][:, 1].real.astype(int),
        ppci["branch"][:, 0].real.astype(int)))

    GEN = np.arange(ppci["gen"][:, 0].size, dtype=int)
    GENBUS = np.unique(ppci["gen"][:, 0]).astype(int)

    TIME = [0]

    # TODO(timestepping): implement for more than one timestep
    BUS_TIME = list(zip(BUS, np.zeros(len(BUS), dtype=int)))
    GEN_TIME = list(zip(GEN, np.zeros(len(GEN), dtype=int)))

    BASE_MVA = ppci["baseMVA"] * 1e+03

    RATE_A = np.array([i[5].real if i[5].real > 0. else 1e+09 for i in ppci["branch"]])

    R = ppci["branch"][:, 2].real
    X = ppci["branch"][:, 3].real
    Z_cmplx = R + 1j*X
    Y_cmplx = Z_cmplx**(-1)
    Y = np.absolute(Y_cmplx)
    SHIFT = np.angle(Y_cmplx)
    G = Y_cmplx.real
    B = Y_cmplx.imag

    model.Bus = Set(within=NonNegativeIntegers, initialize=BUS)
    model.PQ = Set(within=NonNegativeIntegers, initialize=PQ_BUS)
    model.PV = Set(within=NonNegativeIntegers, initialize=PV_BUS)

    model.Slack = Set(within=NonNegativeIntegers, initialize=SLACK_BUS)
    model.Branch = Set(within=model.Bus * model.Bus, initialize=BRANCH)
    model.Branch_inverse = Set(within=model.Bus * model.Bus, initialize=BRANCH_INVERSE)
    model.Gen = Set(within=model.Bus, initialize=GEN, doc="all gens")
    model.Genbus = Set(within=model.Bus, initialize=GENBUS, doc="all buses with gen(s)")
    model.Time = Set(initialize=TIME)

    model.RATE_A = Param(model.Branch, domain=Reals, initialize=cnvrt(BRANCH, RATE_A, BASE_MVA))
    model.G = Param(model.Branch, domain=Reals, initialize=cnvrt(BRANCH, G, BASE_MVA))
    model.B = Param(model.Branch, domain=Reals, initialize=cnvrt(BRANCH, B, BASE_MVA))
    model.R = Param(model.Branch, domain=Reals, initialize=cnvrt(BRANCH, R, BASE_MVA))
    model.Y = Param(model.Branch, domain=Reals, initialize=cnvrt(BRANCH, Y, BASE_MVA))
    model.SHIFT = Param(model.Branch, domain=Reals, initialize=cnvrt(BRANCH, SHIFT, 1))

    model.B_ch = Param(
        model.Branch, domain=Reals, initialize=cnvrt(BRANCH, ppci["branch"][:, 4].real, BASE_MVA))
    model.G_sh = Param(model.Bus, domain=Reals, initialize=cnvrt(BUS, ppci["bus"][:, 4], BASE_MVA))
    model.B_sh = Param(model.Bus, domain=Reals, initialize=cnvrt(BUS, ppci["bus"][:, 5], BASE_MVA))
    model.TAP = Param(
        model.Branch, domain=Reals, initialize=cnvrt(BRANCH, ppci["branch"][:, 8].real, 1))

    model.Gens_buses = Param(
        model.Gen, domain=NonNegativeIntegers,
        initialize=dict(zip(GEN, ppci["gen"][:, 0].astype(int))),
        doc="buses gens are connected to")

    model.P_dem = Param(
        model.Bus, model.Time, domain=Reals,
        initialize=cnvrt(BUS_TIME, ppci["bus"][:, 2], BASE_MVA))
    model.Q_dem = Param(
        model.Bus, model.Time, domain=Reals,
        initialize=cnvrt(BUS_TIME, ppci["bus"][:, 3], BASE_MVA))

    model.V_set = Param(
        model.Genbus, domain=Reals,
        initialize=cnvrt(GENBUS, ppci["gen"][:, 5], 1))

    model.P_gen_max = Param(
        model.Gen, domain=Reals,
        initialize=cnvrt(GEN, ppci["gen"][:, 9], -BASE_MVA))
    model.P_gen_min = Param(
        model.Gen, domain=Reals,
        initialize=cnvrt(GEN, ppci["gen"][:, 8], -BASE_MVA))
    model.Q_gen_max = Param(
        model.Gen, domain=Reals,
        initialize=cnvrt(GEN, ppci["gen"][:, 4], -BASE_MVA))
    model.Q_gen_min = Param(
        model.Gen, domain=Reals,
        initialize=cnvrt(GEN, ppci["gen"][:, 3], -BASE_MVA))
    model.V_max = Param(model.Bus, domain=Reals, initialize=cnvrt(BUS, ppci["bus"][:, 11], 1))
    model.V_min = Param(model.Bus, domain=Reals, initialize=cnvrt(BUS, ppci["bus"][:, 12], 1))

    def read_cost_from_gencost(g):
        """Read cost data from ppci. If no q cost is available it is set to 0.

        Parameters
        ----------
        g : Number of gens in the grid

        Notes
        -----
        n.a.

        """
        # TODO(COST): Only ploynomial (degree 1) cost implemented -> allow other cost models
        # TODO(COST): extent cost model by prioviding polynomial cost with degree > 1
        # TODO(COST): check how the two gens at one bus problem can be solved
        l = len(ppci["gencost"])
        if l > g:
            COST_p = np.array([i[5] for i in ppci["gencost"][:l//2]])
            COST_q = np.array([i[5] for i in ppci["gencost"][l//2:]])
        else:
            COST_p = np.array([i[5] for i in ppci["gencost"]])
            COST_q = np.array([0. for i in ppci["gencost"]])

        return(COST_p, COST_q)

    COST_p, COST_q = read_cost_from_gencost(len(GEN))

    model.cost_p = Param(
        model.Gen, model.Time, domain=Reals,
        initialize=cnvrt(GEN_TIME, COST_p, 1),
        doc="cost for active power (€/kW per time step)")
    model.cost_q = Param(
        model.Gen, model.Time, domain=Reals,
        initialize=cnvrt(GEN_TIME, COST_q, 1),
        doc="cost for active power (€/kW per time step)")

    return(model)


def create_variable(model):
    """Var definitions.

    Define all Variables of the OPF model.
    """
    model.v = Var(model.Bus, model.Time, domain=NonNegativeReals, initialize=1.0, bounds=(0.9, 1.1))
    model.theta = Var(
        model.Bus, model.Time, domain=Reals, initialize=0., bounds=(-2 * np.pi, 2 * np.pi))
    # NOTE(create_variable): use inverse Branch Set to reduce the number of items in p to 2n
    model.p = Var(model.Branch | model.Branch_inverse, model.Time, domain=Reals)
    model.q = Var(model.Branch | model.Branch_inverse, model.Time, domain=Reals)
    model.p_gen = Var(model.Gen, model.Time, domain=Reals)
    model.q_gen = Var(model.Gen, model.Time, domain=Reals)
    model.p_gen_max = Var(model.Gen, model.Time, domain=Reals)

    return(model)
