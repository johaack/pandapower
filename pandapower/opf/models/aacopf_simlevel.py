"""Simlevel algebraic AC OPF impelementation based on GlavitschBacher1991.

Glavitsch, H. & Bacher, R. LEONDES, C. (Ed.)
Optimal Power Flow Algorithms
Analysis and Control System Techniques for Electric Power Systems, Part 1 of 4, Academic Press,
1991, 41, Part 1 , 135 - 205
DOI: http://dx.doi.org/10.1016/B978-0-12-012741-2.50008-7
"""

import numpy as np

from pyomo.environ import Constraint
from pyomo.environ import cos
from pyomo.environ import sin


def create_Constraint(model):
    """Constraint Functions OPF algebraic model.

    Functions representing the constraints in an algebraic optimal power flow model
    """
    def make_constraint(arg1, arg2):
        """Add a constraint to the pyomo model exploting a Constraint rule."""
        def wrapper(func):
            m = model
            name = func.__name__.split('_rule')[0]
            return(setattr(m, name, Constraint(arg1, arg2, rule=func, doc=func.__doc__)))
        return(wrapper)

    @make_constraint(model.Bus, model.Time)
    def p_balance_rule(model, b, t):
        """Real power balance at each node (is valid for PQ and PV buses).

        positive active power is power consumption, negative active power is power generation

        The power flow values for branch elements (lines & transformer) are always defined as the
        power flow into the branch element.
        """
        gen = 0.
        if b in model.Genbus:
            for g in np.where(model.Gens_buses == b)[0]:
                gen += model.p_gen[g, t]

        dem = model.P_dem[b, t]  # positive for consumption
        flow = 0.
        for n in model.Bus:
            if (b, n) in model.Branch:  # flow from bus b
                flow -= model.p[b, n, t]
            if (n, b) in model.Branch:  # flow to bus b
                flow += model.p[n, b, t]
        shunt = -model.v[b, t]**2 * model.G_sh[b]
        return(gen + dem + flow + shunt == 0)

    @make_constraint(model.PQ, model.Time)
    def q_balance_PQ_rule(model, b, t):
        """Reactive power balance at each PQ node.

        positive reactive power is inductive consumption
        negative reactive power is capacitive consumption
        """
        gen = 0.
        if b in model.Genbus:
            for g in np.where(model.Gens_buses == b)[0]:
                gen += model.q_gen[g, t]

        dem = model.Q_dem[b, t]
        flow = 0.
        for n in model.Bus:
            if (b, n) in model.Branch:
                flow += model.q[b, n, t]
            if (n, b) in model.Branch:
                flow -= model.q[n, b, t]
        shunt = -model.v[b, t]**2 * model.B_sh[b]
        return(gen + dem + flow + shunt == 0)

    @make_constraint(model.PV, model.Time)
    def q_balance_PV_rule(model, b, t):
        """Reactive power balance at each PV node."""
        gen = 0.
        if b in model.Genbus:
            for g in np.where(model.Gens_buses == b)[0]:
                gen += model.q_gen[g, t]

        dem = model.Q_dem[b, t]
        flow = 0.
        for n in model.Bus:
            if (b, n) in model.Branch:
                flow += model.q[b, n, t]
            if (n, b) in model.Branch:
                flow -= model.q[n, b, t]
        shunt = -model.v[b, t]**2 * model.B_sh[b]
        return(gen + dem + flow + shunt == 0)

    @make_constraint(model.Slack, model.Time)
    def q_balance_Slack_rule(model, b, t):
        """Reactive power balance at each Slack node."""
        gen = 0.
        if b in model.Genbus:
            for g in np.where(model.Gens_buses == b)[0]:
                gen += model.q_gen[g, t]
        dem = model.Q_dem[b, t]
        flow = 0.
        for n in model.Bus:
            if (b, n) in model.Branch:
                flow += model.q[b, n, t]
            if (n, b) in model.Branch:
                flow -= model.q[n, b, t]
        shunt = -model.v[b, t]**2 * model.B_sh[b]
        return(gen + dem + flow + shunt == 0)

    @make_constraint(model.Branch, model.Time)
    def p_flow_from_rule(model, m, n, t):
        """Active power flow from m to n."""
        return(
            model.p[m, n, t] ==
            1 / model.TAP[m, n]**2 * model.v[m, t]**2 * model.Y[m, n] * cos(model.SHIFT[m, n]) -
            1 / model.TAP[m, n] * model.v[m, t] * model.v[n, t] * model.Y[m, n] * cos(
                model.theta[m, t] - model.theta[n, t] - model.SHIFT[m, n]))

    @make_constraint(model.Branch, model.Time)
    def p_flow_to_rule(model, m, n, t):
        """Active power flow to m from n."""
        return(
            model.p[n, m, t] ==
            model.v[n, t]**2 * model.Y[m, n] * cos(model.SHIFT[m, n]) -
            1 / model.TAP[m, n] * model.v[m, t] * model.v[n, t] * model.Y[m, n] * cos(
                model.theta[n, t] - model.theta[m, t] - model.SHIFT[m, n]))

    @make_constraint(model.Branch, model.Time)
    def q_flow_from_rule(model, m, n, t):
        """Reactive power flow from m to n."""
        return(
            model.q[m, n, t] == 1 / model.TAP[m, n]**2 * model.v[m, t]**2 *
            (-model.B_ch[m, n] / 2 - model.Y[m, n] * sin(model.SHIFT[m, n])) -
            1 / model.TAP[m, n] * model.v[m, t] * model.v[n, t] * model.Y[m, n] * sin(
                model.theta[m, t] - model.theta[n, t] - model.SHIFT[m, n]))

    @make_constraint(model.Branch, model.Time)
    def q_flow_to_rule(model, m, n, t):
        """Reactive power flow to m from n."""
        return(
            model.q[n, m, t] == model.v[n, t]**2 *
            (-model.B_ch[m, n] / 2 - model.Y[m, n] * sin(model.SHIFT[m, n])) -
            1 / model.TAP[m, n] * model.v[m, t] * model.v[n, t] * model.Y[m, n] * sin(
                model.theta[n, t] - model.theta[m, t] - model.SHIFT[m, n]))

# NOTE(v): implicitly enforced (at least for single ext_grid -> check again when extending)
    @make_constraint(model.Slack, model.Time)
    def limit_theta_ext_grid_rule(model, b, t):
        """GlavitschBacher1991 (10) Limit voltage angle at slack to 0."""
        return(model.theta[b, t] == 0)

    @make_constraint(model.Gen, model.Time)
    def limit_p_gen_rule(model, g, t):
        """Limit real power generation (which is negative for all except slack) (25)."""
        return(model.P_gen_min[g] <= model.p_gen[g, t] <= model.P_gen_max[g])

# NOTE(v): replaced by bounds in create_Variable -> relevant for the performance?
    # @make_constraint(model.Bus, model.Time)
    # def limit_v_rule(model, b, t):
    #     """Limit voltage magnitude (26)."""
    #     return(model.V_min[b] <= model.v[b, t] <= model.V_max[b])

# TODO(TAP): TAP, SHIFT and shunt are fixed in this model -> needs to be extended
    # @make_constraint(model.Branch, model.Time)
    # def limit_TAP_rule(model, b, t):
    #     """Limit tap positions of a transformer (27)."""
    #     return(model.TAP_min[b] <= model.TAP[b, t] <= model.TAP_max[b])

# TODO(SHIFT): TAP, SHIFT and shunt are fixed in this model -> needs to be extended
    # @make_constraint(model.Branch, model.Time)
    # def limit_SHIFT_rule(model, b, t):
    #     """Limit phase shift angles of a transformer (28)."""
    #     return(model.SHIFT_min[b] <= model.SHIFT[b, t] <= model.SHIFT_max[b])

# TODO(shunt): TAP, SHIFT and shunt are fixed in this model -> needs to be extended
    # @make_constraint(model.Bus, model.Time)
    # def limit_shunt_rule(model, b, t):
    #     """shunt capacitances or reactances (29)."""
    #     return(model.shunt_min[b] <= model.shunt[b, t] <= model.shunt_max[b])

# FIXME(): this ruins the optimisation and I dont unterstand why.
    # @make_constraint(model.Gen, model.Time)
    # def limit_q_gen_rule(model, g, t):
    #     """Limit reactive power generation (30).
    #
    #     Attention: the reactive limits on a generator are complex and usually state dependent.
    #         this rule is a simplication of the limits.
    #     """
    #     return(model.Q_gen_min[g] <= model.q_gen[g, t] <= model.Q_gen_max[g])

    @make_constraint(model.Branch, model.Time)
    def limit_p_flow_rule(model, m, n, t):
        """Limit real power flow to branch limit (31)."""
        return(-model.RATE_A[m, n] <= model.p[m, n, t] <= model.RATE_A[m, n])

    @make_constraint(model.Branch, model.Time)
    def limit_s_flow_rule(model, m, n, t):
        """Limit apparent power flow to branch limit (32)."""
        return(
            -model.RATE_A[m, n]**2 <=
            model.q[m, n, t]**2 + model.p[m, n, t]**2 <=
            model.RATE_A[m, n]**2)

# NOTE(currents): currents are not included in this model, as they are implicitly contained in power
    # @make_constraint(model.Branch, model.Time)
    # def limit_i_flow_rule(model, m, n, t):
    #     """Upper limits on current magnitudes in transmission lines or transformers (33)."""
    #     return(model.i[m, n, t] <= model.I_max[m, n])

    @make_constraint(model.Branch, model.Time)
    def limit_angle_dif_rule(model, m, n, t):
        """Limit voltage angles between nodes (34)."""
        return(-np.pi / 3 <= (model.theta[m, t] - model.theta[n, t]) <= np.pi / 3)

# NOTE(GlavitschBacher1991(35)-(36)): Constraints limiting flow between areas ignored
# NOTE(): the following constraints are not part of GlavitschBacher1991

    @make_constraint(model.PV | model.Slack, model.Time)
    def PV_node_v_set_rule(model, b, t):
        """Set voltage magnitude to V_set at all PV buses."""
        return(model.v[b, t] == model.V_set[b])

    return(model)


if __name__ == '__main__':
    pass
