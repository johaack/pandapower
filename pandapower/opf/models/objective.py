"""A basic colletion of objective functions for the aacopf."""

from pyomo.environ import cos
from pyomo.environ import maximize
from pyomo.environ import minimize
from pyomo.environ import Objective


class CreateObjective(object):
    """Create Objective function and allocate to the targeted model/submoel."""

    def __init__(self, model, submodel, objective):
        """Initalise CreateObjective.

        Parameters
        ----------
        model : pyomo model object the Set, Parameter and Variable data is taken from_ppci
        submodel : pyomo model the objective is assigned to (pyomo model object)
        objective : chosen objective function (string)

        """
        super(CreateObjective, self).__init__()
        self.model = model
        self.submodel = submodel
        self.objective = objective

    def create_Objective(self):
        """Objective Function OPF algebraic mod.

        """
        # TODO(create_Objective): Simlevel restriction of p_gen_max -> influence on the losses?
        # TODO(create_Objective): include TAP as var
        model = self.model

        def make_objective(arg1):
            """Add a constraint to the pyomo mod.exploting an Objective rule."""
            def wrapper(func):
                m = self.submodel
                name = func.__name__
                if name is self.objective:
                    return(setattr(m, name, Objective(rule=func, sense=arg1)))
            return(wrapper)

        @make_objective(minimize)
        def obj_minimize_losses(self):
            """Objective function - Minimize branch losses."""
            mod = model
            return sum(
                (mod.v[m, t]**2-mod.v[n, t]**2) * mod.Y[m, n] * cos(mod.SHIFT[m, n]) -
                mod.v[m, t] * mod.v[n, t] * (
                    cos(mod.theta[m, t] - mod.theta[n, t] - mod.SHIFT[m, n]) -
                    cos(mod.theta[n, t] - mod.theta[m, t] - mod.SHIFT[m, n]))
                for m, n in mod.Branch for t in mod.Time)

        @make_objective(minimize)
        def obj_min_losses_angle_voltage(self):
            """Objective function - Minimize branch losses."""
            mod = model
            return(
                sum(
                    (mod.v[m, t]**2-mod.v[n, t]**2) * mod.Y[m, n] * cos(mod.SHIFT[m, n]) -
                    mod.v[m, t] * mod.v[n, t] * (
                        cos(mod.theta[m, t] - mod.theta[n, t] - mod.SHIFT[m, n]) -
                        cos(mod.theta[n, t] - mod.theta[m, t] - mod.SHIFT[m, n]))
                    for m, n in mod.Branch for t in mod.Time) +
                sum((1. - model.v[b, t])**2 for b in model.Bus for t in model.Time) +
                sum((model.theta[m, t] - model.theta[n, t])**2
                    for m, n in model.Branch for t in model.Time))

        @make_objective(maximize)
        def obj_maximize_losses(self):
            """Objective function - Maximize branch losses."""
            mod = model
            return sum(
                (mod.v[m, t]**2-mod.v[n, t]**2) * mod.Y[m, n] * cos(mod.SHIFT[m, n]) -
                mod.v[m, t] * mod.v[n, t] * (
                    cos(mod.theta[m, t] - mod.theta[n, t] - mod.SHIFT[m, n]) -
                    cos(mod.theta[n, t] - mod.theta[m, t] - mod.SHIFT[m, n]))
                for m, n in mod.Branch for t in mod.Time)

        @make_objective(minimize)
        def obj_minimize_voltage_deviation(self):
            """Objective function - minimize voltage deviation."""
            mod = model
            return sum(
                (1 - mod.v[b, t])**2
                for b in mod.Bus for t in mod.Time)

        @make_objective(minimize)
        def obj_minimize_angle_diffs(self):
            """Objective function - minimize angle diffs."""
            mod = model
            return sum(
                (mod.theta[m, t] - mod.theta[n, t])**2
                for m, n in mod.Branch for t in mod.Time)

        @make_objective(minimize)
        def obj_minimize_p_cost(self):
            """Objective function - minimize p cost."""
            mod = model
            return sum(
                mod.cost_p[g, t] * mod.p_gen[g, t]
                for g in mod.Gen for t in mod.Time)

        @make_objective(minimize)
        def obj_minimize_z(self):
            """Objective function - minimize cost.

            This is a generic objective function, z has to be defined a cost_rule as a constraint.
            """
            mod = model
            return sum(
                mod.z[i, k]
                for i in mod.I for k in mod.K)

        return(model)
