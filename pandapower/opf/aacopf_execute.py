"""Execute the OPF."""

from pyomo.opt import SolverFactory

from pandapower.optimal_powerflow import warnings

# try:
#     import pplog as logging
# except ImportError:
#     import logging
#
# logger = logging.getLogger(__name__)


def aacopf_execute(om, ppopt, suppress_warnings=True):
    """Execute the OPF specified by an OPF model object.

    Analogously to pandapower.opf.opf_execute.py

    """
    if suppress_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            om = _solve_aacopf(om)
    else:
        om = _solve_aacopf(om)

    # TODO: wirte optimisation status of success
    success = True
    # TODO: wirte optimisation raw output
    raw = None

    return(success, raw)


def _solve_aacopf(om):
    """Create, run and return results of the pyomo algebraic OPF model.

    result = opf(ppci, ppopt)
    """
    opt = SolverFactory(om.solveroptions['solver'], tee=False)

    for key, value in om.solveroptions['solver-options'].items():
        setattr(opt.options, key, value)

    instance = om.model
    om.print_model(label="unsolved")
    status = opt.solve(instance, tee=False, logfile=om.log_file + "_solver.log")

    om.print_model(label="results")

    om.set_instance(instance)
    om.set_status(status)

    return(om)
