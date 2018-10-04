"""Solves an optimal power flow."""

from time import time
import numpy as np

from pyomo.environ import Var

from pandapower.opf.aacopf_execute import aacopf_execute
from pandapower.opf.aacopf_setup import aacopf_setup

try:
    import pplog as logging
except ImportError:
    import logging

logger = logging.getLogger(__name__)


def aacopf(aacopop, ppc, ppci, ppopt, **kwargs):
    """Solves an algebraic AC optimal power flow.

    Analogously to pandapower.opf.opf.py

    """
    t0 = time()  # start timer
    om = aacopf_setup(ppc, ppci, ppopt, aacopop)  # construct OPF model object
    success, raw = aacopf_execute(om, ppopt)  # execute the OPF
    results = _instance2result(om.instance, om.status, om.ppci)  # Extract the results form the om
    et = time() - t0  # compute elapsed time

    results['et'] = et
    results['success'] = success
    results['raw'] = raw

    return results


def _instance2result(instance, status, ppci):
    """Load results from instance and status and write them to result.

    # NOTE(_instance2result): all data from the optimisation handed to result for developement
    """
    def gen_res(result):
        """Workaround to catch non set generator power before scaling back."""
        # TODO(q_gen): calculate q_gen for PQ buses with gens
        if None in result:
            result[:] = np.nan
            logger.warn("Gen result table contains NaN")
        else:
            result *= BASE_MVA
        return(result)

    res = {
        'aacopf': {
            'Problem': {},
            'Solver': {'Status': 'ok', 'Time': 42, 'Id': None},
            'Solution': {}},
        'Var': {}}
    for i in instance.component_objects(Var, active=True):
        varobject = getattr(instance, str(i))
        res['Var'][varobject.getname()] = varobject.get_values()

    for i, j in status.items():
        for k in j.keys():
            res['aacopf'][i][k] = getattr(status[i], k)

    # TODO(results): select data needed (index) an delete all other values -> more explicit
    result = ppci

    BASE_MVA = ppci["baseMVA"] * 1e+03
    # TODO(result): manage the solver output which may differ from solver to solver
    result['success'] = res['aacopf']['Solver']['Status']
    result['f'] = getattr(status.Solution, 'Objective')
    result['et'] = res['aacopf']['Solver']['Time']
    result['internal']['aacopf'] = res['aacopf']

    # TODO(instance2result): refactor to speed up writing of the values to result
    # TODO(infeasibilities): shift this check that keep values from being written to results
    if res['aacopf']['Solver']['Id'] is 220:
        logger.error("{}".format("infesible -> see solver warning"))
        # NOTE(infeasibilities): Workaround to identify that the problem was infeasible
        result['success'] = False
        result["bus"][:] = np.nan
        result["branch"][:] = np.nan
        result["gen"][:] = np.nan
    else:
        # TODO(sign): check if all signs are correctly written to results
        bus_index = result["bus"][:, 0].astype(int)
        result["bus"][:, 7] = np.array([res['Var']['v'][(k, 0)] for k in bus_index])
        angles_in_rad = np.array([res['Var']['theta'][(k, 0)] for k in bus_index])
        if any(angles_in_rad / (2*np.pi) >= 1):
            logger.warn("{}".format("angles are to big!"))
        result["bus"][:, 8] = angles_in_rad * 180 / np.pi

        # NOTE(gen_index): Do not use the gens_buses as the index, but the set of gens
        # gen_index = result["gen"][:, 0].astype(int)
        gen_index = np.arange(result["gen"][:, 0].size, dtype=int)

        result["gen"][:, 1] = gen_res(np.array([res['Var']['p_gen'][(k, 0)] for k in gen_index]))
        result["gen"][:, 2] = gen_res(np.array([res['Var']['q_gen'][(k, 0)] for k in gen_index]))

        branch_index = result["branch"][:, 0:2].astype(int)
        # TODO(BASE_MVA): check if those values need to be multiplied by BASE_MVA
        result["branch"][:, 13] = np.array([res['Var']['p'][(k, l, 0)] for k, l in branch_index]) * BASE_MVA
        result["branch"][:, 14] = np.array([res['Var']['q'][(k, l, 0)] for k, l in branch_index]) * BASE_MVA
        result["branch"][:, 15] = np.array([res['Var']['p'][(l, k, 0)] for k, l in branch_index]) * BASE_MVA
        result["branch"][:, 16] = np.array([res['Var']['q'][(l, k, 0)] for k, l in branch_index]) * BASE_MVA

    return(result)
