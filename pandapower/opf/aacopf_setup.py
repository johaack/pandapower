"""Constructs a pyomo OPF model object from a PYPOWER case dict and the passed options."""

from pandapower.opf.aacopf_model import aacopf_model

# try:
#     import pplog as logging
# except ImportError:
#     import logging
#
# logger = logging.getLogger(__name__)


def aacopf_setup(ppc, ppci, ppopt, aacopop, **kwargs):
    """Construct an AACOPF model object from a PYPOWER case dict.

    Analogously to pandapower.opf.opf_setup.py

    """
    om = aacopf_model(ppc, ppci, aacopop)
    om.add_set_parameter_var()
    om.add_objective()
    om.add_constraints()
    return om
