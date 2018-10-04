"""Implements the OPF model object, used to encapsulate a given OPF problem formulation."""

import importlib
import time

from pyomo.environ import ConcreteModel

from pandapower.opf.models import objective

try:
    import pplog as logging
except ImportError:
    import logging

logger = logging.getLogger(__name__)


class aacopf_model(object):
    """Implement the OPF model object, used to encapsulate a given OPF problem formulation.

    Analogously to pandapower.opf.opf_model

    """

    def __init__(self, ppc, ppci, aacopop, data_file=None):
        """Use a PYPOWER case dict to build the pyomo optimisation model object.

        The constraints are loaded from the file **model_file**, that has to contain the
        **create_Constraint** function and implement a specific simlevel or bilevel AC OPF
        algorithm.

        aacopop
        -------
        # opt_level     : 'simlevel' or 'bilevel'
        model_file    : name of the file in resilience.powersystem.models that
            contains the constraint definitions
        obj_func      : objective function(s) for simlevel (SL), bilevel
            upper level (UL) and bilevel lower level (LL) to load from objective
        solveroptions : solver and options to be passed to the solver

        """
        self.ppc = ppc
        self.ppci = ppci
        self.from_ppci = aacopop['from_ppci']
        self.model_file = 'pandapower.opf.models.' + aacopop['model_file']
        self.definitions_file = 'pandapower.opf.models.' + aacopop['definitions_file']
        self.obj_func = aacopop['obj_func']
        self.solveroptions = aacopop['solveroptions']

        self.log_file = aacopop['model_file']

        # TODO:(aacopf_model): implement loading the optimisation data from a .dat file (speedup)
        self.data_file = data_file

        # NOTE(imports): this imports th .py file containing the sets, parameters and variables
        self.definitions = importlib.import_module(self.definitions_file)
        # NOTE(imports): this imports th .py file containing the constraints
        self.mod = importlib.import_module(self.model_file)

        self.model = ConcreteModel()

    def __repr__(self):
        """Return a string representation of the object."""
        # TODO:(__repr__): add a readable output representation of the om object
        pass

    def print_model(self, label=""):
        """Return the pyomo model to a log file."""
        file_name = "{}_{}.log".format(self.log_file, label)
        with open(file_name, "w") as f:
            f.write("Content of the pyomo model {}\n--{}--\n".format(time.ctime(), file_name))
            self.model.pprint(ostream=f)
            logger.info("Pyomo model written to {}".format(file_name))

    def model_file_consistency_check(self):
        """Make sure, that all sets, parameters and variables a model_file needs are available."""
        # TODO:(checks): add consistency checks for the different inputs (model_file, set_pa.., ...)
        pass

    def add_set_parameter_var(self):
        """Add all sets, parameters and variables from definitions_file to the om."""
        if self.from_ppci:
            self.model = self.definitions.set_and_parameter_from_ppci(self.model, self.ppci)
        else:
            self.model = self.definitions.set_and_parameter_from_dat(self.model)

        self.model = self.definitions.create_variable(self.model)
        logger.info("Sets, Parameters and Variables from {} added.".format(self.log_file))

    def add_constraints(self):
        """Add all contraints from a model_file to the om."""
        self.model = self.mod.create_Constraint(self.model)
        logger.info("Constraints from {} added.".format(self.log_file))

    def add_objective(self):
        """Add an objective function to the om.

        For the SimLevel optimisation this function
            1. adds a SL objective function

        """
        obj_SL = self.obj_func['SL']
        objective.CreateObjective(self.model, self.model, obj_SL).create_Objective()

        logger.info("Simlevel pyomo objective function added.")

    def set_solveroptions(self, solveroptions):
        """Set solver options."""
        self.solveroptions = solveroptions

    def set_result(self, results):
        """Set results."""
        self.results = results

    def set_raw(self, raw):
        """Set raw."""
        self.raw = raw

    def set_et(self, et):
        """Set et."""
        self.et = et

    def set_success(self, success):
        """Set success."""
        self.success = success

    def set_instance(self, instance):
        """Set an instance containing the results of the model object."""
        self.instance = instance

    def set_status(self, status):
        """Set the solver status  of the model object."""
        self.status = status

    def get_ppc(self):
        """Return the PYPOWER case dict."""
        return self.ppc

    def get_ppci(self):
        """Return the ppci."""
        return self.ppci
