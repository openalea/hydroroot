import math

from dataclasses import dataclass

from openalea.mtg import traversal
from openalea.metafspm.component import Model, declare
from openalea.metafspm.component_factory import *

from openalea.hydroroot import radius, conductance, main, flux
from openalea.hydroroot.generator import markov
from openalea.hydroroot.law import length_law
from openalea.hydroroot.length import fit_law
from openalea.hydroroot.water_solute_transport import (pressure_calculation_no_non_permeating_solutes,
                                                       init_some_MTG_properties, pressure_calculation)
from openalea.hydroroot.analysis import intercept

EPS = 1.0e-9 # epsilon value used for convergence of Newton-Raphson loop

@dataclass
class HydroRootModel(Model):
    # Input Parameters related to architecture
    length_data: list = declare(default=None, unit="[adim,m]", unit_comment="relative distance to the tip, meter",
                                description="list of float or list of 2 list of float, lateral length vs relat. dist. to tip",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    primary_length: float = declare(default=0.16, unit="m", unit_comment="", description="length of the primary root",
                                    min_value="0", max_value="", value_comment="", references="", DOI="",
                                    variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    seed: int = declare(default=None, unit="", unit_comment="", description="seed used to generate the architecture",
                        min_value="0", max_value="", value_comment="", references="", DOI="",
                        variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    branching_delay: float = declare(default=2.0e-3, unit="m", unit_comment="",
                                     description="distance between branching laterals",
                                     min_value="0", max_value="", value_comment="", references="", DOI="",
                                     variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    branching_variability: float = declare(default=0.25, unit="", unit_comment="",
                                           description="variability of the branching laterals distance",
                                           min_value="0", max_value="1", value_comment="", references="", DOI="",
                                           variable_type="parameter", by="RsaModel", state_variable_type="",
                                           edit_by="user")

    order_max: int = declare(default=4, unit="", unit_comment="", description="maximum order of laterals",
                             min_value="0", max_value="1", value_comment="", references="", DOI="",
                             variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    length: float = declare(default=1.0e-4, unit="m", unit_comment="", description="vertices lentgh",
                                    min_value="0", max_value="", value_comment="", references="", DOI="",
                                    variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    nude_length: float = declare(default=0.02, unit="m", unit_comment="",
                                 description="part of roots without any lateral root, distance from tip",
                                 min_value="0", max_value="", value_comment="", references="", DOI="",
                                 variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    ref_radius: float = declare(default=7.0e-5, unit="m", unit_comment="",
                                description="reference radius or radius of the primary root",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    order_decrease_factor: float = declare(default=0.7, unit="", unit_comment="",
                                           description="radius decrease factor applied when increasing order",
                                           min_value="0", max_value="1", value_comment="", references="", DOI="",
                                           variable_type="parameter", by="RsaModel", state_variable_type="",
                                           edit_by="user")

    # hydraulic parameters
    axial_conductance_data: list = declare(default=None, unit="[m,10-9 m4.MPa-1.s-1]", unit_comment="distance to the tip, K",
                                description="(2 list of Float) axial conductivity versus dist. to tip",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    k0: list = declare(default=None, unit="[m,10-9 m.MPa-1.s-1]",
                                            unit_comment="distance to the tip, k",
                                            description="(2 list of Float) radial conductivity versus dist. to tip",
                                            min_value="0", max_value="", value_comment="", references="", DOI="",
                                            variable_type="parameter", by="RsaModel", state_variable_type="",
                                            edit_by="user")

    # Soil1DModel input
    psi_e: float = declare(default=0.401325, unit="MPa", unit_comment="",
                           description="Homogeneous hydrostatic potential at the vertex boundary",
                           min_value="0", max_value="", value_comment="", references="", DOI="",
                           variable_type="input", by="Soil1DModel", state_variable_type="", edit_by="user")

    # HydrostaticModel parameter
    psi_base: float = declare(default=0.101325, unit="MPa", unit_comment="",
                              description="Hydrostatic potential at the root base",
                              min_value="0", max_value="", value_comment="", references="", DOI="",
                              variable_type="parameter", by="HydrostaticModel", state_variable_type="", edit_by="user")

    Jv: float = declare(default=None, unit="microL.s-1", unit_comment="",
                        description="water flux at the root base, input if invert_model=False ",
                        min_value="0", max_value="", value_comment="", references="", DOI="",
                        variable_type="parameter", by="HydrostaticModel", state_variable_type="", edit_by="user")

    invert_model: bool = declare(default=True, unit="", unit_comment="",
                                 description="when false, distribute output flux within the root ; "
                                             "when true, compute the output flux for the given root and conditions.",
                                 min_value="", max_value="", value_comment="", references="", DOI="",
                                 variable_type="parameter", by="HydrostaticModel", state_variable_type="",
                                 edit_by="user")

    cut_and_flow: bool = declare(default=False, unit="", unit_comment="",
                                 description="Use specific model to compute conductance at tips with cut & flow.",
                                 min_value="", max_value="", value_comment="", references="", DOI="",
                                 variable_type="parameter", by="HydrostaticModel", state_variable_type="",
                                 edit_by="user")

    # Solute parameters
    J_s: float = declare(default=1.0e-7, unit="mol.m-2.s-1", unit_comment="", description="Active pumping rate",
                         min_value="0", max_value="", value_comment="", references="", DOI="",
                         variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    P_s: float = declare(default=1.0e-9, unit="m.s-1", unit_comment="", description="permeability coefficient",
                         min_value="0", max_value="", value_comment="", references="", DOI="",
                         variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    Cse: float = declare(default=13.96, unit="mol.m-3", unit_comment="",
                         description="concentration of permeating solutes",
                         min_value="0", max_value="", value_comment="", references="", DOI="",
                         variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    Ce: float = declare(default=0.0, unit="mol.m-3", unit_comment="",
                        description="concentration of non-permeating solutes",
                        min_value="0", max_value="", value_comment="", references="", DOI="",
                        variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    sigma: float = declare(default=1.0, unit="", unit_comment="", description="reflexion coefficient",
                           min_value="0", max_value="", value_comment="", references="", DOI="",
                           variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    temp: float = declare(default=298.0, unit="K", unit_comment="", description="sap temperature in degree K",
                          min_value="", max_value="", value_comment="", references="", DOI="",
                          variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    c_base: float = declare(default=None, unit="mol.m-3", unit_comment="",
                            description="solute concentration at the root's base",
                            min_value="", max_value="", value_comment="", references="", DOI="",
                            variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    # Plant scale variable
    keq: float = declare(default=None, unit="10-9 m4.MPa-1.s-1", unit_comment="",
                         description="equivalent conductance of the RSA",
                         min_value="0", max_value="", value_comment="", references="", DOI="",
                         variable_type="plant_scale_state", by="HydrostaticModel", state_variable_type="",
                         edit_by="user")

    surface: float = declare(default=None, unit="m2", unit_comment="", description="total surface of the RSA",
                             min_value="0", max_value="", value_comment="", references="", DOI="",
                             variable_type="plant_scale_state", by="RsaModel", state_variable_type="", edit_by="user")

    volume: float = declare(default=None, unit="m3", unit_comment="", description="total surface of the RSA",
                            min_value="0", max_value="", value_comment="", references="", DOI="",
                            variable_type="plant_scale_state", by="RsaModel", state_variable_type="", edit_by="user")

    def __init__(self, g, time_step, **scenario):
        # 4/21/26 FB: copied from root_cynaps RootWaterModel
        self.g = g
        self.props = self.g.properties()
        self.time_step = time_step
        # self.choregrapher(module_family=self.__class__.__name__)
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)
        self.link_self_to_mtg()

    def generate_mtg(self):
        """
        generate a MTG according to the input parameters using self.props
        """
        nb_vertices = int(self.primary_length / self.length)
        # some distance expressed in term of nb of vertices
        nb_nude_vertices = int(self.nude_length / self.length)
        nb_branching_delay = int(self.branching_delay / self.length)

        if isinstance(self.length_data, list):
            length_max_secondary = self.length_data[0].iloc[:, 0]
            law_order1 = length_law(self.length_data[0], scale_x=self.primary_length / 100., scale=self.length)
            law_order2 = length_law(self.length_data[1], scale_x=length_max_secondary / 100., scale=self.length)
            _length_law = [law_order1, law_order2]

        self.g = markov.markov_binary_tree(nb_vertices=nb_vertices,
                                           branching_variability=self.branching_variability,
                                           branching_delay=nb_branching_delay,
                                           length_law=length_law,
                                           nude_tip_length=nb_nude_vertices,
                                           order_max=self.order_max,
                                           seed=self.seed)

        self.compute_radius()
        self.compute_length()
        self.compute_dist_to_tip()
        self.compute_surface_volume()
        self.compute_radius()

    @actual
    def compute_metrics(self):
        if not self.g.property('radius'):
            # TODO ref_radius etc.
            self.compute_radius()
        if not self.g.property('length'):
            self.compute_length()

        self.compute_dist_to_base()
        self.compute_dist_to_tip()
        self.compute_surface_volume()
        self.g.property('surface')[1] = self.surface
        self.g.property('volume')[1] = self.volume

    @actual
    @state
    def compute_length(self):
        """
        compute vertices length according to self.length
        :return:
        """
        self.g = radius.compute_length(self.g, self.length)

    @actual
    @state
    def compute_dist_to_base(self):
        """
        Calculation of the distance from base of each vertex, used for cut and flow
        :return:
        """
        g = self.g
        _dist_to_base = {}
        segment_length = self.length
        for v in traversal.pre_order2(g, 1):
            pid = g.parent(v)
            _dist_to_base[v] = _dist_to_base[pid] + segment_length if pid else segment_length
        g.properties()['dist_to_base'] = _dist_to_base

    @actual
    @state
    def compute_dist_to_tip(self):
        """
        For each vertex compute the distance from the tip of the axis bearing it in absolute
        and relative to the length of the axis bearing it
        :return:
        """
        self.g = radius.compute_relative_position(self.g)

    @actual
    @state
    def compute_radius(self):
        """
        compute the radius of each vertex from self.ref_radius with fixed decrease between each order.
        :return:
        """
        self.g = radius.ordered_radius(self.g, self.ref_radius, self.order_decrease_factor)

    @actual
    @state
    def compute_surface_volume(self):
        """
        set self.g and compute surface, volume, etc.
        """
        self.g, self.surface = radius.compute_surface(self.g)
        self.g, self.volume = radius.compute_volume(self.g)

    @actual
    @axial
    @state
    def compute_axial_conductance(self, data = None):
        """
        compute axial conductance of each vertex from self.ref_axial_conductivity
        :return:
        """
        if data:
            self.axial_conductance_data = data

        if self.axial_conductance_data:
            # Compute K using axial conductance data
            xa, ya = self.axial_conductance_data
            axial_conductivity_law = fit_law(xa, ya)

            self.g = conductance.fit_property_from_spline(self.g, axial_conductivity_law, 'position', 'K_exp')
            self.g = conductance.compute_K(self.g)

    @actual
    @state
    def compute_radial_conductance(self, data = None):
        """
        compute radial conductance of each vertex from self.ref_radial_conductivity
        :return:
        """
        if data:
            self.k0 = data

        if self.k0:
            xr, yr = self.k0
            radial_conductivity_law = fit_law(xr, yr)

            self.g = conductance.fit_property_from_spline(self.g, radial_conductivity_law, 'position', 'k0')
            self.g = conductance.compute_k(self.g, k0='k0')

    def compute_intercepts(self, dists, dl=1e-4, max_order=None):
        return intercept(self.g, dists, dl=dl, max_order=max_order)


    @actual
    @state
    def hydrostatic_solver_flux(self):
        """
        Compute flux according to psi_e and psi_base
        If self.psi_e is None use self.g.property('psi_e')
        :return:
        """
        _psi_e = self.psi_e
        f = flux.Flux(self.g, self.Jv, _psi_e, self.psi_base, self.invert_model, cut_and_flow=self.cut_and_flow)
        f.run()
        self.g = f.g
        v_base = next(self.g.component_roots_at_scale_iter(self.g.root, scale=1))
        self.keq = self.g.property('Keq')[v_base]
        self.Jv = self.g.property('J_out')[v_base]

    @actual
    @state
    def solute_solver_flux(self):
        """
        Compute flux according to psi_e and psi_base
        If self.psi_e is None use self.g.property('psi_e')
        :return:
        """

        # init the MTG with no solute
        _psi_e = self.psi_e
        f = flux.Flux(self.g, self.Jv, _psi_e, self.psi_base, self.invert_model, cut_and_flow=self.cut_and_flow)
        f.run()
        self.g = f.g

        self.g = init_some_MTG_properties(self.g, tau=self.J_s, Cini=self.Cse*1.0e-9, Cpeg_ini=self.Ce*1.0e-9, t=1, Ps=self.P_s)

        if self.Ce <= 0.0:
            calculation = pressure_calculation_no_non_permeating_solutes
        else:
            calculation = pressure_calculation

        # Newton-Raphson loop
        nb_v = self.g.nb_vertices()
        Fdx = 1.0
        Fdx_old = 1.
        while Fdx > EPS:
            g, dx, data, row, col = calculation(g, Temp=self.temp, sigma=self.sigma, Ce=self.Ce*1.0e-9, Cse=self.Cse*1.0e-9,
                                                Pe=self.psi_e, Pbase=self.psi_base, C_base=self.c_base*1.0e-9)
            Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
            if abs(Fdx - Fdx_old) < EPS: break
            Fdx_old = Fdx

        v_base = next(self.g.component_roots_at_scale_iter(self.g.root, scale=1))
        self.Jv = self.g.property('J_out')[v_base]

@dataclass
class HydrostaticModel(Model):
    # Hydraulic input from RsaModel
    length: float = declare(default=1.0e-4, unit="m", unit_comment="", description="vertices lentgh",
                                    min_value="0", max_value="", value_comment="", references="", DOI="",
                                    variable_type="parameter", by="RsaModel", state_variable_type="", edit_by="user")

    axial_conductance_data: list = declare(default=None, unit="[m,10-9 m4.MPa-1.s-1]", unit_comment="distance to the tip, K",
                                description="(2 list of Float) axial conductivity versus dist. to tip",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="RsaModel", state_variable_type="", edit_by="user")

    k0: list = declare(default=None, unit="[m,10-9 m.MPa-1.s-1]", unit_comment="distance to the tip, k",
                                            description="(2 list of Float) radial conductivity versus dist. to tip",
                                            min_value="0", max_value="", value_comment="", references="", DOI="",
                                            variable_type="input", by="RsaModel", state_variable_type="",
                                            edit_by="user")

    # Soil1DModel input
    psi_e: float = declare(default=0.401325, unit="MPa", unit_comment="",
                                description="Homogeneous hydrostatic potential at the vertex boundary",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="Soil1DModel", state_variable_type="", edit_by="user")

    # HydrostaticModel parameter
    psi_base: float = declare(default=0.101325, unit="MPa", unit_comment="",
                                description="Hydrostatic potential at the root base",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="HydrostaticModel", state_variable_type="", edit_by="user")

    Jv: float = declare(default=None, unit="microL.s-1", unit_comment="",
                                description="water flux at the root base, input if invert_model=False ",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="HydrostaticModel", state_variable_type="", edit_by="user")

    invert_model: bool = declare(default=True, unit="", unit_comment="",
                                description="when false, distribute output flux within the root ; "
                                            "when true, compute the output flux for the given root and conditions.",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="HydrostaticModel", state_variable_type="", edit_by="user")

    cut_and_flow: bool = declare(default=False, unit="", unit_comment="",
                                description="Use specific model to compute conductance at tips with cut & flow.",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="HydrostaticModel", state_variable_type="", edit_by="user")

    keq: float = declare(default=None, unit="10-9 m4.MPa-1.s-1", unit_comment="", description="equivalent conductance of the RSA",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="plant_scale_state", by="HydrostaticModel", state_variable_type="", edit_by="user")


    def __init__(self, g, time_step, **scenario):
        # 4/21/26 FB: copied from root_cynaps RootWaterModel
        self.g = g
        self.props = self.g.properties()
        self.time_step = time_step
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)
        self.link_self_to_mtg()


@dataclass
class SoluteModel(Model):
    # Hydraulic input from RsaModel
    axial_conductance_data: list = declare(default=None, unit="[m,10-9 m4.MPa-1.s-1]", unit_comment="distance to the tip, K",
                                description="(2 list of Float) axial conductivity versus dist. to tip",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="RsaModel", state_variable_type="", edit_by="user")

    k0: list = declare(default=None, unit="[m,10-9 m.MPa-1.s-1]", unit_comment="distance to the tip, k",
                                            description="(2 list of Float) radial conductivity versus dist. to tip",
                                            min_value="0", max_value="", value_comment="", references="", DOI="",
                                            variable_type="input", by="RsaModel", state_variable_type="",
                                            edit_by="user")

    # Soil1DModel input
    psi_e: float = declare(default=0.401325, unit="MPa", unit_comment="",
                                description="Homogeneous hydrostatic potential at the vertex boundary",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="Soil1DModel", state_variable_type="", edit_by="user")

    # SoluteModel parameter
    psi_base: float = declare(default=0.101325, unit="MPa", unit_comment="",
                                description="Hydrostatic potential at the root base",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    Jv: float = declare(default=None, unit="microL.s-1", unit_comment="",
                                description="water flux at the root base, input if invert_model=False ",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    cut_and_flow: bool = declare(default=False, unit="", unit_comment="",
                                description="Use specific model to compute conductance at tips with cut & flow.",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    J_s: float = declare(default=1.0e-7, unit="mol.m-2.s-1", unit_comment="", description="Active pumping rate",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    P_s: float = declare(default=1.0e-9, unit="m.s-1", unit_comment="", description="permeability coefficient",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    Cse: float = declare(default=13.96, unit="mol.m-3", unit_comment="", description="concentration of permeating solutes",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    Ce: float = declare(default=0.0, unit="mol.m-3", unit_comment="", description="concentration of non-permeating solutes",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    sigma: float = declare(default=1.0, unit="", unit_comment="", description="reflexion coefficient",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    temp: float = declare(default=298.0, unit="K", unit_comment="", description="sap temperature in degree K",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    c_base: float = declare(default=None, unit="mol.m-3", unit_comment="", description="solute concentration at the root's base",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="SoluteModel", state_variable_type="", edit_by="user")

    def __init__(self, g, time_step, **scenario):
        # 4/21/26 FB: copied from root_cynaps RootWaterModel
        self.g = g
        self.props = self.g.properties()
        self.time_step = time_step
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)
        self.link_self_to_mtg()

    def flux(self):
        """
        Compute flux according to psi_e and psi_base
        If self.psi_e is None use self.g.property('psi_e')
        :return:
        """

        # init the MTG with no solute
        _psi_e = self.psi_e
        f = flux.Flux(self.g, self.Jv, _psi_e, self.psi_base, self.invert_model, cut_and_flow=self.cut_and_flow)
        f.run()
        self.g = f.g

        self.g = init_some_MTG_properties(self.g, tau=self.J_s, Cini=self.Cse*1.0e-9, Cpeg_ini=self.Ce*1.0e-9, t=1, Ps=self.P_s)

        if self.Ce <= 0.0:
            calculation = pressure_calculation_no_non_permeating_solutes
        else:
            calculation = pressure_calculation

        # Newton-Raphson loop
        nb_v = self.g.nb_vertices()
        Fdx = 1.0
        Fdx_old = 1.
        while Fdx > EPS:
            g, dx, data, row, col = calculation(g, Temp=self.temp, sigma=self.sigma, Ce=self.Ce*1.0e-9, Cse=self.Cse*1.0e-9,
                                                Pe=self.psi_e, Pbase=self.psi_base, C_base=self.c_base*1.0e-9)
            Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
            if abs(Fdx - Fdx_old) < EPS: break
            Fdx_old = Fdx

        v_base = next(self.g.component_roots_at_scale_iter(self.g.root, scale=1))
        self.Jv = self.g.property('J_out')[v_base]

@dataclass
class Soil1DModel(Model):
    # Soil1DModel parameters
    soil_data: tuple = declare(default=None, unit="[m, MPa]", unit_comment="",
                                description="tuple of 2 lists, (z,psi_e) z=depth, psi_e=water potential",
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="Soil1DModel", state_variable_type="", edit_by="user")

    # state variable
    psi_e: float = declare(default=0.401325, unit="MPa", unit_comment="",
                                description="Soil hydrostatic potential at the vertex boundary",
                                min_value="0", max_value="", value_comment="", references="", DOI="",
                                variable_type="state_variable", by="HydrostaticModel", state_variable_type="", edit_by="user")

    def __init__(self, g, time_step, **scenario):
        # 4/21/26 FB: copied from root_cynaps RootWaterModel
        self.g = g
        self.props = self.g.properties()
        self.time_step = time_step
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)
        self.link_self_to_mtg()

    def compute_doil_psi_e(self):
        """
        add a soil as heterogeneous water potential along z
        :return:
        """
        self.g = main.soil_1D(self.g, self.soil_data)
