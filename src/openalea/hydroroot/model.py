from dataclasses import dataclass

from openalea.metafspm.component import Model, declare

from openalea.hydroroot.generator import markov
from openalea.hydroroot import radius
from openalea.hydroroot.law import length_law

@dataclass
class RsaModel(Model):
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

    segment_length: float = declare(default=1.0e-4, unit="m", unit_comment="", description="vertices lentgh",
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

    # Plant scale variable
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
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)
        self.link_self_to_mtg()

    def generate_mtg(self):
        """
        generate a MTG according to the input parameters using self.props
        """
        nb_vertices = int(self.primary_length / self.segment_length)
        # some distance expressed in term of nb of vertices
        nb_nude_vertices = int(self.nude_length / self.segment_length)
        nb_branching_delay = int(self.branching_delay / self.segment_length)

        if isinstance(self.length_data, list):
            length_max_secondary = self.length_data[0].iloc[:, 0]
            law_order1 = length_law(self.length_data[0], scale_x=self.primary_length / 100., scale=self.segment_length)
            law_order2 = length_law(self.length_data[1], scale_x=length_max_secondary / 100., scale=self.segment_length)
            _length_law = [law_order1, law_order2]

        self.g = markov.markov_binary_tree(nb_vertices=nb_vertices,
                                           branching_variability=self.branching_variability,
                                           branching_delay=nb_branching_delay,
                                           length_law=length_law,
                                           nude_tip_length=nb_nude_vertices,
                                           order_max=self.order_max,
                                           seed=self.seed)
        _g, self.surface = radius.compute_surface(self.g)
        _g, self.volume = radius.compute_volume(self.g)

    def set_mtg(self, _g):
        """
        set self.g and compute surface, volume, etc.
        :param _g: (MTG) existing mtg, for instance, comming from file reading
        """
        self.g = _g
        _g, self.surface = radius.compute_surface(self.g)
        _g, self.volume = radius.compute_volume(self.g)

