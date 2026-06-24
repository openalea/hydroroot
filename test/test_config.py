
from openalea.hydroroot.init_parameter import Parameters, hydroroot_parameters_to_config
from openalea.core.config import Config, ModelUnit, Parameter
import yaml


params = Parameters()
cfg1 = hydroroot_parameters_to_config(params)
cfg1.custom_comments["length_file"] = "length file"
    
cfg1.dump("hydro_config4.json")
cfg1.dump("hydro.yml")
cfg2 = Config([]).load("hydro_config4.json")
assert dict(cfg1) == dict(cfg2)




