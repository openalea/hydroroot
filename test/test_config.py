from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.model_unit import HydroRootModelUnit
from openalea.core.config import Config

params = Parameters()
unit = HydroRootModelUnit(params)

cfg = Config([unit])

cfg.dump("hydro_config.json")
