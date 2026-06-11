from openalea.hydroroot.init_parameter import Parameters
from openalea.core.config import Config

class HydroRootModelUnit:
    def __init__(self, parameters: Parameters):
        self.parameters = parameters

    def to_dict(self):
        return {
            "archi": self.parameters.archi,
            "hydro": self.parameters.hydro,
            "solute": self.parameters.solute,
            "exp": self.parameters.exp,
            "output": self.parameters.output,
        }
