""" 
This file registers all allowed algorithms for factorization,
including their parameters and default values.
"""

from pydantic import BaseModel, validator
from typing import Union, Optional
import sys
import inspect

ALLOWED_ALGS = ["LLL"]
DEFAULT_ALG = "LLL"
MAX_DEG = 300
LLL_CAP = 20


class LLLOptions(BaseModel):
    alg = "LLL"
    precision: Optional[int] = 64
    delta: Optional[float] = 0.5

    def to_list(cls, input_polynomial: str):
        """Return the command used to run this algorithm."""
        return [
            "./bin/lll_factor",
            str(input_polynomial),
            "-d",
            str(cls.delta),
            "-p",
            str(cls.precision),
            "-t",
            "-v",
            "-newline",
            "-stop",
            str(LLL_CAP),
        ]

    @validator("precision")
    def precision_in_range(cls, precision):
        assert precision >= 32, "Must have at least 32 bits of precision"
        return precision

    @validator("delta")
    def delta_in_range(cls, delta):
        assert delta > 0.25 and delta < 1, "delta parameter must be in (0.25,1)"
        return delta

    @validator("alg")
    def alg_match(cls, alg, values):
        assert alg in ALLOWED_ALGS, f"Unrecognized algorithm {alg}"
        assert alg == "LLL", "Mismatching options for specified algorithm"
        return alg


"""
Auto-create OptType and DefaultOption
"""

current_module = sys.modules[__name__]
option_types = ()
DefaultOption = None
for name, obj in inspect.getmembers(sys.modules[__name__]):
    if inspect.isclass(obj) and hasattr(obj(), "alg"):
        option_types += (obj,)
        if obj().alg == DEFAULT_ALG:
            DefaultOption = obj
OptsType = Union[option_types]
assert DefaultOption is not None, f"Warning: Default option {DEFAULT_ALG} not in allowed algorithms list: {ALLOWED_ALGS}."
