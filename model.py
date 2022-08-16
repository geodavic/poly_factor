from pydantic import BaseModel, validator
from typing import Union, Optional
from utils.poly_parse import parse_poly

# TODO: validation on first two classes

ALLOWED_ALGS = ["LLL"]
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


class FactorRequest(BaseModel):
    poly: str
    opts: Optional[Union[LLLOptions]] = LLLOptions()

    @validator("poly")
    def polynomial_parses(cls, poly):
        try:
            parsed = parse_poly(poly, max_deg=MAX_DEG)
            return parsed
        except Exception as e:
            raise ValueError(str(e))


class FactorResponse(BaseModel):
    factors: list
    time: str


class HTMLFactorResponse(BaseModel):
    html: str
