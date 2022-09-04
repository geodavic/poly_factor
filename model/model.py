from pydantic import BaseModel, validator, Extra
from typing import Optional
from utils.poly_parse import parse_poly
from fastapi import Form
from .opts_model import OptsType, DefaultOption, LLLOptions, MAX_DEG


class FactorRequest(BaseModel):
    poly: str
    opts: Optional[OptsType] = DefaultOption()

    @validator("poly")
    def polynomial_parses(cls, poly):
        try:
            parsed = parse_poly(poly, max_deg=MAX_DEG)
            return parsed
        except Exception as e:
            raise ValueError(str(e))


class LLLFormData(BaseModel, extra=Extra.allow):
    """Form data request for LLL Factorization. Slight hack to
    run it through the FactorRequest and Options validations.

    I wish there was a way to make a generalized form data class...
    """

    poly: str
    precision: str = None
    delta: str = None
    verbose: str = None

    @classmethod
    def as_form(
        cls,
        poly: str = Form(...),
        precision: str = Form(None),
        delta: str = Form(None),
        verbose: str = Form(None),
    ):
        opt_kw = {"precision": precision, "delta": delta}
        opt_kw = {k: v for k, v in opt_kw.items() if v is not None}
        opts = LLLOptions(**opt_kw)

        rq = FactorRequest(poly=poly, opts=opts)
        obj = cls(
            poly=rq.poly,
            precision=opts.precision,
            delta=opts.delta,
            verbose=verbose,
            opts=opts,
        )
        return obj


class FactorResponse(BaseModel):
    factors: list
    time: str
