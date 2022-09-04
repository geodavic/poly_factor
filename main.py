from fastapi import FastAPI, Depends, HTTPException, Form
from fastapi.responses import HTMLResponse
from model import FactorRequest, FactorResponse, LLLFormData
from opts_model import OptsType
from utils.std_parse import parse_output, parse_output_html
import uvicorn
import subprocess


tags_metadata = [
    {
        "name": "factor",
        "description": "The main endpoint",
        "externalDocs": {
            "description": "Source code",
            "url": "https://github.com/geodavic/poly_factor",
        },
    },
    {"name": "web", "description": "Endpoints for the Web UI"},
]
app = FastAPI(
    title="Polynomial Factorization over Z[x]",
    description="An implementation of LLL for polynomial factorization.",
    version="0.1",
    contact={
        "name": "George D. Torres",
        "url": "http://web.ma.utexas.edu/users/gdavtor",
    },
    openapi_tags=tags_metadata,
)


def base_factor(poly: str, opts: OptsType):
    """Base LLL request processor"""
    # Get data from request
    command = opts.to_list(poly)

    # Execute
    print(command, flush=True)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exit_code = proc.wait()

    # Catch C-level errors
    if exit_code == 1:
        raise HTTPException(status_code=500, detail=stderr.decode("utf-8"))

    return stdout.decode("utf-8")


@app.post(
    "/factor",
    tags=["factor"],
    responses={
        200: {"model": FactorResponse},
    },
)
def factor(request: FactorRequest):
    """Factor a polynomial and return result as json response."""

    out = base_factor(request.poly, request.opts)
    response = parse_output(out)
    return FactorResponse(factors=response["factors"], time=response["time"])


@app.post("/lll_form_data_factor", tags=["web"])
def lll_form_data_factor(request: LLLFormData = Depends(LLLFormData.as_form)):
    """Factor a polynomial using LLL and return result as an html string."""

    failed = False
    try:
        opts = request.opts
        out = base_factor(request.poly, opts)
    except HTTPException as e:
        failed = True
        out = str(e)

    html = parse_output_html(out, verbose=request.verbose, failed=failed)
    return HTMLResponse(content=html, status_code=200)


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5000)
