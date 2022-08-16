from fastapi import FastAPI, HTTPException
from model import FactorRequest, FactorResponse
from utils.poly_parse import parse_output
import uvicorn
import subprocess

app = FastAPI(
    title="Polynomial Factorization over Z[x]",
    description="An implementation of LLL for polynomial factorization",
    version="1.0",
)

tags_metadata = [{"name": "factor", "description": "Factor a monic integer polynomial"}]

# TODO: make html factor endpoint

@app.post(
    "/factor",
    tags=["factor"],
    responses={
        200: {"model": FactorResponse},
    },
)
def factor(request: FactorRequest):
    """Factor a polynomial and return result as json response."""

    # Get data from request
    binput = request.poly
    command = request.opts.to_list(binput)

    # Execute
    print(command, flush=True)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exit_code = proc.wait()

    # Catch C-level errors
    if exit_code == 1:
        raise HTTPException(status_code=500, detail=stderr.decode("utf-8"))

    # Success
    response = parse_output(stdout.decode("utf-8"))
    return FactorResponse(factors=response["factors"], time=response["time"])


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5000)
