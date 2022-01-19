from flask import Flask, request
import werkzeug
from utils.poly_parse import parse_poly, parse_opts
import os
import json
import subprocess

MAX_DEG = 20
ALLOWED_OPTS = {"d","v","t","p","newline"}

app = Flask(__name__)

@app.route("/factor",methods=["POST"])
def factor():
    body = json.loads(request.data.decode('utf-8'))

    # Check for polynomial input
    try:
        poly_str = body['poly']
        binput = parse_poly(poly_str)
    except Exception as e:
        raise werkzeug.exceptions.BadRequest(e)

    # Check for options input
    try:
        opts = parse_opts(body.get('opts'),ALLOWED_OPTS)
    except AssertionError as e:
        raise werkzeug.exceptions.BadRequest(e)

    # Limit degree
    degree = len(binput.split(",")) - 2
    if degree  > MAX_DEG:
        raise werkzeug.exceptions.RequestEntityTooLarge

    # Execute
    print(["./bin/factor_poly",binput]+opts,flush=True)
    proc = subprocess.Popen(["./bin/factor_poly",binput]+opts,stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    return out


@app.errorhandler(werkzeug.exceptions.BadRequest)
def handle_bad_request(e):
    return f"Bad request: {e}",400

@app.errorhandler(werkzeug.exceptions.RequestEntityTooLarge)
def handle_poly_too_large(e):
    return f"Polynomial degree exceeded (maximum degree allowed: {MAX_DEG})",413


if __name__ == "__main__":
    app.run(host='0.0.0.0') 
