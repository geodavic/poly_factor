from flask import Flask, request
import werkzeug
from utils.poly_parse import parse_poly, parse_opts, parse_output, parse_output_html
import os
import json
import subprocess

MAX_DEG = 300
ALLOWED_OPTS = {"delta":"-d","precision":"-p"}
OPTS_DEFAULTS = {"delta":0.5,"precision":64}

app = Flask(__name__)

# TODO : waiting page in iframe
#        OPTIONS method

@app.after_request
def after_request(response):
  response.headers.add('Access-Control-Allow-Origin', '*')
  response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
  response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
  return response

@app.route("/factor",methods=["POST"])
def factor():

    got_form = False
    if request.data:
        body = json.loads(request.data.decode('utf-8'))
    else:       
        body = dict(request.form)
        got_form = True

    # Check if verbose output is needed
    verbose = body.pop("verbose","off")

    # Check for polynomial input
    try:
        poly_str = body['poly']
        binput = parse_poly(poly_str)
    except Exception as e:
        raise werkzeug.exceptions.BadRequest(e)

    # Check for options input
    try:
        opts_dict = body.get('opts') or {k:v for k,v in body.items() if k!= 'poly'}
        opts = parse_opts(opts_dict,ALLOWED_OPTS,OPTS_DEFAULTS)
    except Exception as e:
        raise werkzeug.exceptions.BadRequest(e)

    # Limit degree
    degree = len(binput.split(",")) - 2
    if degree  > MAX_DEG:
        raise werkzeug.exceptions.RequestEntityTooLarge

    # Execute
    command = ["./bin/factor_poly",binput]+opts+["-t","-v","-newline"]
    print(command,flush=True)
    proc = subprocess.Popen(command,stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    exit_code = proc.wait()

    # Catch exit code for any C-level errors
    if exit_code == 1:
        raise werkzeug.exceptions.InternalServerError(out)

    # format output
    if got_form:
        out = parse_output_html(out.decode("utf-8"),verbose=verbose)
    else:
        out = parse_output(out.decode("utf-8"))

    return out


@app.errorhandler(werkzeug.exceptions.BadRequest)
def handle_bad_request(e):
    return f"Bad request: {e}",400

@app.errorhandler(werkzeug.exceptions.RequestEntityTooLarge)
def handle_poly_too_large(e):
    return f"Polynomial degree exceeded (maximum degree allowed: {MAX_DEG})",413


if __name__ == "__main__":
    app.run(host='0.0.0.0') 
