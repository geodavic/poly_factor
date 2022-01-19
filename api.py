from flask import Flask, request
from utils.poly_parse import parse_poly
import os
import json
import subprocess

app = Flask(__name__)

@app.route("/factor",methods=["POST"])
def factor():
    body = json.loads(request.data.decode('utf-8'))

    poly_str = body['poly']
    opts = body.get('opts',[])
    binput = parse_poly(poly_str)

    print(["./bin/factor_poly",binput]+opts,flush=True)
    proc = subprocess.Popen(["./bin/factor_poly",binput]+opts,stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    return out

if __name__ == "__main__":
    app.run(host='0.0.0.0') 
