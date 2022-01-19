import flask

app = flask.Flask(__name__)

@app.route("/factor",methods=["GET"])
def factor():
    return "test"

if __name__ == "__main__":
    app.run() 
