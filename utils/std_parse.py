import re


def parse_output(out):
    """Parse the verbose output of lll_factor.
    Keys in return: factors (list), time (float)
    """
    answers = out.split("Factorization:")[-1].split("\n")
    answers = [s for s in answers if s]

    time = answers[-1].split(":")[-1].strip()
    factors = answers[:-1]

    rval = {"time": time, "factors": factors}
    return rval


def parse_output_html(out, verbose=True, failed=False):
    """Parse the verbose output of lll_factor to an html string."""
    font_family = "Courier New"

    if verbose:
        body = out.replace("\n", "<br>")
        body = re.sub(r"={10,}<br>", "<hr>", body)
        if "Factorization:" in body:
            body = body.replace("Factorization:", "<b>Factorization:")
            body += "</b>"
        body = body.replace("-->", "<b>&gt;&gt;&nbsp;")
        body = body.replace("<--", "&nbsp; &lt;&lt;</b>")
    else:
        body = "".join(parse_output(out)["factors"])

    error_str = ""
    if failed:
        error_str = "<h2> Factorization failed: </h2>"

    rval = (
        f'<html style="font-family: {font_family}; font-size:12">'
        + error_str
        + body
        + "</html>"
    )
    return rval
