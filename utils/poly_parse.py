import sys
import re
import numpy as np


class PolynomialConstraintError(Exception):
    def __init__(self,degree):
        self.degree = degree
        if degree == np.infty:
            self.message = f"Polynomial must be monic."
        else:
            self.message = f"Polynomial degree must not exceed {degree} and must be monic."
        super().__init__(self.message)


class PolynomialFormatError(Exception):
    def __init__(self):
        self.message = "Polynomial have integer coefficients, must be in the variable 'x' with nonnegative integer exponents, and be in standard form."
        super().__init__(self.message)



def parse_poly(polystr,max_deg=None):
    """
    parse polynomial into csv format
    e.g. x^4-1  ->  -1,0,0,1
    """
    if max_deg is None:
        max_deg = np.infty
    polystr = polystr.replace(" ","").lower()
    L=polystr.replace("-","+-").split("+")
    degree=0

    #regularize terms
    for i in range(len(L)):
        term=L[i]
        term=term.replace('-x','-1*x')
        if term.startswith('x'):
            term=term.replace('x','1*x')
        if term.endswith('x'):
            term=term.replace('x','x^1')
        if '*' not in term:
            term=term.replace('x','*x')
        L[i]=term
        #get degree of polynomial
        d=0
        try:
            if("^" in term):
                d=int(term.split('^')[-1])
        except ValueError as e:
            if len(term)>0:
                raise PolynomialFormatError()
            else:
                continue
        if d>degree:
            degree=d

    coefs=[0 for i in range(degree+1)]
    for i in range(len(L)):
        term=L[i]
        try:
            c=int(term)
            coefs[0]+=c
        except ValueError:
            if len(term)>0:
                d=int(term.split('^')[-1])
                try:
                    c=int(term.split('*')[0])
                except ValueError:
                    raise PolynomialFormatError() 
                coefs[d]+=c
            else:
                raise PolynomialFormatError()

    if coefs[-1] != 1 or degree > max_deg:
        raise PolynomialConstraintError(max_deg)
    
    rstring = ",".join(map(str,coefs))
    return rstring

def parse_output(out):
    """ Parse the verbose output of lll_factor.
    Keys in return: factors (list), time (float)
    """
    answers = out.split("Factorization:")[-1].split("\n")
    answers = [s for s in answers if s]

    time = answers[-1].split(":")[-1].strip()
    factors = answers[:-1]

    rval = {"time":time,"factors":factors}
    return rval

def parse_output_html(out,verbose="on",error=False):
    """ Parse the verbose output of lll_factor to an html string.
    """
    font_family = "Courier New"

    if verbose == "on":
        body = out.replace("\n","<br>")
        body = re.sub(r"={10,}<br>","<hr>",body)
        if "Factorization:" in body:
            body = body.replace("Factorization:","<b>Factorization:")
            body += "</b>"
        body = body.replace("-->","<b>&gt;&gt;&nbsp;")
        body = body.replace("<--","&nbsp; &lt;&lt;</b>")
    else:
        body = "".join(parse_output(out)["factors"])

    error_str = ""
    if error:
        error_str = "<h2> Factorization failed: </h2>"

    rval = f'<html style="font-family: {font_family}; font-size:12">'+ error_str + body + "</html>"
    return rval

if __name__=="__main__":

    argc=len(sys.argv)
    if argc>2:
        print("e e e") #three arguments, meant to throw error in c program this feeds
    elif argc<2:
        print("") #no arguments
    else:
        polystr=sys.argv[1]
        print(parse_poly(polystr))
