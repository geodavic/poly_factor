import sys

def parse_poly(polystr):
    """
    parse polynomial into csv format
    e.g. x^4-1  ->  -1,0,0,1
    """
    polystr = polystr.replace(" ","")
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
        except ValueError:
            if len(term)>0:
                return "error"
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
                    return "error"
                coefs[d]+=c
    rstring=''
    for c in coefs:
        rstring+='%d,'%c

    return rstring

def parse_opts(opts_list,allowed_opts):
    """ Parse options from api request
    """
    if opts_list is None:
        return []


    for i in range(len(opts_list)):
        opt = opts_list[i]
        if opt[0] != "-":
            opts_list[i] = "-" + opt
        o1 = opt.split(" ")[0].replace("-","")
        assert o1 in allowed_opts,f"Unrecognized option: {opt}"

    
    rval = " ".join(opts_list)
    rval = rval.split(" ")
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
