def degree(poly):
	d=len(poly);
	i=0;
	while poly[i]==0 and i<d-1:
		i+=1;
		d-=1;
	if(max(poly,key=abs)==0):
		return -1
	return d-1;

def extended_synthetic_division(dividend, divisor):
    out = list(dividend) # Copy the dividend
    offset=len(divisor)-degree(divisor);
    first_nz=0;
    for i in range(len(divisor)):
        if divisor[i] !=0 :
            first_nz=i
            break
    lc=divisor[first_nz]
    for i in xrange(len(dividend)-(degree(divisor))):
        out[i] /= 1.0*lc
        coef = out[i];
        if coef != 0: # useless to multiply if coef is 0
            for j in xrange(offset, len(divisor)): 
                out[i + j -offset+1] += -divisor[j] * coef
    separator = -(degree(divisor))
    if(degree(divisor)==0):
        return out,[0] #edge case
    else:
        return out[:separator], out[separator:] # return quotient, remainder.


def gcd(p1,p2):
	a=p1
	b=p2
	c=0
	while b != [0] and c<10:
		c+=1
		r=extended_synthetic_division(a,b)[1]
		a=b
		b=r
		print [a,b]

		if degree(r)==-1:
			return b
	return -1
 
#p=[1,3,-6,-9,9,0,0,-1,-3,3];
#d=[0,0,0,0,0,0,0,0,0,1,3,-3];
#p=[5,0,0,0,0];
#d=[0,0,0,0,-1];
#print(extended_synthetic_division(p,d))

p1=[1,0,1,0,0,0,1,0,1]
p2=[0,8,0,6,0,0,0,2,0]

gcd(p1,p2)

