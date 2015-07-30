from scipy.optimize import leastsq
from numpy import *

def plane(b,c,d,x0,y0,x,y):
    return b*(x-x0) + c*(y-y0) + d

def fillplane(b,c,d,x0,y0,arr):
    output = zeros(arr.shape)
    for i in range(len(output)):
        for j in range(len(output[0])):
            output[i,j] = plane(b,c,d,x0,y0,j,i)
    return output

def dflyi(params,x0,x,y0,y,Hi):
    a,b,c,d = params
    return a*Hi + plane(b,c,d,x0,y0,x,y)

def residuals(params,Di,x0,x,y0,y,Hi):
    return abs(Di - dflyi(params,x0,x,y0,y,Hi))

def subBGplane(Di,Hi,p0):
    y0 = Di.shape[0]/2
    x0 = Di.shape[1]/2
    inds = where((Di!=0) & (isnan(Hi) == False))
    y = inds[0]
    x = inds[1]
    ps = leastsq(residuals,p0,args = (Di[inds],x0,x,y0,y,Hi[inds]),
                 full_output=1)
    success = ps[-1]
    allowedsuccess = [1,2,3,4]
    if success in allowedsuccess:
        print 'Successful fit'
        a,b,c,d = ps[0]
        Bgplane = fillplane(b,c,d,x0,y0,Di)
        return Di - Bgplane, Bgplane
    if not success in allowedsuccess:
        print 'Failed to fit'
        return 0,0
