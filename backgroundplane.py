"""
backgroundplane - contains functions to subtract a background plane
    from an image

Requires the following modules: scipy, numpy, docopt

Contains the following functions: plane, fillplane, dflyi, residuals
                                  subBGplane

"""

########################## IMPORT PACKAGES ###########################

from scipy.optimize import leastsq
from numpy import *

########################## FUNCTIONS ###########################

def plane(b,c,d,x0,y0,x,y):
    """
    Basic equation for a plane

    b,c,d:      constant parameters of plane
    x0,y0:      coordinates of centre of the plane
    x,y:        x,y arrays, must have the same shape

    Returns an array with the shape of x (or equivalently y)

    """
    return b*(x-x0) + c*(y-y0) + d

def fillplane(b,c,d,x0,y0,arr):
    """
    Fill a 2D array with appropriate plane values.
        Assumes that the x,y coordinates in the plane 
        correspond to the indices of the array

    b,c,d:      constant parameters of plane
    x0,y0:      coordinates of centre of the plane 
    arr:        array to fill with plane values

    Returns a 2D array with shape of arr   

    """
    output = zeros(arr.shape)
    for i in range(len(output)):
        for j in range(len(output[0])):
            output[i,j] = plane(b,c,d,x0,y0,j,i)
    return output

def dflyi(params,x0,x,y0,y,Hi):
    """
    A theoretical dragonfly image constructed from a Herschel image
        and a background plane

    params:     a set parameters including a multiplicative factor
                for the Herschel image, and three constants for 
                the plane
    x0,y0:      coordinates of centre of the plane
    x,y:        x,y arrays, must have the same shape
    Hi:         Herschel image with the same shape as x and y

    Returns an array in shape of Hi

    """
    a,b,c,d = params
    return a*Hi + plane(b,c,d,x0,y0,x,y)

def residuals(params,Di,x0,x,y0,y,Hi):
    """
    Find the absolute difference between a Dragonfly image and the
        theoretical image produced by dflyi()
    
    params:     a set parameters including a multiplicative factor
                for the Herschel image, and three constants for 
                the plane
    x0,y0:      coordinates of centre of the plane
    x,y:        x,y arrays, must have the same shape
    Hi:         Herschel image with the same shape as x and y
    Di:         Dragonfly image with the same shape as x and y


    """
    return abs(Di - dflyi(params,x0,x,y0,y,Hi))

def subBGplane(Di,Hi,p0):
    """
    Fit a background plane to a Dragonfly image

    Di:     Dragonfly image
    Hi:     Herschel image
    p0:     preliminary guess for fit parameters

    Returns background subtracted dragonfly image
    """
    # Find central pixel coordinates
    y0 = Di.shape[0]/2
    x0 = Di.shape[1]/2
    # Remove masked pixels and NAN pixels
    inds = where((Di!=0) & (isnan(Hi) == False))
    y = inds[0]
    x = inds[1]
    # Fit background plane
    ps = leastsq(residuals,p0,args = (Di[inds],x0,x,y0,y,Hi[inds]),
                 full_output=1)
    # Confirm that fit is successful and return results accordingly
    success = ps[-1]
    allowedsuccess = [1,2,3,4]
    if success in allowedsuccess:
        print 'Successful fit'
        a,b,c,d = ps[0]
        Bgplane = fillplane(b,c,d,x0,y0,Di)
        return Di - Bgplane, Bgplane, [b,c,d,x0,y0]
    if not success in allowedsuccess:
        print 'Failed to fit'
        return 0,0,[]
