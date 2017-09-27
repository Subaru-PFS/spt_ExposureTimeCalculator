#
# A couple of useful functions to calculate PFS's SHA-1s
#
import hashlib
import numpy as np

def calculate_pfsVisitHash(visits):
    """Calculate and return the 64-bit SHA-1 from a list of visits"""

    if len(visits) == 1 and visits[0] == 0:
        return 0x0
    
    m = hashlib.sha1()
        
    for visit in visits:
        m.update("%d" % (visit))

    # convert to int and truncate to 8 bytes
    sha1 = int(m.hexdigest(), 16) & 0xffffffff
    
    return sha1

def calculate_pfsConfigId(fiberIds, ras, decs):
    """Calculate and return the 64-bit SHA-1 from a set of lists of
    fiberId, ra, and dec"""

    if fiberIds is None:
        if ras is None and decs is None:
            return 0x0

        raise RuntimeError("Either all or none of fiberId, ra, and dec may be None")

    if (fiberIds == 1 + np.arange(len(fiberIds))).all() and (ras == 0.0).all() and (decs == 0.0).all():
        return 0x0    
    
    m = hashlib.sha1()

    for fiberId, ra, dec in zip(fiberIds, ras, decs):
        m.update("%d %.0f %.0f" % (fiberId, ra, dec))

    # convert to int and truncate to 8 bytes
    sha1 = int(m.hexdigest(), 16) & 0xffffffff

    return sha1

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def makeFullCovariance(covar):
    """Given a matrix of the diagonal part of the covariance matrix return a full matrix

    This is mostly useful for visualising the COVAR arrays in pfsArm/pfsObject files
    
    Specifically,
       covar[0, 0:]    Diagonal
       covar[1, 0:-1]  +-1 off diagonal
       covar[2, 0:-2]  +-2 off diagonal
    """
    nband = covar.shape[0]
    C = np.zeros(covar.shape[1]**2).reshape((covar.shape[1], covar.shape[1]))
    
    i = np.arange(C.shape[0])
    
    C[i, i] = covar[0]                     # diagonal    
    for j in range(1, covar.shape[0]):     # bands near the diagonal
        C[i[j:], i[:-j]] = covar[j][0:-j]         # above the diagonal
        C[i[:-j], i[j:]] = covar[j][0:-j]         # below the diagonal
        
    return C
