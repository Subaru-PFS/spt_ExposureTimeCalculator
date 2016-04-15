#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np

for i in range(2):
    try:
        from pfs.datamodel.pfsConfig import PfsConfig
    except ImportError as e:
        if i == 0:
            sys.path.append("datamodel/python")

from pfs.datamodel.pfsArm import PfsArmSet
from pfs.datamodel.pfsObject import PfsObject, makePfsObject

#######################
offset = 0.01
#######################

def arm_name(arm_num):
    arm_num = np.array(arm_num)
    return np.where(arm_num == 0, 'b',
                    np.where(arm_num == 1, 'r', 
                    np.where(arm_num == 2, 'n',
                    'm')))

def arm_number(armStr):
    return dict(b=0, r=1, n=2, m=3)[armStr]

def makeFakePfsConfig(tract, patch, ra, dec, catId, objId, objectMags, nFiber=1):
    """Make and return a PfsConfig with nFiber entries referring to the same object

    Successive objIds are applied to the successive fibers
    """
    fiberId = np.arange(1, nFiber+1, dtype=np.int32)

    tmp, ra = ra, np.empty(nFiber)
    ra.fill(tmp)

    tmp, dec = dec, np.empty(nFiber);
    dec.fill(tmp)

    catIds = np.empty(nFiber)
    catIds.fill(catId)

    tmp, tract = tract, np.empty(nFiber, dtype=np.int32);
    tract.fill(tmp)

    tmp, patch = patch, np.empty(nFiber, dtype=str);
    patch.fill(tmp)
    patch = nFiber*[tmp]

    objIds = objId + np.arange(nFiber, dtype=np.int64)

    fiberMag = np.empty((nFiber, 5))
    for i in range(nFiber):
        fiberMag[i] = objectMags

    mpsCen = np.zeros((nFiber, 2))

    return PfsConfig(None, tract, patch, fiberId, ra, dec, catId=catIds,
                     objId=objIds, mpsCen=mpsCen, fiberMag=fiberMag)

def calculateFiberMagnitude(wav, mag, filterName):
    """Calculate the average magnitude over the bandpass"""
    #                           50%/nm  50%/nm peak
    filterBandpasses = dict(g=(399.5,  546.5, 0.97),
                            r=(542.5,  696.5, 0.95),
                            i=(698.5,  853.3, 0.90), # i2
                            z=(852.5,  932.0, 0.97),
                            y=(943.0, 1072.0, 0.95)
                            )
    wav0, wav1, peak = filterBandpasses[filterName]

    counts = np.exp(-mag)
    bandpass = np.where(np.logical_and(wav >= wav0, wav <= wav1), peak, 0)
    fiberMag = -np.log(np.trapz(bandpass*counts, wav)/np.trapz(bandpass, wav))
    return fiberMag

def write_ascii(aset, asciiTable, outDir):
    """Write an ascii table"""

    outFile = os.path.join(outDir, asciiTable)
    nFiber = len(aset.data.values()[0].flux)

    for i in range(nFiber):
        with open(outFile if nFiber == 1 else '%s.%d' % (outFile, i), "w") as fd:
            fd.write('''#  1  WAVELENGTH  [nm]
#  2  FLUX        [10^-17 erg/s/cm^2/A]
#  3  ERROR       [10^-17 erg/s/cm^2/A]
#  4  MASK        [1=masked]
#  5  SKY         [10^-17 erg/s/cm^2/A]
#  6  ARM         [0=blue,1=red,2=NIR,3=redMR]
'''
            )

            for armStr, arm in aset.data.items():
                lam = arm.lam[i]
                flux = arm.flux[i]
                sigma = np.sqrt(arm.covar[i][0])
                mask = arm.mask[i]
                sky = arm.sky[i]
                armNum = arm_number(armStr)

                for j in range(len(lam)):
                    fd.write('%8.3f %12.4e %12.4e %2d %12.4e %1d\n' % \
                             (lam[j], flux[j], sigma[j], mask[j], sky[j], armNum))

def main(inFile,
         magFile,
         outDir=".",
         tract=0,
         patch='0,0',
         visit=1,
         spectrograph=1,
         catId=0,
         objId=1,
         nRealize=1,
         nExposure=1,
         countsMin=0.1,
         asciiTable=None,
         writeFits=True, writePfsArm=True,
         plotArmSet=False, plotObject=False,
         ):
    ## read input file ##
    # arm: 0-3 telling us which arm the data comes from (arm_name will convert to b, r, n, m)
    # wav: wavelength
    # nsv: per-pixel (instrumental + sky) variance (counts^2)
    # trn: conversion from I_nu to counts
    # smp: samplingFactor.  A fiddle factor for the Poisson noise in HgCdTe devices
    # skm: sky flux
    arm, wav, nsv, trn, smp, skm = np.loadtxt(inFile, usecols=(0, 2, 5, 8, 9, 10), unpack=True)
    arm = arm.astype(int)
    trn[trn < 1.0e26] = 1.0e26

    ## load magnitude or filename ##
    if os.path.exists(magFile):
        _lam, _mag = np.loadtxt(magFile, usecols=(0, 1), unpack=True)
        mag = np.interp(wav, _lam, _mag)
    else:
        mag = np.empty(len(wav))
        mag.fill(float(magFile))

    ##
    ## Calculate the flux etc. in observed units
    ##
    fnu   = 10**(-0.4*(mag + 48.6))
    flam  = 3.0e18*fnu/(10*wav)**2/1e-17

    counts = trn*fnu
    if (counts == 0).any():
        print >> sys.stderr, "counts == 0 detected in some pixels; setting to %g for variance" % countsMin
        countsp = np.where(counts == 0, countsMin, counts) # version of counts with zero pixels replaced
    else:
        countsp = counts
        
    snr = countsp/np.sqrt(smp*countsp + nsv)*np.sqrt(nExposure)
    sigma = flam/snr
    msk = np.zeros_like(wav, dtype=np.int32)
    sky = 3.0e18*(skm/trn)/(10*wav)**2/1e-17

    arm = arm_name(arm)
    arms = np.array(sorted(set(arm), key=lambda x: dict(b=0, r=1, m=1.5, n=2)[x])) # unique values of arm
    #
    # Create and populate the objects corresponding to the datamodel
    #
    # First the parameters describing the observation, in PfsConfig
    objectMags = [calculateFiberMagnitude(wav, mag, b) for b in "grizy"]
    pfsConfig = makeFakePfsConfig(tract, patch,
                                  150.0, 2.0, catId, objId, objectMags, nFiber=nRealize)
    #
    # Create the PfsArmSet;  we'll put each realisation into a different fibre
    #
    pfsConfigId = pfsConfig.pfsConfigId
    pfsArmSet = PfsArmSet(visit, spectrograph, arms=arms, pfsConfig=pfsConfig)

    for armStr, data in pfsArmSet.data.items():
        thisArm = (arm == armStr)
        nPt = np.sum(thisArm)

        for i in range(nRealize):
            data.lam.append(wav[thisArm])
            data.flux.append(flam[thisArm] + np.random.normal(0.0, sigma[thisArm]))
            data.sky.append(sky[thisArm])
            data.mask.append(msk[thisArm])
            covar = np.zeros(3*nPt).reshape((3, nPt))
            covar[0] = sigma[thisArm]**2
            data.covar.append(covar)

    if plotArmSet:
        pfsArmSet.plot(showFlux=True, showMask=False, showSky=False, showCovar=False)
    #
    # Now make the PfsObject from the PfsArmSet
    #
    pfsObject = makePfsObject(objId, [pfsArmSet])

    if plotObject:
        pfsObject.plot()
    #
    # Time for I/O
    #
    # Ascii
    #
    if asciiTable:
        write_ascii(pfsArmSet, asciiTable, outDir)
        print "ASCII table %s was generated" % asciiTable
    #
    # Fits
    #
    if writeFits:
        pfsConfig.write(outDir)         # pfsConfig file

        if writePfsArm:                 # write pfsArm files
            for data in pfsArmSet.data.values():
                data.write(outDir)      

        pfsObject.write(outDir)         # pfsObject file
    
def strToBool(val):
    if val.lower()   in ("1", "t", "true"):
        return True
    elif val.lower() in ("0", "f", "false",):
        args.writeFits = False
    else:
        exit("Unable to interpret \"%s\" as bool" % val)

if __name__ == '__main__':
    def convert_arg_line_to_args(arg_line):
        """Make argparse handle Hirata-style parameter files"""
        for i, arg in enumerate(arg_line.split()):
            if not arg.strip():
                continue
            if arg[0] == '#':
                return
            
            if i == 0:
                arg = "--%s" % arg
            yield arg
    
    parser = argparse.ArgumentParser(description='PFS Spectral Simulator developed by Kiyoto Yabe',
                                     fromfile_prefix_chars='@')
    parser.convert_arg_line_to_args = convert_arg_line_to_args
    
    parser.add_argument("--EXP_NUM", type=int, help="Number of exposures", default=1)
    parser.add_argument("--MAG_FILE", type=str, help="magnitude input file", default="22.5")
    parser.add_argument("--etcFile", type=str, help="continuum results input file", default="etc-lr.dat")
    parser.add_argument("--nrealize", type=int, help="the number of realization", default=1)
    parser.add_argument("--outDir", type=str, help="Directory for outputs", default=".")
    parser.add_argument("--tract", type=int, help="tract", default=0)
    parser.add_argument("--patch", type=str, help="patch", default='0,0')
    parser.add_argument("--visit", type=int, help="visit number", default=1)
    parser.add_argument("--objId", type=int, help="object id", default=1)
    parser.add_argument("--catId", type=int, help="catalogue id", default=0)
    parser.add_argument("--spectrograph", type=int, help="spectrograph number", default=1)
    parser.add_argument("--countsMin", type=float, help="Minimum counts per pixel for noise", default=0.1)
    parser.add_argument("--asciiTable", type=str, help="simulated spectrum output ASCII file", default=None)
    parser.add_argument("--writeFits", type=str, help="Write FITS files", default="True")
    parser.add_argument("--writePfsArm", type=str, help="Write pfsArm files (writeFits must be set)",
                        default="True")
    parser.add_argument("--plotArmSet", action='store_true', help="Plot the pfsArmSet data", default=False)
    parser.add_argument("--plotObject", action='store_true', help="Plot the pfsObject data", default=False)

    args = parser.parse_args()

    if args.asciiTable == 'None':
        args.asciiTable = None

    args.writeFits = strToBool(args.writeFits)
        
    if not args.writeFits and not args.asciiTable:
        exit("Please specify --asciiTable or omit --writeFits (or say --writeFits true)")

    if not os.path.exists(args.outDir):
        try:
            os.makedirs(args.outDir)
        except OSError as e:
            exit("Unable to create args.OutDir: %s" % e)

    main(inFile=args.etcFile,
         magFile=args.MAG_FILE,
         outDir=args.outDir,
         tract=args.tract,
         patch=args.patch,
         visit=args.visit,
         spectrograph=args.spectrograph,
         catId=args.catId,
         objId=args.objId,
         nRealize=args.nrealize,
         nExposure=args.EXP_NUM,
         countsMin=args.countsMin,
         asciiTable=args.asciiTable,
         writeFits=args.writeFits, writePfsArm=args.writePfsArm,
         plotArmSet=args.plotArmSet, plotObject=args.plotObject,
         )
