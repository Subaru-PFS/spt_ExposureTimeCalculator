# -*- coding: utf-8 -*-

from __future__ import print_function, division

import sys
import os
from os import path
import numpy as np
import scipy as sp


HOME_DIR = path.dirname(path.abspath(__file__))
''' import datamodel module '''
sys.path.append(HOME_DIR + "/datamodel/python")
from pfs.datamodel.pfsConfig import PfsConfig
from pfs.datamodel.pfsArm import PfsArmSet
from pfs.datamodel.pfsObject import PfsObject, makePfsObject


def arm_name(arm_num):
    arm_num = np.array(arm_num)
    return np.where(arm_num == 0, 'b',
                    np.where(arm_num == 1, 'r',
                             np.where(arm_num == 2, 'n',
                                      'm')))


def arm_number(armStr):
    return dict(b=0, r=1, n=2, m=3)[armStr]


def makeFakePfsConfig(tract, patch, ra, dec, catId, startingObjId, objectMags, nFiber=1):
    """Make and return a PfsConfig with nFiber entries referring to the same object

    Successive fibres are given object IDs that increase by one, starting with objId
    """
    fiberId = np.arange(1, nFiber + 1, dtype=np.int32)

    tmp, ra = ra, np.empty(nFiber)
    ra.fill(tmp)

    tmp, dec = dec, np.empty(nFiber)
    dec.fill(tmp)

    catIds = np.empty(nFiber)
    catIds.fill(catId)

    tmp, tract = tract, np.empty(nFiber, dtype=np.int32)
    tract.fill(tmp)

    tmp, patch = patch, np.empty(nFiber, dtype=str)
    patch.fill(tmp)
    patch = nFiber * [tmp]

    objIds = startingObjId + np.arange(nFiber, dtype=np.int64)

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
                            i=(698.5,  853.3, 0.90),  # i2
                            z=(852.5,  932.0, 0.97),
                            y=(943.0, 1072.0, 0.95)
                            )
    wav0, wav1, peak = filterBandpasses[filterName]

    counts = np.exp(-mag)
    bandpass = np.where(np.logical_and(wav >= wav0, wav <= wav1), peak, 0)
    fiberMag = -np.log(np.trapz(bandpass * counts, wav) / np.trapz(bandpass, wav))
    return fiberMag


def write_ascii(aset, asciiTable, outDir):
    """Write an ascii table"""

    outFile = os.path.join(outDir, asciiTable)
    nFiber = len(aset.data.values()[0].flux)

    for i in range(nFiber):
        with open('%s.dat' % (outFile) if nFiber == 1 else '%s.%d.dat' % (outFile, i), "w") as fd:
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
                    fd.write('%8.3f %12.4e %12.4e %2d %12.4e %1d\n' %
                             (lam[j], flux[j], sigma[j], mask[j], sky[j], armNum))


def strToBool(val):
    if val.lower() in ("1", "t", "true"):
        return True
    elif val.lower() in ("0", "f", "false"):
        return False
    else:
        exit("Unable to interpret \"%s\" as bool" % val)


class Pfsspec(object):

    def __init__(self):
        self.params = {'EXP_NUM': '8',
                       'MAG_FILE': '22.5',
                       'countsMin': '0.1',
                       'etcFile': 'out/ref.snc.dat',
                       'nrealize': '1',
                       'outDir': 'out',
                       'asciiTable': 'None',
                       'tract': '0',
                       'patch': '0,0',
                       'visit': '1',
                       'catId': '0',
                       'objId': '1',
                       'spectrograph': '1',
                       'writeFits': 't',
                       'writePfsArm': 't',
                       'plotArmSet': 'f',
                       'plotObject': 'f'
                       }
        return None

    def set_param(self, param_name, param_value):
        if param_name in self.params.keys():
            try:
                self.params[param_name] = param_value
            except:
                print('Error!')
        else:
            print('param_name %s can not be recognized ...' % (param_name))
        return 0

    def load_param_file(self, filename):
        try:
            for line in open(filename, 'r'):
                a = line.split()
                if line[0] != '#' and len(a) > 0:
                    self.params[a[0]] = a[1]
        except:
            exit('Error: maybe file not found')
        return 0

    def make_sim_spec(self):
        self.outdir = self.params['outDir']
        self.tract = int(self.params['tract'])
        self.patch = self.params['patch']
        self.visit = int(self.params['visit'])
        self.catId = int(self.params['catId'])
        self.objId = int(self.params['objId'])
        self.spectrograph = int(self.params['spectrograph'])
        try:
            self.mag_file = '%.4e' % (float(self.params['MAG_FILE']))
        except:
            self.mag_file = self.params['MAG_FILE']
        ''' some checks '''
        self.writeFits = strToBool(self.params['writeFits'])
        self.writePfsArm = strToBool(self.params['writePfsArm'])
        self.plotArmSet = strToBool(self.params['plotArmSet'])
        self.plotObject = strToBool(self.params['plotObject'])
        self.asciiTable = self.params['asciiTable']

        if not self.writeFits and not self.asciiTable:
            exit("Please specify asciiTable or omit writeFits (or say writeFits true)")
        if not os.path.exists(self.outdir):
            try:
                os.makedirs(self.outdir)
            except OSError as e:
                exit("Unable to create outDir: %s" % e)
        if int(self.params['nrealize']) <= 0:
            exit("Please specify at least one realization")
        ''' check magfile '''
        if os.path.exists(self.mag_file):
            dat = np.loadtxt(self.mag_file)
            nobj = dat.shape[1] - 1
        else:
            nobj = 1
        if nobj > 1:
            if int(self.params['nrealize']) > 1:
                exit("The number of realization should be one for multiple input template")
            else:
                objIds = range(self.objId, self.objId + nobj)
        else:
            objIds = range(self.objId, self.objId + int(self.params['nrealize']))
        '''
            ## read input file ##
            # arm: 0-3 telling us which arm the data comes from (arm_name will convert to b, r, n, m)
            # wav: wavelength
            # nsv: per-pixel (instrumental + sky) variance (counts^2)
            # trn: conversion from I_nu to counts
            # smp: samplingFactor.  A fiddle factor for the Poisson noise in HgCdTe devices
            # skm: sky flux
        '''
        arm, wav, nsv, trn, smp, skm = np.loadtxt(self.params['etcFile'], usecols=(0, 2, 5, 8, 9, 10), unpack=True)
        arm = arm.astype(int)
        trn[trn < 1.0e26] = 1.0e26
        ''' load magnitude or filename '''
        if os.path.exists(self.mag_file):
            dat = np.loadtxt(self.mag_file)
            nobj = dat.shape[1] - 1
            _lam = dat[:, 0]
            mag = np.empty((len(wav), nobj))
            for i in range(nobj):
                _mag = dat[:, i + 1]
                mag[:, i] = np.interp(wav, _lam, _mag)
        else:
            nobj = 1
            mag = np.empty((len(wav), nobj))
            mag[:, 0].fill(self.mag_file)
        wav_mtrx = np.empty((len(wav), nobj))
        trn_mtrx = np.empty((len(wav), nobj))
        smp_mtrx = np.empty((len(wav), nobj))
        nsv_mtrx = np.empty((len(wav), nobj))
        for i in range(nobj):
            wav_mtrx[:, i] = wav
            trn_mtrx[:, i] = trn
            smp_mtrx[:, i] = smp
            nsv_mtrx[:, i] = nsv
        ''' calculate the flux etc. in observed units '''
        fnu = 10**(-0.4 * (mag + 48.6))
        flam = 3.0e18 * fnu / (10 * wav_mtrx)**2 / 1e-17
        counts = trn_mtrx * fnu
        if (counts == 0).any():
            print("counts == 0 detected in some pixels; setting to %g for variance" % (float(self.params['countsMin'])), file=sys.stderr)
            countsp = np.where(counts == 0, float(self.params['countsMin']), counts)  # version of counts with zero pixels replaced
        else:
            countsp = counts
        snr = countsp / np.sqrt(smp_mtrx * countsp + nsv_mtrx) * np.sqrt(int(self.params['EXP_NUM']))
        sigma = flam / snr
        msk = np.zeros_like(wav, dtype=np.int32)
        sky = 3.0e18 * (skm / trn) / (10 * wav)**2 / 1e-17
        arm = arm_name(arm)
        arms = np.array(sorted(set(arm), key=lambda x: dict(b=0, r=1, m=1.5, n=2)[x]))  # unique values of arm
        '''
            Create and populate the objects corresponding to the datamodel
    
            First the parameters describing the observation, in PfsConfig
        '''
        objectMags = [calculateFiberMagnitude(wav, mag[:, 0], b) for b in "grizy"]
        if nobj > 1:
            pfsConfig = makeFakePfsConfig(self.tract, self.patch, 150.0, 2.0, self.catId, objIds[0], objectMags, nFiber=nobj)
        else:
            pfsConfig = makeFakePfsConfig(self.tract, self.patch, 150.0, 2.0, self.catId, objIds[0], objectMags, nFiber=int(self.params['nrealize']))
        '''
            Create the PfsArmSet;  we'll put each realisation into a different fibre
        '''
        pfsConfigId = pfsConfig.pfsConfigId
        pfsArmSet = PfsArmSet(self.visit, self.spectrograph, arms=arms, pfsConfig=pfsConfig)

        for armStr, data in pfsArmSet.data.items():
            thisArm = (arm == armStr)
            nPt = np.sum(thisArm)
            if nobj > 1:
                for i in range(nobj):
                    data.lam.append(wav[thisArm])
                    data.flux.append(flam[thisArm, i] + np.random.normal(0.0, sigma[thisArm, i]))
                    data.sky.append(sky[thisArm])
                    data.mask.append(msk[thisArm])
                    covar = np.zeros(3 * nPt).reshape((3, nPt))
                    covar[0] = sigma[thisArm, i]**2
                    data.covar.append(covar)
            else:
                for i in range(int(self.params['nrealize'])):
                    data.lam.append(wav[thisArm])
                    data.flux.append(flam[thisArm, 0] + np.random.normal(0.0, sigma[thisArm, 0]))
                    data.sky.append(sky[thisArm])
                    data.mask.append(msk[thisArm])
                    covar = np.zeros(3 * nPt).reshape((3, nPt))
                    covar[0] = sigma[thisArm, 0]**2
                    data.covar.append(covar)
        if self.plotArmSet:
            pfsArmSet.plot(showFlux=True, showMask=False, showSky=False, showCovar=False)
        '''
            Time for I/O
        '''
        ''' Fits '''
        if self.writeFits:
            pfsConfig.write(self.outdir)         # pfsConfig file
            if self.writePfsArm:                 # write pfsArm files
                for data in pfsArmSet.data.values():
                    data.write(self.outdir)
        ''' Ascii '''
        if self.asciiTable != "None":
            write_ascii(pfsArmSet, self.asciiTable, self.outdir)
            print("ASCII table %s was generated" % self.asciiTable)
        '''
            Now make the PfsObject from the PfsArmSet
        '''
        for objId in objIds:
            pfsObject = makePfsObject(objId, [pfsArmSet])

            if self.plotObject:
                pfsObject.plot()

            if self.writeFits:
                pfsObject.write(self.outdir)         # pfsObject file

        return 0
