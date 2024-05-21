# -*- coding: utf-8 -*-

from __future__ import print_function, division

import sys
import os
from os import path
import numpy as np
import collections
import matplotlib.pyplot as plt

from . import dm_utils

WAV_ERR_SHIFT = 0.002
PSF_VAR_SIGMA = 0.01


def arm_name(arm_num):
    arm_num = np.array(arm_num)
    return np.where(arm_num == 0, 'b',
                    np.where(arm_num == 1, 'r',
                             np.where(arm_num == 2, 'n',
                                      'm')))


def arm_number(armStr):
    return dict(b=0, r=1, n=2, m=3)[armStr]


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
    fiberMag = -np.log(np.trapz(bandpass * counts, wav) /
                       np.trapz(bandpass, wav))
    return fiberMag


def write_ascii(aset, arms, asciiTable, outDir):
    """Write an ascii table"""

    outFile = os.path.join(outDir, asciiTable)
    nFiber = len(aset[0].flux)
    for i in range(nFiber):
        with open('%s.dat' % (outFile) if nFiber == 1 else '%s.%d.dat' % (outFile, i), "w") as fd:
            fd.write('''#  1  WAVELENGTH  [nm]
#  2  FLUX        [nJy]
#  3  ERROR       [nJy]
#  4  MASK        [1=masked]
#  5  SKY         [nJy]
#  6  ARM         [0=blue,1=red,2=NIR,3=redMR]
'''
                     )

            for a, armStr in zip(aset, arms):
                lam = a.wavelength[i]
                flux = a.flux[i]
                sigma = np.sqrt(a.covar[i][0])
                mask = a.mask[i]
                sky = a.sky[i]
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
        sys.exit("Unable to interpret \"%s\" as bool" % val)


def plotPfsObject(pfsObjects, fig=None):
    if fig is None:
        fig = plt.figure(figsize=(7, 4))
    axe = fig.add_subplot()
    axe.set_xlabel('wavelength (nm)')
    axe.set_ylabel('flux (nJy)')
    for pfsObject in pfsObjects:
        axe.plot(pfsObject.wavelength, pfsObject.flux)
    plt.show()


def plotPfsArm(pfsArmSet, fig=None):
    if fig is None:
        fig = plt.figure(figsize=(7, 4))
    axe = fig.add_subplot()
    axe.set_xlabel('wavelength (nm)')
    axe.set_ylabel('flux (electron/pix)')
    for i in range(len(pfsArmSet[0])):
        axe.plot(pfsArmSet[0].wavelength[i],
                 pfsArmSet[0].flux[i], color=f'C{i}')
        axe.plot(pfsArmSet[1].wavelength[i],
                 pfsArmSet[1].flux[i], color=f'C{i}')
        axe.plot(pfsArmSet[2].wavelength[i],
                 pfsArmSet[2].flux[i], color=f'C{i}')
    plt.show()


class Pfsspec(object):

    def __init__(self):
        self.params = {'EXP_NUM': '4',
                       'MAG_FILE': '22.5',
                       'countsMin': '0.1',
                       'etcFile': 'out/ref.snc.dat',
                       'nrealize': '1',
                       'outDir': 'out',
                       'asciiTable': 'None',
                       'ra': 150.0,
                       'dec': 2.0,
                       'tract': 0,
                       'patch': '0,0',
                       'visit0': 1,
                       'catId': 0,
                       'objId': 1,
                       'fiberId': 1,
                       'fiberMag': [22.5, 22.5, 22.5, 22.5, 22.5],
                       'filterName': ['hcs_g', 'hcs_r', 'hcs_i', 'hcs_z', 'hcs_y'],
                       'spectrograph': 1,
                       'pfsConfigFull': 'f',
                       'writeFits': 't',
                       'writePfsArm': 't',
                       'plotArmSet': False,
                       'plotObject': False,
                       'SKY_SUB_FLOOR': '0.01',
                       'SKY_SUB_MODE': 'random',
                       'SKY_SUB_SEED': 0
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
            sys.exit('Error: maybe file not found')
        return 0

    def make_sim_spec(self):
        self.outdir = self.params['outDir']
        self.tract = self.params['tract']
        self.patch = self.params['patch']
        self.visit0 = int(self.params['visit0'])
        self.fiberId = self.params['fiberId']
        self.ra = self.params['ra']
        self.dec = self.params['dec']
        self.catId = self.params['catId']
        self.objId = self.params['objId']
        self.spectrograph = self.params['spectrograph']
        try:
            self.mag_file = '%.4e' % (float(self.params['MAG_FILE']))
        except:
            self.mag_file = self.params['MAG_FILE']
        self.fiberMag = self.params['fiberMag']
        self.filterName = self.params['filterName']

        # do some checks
        self.writeFits = strToBool(self.params['writeFits'])
        self.writePfsArm = strToBool(self.params['writePfsArm'])
        self.plotArmSet = self.params['plotArmSet']
        self.plotObject = self.params['plotObject']
        self.asciiTable = self.params['asciiTable']
        self.pfsConfigFull = strToBool(self.params['pfsConfigFull'])
        self.sky_sub_err = float(self.params['SKY_SUB_FLOOR'])
        self.sky_sub_mode = self.params['SKY_SUB_MODE']
        self.sky_sub_seed = self.params['SKY_SUB_SEED']
        nrealize = int(self.params['nrealize'])
        nexp = int(self.params['EXP_NUM'])

        try:
            if len(self.fiberId) > 0:
                self.multi_info = 1
            else:
                self.multi_info = 0
        except:
            self.multi_info = 0

        if not self.writeFits and not self.asciiTable:
            sys.exit(
                "Please specify asciiTable or omit writeFits (or say writeFits true)")
        if not os.path.exists(self.outdir):
            try:
                os.makedirs(self.outdir)
            except OSError as e:
                sys.exit("Unable to create outDir: %s" % e)
        if nrealize <= 0:
            sys.exit("Please specify at least one realization")

        # check magfile
        if os.path.exists(self.mag_file):
            dat = np.loadtxt(self.mag_file)
            nobj = dat.shape[1] - 1
        else:
            nobj = 1
        if nobj > 1:
            if nrealize > 1:
                sys.exit(
                    "The number of realization should be one for multiple input template")
            else:
                if self.multi_info == 0:
                    objIds = np.arange(self.objId, self.objId + nobj)
                    tmp, catIds = self.catId, np.empty(nobj, dtype=np.int32)
                    catIds.fill(tmp)
                    fiberIds = np.arange(self.fiberId, self.fiberId + nobj)
                    tmp, ras = self.ra, np.empty(nobj, dtype=np.float32)
                    ras.fill(tmp)
                    tmp, decs = self.dec, np.empty(nobj, dtype=np.float32)
                    decs.fill(tmp)
                    tmp, tracts = self.tract, np.empty(nobj, dtype=np.int32)
                    tracts.fill(tmp)
                    tmp, patches = self.patch, np.empty(nobj, dtype='U3')
                    patches.fill(tmp)
                    fiberMags = [self.fiberMag for i in range(nobj)]
                    filterNames = [self.filterName for i in range(nobj)]
                else:
                    objIds = np.array(self.objId)
                    catIds = np.array(self.catId)
                    fiberIds = np.array(self.fiberId)
                    ras = np.array(self.ra)
                    decs = np.array(self.dec)
                    tracts = np.array(self.tract)
                    patches = np.array(self.patch)
                    fiberMags = self.fiberMag
                    filterNames = self.filterName
        else:
            if nrealize > 1:
                try:
                    if len(self.objId) == 1:
                        objIds = np.arange(
                            self.objId[0], self.objId[0] + nrealize)
                        tmp, catIds = self.catId[0], np.empty(
                            nrealize, dtype=np.int32)
                        catIds.fill(tmp)
                        fiberIds = np.arange(
                            self.fiberId[0], self.fiberId[0] + nrealize)
                        tmp, ras = self.ra[0], np.empty(
                            nrealize, dtype=np.float32)
                        ras.fill(tmp)
                        tmp, decs = self.dec[0], np.empty(
                            nrealize, dtype=np.float32)
                        decs.fill(tmp)
                        tmp, tracts = self.tract[0], np.empty(
                            nrealize, dtype=np.int32)
                        tracts.fill(tmp)
                        tmp, patches = self.patch[0], np.empty(
                            nrealize, dtype='U3')
                        patches.fill(tmp)
                        fiberMags = [self.fiberMag[0] for i in range(nrealize)]
                        filterNames = [self.filterName[0]
                                       for i in range(nrealize)]
                except:
                    objIds = np.arange(self.objId, self.objId + nrealize)
                    tmp, catIds = self.catId, np.empty(
                        nrealize, dtype=np.int32)
                    catIds.fill(tmp)
                    fiberIds = np.arange(self.fiberId, self.fiberId + nrealize)
                    tmp, ras = self.ra, np.empty(nrealize, dtype=np.float32)
                    ras.fill(tmp)
                    tmp, decs = self.dec, np.empty(nrealize, dtype=np.float32)
                    decs.fill(tmp)
                    tmp, tracts = self.tract, np.empty(
                        nrealize, dtype=np.int32)
                    tracts.fill(tmp)
                    tmp, patches = self.patch, np.empty(nrealize, dtype='U3')
                    patches.fill(tmp)
                    fiberMags = [self.fiberMag for i in range(nrealize)]
                    filterNames = [self.filterName for i in range(nrealize)]
            else:
                try:
                    if len(self.objId) == 1:
                        objIds = np.array(self.objId)
                        catIds = np.array(self.catId)
                        fiberIds = np.array(self.fiberId)
                        ras = np.array(self.ra)
                        decs = np.array(self.dec)
                        tracts = np.array(self.tract)
                        patches = np.array(self.patch)
                        fiberMags = self.fiberMag
                        filterNames = self.filterName
                except:
                    objIds = np.array([self.objId])
                    catIds = np.array([self.catId])
                    fiberIds = np.array([self.fiberId])
                    ras = np.array([self.ra])
                    decs = np.array([self.dec])
                    tracts = np.array([self.tract])
                    patches = np.array([self.patch])
                    fiberMags = [self.fiberMag]
                    filterNames = [self.filterName]
        '''
            ## read input file ##
            # arm: 0-3 telling us which arm the data comes from (arm_name will convert to b, r, n, m)
            # wav: wavelength
            # nsv: per-pixel (instrumental + sky) variance (counts^2)
            # trn: conversion from I_nu to counts
            # smp: samplingFactor.  A fiddle factor for the Poisson noise in HgCdTe devices
            # skm: sky flux
        '''
        with open(self.params['etcFile'], 'r') as f:
            for line in f.readlines():
                if "EXP_NUM" in line:
                    nexp_etc = int(line.split()[2])
        arm, wav, nsv, trn, smp, skm = np.loadtxt(
            self.params['etcFile'], usecols=(0, 2, 5, 8, 9, 10), unpack=True)

        # remove sky systematics
        skm_sysref = skm.copy()
        skmp = np.roll(skm_sysref, 1)
        skmp[0] = 0.0
        skmm = np.roll(skm_sysref, -1)
        skmm[-1] = 0.0
        skm_sysref = np.amax([skm_sysref, skmm, skmp], axis=0)
        nsv_sys = (self.sky_sub_err * np.sqrt(nexp_etc) * skm_sysref)**2
        nsv_rnd = nsv - nsv_sys
        arm = arm.astype(int)
        trn[trn < 1.0e26] = 1.0e26

        # load magnitude or filename
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
        nsv_rnd_mtrx = np.empty((len(wav), nobj))
        nsv_sys_mtrx = np.empty((len(wav), nobj))
        for i in range(nobj):
            wav_mtrx[:, i] = wav
            trn_mtrx[:, i] = trn
            smp_mtrx[:, i] = smp
            nsv_rnd_mtrx[:, i] = nsv_rnd
            nsv_sys_mtrx[:, i] = nsv_sys * float(nexp) / float(nexp_etc)

        # calculate the flux etc. in observed units
        fnu = 10**(-0.4 * (mag + 48.6))
        fnu_in_njy = fnu / 1e-32
        counts = trn_mtrx * fnu
        if (counts == 0).any():
            print("counts == 0 detected in some pixels; setting to %g for variance" %
                  (float(self.params['countsMin'])), file=sys.stderr)
            # version of counts with zero pixels replaced
            countsp = np.where(counts == 0, float(
                self.params['countsMin']), counts)
        else:
            countsp = counts
        if self.sky_sub_mode == 'random':
            snr1 = countsp / \
                np.sqrt(smp_mtrx * countsp +
                        (nsv_rnd_mtrx + nsv_sys_mtrx)) * np.sqrt(nexp)
            snr2 = countsp / \
                np.sqrt(smp_mtrx * countsp +
                        (nsv_rnd_mtrx + nsv_sys_mtrx)) * np.sqrt(nexp)
        else:
            snr1 = countsp / \
                np.sqrt(smp_mtrx * countsp + nsv_rnd_mtrx) * np.sqrt(nexp)
            snr2 = countsp / \
                np.sqrt(smp_mtrx * countsp +
                        (nsv_rnd_mtrx + nsv_sys_mtrx)) * np.sqrt(nexp)
        sigma1 = fnu_in_njy / snr1
        sigma2 = fnu_in_njy / snr2

        msk = np.zeros_like(wav, dtype=np.int32)
        sky = (skm / trn) / 1e-32
        skm_sysref = sky.copy()
        skmp = np.roll(skm_sysref, 1)
        skmp[0] = 0.0
        skmm = np.roll(skm_sysref, -1)
        skmm[-1] = 0.0
        skm_sysref = np.amax([skm_sysref, skmm, skmp], axis=0)
        arm = arm_name(arm)
        arms = np.array(sorted(set(arm), key=lambda x: dict(
            b=0, r=1, m=1.5, n=2)[x]))  # unique values of arm
        '''
            Create and populate the objects corresponding to the datamodel

            First the parameters describing the observation, in PfsDesign and PfsConfig
        '''
        objectMags = []
        if nobj > 1:
            for i in range(nobj):
                objectMags.append([calculateFiberMagnitude(
                    wav, mag[:, i], b) for b in "grizy"])
        else:
            for i in range(nrealize):
                objectMags.append([calculateFiberMagnitude(
                    wav, mag[:, 0], b) for b in "grizy"])

        pfsDesign = dm_utils.makePfsDesign(tracts, patches, fiberIds, ras,
                                           decs, catIds, objIds, fiberMags, filterNames)

        pfsConfig = dm_utils.makePfsConfig(pfsDesign.pfsDesignId, self.visit0, tracts,
                                           patches, fiberIds, ras, decs, catIds, objIds, fiberMags, filterNames)

        '''
            Create the PfsArm;  we'll put each realisation into a different fibre
            (currently the content is not the one expected in the real pfsArm files so FIXME)
        '''
        metadata = {}
        mapper = {"NO_DATA": 1}
        flags = dm_utils.MaskHelper(**mapper)
        pfsArmSet = []
        for armStr in arms:
            thisArm = (arm == armStr)
            identity = dm_utils.Identity(visit=self.visit0, arm=armStr,
                                         spectrograph=self.spectrograph, pfsDesignId=pfsDesign.pfsDesignId)
            sky_res_fac = 0.0
            # np.random.seed(self.sky_sub_seed)
            nPt = np.sum(thisArm)
            datalam = []
            dataflux = []
            datasky = []
            datanorm = []
            datamask = []
            datacovar = []
            if nobj > 1:
                for i in range(nobj):
                    datalam.append(wav[thisArm])
                    if self.sky_sub_mode == 'systematic':
                        flux = []
                        for j in range(nexp):
                            sky_res_fac = np.random.normal(
                                0.0, self.sky_sub_err)
                            skyres = skm_sysref[thisArm] * sky_res_fac
                            flux.append(fnu_in_njy[thisArm, i] + np.random.normal(0.0,
                                        abs(sigma1[thisArm, i]) * np.sqrt(nexp)) + skyres)
                        dataflux.append(np.nanmean(flux, axis=0))
                    elif self.sky_sub_mode == 'wavecalib':
                        flux = []
                        for j in range(nexp):
                            wav_err_shift = np.random.normal(
                                0.0, WAV_ERR_SHIFT)
                            # skyres = (skm_sysref[thisArm] - np.roll(skm_sysref[thisArm], wav_err_shift))
                            skyres = skm_sysref[thisArm] - \
                                np.interp(wav[thisArm] + wav_err_shift,
                                          wav[thisArm], skm_sysref[thisArm])
                            flux.append(fnu_in_njy[thisArm, i] + np.random.normal(0.0,
                                        abs(sigma1[thisArm, i]) * np.sqrt(nexp)) + skyres)
                        dataflux.append(np.nanmean(flux, axis=0))
                    elif self.sky_sub_mode == 'psfvar':
                        flux = []
                        wav_fine = np.linspace(min(wav[thisArm]), max(
                            wav[thisArm]), len(wav[thisArm]) * 10)
                        skm_sysref_fine = np.interp(
                            wav_fine, wav[thisArm], skm_sysref[thisArm])
                        for j in range(nexp):
                            psf_var_sigma_pix = PSF_VAR_SIGMA / 0.07 * 10
                            gauss_kernel_sigma1 = np.random.uniform(
                                psf_var_sigma_pix * 0.5, psf_var_sigma_pix * 2.0)
                            gauss_kernel_sigma2 = np.random.uniform(
                                psf_var_sigma_pix * 0.5, psf_var_sigma_pix * 2.0)
                            x = np.arange(-100, 101, 1)
                            gauss_kernel1 = (1.0 / np.sqrt(2 * np.pi * gauss_kernel_sigma1**2)
                                             ) * np.exp(-1.0 * x**2 / (2 * gauss_kernel_sigma1**2))
                            gauss_kernel2 = (1.0 / np.sqrt(2 * np.pi * gauss_kernel_sigma2**2)
                                             ) * np.exp(-1.0 * x**2 / (2 * gauss_kernel_sigma2**2))
                            skm_sysref_fine_conv1 = np.convolve(
                                skm_sysref_fine, gauss_kernel1, mode='same')
                            skm_sysref_fine_conv2 = np.convolve(
                                skm_sysref_fine, gauss_kernel2, mode='same')
                            skyres_fine = skm_sysref_fine_conv1 - skm_sysref_fine_conv2
                            skyres = np.interp(
                                wav[thisArm], wav_fine, skyres_fine)
                            flux.append(fnu_in_njy[thisArm, i] + np.random.normal(0.0,
                                        abs(sigma1[thisArm, i]) * np.sqrt(nexp)) + skyres)
                        dataflux.append(np.nanmean(flux, axis=0))
                    else:
                        dataflux.append(fnu_in_njy[thisArm, i] +
                                        np.random.normal(0.0, abs(sigma1[thisArm, i])))
                    datasky.append(sky[thisArm])
                    datamask.append(msk[thisArm])
                    covar = np.zeros(3 * nPt).reshape((3, nPt))
                    covar[0] = sigma2[thisArm, i]**2
                    datacovar.append(covar)
            else:
                for i in range(nrealize):
                    datalam.append(wav[thisArm])
                    if self.sky_sub_mode == 'systematic':
                        flux = []
                        for j in range(nexp):
                            sky_res_fac = np.random.normal(
                                0.0, self.sky_sub_err)
                            skyres = skm_sysref[thisArm] * sky_res_fac
                            flux.append(fnu_in_njy[thisArm, 0] + np.random.normal(0.0,
                                        abs(sigma1[thisArm, 0]) * np.sqrt(nexp)) + skyres)
                        dataflux.append(np.nanmean(flux, axis=0))
                    elif self.sky_sub_mode == 'wavecalib':
                        flux = []
                        wav_fine = np.linspace(min(wav[thisArm]), max(
                            wav[thisArm]), len(wav[thisArm]) * 10)
                        skm_sysref_fine = np.interp(
                            wav_fine, wav[thisArm], skm_sysref[thisArm])
                        for j in range(nexp):
                            wav_err_shift = np.random.normal(
                                0.0, WAV_ERR_SHIFT)
                            # skyres = (skm_sysref[thisArm] - np.roll(skm_sysref[thisArm], wav_err_shift))
                            skyres_fine = skm_sysref_fine - \
                                np.interp(wav_fine + wav_err_shift,
                                          wav_fine, skm_sysref_fine)
                            skyres = np.interp(
                                wav[thisArm], wav_fine, skyres_fine)
                            flux.append(fnu_in_njy[thisArm, 0] + np.random.normal(0.0,
                                        abs(sigma1[thisArm, 0]) * np.sqrt(nexp)) + skyres)
                        dataflux.append(np.nanmean(flux, axis=0))
                    elif self.sky_sub_mode == 'psfvar':
                        flux = []
                        wav_fine = np.linspace(min(wav[thisArm]), max(
                            wav[thisArm]), len(wav[thisArm]) * 10)
                        skm_sysref_fine = np.interp(
                            wav_fine, wav[thisArm], skm_sysref[thisArm])
                        for j in range(nexp):
                            psf_var_sigma_pix = PSF_VAR_SIGMA / 0.07 * 10
                            gauss_kernel_sigma1 = np.random.uniform(
                                psf_var_sigma_pix * 0.3, psf_var_sigma_pix * 3.0)
                            gauss_kernel_sigma2 = np.random.uniform(
                                psf_var_sigma_pix * 0.3, psf_var_sigma_pix * 3.0)
                            x = np.arange(-100, 101, 1)
                            gauss_kernel1 = (1.0 / np.sqrt(2 * np.pi * gauss_kernel_sigma1**2)
                                             ) * np.exp(-1.0 * x**2 / (2 * gauss_kernel_sigma1**2))
                            gauss_kernel2 = (1.0 / np.sqrt(2 * np.pi * gauss_kernel_sigma2**2)
                                             ) * np.exp(-1.0 * x**2 / (2 * gauss_kernel_sigma2**2))
                            skm_sysref_fine_conv1 = np.convolve(
                                skm_sysref_fine, gauss_kernel1, mode='same')
                            skm_sysref_fine_conv2 = np.convolve(
                                skm_sysref_fine, gauss_kernel2, mode='same')
                            skyres_fine = skm_sysref_fine_conv1 - skm_sysref_fine_conv2
                            skyres = np.interp(
                                wav[thisArm], wav_fine, skyres_fine)
                            flux.append(fnu_in_njy[thisArm, 0] + np.random.normal(0.0,
                                        abs(sigma1[thisArm, 0]) * np.sqrt(nexp)) + skyres)
                        dataflux.append(np.nanmean(flux, axis=0))
                    else:
                        dataflux.append(fnu_in_njy[thisArm, 0] +
                                        np.random.normal(0.0, abs(sigma1[thisArm, 0])))
                    datasky.append(sky[thisArm])
                    datanorm.append([1.0 for _ in sky[thisArm]])
                    datamask.append(msk[thisArm])
                    covar = np.zeros(3 * nPt).reshape((3, nPt))
                    covar[0] = sigma2[thisArm, 0]**2
                    datacovar.append(covar)
            pfsArm = dm_utils.PfsArm(identity=identity,
                                     fiberId=fiberIds,
                                     wavelength=np.array(datalam, dtype='f4'),
                                     flux=np.array(dataflux, dtype='f4'),
                                     mask=np.array(datamask, dtype='i4'),
                                     sky=np.array(datasky, dtype='f4'),
                                     norm=np.array(datanorm, dtype='f4'),
                                     covar=np.array(datacovar, dtype='f4'),
                                     flags=flags,
                                     metadata=metadata
                                     )
            pfsArmSet.append(pfsArm)
        if self.plotArmSet:
            plotPfsArm(pfsArmSet)

        '''
            Time for I/O
        '''
        # writ to Fits
        if self.writeFits:
            print(f'pfsDesignId: {pfsDesign.pfsDesignId:#013x}')
            pfsDesign.write(dirName=self.outdir)  # pfsDesign file
            pfsConfig.write(dirName=self.outdir)  # pfsConfig file
            if self.writePfsArm:                  # write pfsArm files
                for pfsArm in pfsArmSet:
                    pfsArm.write(self.outdir)
        # write to Ascii
        if self.asciiTable != "None":
            write_ascii(pfsArmSet, arms, self.asciiTable, self.outdir)
            # print("ASCII table %s was generated" % self.asciiTable)

        # Now make the PfsObject from the PfsArmSet

        visits = self.visit0 + np.arange(nexp)
        # identityList = [dm_utils.Identity(visit=v, arm=arm[0], spectrograph=self.spectrograph, pfsDesignId=pfsDesign.pfsDesignId).getDict() for v in visits]
        pfsObjects, pfsVisitHashes = dm_utils.makePfsObject(
            pfsConfig=pfsConfig, pfsArmSet=pfsArmSet, visits=visits, minWavelength=350., maxWavelength=1260., dWavelength=0.08)
        for pfsObject in pfsObjects:
            if self.writeFits:                   # write FITS file
                pfsObject.write(self.outdir)
        if self.plotObject:                  # plot pfsObject
            plotPfsObject(pfsObjects)
        self.pfsObjects = pfsObjects
        self.pfsVisitHashes = pfsVisitHashes
        return 0

    def proc_multi(self, inputs):
        for k, v in inputs.items():
            print(k, v)
            self.params[k] = v
        self.make_sim_spec()
        return 0

    def make_sim_spec_multi(self, nproc, params):
        from multiprocessing import Pool
        p = Pool(nproc)
        result = p.map(self.proc_multi, params)
        return 0
