# -*- coding: utf-8 -*-

from __future__ import print_function, division

import sys
import os
from os import path
import numpy as np
import collections

HOME_DIR = path.dirname(path.abspath(__file__))
''' import datamodel module '''
sys.path.append(HOME_DIR + "/datamodel/python")
from pfs.datamodel.pfsConfig import PfsDesign, PfsConfig
from pfs.datamodel.drp import PfsArm, PfsMerged, PfsObject
from pfs.datamodel.masks import MaskHelper
from pfs.datamodel.target import TargetData, TargetObservations
from pfs.datamodel import utils
from pfs.datamodel import FluxTable
from pfs.datamodel.wavelengthArray import WavelengthArray


def makePfsDesign(tracts, patches, fiberIds, ras, decs, catIds, objIds, objectMags):
    """
        Make and return a PfsDesign with real information
    """
    nFiber = len(fiberIds)
    fiberMag = np.empty((nFiber, 5))
    for i in range(nFiber):
        fiberMag[i] = objectMags[i]
    raBoresight = np.median(ras)
    decBoresight = np.median(decs)
    targetTypes = np.array([1 for i in range(nFiber)], dtype='i4')
    fiberMags = [[mag] for mag in fiberMag[:, 0]]
    filterNames = [['g'] for i in range(nFiber)]
    pfiNominals = np.zeros((nFiber, 2))
    pfsDesignId = utils.calculate_pfsDesignId(fiberIds, ras, decs)
    return PfsDesign(pfsDesignId=pfsDesignId, raBoresight=raBoresight, decBoresight=decBoresight,
                     fiberId=fiberIds, tract=tracts, patch=patches, ra=ras, dec=decs,
                     catId=catIds, objId=objIds, targetType=targetTypes,
                     fiberMag=fiberMags, filterNames=filterNames, pfiNominal=pfiNominals)


def makePfsConfig(pfsDesignId, visit0, tracts, patches, fiberIds, ras, decs, catIds, objIds, objectMags):
    """
        Make and return a PfsConfig with real information
    """
    nFiber = len(fiberIds)
    fiberMag = np.empty((nFiber, 5))
    for i in range(nFiber):
        fiberMag[i] = objectMags[i]
    raBoresight = np.median(ras)
    decBoresight = np.median(decs)
    targetTypes = np.array([1 for i in range(nFiber)], dtype='i4')
    fiberMags = [[mag] for mag in fiberMag[:, 0]]
    filterNames = [['g'] for i in range(nFiber)]
    pfiNominals = np.zeros((nFiber, 2))
    pfiCenters = np.zeros((nFiber, 2))
    return PfsConfig(pfsDesignId=pfsDesignId, visit0=visit0, raBoresight=raBoresight, decBoresight=decBoresight,
                     fiberId=fiberIds, tract=tracts, patch=patches, ra=ras, dec=decs,
                     catId=catIds, objId=objIds, targetType=targetTypes,
                     fiberMag=fiberMags, filterNames=filterNames,
                     pfiCenter=pfiCenters, pfiNominal=pfiNominals)


def makePfsObjects(pfsConfig, visit0, pfsArmSet, minWavelength, maxWavelength, dWavelength):
    minWl = minWavelength
    maxWl = maxWavelength
    dWl = dWavelength
    #wavelength = minWl + dWl * np.arange(int((maxWl - minWl) / dWl) + 1, dtype=float)
    length = int((maxWl - minWl) / dWl) + 1
    wavelength = WavelengthArray(minWl, maxWl, length)

    def combine(spectra, flags):
        """Combine spectra

        Parameters
        ----------
        spectra : iterable of `pfs.datamodel.PfsSpectra`
            List of spectra to combine. These should already have been
            resampled to a common wavelength representation.
        flags : `pfs.datamodel.MaskHelper`
            Mask interpreter, for identifying bad pixels.

        Returns
        -------
        wavelength : `numpy.ndarray` of `float`
            Wavelengths for combined spectrum.
        flux : `numpy.ndarray` of `float`
            Flux measurements for combined spectrum.
        sky : `numpy.ndarray` of `float`
            Sky measurements for combined spectrum.
        covar : `numpy.ndarray` of `float`
            Covariance matrix for combined spectrum.
        mask : `numpy.ndarray` of `int`
            Mask for combined spectrum.
        """
        archetype = spectra[0]
        mask = np.zeros_like(archetype.mask)
        flux = np.zeros_like(archetype.flux)
        sky = np.zeros_like(archetype.sky)
        covar = np.zeros_like(archetype.covar)
        sumWeights = np.zeros_like(archetype.flux)

        for ss in spectra:
            good = ((ss.mask & ss.flags.get(*["NO_DATA"])) == 0) & (ss.covar[:, 0] > 0)
            weight = np.zeros_like(ss.flux)
            weight[good] = 1.0 / ss.covar[:, 0][good]
            flux += ss.flux * weight
            sky += ss.sky * weight
            mask[good] |= ss.mask[good]
            sumWeights += weight

        good = sumWeights > 0
        flux[good] /= sumWeights[good]
        sky[good] /= sumWeights[good]
        covar[:, 0][good] = 1.0 / sumWeights[good]
        covar[:, 0][~good] = np.inf
        covar[:, 1:2] = np.where(good, 0.0, np.inf)[:, np.newaxis]
        mask[~good] = flags["NO_DATA"]
        covar2 = np.zeros((1, 1), dtype=archetype.covar.dtype)
        Struct = collections.namedtuple('Struct', 'wavelength flux sky covar mask covar2')
        return Struct(archetype.wavelength, flux, sky, covar, mask, covar2)
#        return Struct(wavelength=archetype.wavelength, flux=flux, sky=sky, covar=covar, mask = mask, covar2 = covar2)

    def mergeSpectra(spectraList, identityKeys):
        """Combine all spectra from the same exposure

        All spectra should have the same fibers, so we simply iterate over the
        fibers, combining each spectrum from that fiber.

        Parameters
        ----------
        spectraList : iterable of `pfs.datamodel.PfsSpectra`
            List of spectra to coadd.
        identityKeys : iterable of `str`
            Keys to select from the input spectra's ``identity`` for the
            merged spectra's ``identity``.

        Returns
        -------
        result : `pfs.datamodel.PfsMerged`
            Merged spectra.
        """
        archetype = spectraList[0]
        identity = {key: archetype.identity[key] for key in identityKeys}
        fiberId = archetype.fiberId
        if any(np.any(ss.fiberId != fiberId) for ss in spectraList):
            raise RuntimeError("Selection of fibers differs")
        resampled = [ss.resample(wavelength) for ss in spectraList]
        flags = MaskHelper.fromMerge([ss.flags for ss in spectraList])
        combination = combine(resampled, flags)
        return PfsMerged(identity, fiberId, combination.wavelength, combination.flux, combination.mask, combination.sky, combination.covar, flags, archetype.metadata), combination.covar2
    """ make arm merged spectra """
    sm, covar2 = mergeSpectra(pfsArmSet, ["visit", "spectrograph"])
    # print(sm.flux)
    # import matplotlib.pyplot as plt
    # plt.plot(sm.wavelength[1], sm.flux[1])

    """ make pfsObject """
    pfsObjects = []
    pfsVisitHashes = []
    for i in range(len(sm.fiberId)):
        fiberId = pfsConfig.fiberId[i]
        catId = pfsConfig.catId[i]
        objId = pfsConfig.objId[i]
        tract = pfsConfig.tract[i]
        patch = pfsConfig.patch[i]
        ra = pfsConfig.ra[i]
        dec = pfsConfig.dec[i]
        targetType = pfsConfig.targetType[i]
        fiberMags = collections.defaultdict(list)
        for ff, mag in zip(pfsConfig.filterNames[i], pfsConfig.fiberMag[i]):
            fiberMags[ff].append(mag)
        targetData = TargetData(catId, tract, patch, objId, ra, dec, targetType, dict(**fiberMags))
        identityList = [{'visit': visit0}]
        pfiNominal = pfsConfig.pfiNominal[i]
        pfiCenter = pfsConfig.pfiCenter[i]
        observations = TargetObservations(identityList,
                                          np.array([fiberId]),
                                          np.array([pfiNominal]),
                                          np.array([pfiCenter])
                                          )
        fluxTable = FluxTable(np.array([d.wavelength[0] for d in pfsArmSet]).flatten(),
                              np.array([d.flux[0] for d in pfsArmSet]).flatten(),
                              np.array([d.covar[0][0] for d in pfsArmSet]).flatten(),
                              np.array([d.mask[0] for d in pfsArmSet]).flatten(),
                              MaskHelper(missing=1)
                              )
        pfsObject = PfsObject(targetData,
                              observations,
                              wavelength,
                              sm.flux[i],
                              sm.mask[i],
                              sm.sky[i],
                              sm.covar[i],
                              covar2,
                              sm.flags,
                              fluxTable
                              )
        pfsObjects.append(pfsObject)
        pfsVisitHash = observations.getIdentity()['pfsVisitHash']
        pfsVisitHashes.append(pfsVisitHash)
    return pfsObjects, pfsVisitHashes
