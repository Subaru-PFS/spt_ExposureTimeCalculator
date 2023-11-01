# -*- coding: utf-8 -*-

from __future__ import print_function, division

import sys
import os
from os import path
import numpy as np
from scipy.interpolate import interp1d
import collections
import importlib

HOME_DIR = path.dirname(path.abspath(__file__))
''' import datamodel module '''
sys.path.append(HOME_DIR + "/datamodel/python")
from pfs.datamodel.pfsConfig import PfsDesign, PfsConfig
from pfs.datamodel.drp import PfsArm, PfsMerged, PfsObject, PfsFiberArray
from pfs.datamodel.masks import MaskHelper
from pfs.datamodel.target import Target
from pfs.datamodel.observations import Observations
from pfs.datamodel.identity import Identity
from pfs.datamodel import utils
from pfs.datamodel import FluxTable
from pfs.datamodel.wavelengthArray import WavelengthArray


def getPfsVersions(prefix="VERSION_"):
    """Retrieve the software versions of PFS DRP-2D software

    The version values come from the ``<module>.version.__version__``
    attribute, which gets set at build time.

    The versions are in a format suitable for inclusion in a FITS header.

    Parameters
    ----------
    prefix : `str`, optional
        Prefix to add to software product name, to distinguish in the FITS
        header.

    Returns
    -------
    versions : `dict` (`str`: `str`)
        Versions of software products.
    """
    versions = {prefix + 'datamodel': '5.2'}
    # for name, module in (("datamodel", "pfs.datamodel"),
    #                     ("obs_pfs", "lsst.obs.pfs"),
    #                     ("drp_stella", "pfs.drp.stella"),
    #                     ):
    #    importlib.import_module(module + ".version")
    #    versions[prefix + name] = sys.modules[module + ".version"].__version__
    return versions


def makePfsDesign(tracts, patches, fiberIds, ras, decs, catIds, objIds, fiberMags, filterNames, arms='br'):
    """
        Make and return a PfsDesign with real information
    """
    nFiber = len(fiberIds)
    #fiberMag = np.empty((nFiber, 5))
    # for i in range(nFiber):
    #    fiberMag[i] = objectMags[i]
    raBoresight = np.median(ras)
    decBoresight = np.median(decs)
    posAng = 0.0
    targetTypes = np.array([1 for i in range(nFiber)], dtype='i4') # SCIENCE
    fiberStatus = np.array([1 for i in range(nFiber)], dtype='i4') # GOOD
    epoch = np.array(['J2000.0' for _ in range(nFiber)], dtype=object)
    pmRa = np.array([0.0 for _ in range(nFiber)], dtype='f4')
    pmDec = np.array([0.0 for _ in range(nFiber)], dtype='f4')
    parallax = np.array([0.0 for _ in range(nFiber)], dtype='f4')
    proposalId = np.array(['S24B-001QN' for _ in range(nFiber)], dtype=object)
    obCode = np.array([f'{oid}' for oid in objIds], dtype=object)
    fiberFlux = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)], dtype='f4')
    fiberFluxErr = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)], dtype='f4')
    psfFlux = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)], dtype='f4')
    psfFluxErr = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)], dtype='f4')
    totalFlux = np.array([[10**(-0.4*(m-8.9))*1e+09 for m in fiberMags[i]] for i in range(nFiber)], dtype='f4')
    totalFluxErr = np.array([[0.0 for _ in fiberMags[i]] for i in range(nFiber)], dtype='f4')
    pfiNominals = np.zeros((nFiber, 2))
    pfsDesignId = utils.calculate_pfsDesignId(fiberIds, ras, decs)
    return PfsDesign(pfsDesignId=pfsDesignId, 
                     raBoresight=raBoresight, decBoresight=decBoresight,
                     posAng=posAng, arms=arms, fiberId=fiberIds, 
                     tract=tracts, patch=patches, ra=ras, dec=decs,
                     catId=catIds, objId=objIds, 
                     targetType=targetTypes, fiberStatus=fiberStatus,
                     epoch=epoch, pmRa=pmRa, pmDec=pmDec, parallax=parallax,
                     proposalId=proposalId, obCode=obCode,
                     fiberFlux=fiberFlux, psfFlux=psfFlux, totalFlux=totalFlux,
                     fiberFluxErr=fiberFluxErr, psfFluxErr=psfFluxErr, totalFluxErr=totalFluxErr,
                     filterNames=filterNames, pfiNominal=pfiNominals,
                     guideStars=None, designName='ETC', 
                     variant=None, designId0=pfsDesignId)


def makePfsConfig(pfsDesignId, visit0, tracts, patches, fiberIds, ras, decs, catIds, objIds, fiberMags, filterNames, arms='br'):
    """
        Make and return a PfsConfig with real information
    """
    nFiber = len(fiberIds)
    raBoresight = np.median(ras)
    decBoresight = np.median(decs)
    posAng = 0.0
    targetTypes = np.array([1 for i in range(nFiber)], dtype='i4') # SCIENCE
    fiberStatus = np.array([1 for i in range(nFiber)], dtype='i4') # GOOD
    epoch = np.array(['J2000.0' for _ in range(nFiber)])
    pmRa = np.array([0.0 for _ in range(nFiber)])
    pmDec = np.array([0.0 for _ in range(nFiber)])
    parallax = np.array([0.0 for _ in range(nFiber)])
    proposalId = np.array(['S24B-001QN' for _ in range(nFiber)])
    obCode = np.array([f'{oid}' for oid in objIds])
    fiberFlux = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)])
    fiberFluxErr = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)])
    psfFlux = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)])
    psfFluxErr = np.array([[np.nan for _ in fiberMags[i]] for i in range(nFiber)])
    totalFlux = np.array([[10**(-0.4*(m-8.9))*1e+09 for m in fiberMags[i]] for i in range(nFiber)])
    totalFluxErr = np.array([[0.0 for _ in fiberMags[i]] for i in range(nFiber)])
    pfiNominals = np.zeros((nFiber, 2))
    pfiCenters = np.zeros((nFiber, 2))
    return PfsConfig(pfsDesignId=pfsDesignId, visit=visit0,
                     raBoresight=raBoresight, decBoresight=decBoresight, 
                     posAng=posAng, arms=arms, fiberId=fiberIds, 
                     tract=tracts, patch=patches, ra=ras, dec=decs,
                     catId=catIds, objId=objIds, 
                     targetType=targetTypes, fiberStatus=fiberStatus,
                     epoch=epoch, pmRa=pmRa, pmDec=pmDec, parallax=parallax,
                     proposalId=proposalId, obCode=obCode,
                     fiberFlux=fiberFlux, psfFlux=psfFlux, totalFlux=totalFlux,
                     fiberFluxErr=fiberFluxErr, psfFluxErr=psfFluxErr, totalFluxErr=totalFluxErr,
                     filterNames=filterNames, 
                     pfiCenter=pfiCenters, pfiNominal=pfiNominals, 
                     guideStars=None, designName='ETC', 
                     variant=None, designId0=pfsDesignId)


def makePfsObject(pfsConfig, pfsArmSet, visits, minWavelength, maxWavelength, dWavelength):

    def interpolateFlux(fromWavelength, fromFlux, toWavelength, fill=0.0):
        """Interpolate a flux-like spectrum

        Basic linear interpolation, suitable for fluxes and flux-like (e.g., maybe
        variances) quantities.

        Parameters
        ----------
        fromWavelength : array-like of `float`
            Source wavelength array.
        fromFlux : array-like of `float`
            Source flux(-like) array.
        toWavelength : array-like of `float`
            Target wavelength array.
        fill : `float`, optional
            Fill value.

        Returns
        -------
        toFlux : `numpy.ndarray` of `float`
            Target flux-(like) array.
        """
        with np.errstate(invalid="ignore"):
            return interp1d(fromWavelength, fromFlux, kind="linear", bounds_error=False,
                            fill_value=fill, copy=True, assume_sorted=True)(toWavelength)

    def interpolateMask(fromWavelength, fromMask, toWavelength, fill=0):
        """Interpolate a mask spectrum

        Linear interpolation for masks.

        Parameters
        ----------
        fromWavelength : array-like of `float`
            Source wavelength array.
        fromMask : array-like of `float`
            Source mask array.
        toWavelength : array-like of `float`
            Target wavelength array.
        fill : `float`, optional
            Fill value.

        Returns
        -------
        toMask : `numpy.ndarray` of `float`
            Target mask array.
        """
        length = len(fromWavelength)
        with np.errstate(invalid="ignore"):
            index = interp1d(fromWavelength, np.arange(length), kind="linear", bounds_error=False,
                             fill_value=fill, copy=True, assume_sorted=True)(toWavelength)
        intIndex = index.astype(int)
        result = np.full(toWavelength.shape, fill, dtype=fromMask.dtype)
        intIndex[(intIndex == index) & (index > 0)] -= 1  # Linear interpolation takes the index before
        select = (intIndex >= 0) & (intIndex < length - 1)
        result[select] = fromMask[intIndex[select]] | fromMask[intIndex[select] + 1]
        return result

    def resample(spectra, wavelength, fiberId=None):
        """Construct a new PfsSpectra resampled to a common wavelength vector
        Parameters
        ----------
        wavelength : `numpy.ndarray` of `float`
            New wavelength values (nm).
        fiberId : `numpy.ndarray` of int, optional
            Fibers to resample. If ``None``, resample all fibers.
        Returns
        -------
        result : `PfsSpectra`
            Resampled spectra.
        """
        if fiberId is None:
            fiberId = spectra.fiberId

        numSpectra = len(fiberId)
        numSamples = len(wavelength)
        flux = np.empty((numSpectra, numSamples), dtype=spectra.flux.dtype)
        mask = np.empty((numSpectra, numSamples), dtype=spectra.mask.dtype)
        sky = np.empty((numSpectra, numSamples), dtype=spectra.sky.dtype)
        norm = np.ones((numSpectra, numSamples), dtype=spectra.sky.dtype)
        covar = np.zeros((numSpectra, 3, numSamples), dtype=spectra.covar.dtype)

        for ii, ff in enumerate(fiberId):
            jj = np.argwhere(spectra.fiberId == ff)[0][0]
            flux[ii] = interpolateFlux(spectra.wavelength[jj], spectra.flux[jj], wavelength)
            sky[ii] = interpolateFlux(spectra.wavelength[jj], spectra.sky[jj], wavelength)
            # XXX dropping covariance on the floor: just doing the variance for now
            covar[ii][0] = interpolateFlux(spectra.wavelength[jj], spectra.covar[jj][0], wavelength, fill=np.inf)
            mask[ii] = interpolateMask(spectra.wavelength[jj], spectra.mask[jj], wavelength,
                                       fill=spectra.flags["NO_DATA"]).astype(spectra.mask.dtype)

        return type(spectra)(spectra.identity, fiberId, np.concatenate([[wavelength]] * numSpectra),
                             flux, mask, sky, norm, covar, spectra.flags, spectra.metadata)

    def combine(spectra, flags):
        """Combine spectra

        Parameters
        ----------
        spectra : iterable of `pfs.datamodel.PfsFiberArraySet`
            List of spectra to combine. These should already have been
            resampled to a common wavelength representation.
        flags : `pfs.datamodel.MaskHelper`
            Mask interpreter, for identifying bad pixels.

        Returns
        -------
        wavelength : `numpy.ndarray` of `float`
            Wavelengths for combined spectrum.
        flux : `numpy.ndarray` of `float`
            Normalised flux measurements for combined spectrum.
        sky : `numpy.ndarray` of `float`
            Sky measurements for combined spectrum.
        norm : `numpy.ndarray` of `float`
            Normalisation of combined spectrum.
        covar : `numpy.ndarray` of `float`
            Covariance matrix for combined spectrum.
        mask : `numpy.ndarray` of `int`
            Mask for combined spectrum.
        """
        archetype = spectra[0]
        mask = np.zeros_like(archetype.mask)
        flux = np.zeros_like(archetype.flux)
        sky = np.zeros_like(archetype.sky)
        norm = np.zeros_like(archetype.norm)
        covar = np.zeros_like(archetype.covar)
        sumWeights = np.zeros_like(archetype.flux)

        for ss in spectra:
            with np.errstate(invalid="ignore", divide="ignore"):
                variance = ss.variance/ss.norm**2
                good = ((ss.mask & ss.flags.get(*["NO_DATA"])) == 0) & (variance > 0)

            weight = np.zeros_like(ss.flux)
            weight[good] = 1.0/variance[good]
            with np.errstate(invalid="ignore"):
                flux[good] += ss.flux[good]*weight[good]/ss.norm[good]
                sky[good] += ss.sky[good]*weight[good]/ss.norm[good]
                norm[good] += ss.norm[good]*weight[good]
            mask[good] |= ss.mask[good]
            sumWeights += weight

        good = sumWeights > 0
        flux[good] /= sumWeights[good]
        sky[good] /= sumWeights[good]
        norm[good] /= sumWeights[good]
        covar[:, 0][good] = 1.0/sumWeights[good]
        covar[:, 0][~good] = np.inf
        covar[:, 1:] = np.where(good, 0.0, np.inf)[:, np.newaxis]

        for ss in spectra:
            mask[~good] |= ss.mask[~good]
        mask[~good] |= flags["NO_DATA"]
        covar2 = np.zeros((1, 1), dtype=archetype.covar.dtype)
        Struct = collections.namedtuple('Struct', 'wavelength flux sky norm covar mask covar2')
        with np.errstate(invalid="ignore"):
            return Struct(wavelength=archetype.wavelength, flux=flux*norm, sky=sky*norm, norm=norm,
                          covar=covar*norm[:, np.newaxis, :]**2, mask=mask, covar2=covar2)

    def mergeSpectra(spectraList):
        """Combine all spectra from the same exposure

        All spectra should have the same fibers, so we simply iterate over the
        fibers, combining each spectrum from that fiber.

        Parameters
        ----------
        spectraList : iterable of `pfs.datamodel.PfsFiberArraySet`
            List of spectra to coadd.

        Returns
        -------
        result : `pfs.datamodel.PfsMerged`
            Merged spectra.
        """
        archetype = spectraList[0]
        identity = Identity.fromMerge([ss.identity for ss in spectraList])
        fiberId = archetype.fiberId
        if any(np.any(ss.fiberId != fiberId) for ss in spectraList):
            raise RuntimeError("Selection of fibers differs")
        #resampled = [ss.resample(wavelength) for ss in spectraList]
        resampled = [resample(ss, wavelength) for ss in spectraList]
        flags = MaskHelper.fromMerge([ss.flags for ss in spectraList])
        combination = combine(resampled, flags)
        notes = PfsMerged.NotesClass.empty(len(archetype))
        #notes = [None for _ in range(len(archetype))]
        for name in ["blackSpotId", "blackSpotDistance", "blackSpotCorrection"]:
            getattr(notes, name)[:] = getattr(archetype.notes, name)

        return PfsMerged(identity, fiberId, combination.wavelength, combination.flux, combination.mask,
                         combination.sky, combination.norm, combination.covar, flags, archetype.metadata, notes), combination.covar2


    def readVisit(spectraList, pfsConfig):
        """Read a single visit

        The input ``pfsArm`` files are read, and the 1D sky subtraction from
        ``sky1d`` is applied.

        Parameters
        ----------
        dataRefList : iterable of `lsst.daf.persistence.ButlerDataRef`
            Butler data references.

        Returns
        -------
        spectra : `list` of `pfs.datamodel.PfsSingle`
            Sky-subtracted, flux-calibrated arm spectra.
        pfsConfig : `pfs.datamodel.PfsConfig`
            Top-end configuration.
        """
        result = []
        for ss in spectraList:
            #result += [ss.extractFiber(PfsFiberArray, pfsConfig, fiberId) for fiberId in ss.fiberId]
            result += [ss[ss.fiberId==fiberId] for fiberId in ss.fiberId]
        Struct = collections.namedtuple('Struct', 'spectra pfsConfig')
        return Struct(spectra=result, pfsConfig=pfsConfig)

    def getTargetData(target, pfsConfigList, indices):
        """Generate a ``TargetData`` for this target

        We combine the various declarations about the target in the
        ``PfsConfig``s.

        Parameters
        ----------
        target : `types.SimpleNamespace`
            Struct with target identity (with ``catId``, ``tract``, ``patch``
            and ``objId``).
        pfsConfigList : iterable of `pfs.datamodel.PfsConfig`
            List of top-end configurations.
        indices : `numpy.ndarray` of `int`
            Indices for the fiber of interest in each of the ``pfsConfigList``.

        Returns
        -------
        result : `pfs.datamodel.TargetData`
            ``TargetData`` for this target
        """
        ra = np.array([pfsConfig.ra[ii] for pfsConfig, ii in zip(pfsConfigList, indices)], dtype='f8')
        dec = np.array([pfsConfig.dec[ii] for pfsConfig, ii in zip(pfsConfigList, indices)], dtype='f8')
        ra = ra.mean()
        dec = dec.mean()
        # radec = averageSpherePoint(radec)

        targetType = collections.Counter([pfsConfig.targetType[ii] for pfsConfig, ii in zip(pfsConfigList, indices)])
        if len(targetType) > 1:
            print("Multiple targetType for target %s (%s); using most common" % (target, targetType))
        targetType = targetType.most_common(1)[0][0]

        fiberMags = collections.defaultdict(list)
        for pfsConfig, ii in zip(pfsConfigList, indices):
            for ff, flx in zip(pfsConfig.filterNames[ii], pfsConfig.fiberFlux[ii]):
                fiberMags[ff].append(-2.5*np.log10(flx*1e-09)+8.9)
        for ff in fiberMags:
            mag = set(fiberMags[ff])
            if len(mag) > 1:
                #print("Multiple %s mag for target %s (%s); using average" % (ff, target, mag))
                mag = np.average(np.array(fiberMags[ff]))
            else:
                mag = mag.pop()
            fiberMags[ff] = mag

        return Target(target.catId, target.tract, target.patch, target.objId,
                      ra, dec, targetType, dict(**fiberMags))

    def getObservations(identityList, pfsConfigList, indices):
        """Construct a list of observations of the target

        Parameters
        ----------
        identityList : iterable of `dict`
            List of sets of keyword-value pairs that identify the observation.
        pfsConfigList : iterable of `pfs.datamodel.PfsConfig`
            List of top-end configurations.
        indices : `numpy.ndarray` of `int`
            Indices for the fiber of interest in each of the ``pfsConfigList``.

        Returns
        -------
        observations : `pfs.datamodel.TargetObservations`
            Observations of the target.
        """
        visit = np.array([ident["visit"] for ident in identityList])
        arm = [ident["arm"] for ident in identityList]
        spectrograph = np.array([ident["spectrograph"] for ident in identityList])
        pfsDesignId = np.array([pfsConfig.pfsDesignId for pfsConfig in pfsConfigList])
        fiberId = np.array([pfsConfig.fiberId[ii] for pfsConfig, ii in zip(pfsConfigList, indices)])
        pfiNominal = np.array([pfsConfig.pfiNominal[ii] for pfsConfig, ii in zip(pfsConfigList, indices)])
        pfiCenter = np.array([pfsConfig.pfiCenter[ii] for pfsConfig, ii in zip(pfsConfigList, indices)])
        return Observations(visit, arm, spectrograph, pfsDesignId, fiberId, pfiNominal, pfiCenter)

    ''' Main part '''
    minWl = minWavelength
    maxWl = maxWavelength
    dWl = dWavelength
    # wavelength = minWl + dWl * np.arange(int((maxWl - minWl) / dWl) + 1, dtype=float)
    length = int((maxWl - minWl) / dWl) + 1
    wavelength = WavelengthArray(minWl, maxWl, length)

    """ make arm merged spectra """
    merged, covar2 = mergeSpectra(pfsArmSet)

    #import matplotlib.pyplot as plt
    #plt.plot(merged.wavelength[0], merged.flux[0])
    #plt.plot(merged.wavelength[1], merged.flux[1])
    #plt.plot(merged.wavelength[2], merged.flux[2])

    """ make pfsObject """
    n_visit = len(visits)
    n_fiber = len(merged.fiberId)
    n_filter = len(pfsConfig.filterNames[0])
    arms = np.array(['b' for v in visits], dtype='U')
    spectrographs = np.array([1 for v in visits], dtype='i4')
    #identityList = [pfsArmSet[0].identity.getDict()]
    spectraList = [merged]
    pfsConfigList = [pfsConfig for v in visits]

    data = [readVisit(spectraList, pfsConfig)]
    targetList = []
    Struct = collections.namedtuple('Struct', 'catId tract patch objId')
    spectra = collections.defaultdict(list)
    for dd in data:
        for ss, config in zip(dd.spectra, dd.pfsConfig):
            spectra[(config.catId, config.tract, config.patch, config.objId)].append(ss)
            targetList.append(Struct(config.catId, config.tract, config.patch, config.objId))
    pfsObjects = []
    pfsVisitHashes = []
    for i, target in enumerate(targetList):
        indices = [pfsConfig.selectTarget(target.catId, target.tract, target.patch, target.objId) for
                   pfsConfig in pfsConfigList]
        targetData = getTargetData(target, pfsConfigList, indices)
        identityList = [{'visit': v, 'arm': 'merged', 'spectrograph': 1} for v in visits]
        observations = getObservations(identityList, pfsConfigList, indices)
        ss = spectra[(target.catId, target.tract, target.patch, target.objId)]
        spectrumList = spectra[(target.catId, target.tract, target.patch, target.objId)]
        flags = MaskHelper.fromMerge([ss.flags for ss in spectrumList])
        combination = combine(spectrumList, flags)
        #fluxTable = fluxTable.run(identityList, spectrumList, flags)
        fluxTable = FluxTable(np.array([d.wavelength[i] for d in pfsArmSet]).flatten(),
                              np.array([d.flux[i] for d in pfsArmSet]).flatten(),
                              np.array([d.covar[i][0] for d in pfsArmSet]).flatten(),
                              np.array([d.mask[i] for d in pfsArmSet]).flatten(),
                              MaskHelper(missing=1)
                              )
        #print(combination.covar2.shape)
        coadd = PfsObject(targetData, observations, wavelength, combination.flux[0],
                          combination.mask[0], combination.sky[0], combination.covar[0], combination.covar2, flags,
                          getPfsVersions(), fluxTable)
        pfsObjects.append(coadd)
        pfsVisitHash = utils.calculatePfsVisitHash(visits)
        pfsVisitHashes.append(pfsVisitHash)
    return pfsObjects, pfsVisitHashes
