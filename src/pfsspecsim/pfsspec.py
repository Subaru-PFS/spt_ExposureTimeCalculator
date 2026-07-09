# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from . import dm_utils

WAV_ERR_SHIFT = 0.002
PSF_VAR_SIGMA = 0.01
# F2: the psfvar Gaussian-kernel abscissa never varies (same range/step
# every exposure, every realization); compute it once.
_PSF_VAR_KERNEL_X = np.arange(-100, 101, 1)


def arm_name(arm_num):
    arm_num = np.array(arm_num)
    return np.where(
        arm_num == 0, "b", np.where(arm_num == 1, "r", np.where(arm_num == 2, "n", "m"))
    )


def arm_number(armStr):
    return dict(b=0, r=1, n=2, m=3)[armStr]


def _replicate_ids(
    n, objId, catId, fiberId, ra, dec, tract, patch, fiberMag, filterName
):
    """Duplicate a single object's identity fields into arrays of length
    `n`: `objId`/`fiberId` become `n` successive integers (`arange`), the
    rest are constant-filled (`full`); `fiberMag`/`filterName` become `n`
    aliased copies of the same list.
    """
    objIds = np.arange(objId, objId + n)
    catIds = np.full(n, catId, dtype=np.int32)
    fiberIds = np.arange(fiberId, fiberId + n)
    ras = np.full(n, ra, dtype=np.float32)
    decs = np.full(n, dec, dtype=np.float32)
    tracts = np.full(n, tract, dtype=np.int32)
    patches = np.full(n, patch, dtype="U3")
    fiberMags = [fiberMag for _ in range(n)]
    filterNames = [filterName for _ in range(n)]
    return objIds, catIds, fiberIds, ras, decs, tracts, patches, fiberMags, filterNames


def calculateFiberMagnitude(wav, mag, filterName):
    """Calculate the average magnitude over the bandpass"""
    #                           50%/nm  50%/nm peak
    filterBandpasses = dict(
        g=(399.5, 546.5, 0.97),
        r=(542.5, 696.5, 0.95),
        i=(698.5, 853.3, 0.90),  # i2
        z=(852.5, 932.0, 0.97),
        y=(943.0, 1072.0, 0.95),
    )
    wav0, wav1, peak = filterBandpasses[filterName]

    counts = np.exp(-mag)
    bandpass = np.where(np.logical_and(wav >= wav0, wav <= wav1), peak, 0)
    fiberMag = -np.log(
        np.trapezoid(bandpass * counts, wav) / np.trapezoid(bandpass, wav)
    )
    return fiberMag


def write_ascii(aset, arms, asciiTable, outDir):
    """Write an ascii table"""

    outFile = Path(outDir) / asciiTable
    nFiber = len(aset[0].flux)
    for i in range(nFiber):
        outPath = f"{outFile}.dat" if nFiber == 1 else f"{outFile}.{i}.dat"
        with open(outPath, "w") as fd:
            fd.write("""#  1  WAVELENGTH  [nm]
#  2  FLUX        [nJy]
#  3  ERROR       [nJy]
#  4  MASK        [1=masked]
#  5  SKY         [nJy]
#  6  ARM         [0=blue,1=red,2=NIR,3=redMR]
""")

            for a, armStr in zip(aset, arms):
                lam = a.wavelength[i]
                flux = a.flux[i]
                sigma = np.sqrt(a.covar[i][0])
                mask = a.mask[i]
                sky = a.sky[i]
                armNum = np.full(len(lam), arm_number(armStr))
                np.savetxt(
                    fd,
                    np.column_stack([lam, flux, sigma, mask, sky, armNum]),
                    fmt="%8.3f %12.4e %12.4e %2d %12.4e %1d",
                )


def strToBool(val):
    if val.lower() in ("1", "t", "true"):
        return True
    elif val.lower() in ("0", "f", "false"):
        return False
    else:
        raise ValueError(f'Unable to interpret "{val}" as bool')


def plotPfsObject(pfsObjects, fig=None):
    if fig is None:
        fig = plt.figure(figsize=(7, 4))
    axe = fig.add_subplot()
    axe.set_xlabel("wavelength (nm)")
    axe.set_ylabel("flux (nJy)")
    for pfsObject in pfsObjects:
        axe.plot(pfsObject.wavelength, pfsObject.flux)
    plt.show()


def plotPfsArm(pfsArmSet, fig=None):
    if fig is None:
        fig = plt.figure(figsize=(7, 4))
    axe = fig.add_subplot()
    axe.set_xlabel("wavelength (nm)")
    axe.set_ylabel("flux (electron/pix)")
    for i in range(len(pfsArmSet[0])):
        axe.plot(pfsArmSet[0].wavelength[i], pfsArmSet[0].flux[i], color=f"C{i}")
        axe.plot(pfsArmSet[1].wavelength[i], pfsArmSet[1].flux[i], color=f"C{i}")
        axe.plot(pfsArmSet[2].wavelength[i], pfsArmSet[2].flux[i], color=f"C{i}")
    plt.show()


class Pfsspec:

    def __init__(self):
        self.params = {
            "EXP_NUM": "4",
            "MAG_FILE": "22.5",
            "countsMin": "0.1",
            "etcFile": "out/ref.snc.dat",
            "nrealize": "1",
            "outDir": "out",
            "asciiTable": "None",
            "ra": 150.0,
            "dec": 2.0,
            "tract": 0,
            "patch": "0,0",
            "visit0": 1,
            "catId": 0,
            "objId": 1,
            "fiberId": 1,
            "fiberMag": [22.5, 22.5, 22.5, 22.5, 22.5],
            "filterName": ["hcs_g", "hcs_r", "hcs_i", "hcs_z", "hcs_y"],
            "spectrograph": 1,
            "pfsConfigFull": "f",
            "writeFits": "t",
            "writePfsArm": "t",
            "plotArmSet": False,
            "plotObject": False,
            "SKY_SUB_FLOOR": "0.01",
            "SKY_SUB_MODE": "random",
            "SKY_SUB_SEED": 0,
        }

    def set_param(self, param_name, param_value):
        if param_name in self.params:
            try:
                self.params[param_name] = param_value
            except Exception:
                print("Error!")
        else:
            print(f"param_name {param_name} can not be recognized ...")
        return 0

    def load_param_file(self, filename):
        try:
            with open(filename, "r") as f:
                for line in f:
                    a = line.split()
                    if line[0] != "#" and len(a) > 0:
                        self.params[a[0]] = a[1]
        except OSError as e:
            raise OSError(f"Error: maybe file not found: {e}") from e
        return 0

    def make_sim_spec(self):
        self.outdir = self.params["outDir"]
        self.tract = self.params["tract"]
        self.patch = self.params["patch"]
        self.visit0 = int(self.params["visit0"])
        self.fiberId = self.params["fiberId"]
        self.ra = self.params["ra"]
        self.dec = self.params["dec"]
        self.catId = self.params["catId"]
        self.objId = self.params["objId"]
        self.spectrograph = self.params["spectrograph"]
        try:
            self.mag_file = f"{float(self.params['MAG_FILE']):.4e}"
        except (ValueError, TypeError):
            self.mag_file = self.params["MAG_FILE"]
        self.fiberMag = self.params["fiberMag"]
        self.filterName = self.params["filterName"]

        # do some checks
        self.writeFits = strToBool(self.params["writeFits"])
        self.writePfsArm = strToBool(self.params["writePfsArm"])
        self.plotArmSet = self.params["plotArmSet"]
        self.plotObject = self.params["plotObject"]
        self.asciiTable = self.params["asciiTable"]
        self.pfsConfigFull = strToBool(self.params["pfsConfigFull"])
        self.sky_sub_err = float(self.params["SKY_SUB_FLOOR"])
        self.sky_sub_mode = self.params["SKY_SUB_MODE"]
        self.sky_sub_seed = self.params["SKY_SUB_SEED"]
        nrealize = int(self.params["nrealize"])
        nexp = int(self.params["EXP_NUM"])

        try:
            if len(self.fiberId) > 0:
                self.multi_info = 1
            else:
                self.multi_info = 0
        except TypeError:
            self.multi_info = 0

        if not self.writeFits and not self.asciiTable:
            raise ValueError(
                "Please specify asciiTable or omit writeFits (or say writeFits true)"
            )
        if not Path(self.outdir).exists():
            try:
                Path(self.outdir).mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise OSError(f"Unable to create outDir: {e}") from e
        if nrealize <= 0:
            raise ValueError("Please specify at least one realization")

        # check magfile
        if Path(self.mag_file).exists():
            dat = np.loadtxt(self.mag_file)
            nobj = dat.shape[1] - 1
        else:
            nobj = 1
        if nobj > 1:
            if nrealize > 1:
                raise ValueError(
                    "The number of realization should be one for multiple input template"
                )
            else:
                if self.multi_info == 0:
                    (
                        objIds,
                        catIds,
                        fiberIds,
                        ras,
                        decs,
                        tracts,
                        patches,
                        fiberMags,
                        filterNames,
                    ) = _replicate_ids(
                        nobj,
                        self.objId,
                        self.catId,
                        self.fiberId,
                        self.ra,
                        self.dec,
                        self.tract,
                        self.patch,
                        self.fiberMag,
                        self.filterName,
                    )
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
                        (
                            objIds,
                            catIds,
                            fiberIds,
                            ras,
                            decs,
                            tracts,
                            patches,
                            fiberMags,
                            filterNames,
                        ) = _replicate_ids(
                            nrealize,
                            self.objId[0],
                            self.catId[0],
                            self.fiberId[0],
                            self.ra[0],
                            self.dec[0],
                            self.tract[0],
                            self.patch[0],
                            self.fiberMag[0],
                            self.filterName[0],
                        )
                except TypeError:
                    (
                        objIds,
                        catIds,
                        fiberIds,
                        ras,
                        decs,
                        tracts,
                        patches,
                        fiberMags,
                        filterNames,
                    ) = _replicate_ids(
                        nrealize,
                        self.objId,
                        self.catId,
                        self.fiberId,
                        self.ra,
                        self.dec,
                        self.tract,
                        self.patch,
                        self.fiberMag,
                        self.filterName,
                    )
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
                except TypeError:
                    objIds = np.array([self.objId])
                    catIds = np.array([self.catId])
                    fiberIds = np.array([self.fiberId])
                    ras = np.array([self.ra])
                    decs = np.array([self.dec])
                    tracts = np.array([self.tract])
                    patches = np.array([self.patch])
                    fiberMags = [self.fiberMag]
                    filterNames = [self.filterName]
        """
            ## read input file ##
            # arm: 0-3 telling us which arm the data comes from (arm_name will convert to b, r, n, m)
            # wav: wavelength
            # nsv: per-pixel (instrumental + sky) variance (counts^2)
            # trn: conversion from I_nu to counts
            # smp: samplingFactor.  A fiddle factor for the Poisson noise in HgCdTe devices
            # skm: sky flux
        """
        try:
            # New (T14): the ETC now writes its SNC output as an Astropy
            # ECSV table with every resolved input parameter (incl.
            # EXP_NUM) embedded in table.meta["params"]. Column names
            # correspond 1:1 to the old whitespace-delimited usecols
            # (0, 2, 5, 8, 9, 10) below.
            etcTable = Table.read(self.params["etcFile"], format="ascii.ecsv")
            nexp_etc = etcTable.meta["params"]["exp_num"]
            arm = np.asarray(etcTable["arm"], dtype=float)
            wav = np.asarray(etcTable["wavelength"], dtype=float)
            nsv = np.asarray(etcTable["noise_variance"], dtype=float)
            trn = np.asarray(etcTable["conversion_factor"], dtype=float)
            smp = np.asarray(etcTable["sampling_factor"], dtype=float)
            skm = np.asarray(etcTable["sky"], dtype=float)
        except Exception:
            # Fall back to the legacy plain-text ETC output format
            # (insurance for old-format files predating the ECSV switch).
            nexp_etc = None
            with open(self.params["etcFile"], "r") as f:
                for line in f.readlines():
                    if "EXP_NUM" in line:
                        nexp_etc = int(line.split()[2])
            if nexp_etc is None:
                raise ValueError(
                    f"Could not read {self.params['etcFile']!r} as either "
                    "an Astropy ECSV table (missing/unreadable "
                    "table.meta['params']['exp_num']) or a legacy "
                    "plain-text ETC output file (no 'EXP_NUM' header line "
                    "found); both parsing attempts failed."
                )
            arm, wav, nsv, trn, smp, skm = np.loadtxt(
                self.params["etcFile"], usecols=(0, 2, 5, 8, 9, 10), unpack=True
            )

        # remove sky systematics
        skm_sysref = skm.copy()
        skmp = np.roll(skm_sysref, 1)
        skmp[0] = 0.0
        skmm = np.roll(skm_sysref, -1)
        skmm[-1] = 0.0
        skm_sysref = np.amax([skm_sysref, skmm, skmp], axis=0)
        nsv_sys = (self.sky_sub_err * np.sqrt(nexp_etc) * skm_sysref) ** 2
        nsv_rnd = nsv - nsv_sys
        arm = arm.astype(int)
        trn[trn < 1.0e26] = 1.0e26

        # load magnitude or filename
        if Path(self.mag_file).exists():
            dat = np.loadtxt(self.mag_file)
            nobj = dat.shape[1] - 1
            _lam = dat[:, 0]
            mag = np.empty((len(wav), nobj))
            for i in range(nobj):
                _mag = dat[:, i + 1]
                mag[:, i] = np.interp(wav, _lam, _mag)
        else:
            nobj = 1
            mag = np.full((len(wav), nobj), float(self.mag_file))
        # F5: broadcast instead of copying the same column `nobj` times
        # (`trn`/`smp`/`nsv_rnd`/the rescaled `nsv_sys` don't vary with i).
        trn_mtrx = np.broadcast_to(trn[:, None], (len(wav), nobj))
        smp_mtrx = np.broadcast_to(smp[:, None], (len(wav), nobj))
        nsv_rnd_mtrx = np.broadcast_to(nsv_rnd[:, None], (len(wav), nobj))
        nsv_sys_mtrx = np.broadcast_to(
            (nsv_sys * float(nexp) / float(nexp_etc))[:, None], (len(wav), nobj)
        )

        # calculate the flux etc. in observed units
        fnu = 10 ** (-0.4 * (mag + 48.6))
        fnu_in_njy = fnu / 1e-32
        counts = trn_mtrx * fnu
        if (counts == 0).any():
            print(
                "counts == 0 detected in some pixels; setting to "
                f"{float(self.params['countsMin']):g} for variance",
                file=sys.stderr,
            )
            # version of counts with zero pixels replaced
            countsp = np.where(counts == 0, float(self.params["countsMin"]), counts)
        else:
            countsp = counts
        if self.sky_sub_mode == "random":
            snr1 = (
                countsp
                / np.sqrt(smp_mtrx * countsp + (nsv_rnd_mtrx + nsv_sys_mtrx))
                * np.sqrt(nexp)
            )
            snr2 = (
                countsp
                / np.sqrt(smp_mtrx * countsp + (nsv_rnd_mtrx + nsv_sys_mtrx))
                * np.sqrt(nexp)
            )
        else:
            snr1 = countsp / np.sqrt(smp_mtrx * countsp + nsv_rnd_mtrx) * np.sqrt(nexp)
            snr2 = (
                countsp
                / np.sqrt(smp_mtrx * countsp + (nsv_rnd_mtrx + nsv_sys_mtrx))
                * np.sqrt(nexp)
            )
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
        arms = np.array(
            sorted(set(arm), key=lambda x: dict(b=0, r=1, m=1.5, n=2)[x])
        )  # unique values of arm
        """
            Create and populate the objects corresponding to the datamodel

            First the parameters describing the observation, in PfsDesign and PfsConfig
        """
        # `nobj > 1` (multi-object, one realization each) and `nobj == 1`
        # (single object, `nrealize` realizations) share the same
        # per-spectrum generation logic below, indexed by `n_spec`/`col`.
        is_multi_obj = nobj > 1
        n_spec = nobj if is_multi_obj else nrealize

        if is_multi_obj:
            objectMags = [
                [calculateFiberMagnitude(wav, mag[:, i], b) for b in "grizy"]
                for i in range(n_spec)
            ]
        else:
            # F3: `col` is always 0 here, so every entry is identical --
            # compute it once instead of `n_spec` times.
            single_mag = [calculateFiberMagnitude(wav, mag[:, 0], b) for b in "grizy"]
            objectMags = [single_mag for _ in range(n_spec)]

        pfsDesign = dm_utils.makePfsDesign(
            tracts, patches, fiberIds, ras, decs, catIds, objIds, fiberMags, filterNames
        )

        pfsConfig = dm_utils.makePfsConfig(
            pfsDesign.pfsDesignId,
            self.visit0,
            tracts,
            patches,
            fiberIds,
            ras,
            decs,
            catIds,
            objIds,
            fiberMags,
            filterNames,
        )

        """
            Create the PfsArm;  we'll put each realisation into a different fibre
            (currently the content is not the one expected in the real pfsArm files so FIXME)
        """
        metadata = {}
        mapper = {"NO_DATA": 1}
        flags = dm_utils.MaskHelper(**mapper)
        pfsArmSet = []
        for armStr in arms:
            thisArm = arm == armStr
            identity = dm_utils.Identity(
                visit=self.visit0,
                arm=armStr,
                spectrograph=self.spectrograph,
                pfsDesignId=pfsDesign.pfsDesignId,
            )
            # np.random.seed(self.sky_sub_seed)
            nPt = np.sum(thisArm)
            datalam = []
            dataflux = []
            datasky = []
            datanorm = []
            datamask = []
            datacovar = []
            # QUIRK, preserved verbatim: the `nobj > 1` (multi-object) and
            # `nobj == 1` (multi-realization) systems historically diverged
            # in three places, kept here as `is_multi_obj`-gated branches
            # rather than silently unified:
            #  - wavecalib: `nobj > 1` interpolates the wavelength-shifted
            #    sky residual directly on the coarse `wav[thisArm]` grid;
            #    the `nrealize` system resamples onto a 10x finer grid first.
            #  - psfvar: the PSF-variation kernel sigma is drawn from
            #    `psf_var_sigma_pix * [0.5, 2.0]` for `nobj > 1` vs
            #    `* [0.3, 3.0]` for `nrealize`.
            #  - datanorm: only populated by the `nrealize` system; `nobj > 1`
            #    leaves it empty (an all-ones norm array is presumably
            #    assumed downstream). Whether one system is simply behind
            #    the other has not been established; report to the maintainer.
            psfvar_sigma_lo, psfvar_sigma_hi = (
                (0.5, 2.0) if is_multi_obj else (0.3, 3.0)
            )

            # F2: `wav_fine`/`skm_sysref_fine` depend only on `thisArm`, not
            # on `i`/`col` -- hoisted out of the spec loop below (previously
            # recomputed identically `n_spec` times for psfvar, and for
            # wavecalib's `nrealize` branch).
            wav_fine = skm_sysref_fine = None
            if self.sky_sub_mode == "psfvar" or (
                self.sky_sub_mode == "wavecalib" and not is_multi_obj
            ):
                wav_fine = np.linspace(
                    min(wav[thisArm]), max(wav[thisArm]), len(wav[thisArm]) * 10
                )
                skm_sysref_fine = np.interp(wav_fine, wav[thisArm], skm_sysref[thisArm])

            # F1: 'random' mode needs no per-exposure loop at all; draw the
            # whole (n_spec, nPt) noise batch for this arm in one shot
            # instead of looping over `i` with one 1D draw + list append
            # each. Note np.random's draw order (and hence, for a fixed
            # seed, the exact realized values) differs from the pre-Stage-4
            # per-`i`/per-exposure loop; only the distribution is preserved.
            random_mode_flux = None
            if self.sky_sub_mode not in ("systematic", "wavecalib", "psfvar"):
                if is_multi_obj:
                    base = fnu_in_njy[thisArm, :].T
                    std = np.abs(sigma1[thisArm, :]).T
                else:
                    base = np.broadcast_to(fnu_in_njy[thisArm, 0], (n_spec, nPt))
                    std = np.broadcast_to(np.abs(sigma1[thisArm, 0]), (n_spec, nPt))
                random_mode_flux = base + np.random.normal(0.0, std)

            for i in range(n_spec):
                col = i if is_multi_obj else 0
                datalam.append(wav[thisArm])
                if self.sky_sub_mode == "systematic":
                    # F1: batch the `nexp` scalar `sky_res_fac` draws and the
                    # `nexp` x `nPt` flux-noise draws instead of looping
                    # python-side and averaging a list.
                    sky_res_fac = np.random.normal(0.0, self.sky_sub_err, size=nexp)
                    skyres = skm_sysref[thisArm][None, :] * sky_res_fac[:, None]
                    noise = np.random.normal(
                        0.0,
                        abs(sigma1[thisArm, col]) * np.sqrt(nexp),
                        size=(nexp, nPt),
                    )
                    flux = fnu_in_njy[thisArm, col][None, :] + noise + skyres
                    dataflux.append(np.nanmean(flux, axis=0))
                elif self.sky_sub_mode == "wavecalib":
                    wav_err_shift = np.random.normal(0.0, WAV_ERR_SHIFT, size=nexp)
                    if is_multi_obj:
                        skyres = np.array(
                            [
                                skm_sysref[thisArm]
                                - np.interp(
                                    wav[thisArm] + shift,
                                    wav[thisArm],
                                    skm_sysref[thisArm],
                                )
                                for shift in wav_err_shift
                            ]
                        )
                    else:
                        skyres = np.array(
                            [
                                np.interp(
                                    wav[thisArm],
                                    wav_fine,
                                    skm_sysref_fine
                                    - np.interp(
                                        wav_fine + shift, wav_fine, skm_sysref_fine
                                    ),
                                )
                                for shift in wav_err_shift
                            ]
                        )
                    noise = np.random.normal(
                        0.0,
                        abs(sigma1[thisArm, col]) * np.sqrt(nexp),
                        size=(nexp, nPt),
                    )
                    flux = fnu_in_njy[thisArm, col][None, :] + noise + skyres
                    dataflux.append(np.nanmean(flux, axis=0))
                elif self.sky_sub_mode == "psfvar":
                    psf_var_sigma_pix = PSF_VAR_SIGMA / 0.07 * 10
                    gauss_kernel_sigma1 = np.random.uniform(
                        psf_var_sigma_pix * psfvar_sigma_lo,
                        psf_var_sigma_pix * psfvar_sigma_hi,
                        size=nexp,
                    )
                    gauss_kernel_sigma2 = np.random.uniform(
                        psf_var_sigma_pix * psfvar_sigma_lo,
                        psf_var_sigma_pix * psfvar_sigma_hi,
                        size=nexp,
                    )
                    skyres = np.empty((nexp, nPt))
                    for j in range(nexp):
                        gauss_kernel1 = (
                            1.0 / np.sqrt(2 * np.pi * gauss_kernel_sigma1[j] ** 2)
                        ) * np.exp(
                            -1.0
                            * _PSF_VAR_KERNEL_X**2
                            / (2 * gauss_kernel_sigma1[j] ** 2)
                        )
                        gauss_kernel2 = (
                            1.0 / np.sqrt(2 * np.pi * gauss_kernel_sigma2[j] ** 2)
                        ) * np.exp(
                            -1.0
                            * _PSF_VAR_KERNEL_X**2
                            / (2 * gauss_kernel_sigma2[j] ** 2)
                        )
                        skm_sysref_fine_conv1 = np.convolve(
                            skm_sysref_fine, gauss_kernel1, mode="same"
                        )
                        skm_sysref_fine_conv2 = np.convolve(
                            skm_sysref_fine, gauss_kernel2, mode="same"
                        )
                        skyres_fine = skm_sysref_fine_conv1 - skm_sysref_fine_conv2
                        skyres[j] = np.interp(wav[thisArm], wav_fine, skyres_fine)
                    noise = np.random.normal(
                        0.0,
                        abs(sigma1[thisArm, col]) * np.sqrt(nexp),
                        size=(nexp, nPt),
                    )
                    flux = fnu_in_njy[thisArm, col][None, :] + noise + skyres
                    dataflux.append(np.nanmean(flux, axis=0))
                else:
                    dataflux.append(random_mode_flux[i])
                datasky.append(sky[thisArm])
                if not is_multi_obj:
                    datanorm.append([1.0 for _ in sky[thisArm]])
                datamask.append(msk[thisArm])
                covar = np.zeros(3 * nPt).reshape((3, nPt))
                covar[0] = sigma2[thisArm, col] ** 2
                datacovar.append(covar)
            pfsArm = dm_utils.PfsArm(
                identity=identity,
                fiberId=fiberIds,
                wavelength=np.array(datalam, dtype="f4"),
                flux=np.array(dataflux, dtype="f4"),
                mask=np.array(datamask, dtype="i4"),
                sky=np.array(datasky, dtype="f4"),
                norm=np.array(datanorm, dtype="f4"),
                covar=np.array(datacovar, dtype="f4"),
                flags=flags,
                metadata=metadata,
            )
            pfsArmSet.append(pfsArm)
        if self.plotArmSet:
            plotPfsArm(pfsArmSet)

        """
            Time for I/O
        """
        # writ to Fits
        if self.writeFits:
            print(f"pfsDesignId: {pfsDesign.pfsDesignId:#013x}")
            pfsDesign.write(dirName=self.outdir)  # pfsDesign file
            pfsConfig.write(dirName=self.outdir)  # pfsConfig file
            if self.writePfsArm:  # write pfsArm files
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
            pfsConfig=pfsConfig,
            pfsArmSet=pfsArmSet,
            visits=visits,
            minWavelength=350.0,
            maxWavelength=1260.0,
            dWavelength=0.08,
        )
        for pfsObject in pfsObjects:
            if self.writeFits:  # write FITS file
                pfsObject.write(self.outdir)
        if self.plotObject:  # plot pfsObject
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
