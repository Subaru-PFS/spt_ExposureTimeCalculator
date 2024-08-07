# -*- coding: utf-8 -*-

from __future__ import print_function, division

import sys
import os
from os import path
import numpy as np
import time
import subprocess
import multiprocessing

### Some parameters ###########################
### CAUTION: ##################################
### CHANGE BELOW ON YOUR OWN RESPONSIBILITY ###
###############################################
SKYMODELS = '11006'
OFFSET_FIB = '0.10'
offset = 0.01


def add_header(filename, seeing, zenith_ang, galactic_ext, moon_zenith_ang, moon_target_ang, moon_phase, texp, nexp, field_ang, mag_file, reff, line_flux, line_width):
    hdr = ''
    hdr += '#  SEEING: %s\n' % (seeing)
    hdr += '#  ZENITH_ANG: %s\n' % (zenith_ang)
    hdr += '#  GALACTIC_EXT: %s\n' % (galactic_ext)
    hdr += '#  MOON_ZENITH_ANG: %s\n' % (moon_zenith_ang)
    hdr += '#  MOON_TARGET_ANG: %s\n' % (moon_target_ang)
    hdr += '#  MOON_PHASE: %s\n' % (moon_phase)
    hdr += '#  EXP_TIME: %s\n' % (texp)
    hdr += '#  EXP_NUM: %s\n' % (nexp)
    hdr += '#  FIELD_ANG: %s\n' % (field_ang)
    hdr += '#  MAG_FILE: %s\n' % (mag_file)
    hdr += '#  REFF: %s\n' % (reff)
    hdr += '#  LINE_FLUX: %s\n' % (line_flux)
    hdr += '#  LINE_WIDTH: %s\n' % (line_width)
    with open(filename, 'r') as f:
        dat = f.read()
    with open(filename, 'w') as f:
        f.write(hdr)
        f.write(dat)
    return 0


def calc_obscuration(dist):
    """Calculate the obscuration by the detector and structures around it

    The value depnds on the distance from the center on the telescope focal plane. See  https://pfspipe.ipmu.jp/jira/browse/PIPE2D-980 for more deitals.

    Parameters
    ----------
    dist : `float`
        The distance from the center on the focal plane.

    Returns
    -------
    obscuration : `float'
        The fraction of light obscuration by the detector anad it's strucuture.
    """
    obsc_etc = 0.19    # obscuration assumed in ETC (independent of PFI location)
    obsc_inner = 0.27  # obscuration at center (PIPE2D-980)
    obsc_outer = 0.36  # obscuration at edge   (PIPE2D-980)
    dist_edge = 0.675  # radius of FoV in deg. (FIELD_ANG)
    obsc = (obsc_outer - obsc_inner) / dist_edge * dist + obsc_inner
    corr = (1 - obsc) / (1 - obsc_etc)
    return obsc, corr


class Etc(object):
    """ PFS ETC

    See https://github.com/Subaru-PFS/spt_ExposureTimeCalculator for details

    Parameters
    ----------
    omp_num_threads : `int` (default:16)
        The number of threads to use for parallelization.

    Returns
    -------

    Examples
    ----------    
    """

    def __init__(self, omp_num_threads=16):
        self.params = {'SEEING': '0.80',
                       'ZENITH_ANG': '35.0',
                       'GALACTIC_EXT': '0.00',
                       'MOON_ZENITH_ANG': '30.0',
                       'MOON_TARGET_ANG': '60.0',
                       'MOON_PHASE': '0.125',
                       'EXP_TIME': '900',
                       'EXP_NUM': '4',
                       'FIELD_ANG': '0.45',
                       'MAG_FILE': '22.5',
                       'REFF': '0.3',
                       'LINE_FLUX': '1.0e-17',
                       'LINE_WIDTH': '70',
                       'NOISE_REUSED': 'N',
                       'MR_MODE': 'N',
                       'OVERWRITE': 'Y',
                       'INFILE_OIICat': '-',
                       'OUTFILE_OIICat': '-',
                       'minSNR': '9.0',
                       'degrade': '1.0',
                       'SKY_SUB_FLOOR': '0.01',
                       'DIFFUSE_STRAY': '0.02',
                       'throughput_model': '20240714',
                       'spectrograph': 'ave',
                       'OUTDIR': 'out',
                       'TMPDIR': 'tmp',
                       'BINDIR': 'bin',
                       'obscFoVDep': 'Y',
                       }
        self.params['OUTFILE_NOISE'] = os.path.join(
            self.params['OUTDIR'], 'ref.noise.dat')
        self.params['OUTFILE_SNC'] = os.path.join(
            self.params['OUTDIR'], 'ref.snc.dat')
        self.params['OUTFILE_SNL'] = os.path.join(
            self.params['OUTDIR'], 'ref.snl.dat')
        self.params['OUTFILE_OII'] = os.path.join(
            self.params['OUTDIR'], 'ref.sno2.dat')

        self.HOME_DIR = path.dirname(path.abspath(__file__))
        if os.path.exists(os.path.join(self.HOME_DIR, self.params['BINDIR'], "gsetc_omp.x")):
            self.ETC_SRC = os.path.join(
                self.HOME_DIR, self.params['BINDIR'], "gsetc_omp.x")
            n_cpu = multiprocessing.cpu_count()
            OMP_MAX_THREADS = 32 if int(
                n_cpu * 1.5) >= 32 else int(n_cpu * 1.5)
            self.omp_num_threads = omp_num_threads if omp_num_threads <= OMP_MAX_THREADS else OMP_MAX_THREADS
            print(
                f"Use OpenMP version of gsetc with {self.omp_num_threads} threads")
        else:
            self.omp_num_threads = 1
            self.ETC_SRC = os.path.join(self.HOME_DIR, self.params['BINDIR'], "gsetc.x")
        #        if not os.path.exists(self.HOME_DIR + '/bin'):
        #            os.mkdir(self.HOME_DIR + '/bin')
        if not os.path.exists(self.ETC_SRC):
            exit("Unable to find ETC engine; please run make first and try again")
        return None

    def set_param(self, param_name, param_value):
        """ Set ETC parameters

        Parameters
        ----------
        param_name : `str` (the parameter name)
        param_value:       (the parameter value)

        Returns
        -------

        Examples
        ----------    
        """
        if param_name in self.params.keys():
            try:
                self.params[param_name] = str(param_value)
            except:
                print('Error!')
        else:
            print('param_name %s can not be recognized ...' % (param_name))
        return 0

    def load_param_file(self, filename):
        for line in open(filename, 'r'):
            a = line.split()
            if line[0] != '#' and len(a) > 0:
                self.params[a[0]] = a[1]
        return 0

    def run(self):
        start = time.time()

        # create directories for output files
        if not os.path.exists(self.params['OUTDIR']):
            os.mkdir(self.params['OUTDIR'])
        if not os.path.exists(self.params['TMPDIR']):
            os.mkdir(self.params['TMPDIR'])

        # define spectrograph
        self.spectrograph = self.params['spectrograph'].lower()

        # select throughput model
        self.throughput_model = self.params['throughput_model']
        if self.spectrograph in ['sm1', 'sm2', 'sm3', 'sm4']:
            self.INSTR_SETUP = self.HOME_DIR + \
                '/config/PFS.%s.%s.dat' % (self.throughput_model,
                                           self.spectrograph)
            self.INSTR_SETUP_MR = self.HOME_DIR + \
                '/config/PFS.redMR.%s.%s.dat' % (
                    self.throughput_model, self.spectrograph)
        else:
            self.INSTR_SETUP = self.HOME_DIR + \
                '/config/PFS.%s.dat' % (self.throughput_model)
            self.INSTR_SETUP_MR = self.HOME_DIR + \
                '/config/PFS.redMR.%s.dat' % (self.throughput_model)

        # Noise reuse flag
        flag = '0'
        if self.params['NOISE_REUSED'].lower() == 'y':
            flag = '1'
        elif self.params['NOISE_REUSED'].lower() == 'n':
            flag = '0'
        self.NOISE_REUSED = flag

        # Medium Resolution Mode ?
        if self.params['MR_MODE'].lower() == 'yes' or self.params['MR_MODE'].lower() == 'y':
            self.INSTR_SETUP = self.INSTR_SETUP_MR
        else:
            self.INSTR_SETUP = self.INSTR_SETUP
        # make continuum magnitude file
        if os.path.exists(self.params['TMPDIR']) is False:
            os.mkdir(self.params['TMPDIR'])
        try:
            _mag = float(self.params['MAG_FILE'])
            file = open(os.path.join(
                self.params['TMPDIR'], 'mag_%s.dat' % (self.params['MAG_FILE'])), 'w')
            file.write('300.0 %.2f\n 1300. %.2f\n' % (_mag, _mag))
            file.close()
            self.mag_file = os.path.join(
                self.params['TMPDIR'], 'mag_%s.dat' % (self.params['MAG_FILE']))
        except:
            ## interpolation so that the resolution is slightly higher than that of instrument ##
            _lam, _mag = np.loadtxt(
                self.params["MAG_FILE"], usecols=(0, 1), unpack=True
            )
            lam = np.arange(300.0, 1300.0, 0.05)
            mag = np.interp(lam, _lam, _mag)
            self.mag_file = os.path.join(self.params['TMPDIR'], '%s' % (self.params['MAG_FILE'].split('/')[-1]))
            np.savetxt(self.mag_file, np.array([lam, mag]).transpose(), fmt="%.4e")
        ''' check file overwritten '''
        C = 0
        if self.params['OVERWRITE'].lower() == 'no' or self.params['OVERWRITE'].lower() == 'n':
            if os.path.exists(self.params['OUTFILE_NOISE']):
                print("Error: %s already exists... " % (args.OUTFILE_NOISE))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNC']):
                print("Error: %s already exists... " %
                      (self.params['OUTFILE_SNC']))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNL']):
                print("Error: %s already exists... " %
                      (self.params['OUTFILE_SNL']))
                C += 1
        if self.params['OVERWRITE'].lower() == 'yes' or self.params['OVERWRITE'].lower() == 'y':
            C = 0

        # apply camera obscuration depending of FoV location
        if self.params['obscFoVDep'].lower() == 'yes' or self.params['obscFoVDep'].lower() == 'y':
            obsc, corr = calc_obscuration(float(self.params['FIELD_ANG']))
            degrade = float(self.params['degrade']) * corr
            self.params['degrade'] = f'{degrade}'
        #print('The throughput degrade: %s' % (self.params['degrade']))

        # run ETC
        # print("Spectrograph setup : %s" % (self.INSTR_SETUP))
        if C != 0:
            exit('No execution of ETC')
        try:
            print('##### starting to run ETC ... (it takes a few min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE, env={
                                    "OMP_NUM_THREADS": f"{self.omp_num_threads}"})
            proc.communicate("\n".join([self.INSTR_SETUP,
                                        self.params['degrade'],
                                        SKYMODELS,
                                        self.params['SEEING'],
                                        self.params['ZENITH_ANG'],
                                        self.params['GALACTIC_EXT'],
                                        self.params['FIELD_ANG'],
                                        OFFSET_FIB,
                                        self.params['MOON_ZENITH_ANG'],
                                        self.params['MOON_TARGET_ANG'],
                                        self.params['MOON_PHASE'],
                                        self.params['EXP_TIME'],
                                        self.params['EXP_NUM'],
                                        self.params['SKY_SUB_FLOOR'],
                                        self.params['DIFFUSE_STRAY'],
                                        self.NOISE_REUSED,
                                        self.params['OUTFILE_NOISE'],
                                        self.params['OUTFILE_OII'],
                                        self.params['OUTFILE_SNL'],
                                        self.params['LINE_FLUX'],
                                        self.params['LINE_WIDTH'],
                                        self.params['OUTFILE_SNC'],
                                        self.params['INFILE_OIICat'],
                                        self.params['OUTFILE_OIICat'],
                                        self.params['minSNR'],
                                        self.mag_file,
                                        self.params['REFF']
                                        ]).encode())
        except OSError as e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)

        # load OUTFILE_NOISE
        if self.params['OUTFILE_NOISE'] != '-':
            add_header(self.params['OUTFILE_NOISE'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.nsm_arms,
                    self.nsm_pixs,
                    self.nsm_lams,
                    self.nsm_nois,
                    self.nsm_skys,
                ) = np.genfromtxt(
                    self.params["OUTFILE_NOISE"], unpack=True, usecols=(0, 1, 2, 3, 4)
                )
            except:
                print('OUTFILE_NOISE is not found ...')
                pass
        # load OUTFILE_SNC
        if self.params['OUTFILE_SNC'] != '-':
            add_header(self.params['OUTFILE_SNC'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.snc_arms,
                    self.snc_pixs,
                    self.snc_lams,
                    self.snc_sncs,
                    self.snc_sigs,
                    self.snc_nois_mobj,
                    self.snc_nois,
                    self.snc_spin,
                    self.snc_conv,
                    self.snc_samp,
                    self.snc_skys,
                ) = np.genfromtxt(
                    self.params["OUTFILE_SNC"],
                    unpack=True,
                    usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                )
            except:
                print('OUTFILE_SNC is not found ...')
                pass

        #  load OUTFILE_SNL
        if self.params['OUTFILE_SNL'] != '-':
            add_header(self.params['OUTFILE_SNL'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.snl_lams,
                    self.snl_fcov,
                    self.snl_effa,
                    self.snl_sna0,
                    self.snl_sna1,
                    self.snl_sna2,
                    self.snl_snls,
                ) = np.genfromtxt(
                    self.params["OUTFILE_SNL"],
                    unpack=True,
                    usecols=(0, 1, 2, 3, 4, 5, 6),
                )
            except:
                print('OUTFILE_SNL is not found ...')
                pass

        # load OUTFILE_OII
        if self.params['OUTFILE_OII'] != '-':
            add_header(self.params['OUTFILE_OII'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.sno2_zsps,
                    self.sno2_lam1,
                    self.sno2_lam2,
                    self.sno2_fcov,
                    self.sno2_effa,
                    self.sno2_sna0,
                    self.sno2_sna1,
                    self.sno2_sna2,
                    self.sno2_sno2,
                ) = np.genfromtxt(
                    self.params["OUTFILE_OII"],
                    unpack=True,
                    usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8),
                )
            except:
                print('OUTFILE_OII is not found ...')
                pass

        # end of the process
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (
            elapsed_time))

        return 0

    def make_noise_model(self):
        start = time.time()

        # set Noise reuse flag
        flag = '0'
        if self.params['NOISE_REUSED'].lower() == 'y':
            flag = '1'
        elif self.params['NOISE_REUSED'].lower() == 'n':
            flag = '0'
        self.NOISE_REUSED = flag

        # Medium Resolution Mode ?
        if self.params['MR_MODE'].lower() == 'yes' or self.params['MR_MODE'].lower() == 'y':
            self.INSTR_SETUP = self.INSTR_SETUP_MR
        else:
            self.INSTR_SETUP = self.INSTR_SETUP

        # make continuum magnitude file
        try:
            _mag = float(self.params['MAG_FILE'])
            if os.path.exists(self.params['TMPDIR']) is False:
                os.mkdir(self.params['TMPDIR'])
            file = open(os.path.join(
                self.params['TMPDIR'], 'mag_%s.dat' % (self.params['MAG_FILE'])), 'w')
            file.write('300.0 %.2f\n 1300. %.2f\n' % (_mag, _mag))
            file.close()
            self.mag_file = os.path.join(
                self.params['TMPDIR'], 'mag_%s.dat' % (self.params['MAG_FILE']))
        except:
            ## interpolation so that the resolution is slightly higher than that of instrument ##
            _lam, _mag = np.loadtxt(
                self.params["MAG_FILE"], usecols=(0, 1), unpack=True
            )
            lam = np.arange(300.0, 1300.0, 0.05)
            mag = np.interp(lam, _lam, _mag)
            self.mag_file = os.path.join(self.params['TMPDIR'], '%s' % (self.params['MAG_FILE'].split('/')[-1]))
            np.savetxt(self.mag_file, np.array([lam, mag]).transpose(), fmt="%.4e")
        ''' check file overwritten '''
        C = 0
        if self.params['OVERWRITE'].lower() == 'no' or self.params['OVERWRITE'].lower() == 'n':
            if os.path.exists(self.params['OUTFILE_NOISE']):
                print("Error: %s already exists... " % (args.OUTFILE_NOISE))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNC']):
                print("Error: %s already exists... " %
                      (self.params['OUTFILE_SNC']))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNL']):
                print("Error: %s already exists... " %
                      (self.params['OUTFILE_SNL']))
                C += 1
        if self.params['OVERWRITE'].lower() == 'yes' or self.params['OVERWRITE'].lower() == 'y':
            C = 0

        # apply camera obscuration depending of FoV location
        if self.params['obscFoVDep'].lower() == 'yes' or self.params['obscFoVDep'].lower() == 'y':
            obsc, corr = calc_obscuration(float(self.params['FIELD_ANG']))
            degrade = float(self.params['degrade']) * corr
            self.params['degrade'] = f'{degrade}'
        print('The throughput degrade: %s' % (self.params['degrade']))

        # run ETC
        # print("Spectrograph setup : %s" % (self.INSTR_SETUP))
        if C != 0:
            exit('No execution of ETC')
        try:
            print(
                '##### starting to make a noise model ... (it takes about 2 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE, env={
                                    "OMP_NUM_THREADS": f"{self.omp_num_threads}"})
            proc.communicate("\n".join([self.INSTR_SETUP,
                                        self.params['degrade'],
                                        SKYMODELS,
                                        self.params['SEEING'],
                                        self.params['ZENITH_ANG'],
                                        self.params['GALACTIC_EXT'],
                                        self.params['FIELD_ANG'],
                                        OFFSET_FIB,
                                        self.params['MOON_ZENITH_ANG'],
                                        self.params['MOON_TARGET_ANG'],
                                        self.params['MOON_PHASE'],
                                        self.params['EXP_TIME'],
                                        self.params['EXP_NUM'],
                                        self.params['SKY_SUB_FLOOR'],
                                        self.params['DIFFUSE_STRAY'],
                                        '0',
                                        self.params['OUTFILE_NOISE'],
                                        '-',
                                        '-',
                                        self.params['LINE_FLUX'],
                                        self.params['LINE_WIDTH'],
                                        '-',
                                        '-',
                                        '-',
                                        self.params['minSNR'],
                                        self.mag_file,
                                        self.params['REFF']
                                        ]).encode())
        except OSError as e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)

        # load OUTFILE_NOISE
        if self.params['OUTFILE_NOISE'] != '-':
            add_header(self.params['OUTFILE_NOISE'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.nsm_arms,
                    self.nsm_pixs,
                    self.nsm_lams,
                    self.nsm_nois,
                    self.nsm_skys,
                ) = np.genfromtxt(
                    self.params["OUTFILE_NOISE"], unpack=True, usecols=(0, 1, 2, 3, 4)
                )
            except:
                print('OUTFILE_NOISE is not found ...')
                pass

        # end of process
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (
            elapsed_time))
        return 0

    def proc_multi(self, inputs):
        self.params[inputs[0]] = str(inputs[1])
        for outFileName in ['OUTFILE_NOISE', 'OUTFILE_SNC', 'OUTFILE_SNL', 'OUTFILE_OII']:
            if self.params[outFileName] != "-":
                self.params[outFileName] += '.'.join(
                    ['', inputs[0], inputs[1]])
        self.run()
        return 0

    def run_multi(self, nproc, param_name, param_values):
        from multiprocessing import Pool
        p = Pool(nproc)
        result = p.map(self.proc_multi, [(param_name, v)
                       for v in param_values])
        return 0

    def get_noise(self):
        return self.nsm_lams, self.nsm_nois

    def make_snc(self):
        # apply camera obscuration depending of FoV location
        if self.params['obscFoVDep'].lower() == 'yes' or self.params['obscFoVDep'].lower() == 'y':
            obsc, corr = calc_obscuration(float(self.params['FIELD_ANG']))
            degrade = float(self.params['degrade']) * corr
            self.params['degrade'] = f'{degrade}'
        print('The throughput degrade: %s' % (self.params['degrade']))

        # run ETC
        start = time.time()
        try:
            print('##### starting to make an SNC model ... (it takes about 1 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE, env={
                                    "OMP_NUM_THREADS": f"{self.omp_num_threads}"})
            proc.communicate("\n".join([self.INSTR_SETUP,
                                        self.params['degrade'],
                                        SKYMODELS,
                                        self.params['SEEING'],
                                        self.params['ZENITH_ANG'],
                                        self.params['GALACTIC_EXT'],
                                        self.params['FIELD_ANG'],
                                        OFFSET_FIB,
                                        self.params['MOON_ZENITH_ANG'],
                                        self.params['MOON_TARGET_ANG'],
                                        self.params['MOON_PHASE'],
                                        self.params['EXP_TIME'],
                                        self.params['EXP_NUM'],
                                        self.params['SKY_SUB_FLOOR'],
                                        self.params['DIFFUSE_STRAY'],
                                        '1',
                                        self.params['OUTFILE_NOISE'],
                                        '-',
                                        '-',
                                        self.params['LINE_FLUX'],
                                        self.params['LINE_WIDTH'],
                                        self.params['OUTFILE_SNC'],
                                        '-',
                                        '-',
                                        self.params['minSNR'],
                                        self.mag_file,
                                        self.params['REFF']
                                        ]).encode())
        except OSError as e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)
        # load OUTFILE_SNC
        if self.params['OUTFILE_SNC'] != '-':
            add_header(self.params['OUTFILE_SNC'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.snc_arms,
                    self.snc_pixs,
                    self.snc_lams,
                    self.snc_sncs,
                    self.snc_sigs,
                    self.snc_nois_mobj,
                    self.snc_nois,
                    self.snc_spin,
                    self.snc_conv,
                    self.snc_samp,
                    self.snc_skys,
                ) = np.genfromtxt(
                    self.params["OUTFILE_SNC"],
                    unpack=True,
                    usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                )
            except:
                print('OUTFILE_SNC is not found ...')
                pass
        # end of process
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (
            elapsed_time))

        return 0

    def get_snc(self):
        return self.snc_lams, self.snc_sncs

    def make_snl(self):
        # apply camera obscuration depending of FoV location
        if self.params['obscFoVDep'].lower() == 'yes' or self.params['obscFoVDep'].lower() == 'y':
            obsc, corr = calc_obscuration(float(self.params['FIELD_ANG']))
            degrade = float(self.params['degrade']) * corr
            self.params['degrade'] = f'{degrade}'
        print('The throughput degrade: %s' % (self.params['degrade']))

        # run ETC
        start = time.time()
        try:
            print('##### starting to make an SNL model ... (it takes about 1 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE, env={
                                    "OMP_NUM_THREADS": f"{self.omp_num_threads}"})
            proc.communicate("\n".join([self.INSTR_SETUP,
                                        self.params['degrade'],
                                        SKYMODELS,
                                        self.params['SEEING'],
                                        self.params['ZENITH_ANG'],
                                        self.params['GALACTIC_EXT'],
                                        self.params['FIELD_ANG'],
                                        OFFSET_FIB,
                                        self.params['MOON_ZENITH_ANG'],
                                        self.params['MOON_TARGET_ANG'],
                                        self.params['MOON_PHASE'],
                                        self.params['EXP_TIME'],
                                        self.params['EXP_NUM'],
                                        self.params['SKY_SUB_FLOOR'],
                                        self.params['DIFFUSE_STRAY'],
                                        '1',
                                        self.params['OUTFILE_NOISE'],
                                        '-',
                                        self.params['OUTFILE_SNL'],
                                        self.params['LINE_FLUX'],
                                        self.params['LINE_WIDTH'],
                                        '-',
                                        '-',
                                        '-',
                                        self.params['minSNR'],
                                        self.mag_file,
                                        self.params['REFF']
                                        ]).encode())
        except OSError as e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)

        # load OUTFILE_SNL
        if self.params['OUTFILE_SNL'] != '-':
            add_header(self.params['OUTFILE_SNL'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.snl_lams,
                    self.snl_fcov,
                    self.snl_effa,
                    self.snl_sna0,
                    self.snl_sna1,
                    self.snl_sna2,
                    self.snl_snls,
                ) = np.genfromtxt(
                    self.params["OUTFILE_SNL"],
                    unpack=True,
                    usecols=(0, 1, 2, 3, 4, 5, 6),
                )
            except:
                print('OUTFILE_SNL is not found ...')
                pass

        # end of process
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (
            elapsed_time))

        return 0

    def get_snl(self):
        return self.snl_lams, self.snl_snls

    def make_sno2(self):
        # apply camera obscuration depending of FoV location
        if self.params['obscFoVDep'].lower() == 'yes' or self.params['obscFoVDep'].lower() == 'y':
            obsc, corr = calc_obscuration(float(self.params['FIELD_ANG']))
            degrade = float(self.params['degrade']) * corr
            self.params['degrade'] = f'{degrade}'
        print('The throughput degrade: %s' % (self.params['degrade']))

        # run ETC
        start = time.time()
        try:
            print('##### starting to make an OII model ... (it takes about 2 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE, env={
                                    "OMP_NUM_THREADS": f"{self.omp_num_threads}"})
            proc.communicate("\n".join([self.INSTR_SETUP,
                                        self.params['degrade'],
                                        SKYMODELS,
                                        self.params['SEEING'],
                                        self.params['ZENITH_ANG'],
                                        self.params['GALACTIC_EXT'],
                                        self.params['FIELD_ANG'],
                                        OFFSET_FIB,
                                        self.params['MOON_ZENITH_ANG'],
                                        self.params['MOON_TARGET_ANG'],
                                        self.params['MOON_PHASE'],
                                        self.params['EXP_TIME'],
                                        self.params['EXP_NUM'],
                                        self.params['SKY_SUB_FLOOR'],
                                        self.params['DIFFUSE_STRAY'],
                                        '1',
                                        self.params['OUTFILE_NOISE'],
                                        self.params['OUTFILE_OII'],
                                        '-',
                                        self.params['LINE_FLUX'],
                                        self.params['LINE_WIDTH'],
                                        '-',
                                        '-',
                                        '-',
                                        self.params['minSNR'],
                                        self.mag_file,
                                        self.params['REFF'],
                                        ]).encode())
        except OSError as e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)

        # load OUTFILE_OII
        if self.params['OUTFILE_OII'] != '-':
            add_header(self.params['OUTFILE_OII'], self.params['SEEING'], self.params['ZENITH_ANG'], self.params['GALACTIC_EXT'], self.params['MOON_ZENITH_ANG'], self.params['MOON_TARGET_ANG'],
                       self.params['MOON_PHASE'], self.params['EXP_TIME'], self.params['EXP_NUM'], self.params['FIELD_ANG'], self.params['MAG_FILE'], self.params['REFF'], self.params['LINE_FLUX'], self.params['LINE_WIDTH'])
            try:
                (
                    self.sno2_zsps,
                    self.sno2_lam1,
                    self.sno2_lam2,
                    self.sno2_fcov,
                    self.sno2_effa,
                    self.sno2_sna0,
                    self.sno2_sna1,
                    self.sno2_sna2,
                    self.sno2_sno2,
                ) = np.genfromtxt(
                    self.params["OUTFILE_OII"],
                    unpack=True,
                    usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8),
                )
            except:
                print('OUTFILE_OII is not found ...')
                pass

        # end of process
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (
            elapsed_time))

        return 0

    def get_sno2(self):
        return self.sno2_zsps, self.sno2_sno2
