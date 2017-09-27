# -*- coding: utf-8 -*-

from __future__ import print_function, division

import os
from os import path
import scipy as sp
import time
import subprocess

### Some parameters ###########################
### CAUTION: ##################################
### CHANGE BELOW ON YOUR OWN RESPONSIBILITY ###
###############################################
SKYMODELS = '11005'
OFFSET_FIB = '0.00'
SKY_SUB_FLOOR = '0.01'
DIFFUSE_STRAY = '0.02'
offset = 0.01


class Etc(object):

    def __init__(self):
        self.params = {'SEEING': '0.80',
                       'ZENITH_ANG': '45.00',
                       'GALACTIC_EXT': '0.00',
                       'MOON_ZENITH_ANG': '30.0',
                       'MOON_TARGET_ANG': '60.0',
                       'MOON_PHASE': '0',
                       'EXP_TIME': '450',
                       'EXP_NUM': '8',
                       'FIELD_ANG': '0.675',
                       'MAG_FILE': '22.5',
                       'REFF': '0.3',
                       'LINE_FLUX': '1.0e-17',
                       'LINE_WIDTH': '70',
                       'NOISE_REUSED': 'N',
                       'OUTFILE_NOISE': 'out/ref.noise.dat',
                       'OUTFILE_SNC': 'out/ref.snc.dat',
                       'OUTFILE_SNL': 'out/ref.snl.dat',
                       'OUTFILE_OII': '-',
                       'MR_MODE': 'N',
                       'OVERWRITE': 'Y',
                       'INFILE_OIICat': '-',
                       'OUTFILE_OIICat': '-',
                       'minSNR': '9.0'
                       }
        self.HOME_DIR = path.dirname(path.abspath(__file__))
        self.ETC_SRC = self.HOME_DIR + '/bin/gsetc.x'
        self.INSTR_SETUP = self.HOME_DIR + '/config/PFS.dat'
        self.INSTR_SETUP_MR = self.HOME_DIR + '/config/PFS.redMR.dat'
#        if not os.path.exists(self.HOME_DIR + '/bin'):
#            os.mkdir(self.HOME_DIR + '/bin')
        if not os.path.exists('out'):
            os.mkdir('out')
        if not os.path.exists(self.ETC_SRC):
            exit("Unable to find ETC engine; please run make first and try again")
        return None

    def set_param(self, param_name, param_value):
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
        ''' Noise reuse flag '''
        flag = '0'
        if self.params['NOISE_REUSED'].lower() == 'y':
            flag = '1'
        elif self.params['NOISE_REUSED'].lower() == 'n':
            flag = '0'
        self.NOISE_REUSED = flag
        ''' Medium Resolution Mode ? '''
        if self.params['MR_MODE'].lower() == 'yes' or self.params['MR_MODE'].lower() == 'y':
            self.INSTR_SETUP = self.INSTR_SETUP_MR
        else:
            self.INSTR_SETUP = self.INSTR_SETUP
        ''' make continuum magnitude file '''
        try:
            _mag = float(self.params['MAG_FILE'])
            if os.path.exists('tmp') == False:
                os.mkdir('tmp')
            file = open('tmp/mag_%s.dat' % (self.params['MAG_FILE']), 'w')
            file.write('300.0 %.2f\n 1300. %.2f\n' % (_mag, _mag))
            file.close()
            self.mag_file = 'tmp/mag_%s.dat' % (self.params['MAG_FILE'])
        except:
            ## interpolation so that the resolution is slightly higher than that of instrument ##
            _lam, _mag = sp.loadtxt(self.params['MAG_FILE'], usecols=(0, 1), unpack=True)
            lam = sp.arange(300., 1300., 0.05)
            mag = sp.interp(lam, _lam, _mag)
            self.mag_file = 'tmp/%s' % (self.params['MAG_FILE'].split('/')[-1])
            sp.savetxt(self.mag_file, sp.array([lam, mag]).transpose(), fmt='%.4e')
        ''' check file overwritten '''
        C = 0
        if self.params['OVERWRITE'].lower() == 'no' or self.params['OVERWRITE'].lower() == 'n':
            if os.path.exists(self.params['OUTFILE_NOISE']):
                print("Error: %s already exists... " % (args.OUTFILE_NOISE))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNC']):
                print("Error: %s already exists... " % (self.params['OUTFILE_SNC']))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNL']):
                print("Error: %s already exists... " % (self.params['OUTFILE_SNL']))
                C += 1
        if self.params['OVERWRITE'].lower() == 'yes' or self.params['OVERWRITE'].lower() == 'y':
            C = 0
        ''' run ETC '''
        #print("Spectrograph setup : %s" % (self.INSTR_SETUP))
        if C != 0:
            exit('No execution of ETC')
        try:
            print('##### starting to run ETC ... (it takes a few min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE)
            proc.communicate("\n".join([self.INSTR_SETUP,
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
                                        SKY_SUB_FLOOR,
                                        DIFFUSE_STRAY,
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
                                        ]))
        except OSError, e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)
        ''' load OUTFILE_NOISE '''
        try:
            self.nsm_arms, self.nsm_pixs, self.nsm_lams, self.nsm_nois, self.nsm_skys = sp.genfromtxt(self.params['OUTFILE_NOISE'], unpack=True, usecols=(0, 1, 2, 3, 4))
        except:
            exit('OUTFILE_NOISE is not found ...')
        ''' load OUTFILE_SNC '''
        try:
            self.snc_arms, self.snc_pixs, self.snc_lams, self.snc_sncs, self.snc_sigs, self.snc_nois_mobj, self.snc_nois, self.snc_spin, self.snc_conv, self.snc_samp, self.snc_skys = sp.genfromtxt(self.params['OUTFILE_SNC'], unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
        except:
            exit('OUTFILE_SNC is not found ...')
        ''' load OUTFILE_SNL '''
        try:
            self.snl_lams, self.snl_fcov, self.snl_effa, self.snl_sna0, self.snl_sna1, self.snl_sna2, self.snl_snls = sp.genfromtxt(self.params['OUTFILE_SNL'], unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6))
        except:
            exit('OUTFILE_SNL is not found ...')
        ''' load OUTFILE_OII '''
        try:
            self.sno2_zsps, self.sno2_lam1, self.sno2_lam2, self.sno2_fcov, self.sno2_effa, self.sno2_sna0, self.sno2_sna1, self.sno2_sna2, self.sno2_sno2 = sp.genfromtxt(self.params['OUTFILE_OII'], unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8))
        except:
            exit('OUTFILE_OII is not found ...')

        ''' end of the process '''
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (elapsed_time))

        return 0

    def make_noise_model(self):
        start = time.time()
        ''' Noise reuse flag '''
        flag = '0'
        if self.params['NOISE_REUSED'].lower() == 'y':
            flag = '1'
        elif self.params['NOISE_REUSED'].lower() == 'n':
            flag = '0'
        self.NOISE_REUSED = flag
        ''' Medium Resolution Mode ? '''
        if self.params['MR_MODE'].lower() == 'yes' or self.params['MR_MODE'].lower() == 'y':
            self.INSTR_SETUP = self.INSTR_SETUP_MR
        else:
            self.INSTR_SETUP = self.INSTR_SETUP
        ''' make continuum magnitude file '''
        try:
            _mag = float(self.params['MAG_FILE'])
            if os.path.exists('tmp') == False:
                os.mkdir('tmp')
            file = open('tmp/mag_%s.dat' % (self.params['MAG_FILE']), 'w')
            file.write('300.0 %.2f\n 1300. %.2f\n' % (_mag, _mag))
            file.close()
            self.mag_file = 'tmp/mag_%s.dat' % (self.params['MAG_FILE'])
        except:
            ## interpolation so that the resolution is slightly higher than that of instrument ##
            _lam, _mag = sp.loadtxt(self.params['MAG_FILE'], usecols=(0, 1), unpack=True)
            lam = sp.arange(300., 1300., 0.05)
            mag = sp.interp(lam, _lam, _mag)
            self.mag_file = 'tmp/%s' % (self.params['MAG_FILE'].split('/')[-1])
            sp.savetxt(self.mag_file, sp.array([lam, mag]).transpose(), fmt='%.4e')
        ''' check file overwritten '''
        C = 0
        if self.params['OVERWRITE'].lower() == 'no' or self.params['OVERWRITE'].lower() == 'n':
            if os.path.exists(self.params['OUTFILE_NOISE']):
                print("Error: %s already exists... " % (args.OUTFILE_NOISE))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNC']):
                print("Error: %s already exists... " % (self.params['OUTFILE_SNC']))
                C += 1
            if os.path.exists(self.params['OUTFILE_SNL']):
                print("Error: %s already exists... " % (self.params['OUTFILE_SNL']))
                C += 1
        if self.params['OVERWRITE'].lower() == 'yes' or self.params['OVERWRITE'].lower() == 'y':
            C = 0
        ''' run ETC '''
        #print("Spectrograph setup : %s" % (self.INSTR_SETUP))
        if C != 0:
            exit('No execution of ETC')
        try:
            print('##### starting to make a noise model ... (it takes about 2 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE)
            proc.communicate("\n".join([self.INSTR_SETUP,
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
                                        SKY_SUB_FLOOR,
                                        DIFFUSE_STRAY,
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
                                        ]))
        except OSError, e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)
        ''' load OUTFILE_NOISE '''
        try:
            self.nsm_arms, self.nsm_pixs, self.nsm_lams, self.nsm_nois, self.nsm_skys = sp.genfromtxt(self.params['OUTFILE_NOISE'], unpack=True, usecols=(0, 1, 2, 3, 4))
        except:
            exit('OUTFILE_NOISE is not found ...')
        ''' end of process '''
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (elapsed_time))

        return 0

    def get_noise(self):
        return self.nsm_lams, self.nsm_nois

    def make_snc(self):
        ''' run ETC '''
        start = time.time()
        try:
            print('##### starting to make an SNC model ... (it takes about 1 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE)
            proc.communicate("\n".join([self.INSTR_SETUP,
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
                                        SKY_SUB_FLOOR,
                                        DIFFUSE_STRAY,
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
                                        ]))
        except OSError, e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)
        ''' load OUTFILE_SNC '''
        try:
            self.snc_arms, self.snc_pixs, self.snc_lams, self.snc_sncs, self.snc_sigs, self.snc_nois_mobj, self.snc_nois, self.snc_spin, self.snc_conv, self.snc_samp, self.snc_skys = sp.genfromtxt(self.params['OUTFILE_SNC'], unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
        except:
            exit('OUTFILE_SNC is not found ...')
        ''' end of process '''
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (elapsed_time))

        return 0

    def get_snc(self):
        return self.snc_lams, self.snc_sncs

    def make_snl(self):
        ''' run ETC '''
        start = time.time()
        try:
            print('##### starting to make an SNL model ... (it takes about 1 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE)
            proc.communicate("\n".join([self.INSTR_SETUP,
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
                                        SKY_SUB_FLOOR,
                                        DIFFUSE_STRAY,
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
                                        ]))
        except OSError, e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)
        ''' load OUTFILE_SNL '''
        try:
            self.snl_lams, self.snl_fcov, self.snl_effa, self.snl_sna0, self.snl_sna1, self.snl_sna2, self.snl_snls = sp.genfromtxt(self.params['OUTFILE_SNL'], unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6))
        except:
            exit('OUTFILE_SNL is not found ...')
        ''' end of process '''
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (elapsed_time))

        return 0

    def get_snl(self):
        return self.snl_lams, self.snl_snls

    def make_sno2(self):
        ''' run ETC '''
        start = time.time()
        try:
            print('##### starting to make an OII model ... (it takes about 2 min.) #####')
            proc = subprocess.Popen([self.ETC_SRC], stdin=subprocess.PIPE)
            proc.communicate("\n".join([self.INSTR_SETUP,
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
                                        SKY_SUB_FLOOR,
                                        DIFFUSE_STRAY,
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
                                        self.params['REFF']
                                        ]))
        except OSError, e:
            exit('Execution error of "%s" (%s)' % self.ETC_SRC, e)
        ''' load OUTFILE_OII '''
        try:
            self.sno2_zsps, self.sno2_lam1, self.sno2_lam2, self.sno2_fcov, self.sno2_effa, self.sno2_sna0, self.sno2_sna1, self.sno2_sna2, self.sno2_sno2 = sp.genfromtxt(self.params['OUTFILE_OII'], unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8))
        except:
            exit('OUTFILE_OII is not found ...')
        ''' end of process '''
        elapsed_time = time.time() - start
        print("##### finished (elapsed_time: %.1f[sec]) #####" % (elapsed_time))

        return 0

    def get_sno2(self):
        return self.sno2_zsps, self.sno2_sno2
