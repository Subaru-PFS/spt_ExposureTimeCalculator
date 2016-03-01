#!/usr/bin/env python
import math,os,sys
import argparse
import random as rn
import time

try:
    import pyfits as pf
    PYFITS_FLG = 1
except:
    PYFITS_FLG = 0
try:
    import numpy as np
    NUMPY_FLG = 1
except:
    NUMPY_FLG = 0
# print PYFITS_FLG,NUMPY_FLG
#######################
offset = 0.01
R0 = 2300
R1 = 3000
R2 = 4300
R3 = 5000
#######################

if __name__ == '__main__':
    start = time.time()
    param_name = []
    param_value = {}
    parser = argparse.ArgumentParser(description='PFS Spectral Simulator developed by Kiyoto Yabe')
    parser.add_argument("params", type=str, help="parameter file")
    parser.add_argument("--SEEING", type=str, help="seeing")
    parser.add_argument("--ZENITH_ANG", type=str, help="zenith angle")
    parser.add_argument("--GALACTIC_EXT", type=str, help="galactic extinction")
    parser.add_argument("--MOON_ZENITH_ANG", type=str, help="moon-zenith angle")
    parser.add_argument("--MOON_TARGET_ANG", type=str, help="moon-target angle")
    parser.add_argument("--MOON_PHASE", type=str, help="moon phase")
    parser.add_argument("--EXP_TIME", type=str, help="exposure time")
    parser.add_argument("--EXP_NUM", type=str, help="exposure number")
    parser.add_argument("--FIELD_ANG", type=str, help="field angle")
    parser.add_argument("--MAG_FILE", type=str, help="magnitude input file")
    parser.add_argument("--REFF", type=str, help="effective radius")
    parser.add_argument("--LINE_FLUX", type=str, help="emission line flux")
    parser.add_argument("--LINE_WIDTH", type=str, help="emission line width (sigma)")
    parser.add_argument("--NOISE_REUSED", type=str, help="noise vector reused flag")
    parser.add_argument("--OUTFILE_NOISE", type=str, help="noise vector output file")
    parser.add_argument("--OUTFILE_SNC", type=str, help="continuum results output file")
    parser.add_argument("--OUTFILE_SNL", type=str, help="emission line results output file")
    parser.add_argument("--OUTFILE_OII", type=str, help="[OII] emission line results output file")
    parser.add_argument("--MR_MODE", type=str, help="medium resolution mode on-off")
    parser.add_argument("--OVERWRITE", type=str, help="overwrite on-off")
    parser.add_argument("-s","--show", help="Show parameter set")
    parser.add_argument("--INFILE_SNC", type=str, help="continuum results input file")
    parser.add_argument("--OUTFILE_SIM", type=str, help="simulated spectrum output file")
#    parser.add_argument("--FILE_TYPE_SIM", type=str, help="simulated spectrum output file type")
    args = parser.parse_args()
    ## read parameter file ##   
    if os.path.exists(args.params):
        for line in open(args.params,'r'):
            line.replace('\r','\n')
            if line[0]!="#" and line!="\n":
                a=line.split()
                if len(a)>0:
                    param_name.append(a[0])
                    if vars(args)[a[0]] is None:
                        param_value[a[0]] = a[1]
                    else:
                        param_value[a[0]] = vars(args)[a[0]]
    ## read input file ##
    arm_list = []
    arm = []
    wav = []
    dsp = []
    mag = []
    snc = []
    for line in open(param_value['INFILE_SNC'],'r'):
        a=line.split()
        if a[0][0] != "#":
            if int(a[0]) not in arm_list:
                arm_list.append(int(a[0]))
            arm.append(int(a[0]))
            wav.append(float(a[2]))
            mag.append(float(a[6]))
            snc.append(float(a[3]))
            if param_value['MR_MODE'].lower() == 'yes':
                if int(a[0]) == 0:
                    dsp.append(float(a[2])/R0)
                elif int(a[0]) == 1:
                    dsp.append(float(a[2])/R3)
                else:
                    dsp.append(float(a[2])/R2)
            else:
                if int(a[0]) == 0:
                    dsp.append(float(a[2])/R0)
                elif int(a[0]) == 1:
                    dsp.append(float(a[2])/R1)
                else:
                    dsp.append(float(a[2])/R2)
## data output ##
## FITS mode ##
    if PYFITS_FLG == 1 and NUMPY_FLG == 1:
        for arm_num in arm_list:
            wav_a = []
            dsp_a = []
            flx = []
            var  = []
            for i in range(len(wav)):
                if arm[i] == arm_num:
                    fnu   = 10**(-0.4*(mag[i]+48.6))
                    flam  = 3.0e18*fnu/(10*wav[i])**2/1e-17
                    sigma = flam/(snc[i]+offset)
                    wav_a.append(wav[i])
                    dsp_a.append(dsp[i])
                    flx.append(flam+rn.gauss(0.0,sigma))
                    var.append(sigma**2)
            Npix = len(wav_a)
            msk = np.zeros(len(wav_a))
            sky = np.zeros(len(wav_a))
            xps = np.zeros(len(wav_a))
            cfg = 'TBD'
            col1 = pf.Column(name = 'FLUX', format = '%dE' % (Npix), array = np.array([np.array(flx)]))
            col2 = pf.Column(name = 'COVAR', format = '%dE' % (Npix), array = np.array([np.array(var)]))
            col3 = pf.Column(name = 'COVAR2', format = '3A', array = np.array(['N/A']))
            col4 = pf.Column(name = 'MASK', format = '%dJ' % (Npix), array = np.array([msk]))
            col5 = pf.Column(name = 'WAVEDISP', format = '%dE' % (Npix), array = np.array([np.array(dsp_a)]))
            col6 = pf.Column(name = 'SKY', format = '%dE' % (Npix), array = np.array([sky]))
            col7 = pf.Column(name = 'CONFIG', format = '3A', array = np.array(['TBD']))
            tbhdu1 = pf.BinTableHDU.from_columns([col1])
            tbhdu2 = pf.BinTableHDU.from_columns([col2])
            tbhdu3 = pf.BinTableHDU.from_columns([col3])
            tbhdu4 = pf.BinTableHDU.from_columns([col4])
            tbhdu5 = pf.BinTableHDU.from_columns([col5])
            tbhdu6 = pf.BinTableHDU.from_columns([col6])
            tbhdu7 = pf.BinTableHDU.from_columns([col7])
            hdu = pf.PrimaryHDU()
            hdulist = pf.HDUList([hdu])
            hdulist.append(tbhdu1)
            hdulist.append(tbhdu2)
            hdulist.append(tbhdu3)
            hdulist.append(tbhdu4)
            hdulist.append(tbhdu5)
            hdulist.append(tbhdu6)
            hdulist.append(tbhdu7)
            hdulist[1].header['TUNIT1'] = '1e^-17 erg/s/cm^2/A'
            hdulist[1].header['EXTNAME'] = 'FLUX'
            hdulist[2].header['TUNIT1'] = '(1e^-17 erg/s/cm^2/A)^2'
            hdulist[2].header['EXTNAME'] = 'COVAR'
            hdulist[3].header['TUNIT1'] = ''
            hdulist[3].header['EXTNAME'] = 'COVAR2'
            hdulist[4].header['TUNIT1'] = ''
            hdulist[4].header['EXTNAME'] = 'MASK'
            hdulist[5].header['TUNIT1'] = 'nm'
            hdulist[5].header['EXTNAME'] = 'WAVEDISP'
            hdulist[6].header['TUNIT1'] = '1e^-17 erg/s/cm^2/A'
            hdulist[6].header['EXTNAME'] = 'SKY'
            hdulist[7].header['TUNIT1'] = ''
            hdulist[7].header['EXTNAME'] = 'CONFIG'
            for i in range(7):
                hdulist[i+1].header['CRPIX'] = 1
                hdulist[i+1].header['CRVAL'] = min(wav_a)
                hdulist[i+1].header['CDELT'] = (max(wav_a)-min(wav_a))/(len(wav_a)*1.0)
            os.system('rm -f %s' % (param_value['OUTFILE_SIM']+'_arm%1d.fits' % (arm_num)))
            hdulist.writeto(param_value['OUTFILE_SIM']+'_arm%1d.fits' % (arm_num))
        print "FITS table generated"
    else:
        print "Error: Install pyfits and numpy to use this mode"    
## ASCII mode ##
    file=open(param_value['OUTFILE_SIM']+'.dat','w')
    file.write('''#  1  WAVELENGTH  [nm]
#  2  FLUX        [10^-17 erg/s/cm^2/A]
#  3  ERROR       [10^-17 erg/s/cm^2/A]
#  4  MASK        
#  5  WAVEDISP    [nm]^2
#  6  SKY         [10^-17 erg/s/cm^2/A]
''')
    for i in range(len(wav)):
        fnu   = 10**(-0.4*(mag[i]+48.6))
        flam  = 3.0e18*fnu/(10*wav[i])**2/1e-17
        sigma = flam/(snc[i]+offset)
        flx = flam+rn.gauss(0.0,sigma)
        var = sigma**2
        msk = 0.0
        sky = 0.0
        xps = 0.0
        cfg = 0.0
        file.write('%8.3f %12.4e %12.4e %2d %6.3f %12.4e\n' % (wav[i],flx,sigma,msk,dsp[i],sky))
    file.close()
    print "ASCII table generated"

