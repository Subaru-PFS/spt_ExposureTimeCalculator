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
    parser.add_argument("--OUTFILE_SIM", type=str, help="simulated spectrum output ASCII file")
    parser.add_argument("--OUTFILE_TRACT", type=str, help="tract")
    parser.add_argument("--OUTFILE_PATCH", type=str, help="patch")
    parser.add_argument("--OUTFILE_CATID", type=str, help="catalogue id")
    parser.add_argument("--OUTFILE_OBJID", type=str, help="object id in hexadecimal")
    parser.add_argument("--OUTFILE_NVISIT", type=str, help="visit number")
    parser.add_argument("--OUTFILE_CONFIG", type=str, help="pfsConfigId in hexadecimal")
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
## ASCII mode ##
    file=open(param_value['OUTFILE_SIM'],'w')
    file.write('''#  1  WAVELENGTH  [nm]
#  2  FLUX        [10^-17 erg/s/cm^2/A]
#  3  ERROR       [10^-17 erg/s/cm^2/A]
#  4  MASK        [1=masked]
#  5  SKY         [10^-17 erg/s/cm^2/A]
#  6  ARM         [0=blue,1=red,2=NIR,3=redMR]
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
        file.write('%8.3f %12.4e %12.4e %2d %12.4e %1d\n' % (wav[i],flx,sigma,msk,sky,arm[i]))
    file.close()
    print "ASCII table was generated"
## FITS mode ##
    if PYFITS_FLG == 1 and NUMPY_FLG == 1:
        for arm_num in arm_list:
            tract = int(param_value['OUTFILE_TRACT'])
            patch = param_value['OUTFILE_PATCH']
            catid = int(param_value['OUTFILE_CATID'])
            objid = param_value['OUTFILE_OBJID']
            nvisit = int(param_value['OUTFILE_NVISIT'])
            config = param_value['OUTFILE_CONFIG']
            output_file = './out/pfsObject-%04d-%3s-%03d-%8s-%02d-0x%8s_arm%1d.fits' % (tract,patch,catid,objid,nvisit,config,arm_num)
            config_file_name = 'pfsConfig-0x%8s.fits' % (config)
            config_file = './out/%s' % (config_file_name)
            wav_a = []
            dsp_a = []
            flx = []
            err = []
            var  = []
            for i in range(len(wav)):
                if arm[i] == arm_num:
                    fnu   = 10**(-0.4*(mag[i]+48.6))
                    flam  = 3.0e18*fnu/(10*wav[i])**2/1e-17
                    sigma = flam/(snc[i]+offset)
                    wav_a.append(wav[i])
                    dsp_a.append(dsp[i])
                    flx.append(flam+rn.gauss(0.0,sigma))
                    err.append(sigma)
                    var.append(sigma**2)
            Npix = len(wav_a)
            msk = np.zeros(Npix)
            sky = np.zeros(Npix)
            xps = np.zeros(Npix)
            cfg = 'TBD'
            col1 = np.array(flx)
            col2w = pf.Column(name = 'WAVELENGTH', format = 'E', array = np.array(wav))
            col2f = pf.Column(name = 'FLUX', format = 'E', array = np.array(flx))
            col2e = pf.Column(name = 'ERROR', format = 'E', array = np.array(err))
            col2m = pf.Column(name = 'MASK', format = 'J', array = msk)
            col3 = np.array([np.zeros(Npix),np.zeros(Npix),np.array(var),np.zeros(Npix),np.zeros(Npix)])
            col4 = np.zeros((10,10))
            col5 = np.array(msk)
            col6 = np.array(sky)
            col71 = pf.Column(name = 'pfsConfigId', format = '10A', array = np.array(['0x%s' % (config)]))
            col72 = pf.Column(name = 'visit', format = 'J', array = np.array([nvisit]))
            hdu1 = pf.ImageHDU(data = col1, header=None)
            hdu2 = pf.BinTableHDU.from_columns([col2w,col2f,col2e,col2m])
            hdu3 = pf.ImageHDU(data = col3, header=None)
            hdu4 = pf.ImageHDU(data = col4, header=None)
            hdu5 = pf.ImageHDU(data = col5, header=None)
            hdu6 = pf.ImageHDU(data = col6, header=None)
            hdu7 = pf.BinTableHDU.from_columns([col71,col72])
            hdu = pf.PrimaryHDU()
            hdulist = pf.HDUList([hdu])
            hdulist.append(hdu1)
            hdulist.append(hdu2)
            hdulist.append(hdu3)
            hdulist.append(hdu4)
            hdulist.append(hdu5)
            hdulist.append(hdu6)
            hdulist.append(hdu7)
            hdulist[1].header['EXTNAME'] = 'FLUX'
            hdulist[2].header['EXTNAME'] = 'FLUXTBL'
            hdulist[3].header['EXTNAME'] = 'COVAR'
            hdulist[4].header['EXTNAME'] = 'COVAR2'
            hdulist[5].header['EXTNAME'] = 'MASK'
            hdulist[6].header['EXTNAME'] = 'SKY'
            hdulist[7].header['EXTNAME'] = 'CONFIG'
            for i in range(6):
                hdulist[i+1].header['CRPIX'] = 1
                hdulist[i+1].header['CRVAL'] = min(wav_a)
                hdulist[i+1].header['CDELT'] = (max(wav_a)-min(wav_a))/(len(wav_a)*1.0)
            os.system('rm -f %s' % (output_file))
            hdulist.writeto(output_file)
## pfsConfig ##
# catId, objId, ra, dec, fiber flux, MPS centroid
        col1 = pf.Column(name = 'catId', format = 'J', array = np.array([catid]))
        col2 = pf.Column(name = 'objId', format = '8A', array = np.array([objid]))
        col3 = pf.Column(name = 'ra', format = 'E', array = np.array([0.0]))
        col4 = pf.Column(name = 'dec', format = 'E', array = np.array([0.0]))
        try:
            mag = float(param_value['MAG_FILE'])
            fnu = 10 ** (-0.4 * (mag + 48.6))
        except:
            mag = 99.0000
            fnu = 0.0
        col5 = pf.Column(name = 'fiber flux', format = 'E', array = np.array([fnu]))
        col6 = pf.Column(name = 'MPS centroid', format = 'E', array = np.array([0.0]))
        hdu1 = pf.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6])
        hdu = pf.PrimaryHDU()
        hdulist = pf.HDUList([hdu,hdu1])
        hdulist[1].header['EXTNAME'] = 'CONFIG'
        os.system('rm -f %s' % (config_file))
        hdulist.writeto(config_file)
        print "FITS file was generated"
    else:
        print "FITS file was NOT generated"
        print "Please install pyfits and numpy to generate FITS format"    

