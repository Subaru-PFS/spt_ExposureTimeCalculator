#!/usr/bin/env python
import math,os,sys
import argparse
import random as rn
import time

try:
    import pyfits as pf
    PYFITS_FLG=1
except:
    PYFITS_FLG=0
try:
    import numpy as np
    NUMPY_FLG=1
except:
    NUMPY_FLG=0
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
    parser.add_argument("--FILE_TYPE_SIM", type=str, help="simulated spectrum output file type")
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
    arm = []
    wav = []
    dsp = []
    mag = []
    snc = []
    for line in open(param_value['INFILE_SNC'],'r'):
        a=line.split()
        if a[0][0] != "#":
            arm.append(int(a[0]))
            wav.append(float(a[2]))
            mag.append(float(a[6]))
            snc.append(float(a[3]))
            if param_value['MR_MODE'].lower() == 'yes':
                if int(a[0]) == 0:
                    dsp.append((float(a[2])/R0)**2)
                elif int(a[0]) == 1:
                    dsp.append((float(a[2])/R3)**2)
                else:
                    dsp.append((float(a[2])/R2)**2)
            else:
                if int(a[0]) == 0:
                    dsp.append((float(a[2])/R0)**2)
                elif int(a[0]) == 1:
                    dsp.append((float(a[2])/R1)**2)
                else:
                    dsp.append((float(a[2])/R2)**2)
    ## data output ##
    if PYFITS_FLG == 1 and NUMPY_FLG == 1:
        if param_value['FILE_TYPE_SIM'].lower()=='fits':
            print "FITS table mode"
            flx = []
            var  = []
            for i in range(len(wav)):
                fnu   = 10**(-0.4*(mag[i]+48.6))
                flam  = 3.0e18*fnu/(10*wav[i])**2/1e-17
                sigma = flam/(snc[i]+offset)
                flx.append(flam+rn.gauss(0.0,sigma))
                var.append(sigma**2)
            msk = np.zeros(len(wav))
            sky = np.zeros(len(wav))
            xps = np.zeros(len(wav))
            cfg = np.zeros(len(wav))
            col1 = pf.Column(name='FLUX', format='E', array=np.array(flx)) 
            col2 = pf.Column(name='COVAR', format='E', array=np.array(var)) 
            col3 = pf.Column(name='MASK', format='J', array=msk) 
            col4 = pf.Column(name='WAVELENGTH', format='E', array=np.array(wav)) 
            col5 = pf.Column(name='WAVEDISP', format='E', array=np.array(dsp)) 
            col6 = pf.Column(name='SKY', format='E', array=sky)
            col7 = pf.Column(name='XPOSITION', format='E', array=xps)
            col8 = pf.Column(name='CONFIG', format='E', array=cfg)
            cols = pf.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])
            tbhdu = pf.BinTableHDU.from_columns(cols)
            os.system('rm -f %s'%(param_value['OUTFILE_SIM']))
            tbhdu.writeto(param_value['OUTFILE_SIM'])
        else:    
            print "ASCII table mode"
            file=open(param_value['OUTFILE_SIM'],'w')
            file.write('''#  NUM  CONTENT  [TYPE]  [UNIT]
#  1  FLUX       [FLOAT] [10^-17 erg/s/cm^2/A]
#  2  COVAR      [FLOAT] [10^-17 erg/s/cm^2/A]^2
#  3  MASK       [INT32] [---]
#  4  WAVELENGTH [FLOAT] [nm]
#  5  WAVEDISP   [FLOAT] [nm]
#  6  SKY        [FLOAT] [10^-17 erg/s/cm^2/A]
#  7  XPOSITION  [FLOAT] [---]
#  8  CONFIG     [---]   [---]
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
                file.write('%12.4e %12.4e %2d %8.3f %6.3f %12.4e %.2f %d\n'%(flx,var,msk,wav[i],dsp[i],sky,xps,cfg))
            file.close()
    else:    
        print "ASCII table mode"
        file=open(param_value['OUTFILE_SIM'],'w')
        file.write('''#  NUM  CONTENT  [TYPE]  [UNIT]
#  1  FLUX       [FLOAT] [10^-17 erg/s/cm^2/A]
#  2  COVAR      [FLOAT] [10^-17 erg/s/cm^2/A]^2
#  3  MASK       [INT32] [---]
#  4  WAVELENGTH [FLOAT] [nm]
#  5  WAVEDISP   [FLOAT] [nm]
#  6  SKY        [FLOAT] [10^-17 erg/s/cm^2/A]
#  7  XPOSITION  [FLOAT] [---]
#  8  CONFIG     [---]   [---]
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
            file.write('%12.4e %12.4e %2d %8.3f %6.3f %12.4e %.2f %d\n'%(flx,var,msk,wav[i],dsp[i],sky,xps,cfg))
        file.close()
