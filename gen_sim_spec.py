#!/usr/bin/env python
import math,os,sys
import argparse
import hashlib
import random as rn
import time

try:
    import pyfits as pf
    PYFITS_FLG = 1
except:
    PYFITS_FLG = 0
try:
    import numpy as np
    import scipy as sp
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

def arm_name(arm_num,res_mode):
    if arm_num == 0:
        return 'b'
    elif arm_num == 1:
        if res_mode == 'MR':
            return 'm'
        else:
            return 'r'
    else:
        return 'n'

def calculate_pfsVisitHash(visits):
    """Calculate and return the 64-bit SHA-1 from a list of visits"""
    if len(visits) == 1 and visits[0] == 0:
        return 0
    m = hashlib.sha1()
    for visit in visits:
        m.update("%d" % (visit))
    # convert to int and truncate to 8 bytes
    sha1 = int(m.hexdigest(), 16) & 0xffffffffffffffff
    return sha1

def calculate_pfsConfigId(fiberIds, ras, decs):
    """Calculate and return the 64-bit SHA-1 from a set of lists of
    fiberId, ra, and dec"""
    m = hashlib.sha1()
    for fiberId, ra, dec in zip(fiberIds, ras, decs):
        m.update("%d %.0f %.0f" % (fiberId, ra, dec))
    # convert to int and truncate to 8 bytes
    sha1 = int(m.hexdigest(), 16) & 0xffffffffffffffff
    return sha1

def driver(pfsConfigId, tract=0, patch='0,0', fiberId=None, dataDir="./out/", objId=None):
    pfsConfig = PfsConfig(pfsConfigId)
    pfsConfig.read(dataDir)

    if objId is None:
        if fiberId is None:
            fiberId = 1
    else:
        try:
            _fiberId = np.where(pfsConfig.objId == objId)[0][0]
        except IndexError:
            print >> sys.stderr, "Unable to find objId %08x in configuration with pfsConfigId 0x%08x" % \
                (objId, pfsConfigId)
            return
        _fiberId += 1                   # 1-indexed
        if fiberId is not None and fiberId != _fiberId:
            print >> sys.stderr, "fiberId %d doesn't match objId %08x's fiber %d" % \
                (fiberId, objId, _fiberId)
            return
        fiberId = _fiberId
    objId = pfsConfig.objId[fiberId - 1]
#    print "fiberId = %d,  objId = 0x%x" % (fiberId, objId)
    pfsArms = PfsArmSet(visits=1, spectrograph=1)
    pfsArms.read(dataDir)
    pfsObject = makePfsObject(tract, patch, objId, pfsArms)
    pfsObject.write(dataDir)
    npfs = PfsObject(tract, patch, objId, visits=[1])
    npfs.read(dataDir)

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
    parser.add_argument("--NREALIZE", type=str, help="the number of realization")
    parser.add_argument("--OUTFILE_SIM", type=str, help="simulated spectrum output ASCII file")
    parser.add_argument("--OUTFILE_TRACT", type=str, help="tract")
    parser.add_argument("--OUTFILE_PATCH", type=str, help="patch")
    parser.add_argument("--OUTFILE_CATID", type=str, help="catalogue id")
    parser.add_argument("--OUTFILE_VISIT", type=str, help="visit number")
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
## The number of simulation ##
    nreal = int(param_value['NREALIZE'])
## load exposure number ##
    nexp = float(param_value['EXP_NUM'])
## resolution mode ##
    if param_value['MR_MODE'].lower() == 'yes' or param_value['MR_MODE'].lower() == 'y':
        res_mode = 'MR'
    else:
        res_mode = 'LR'
    ## read input file ##
    arm_list = []
    arm = []
    wav = []
    dsp = []
    nsv = []
    snc = []
    trn = []
    smp = []
    skm = []
    for line in open(param_value['INFILE_SNC'],'r'):
        a=line.split()
        if a[0][0] != "#":
            if int(a[0]) not in arm_list:
                arm_list.append(int(a[0]))
            arm.append(int(a[0]))
            wav.append(float(a[2]))
            nsv.append(float(a[5]))
            snc.append(float(a[3]))
            if float(a[8])>1.0e+26:
                trn.append(float(a[8]))
            else:
                trn.append(1.0e+26)
            smp.append(float(a[9]))
            skm.append(float(a[10]))
            if res_mode == 'MR':
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
## load magnitude or filename ##
    try:
        input_mag = float(param_value['MAG_FILE'])
        mag = np.ones(len(wav))*input_mag
    except:
        tmp1,tmp2 = np.loadtxt(param_value['MAG_FILE'],usecols=(0,1),unpack=True)
        mag = np.interp(wav,tmp1,tmp2)
## data output ##
    obj_o = np.zeros((nreal,len(wav)))
    flx_o = np.zeros((nreal,len(wav)))
    err_o = np.zeros((nreal,len(wav)))
    var_o = np.zeros((nreal,len(wav)))
    msk_o = np.zeros((nreal,len(wav)))
    sky_o = np.zeros((nreal,len(wav)))
    npix_arm = np.zeros(len(arm_list))
    for Creal in range(nreal):
## ASCII mode ##
        if nreal == 1:
            file = open(param_value['OUTFILE_SIM'],'w')
        else:
            file = open(param_value['OUTFILE_SIM']+'.%d' % (Creal),'w')
        file.write('''#  1  WAVELENGTH  [nm]
    #  2  FLUX        [10^-17 erg/s/cm^2/A]
    #  3  ERROR       [10^-17 erg/s/cm^2/A]
    #  4  MASK        [1=masked]
    #  5  SKY         [10^-17 erg/s/cm^2/A]
    #  6  ARM         [0=blue,1=red,2=NIR,3=redMR]
    #  7  Signal-to-Noise Ratio
    ''')
        for i in range(len(wav)):
            if Creal == 0:
                for j in arm_list:
                    if arm[i] == arm_list[j]:
                        npix_arm[j] += 1
            counts = 3.631e-20 * 10**(-0.4*mag[i]) * trn[i]
            snr = counts/np.sqrt(smp[i]*counts + nsv[i])*np.sqrt(nexp)
            fnu   = 10**(-0.4*(mag[i]+48.6))
            flam  = 3.0e18*fnu/(10*wav[i])**2/1e-17
            sigma = flam/(snr+offset)
            flx = flam+rn.gauss(0.0,sigma)
            var = sigma**2
            msk = 0.0
            sky = 3.0e18 * (skm[i] / trn[i]) / (10*wav[i])**2/1e-17
            obj_o[Creal,i] = flam
            flx_o[Creal,i] = flx
            err_o[Creal,i] = sigma
            var_o[Creal,i] = var
            msk_o[Creal,i] = msk
            sky_o[Creal,i] = sky
            file.write('%8.3f %12.4e %12.4e %2d %12.4e %1d %12.4e\n' % (wav[i],flx,sigma,msk,sky,arm[i],snr))
        file.close()
    print "ASCII tables was generated"
## FITS mode ##
    if res_mode == 'MR':
        resampling = 0.04 # resampled pix scale [nm/pix]
    else:
        resampling = 0.08 # resampled pix scale [nm/pix]
    if PYFITS_FLG == 1 and NUMPY_FLG == 1:
## some parameters ##
        tract = int(param_value['OUTFILE_TRACT'])
        patch = param_value['OUTFILE_PATCH']
        catid = int(param_value['OUTFILE_CATID'])
        visit = int(param_value['OUTFILE_VISIT'])
## check wavelength boundary of arms ##
        arm_wav_min = []
        arm_wav_max = []
        for arm_num in arm_list:
            wav_a = []
            for i in range(len(wav)):
                if arm[i] == arm_num:
                    wav_a.append(wav[i])
            arm_wav_min.append(min(wav_a))
            arm_wav_max.append(max(wav_a))
## here assumed that arm = 0, 1, 2 (1=MR if res_mode == 'MR') ##
        if res_mode == 'MR':
            arm_wav_rng_min = [arm_wav_min[0],arm_wav_min[1],arm_wav_min[2]]
            arm_wav_rng_max = [arm_wav_max[0],arm_wav_max[1],arm_wav_max[2]]
        else:
            arm_wav_rng_min = [arm_wav_min[0],0.5*(arm_wav_min[1]+arm_wav_max[0]),0.5*(arm_wav_min[2]+arm_wav_max[1])]
            arm_wav_rng_max = [0.5*(arm_wav_min[1]+arm_wav_max[0]),0.5*(arm_wav_min[2]+arm_wav_max[1]),arm_wav_max[2]]
## pfsConfig ##
        mpscen = np.zeros((1,2))
        objra = np.zeros(nreal) + 150.0
        objdec = np.zeros(nreal) +2.0
        fiberid = np.arange(nreal) + 1
        objid = np.arange(nreal)
        pfsConfigId = calculate_pfsConfigId(fiberIds=fiberid, ras=objra, decs=objdec)
        pfsVisitHash = calculate_pfsVisitHash(visits=[visit])
#        print pfsConfigId
#        print pfsVisitHash
        config_file_name = 'pfsConfig-0x%08x.fits' % (pfsConfigId)
        config_file = './out/%s' % (config_file_name)
# catId, objId, ra, dec, fiber flux, MPS centroid
        col1 = pf.Column(name = 'fiberId', format = 'J', array = fiberid)
        col2 = pf.Column(name = 'catId', format = 'J', array = np.zeros(nreal) + catid)
        col3 = pf.Column(name = 'objId', format = 'K', array = objid)
        col4 = pf.Column(name = 'ra', format = 'E', array = objra)
        col5 = pf.Column(name = 'dec', format = 'E', array = objdec)
        try:
            mag = float(param_value['MAG_FILE'])
            fnu = 10 ** (-0.4 * (mag + 48.6))
        except:
            mag = 99.0000
            fnu = 0.0
        col6 = pf.Column(name = 'fiber flux', format = 'E', array = np.zeros(nreal) + fnu)
        col7 = pf.Column(name = 'MPS centroid', format = '2E', array = mpscen)
        hdu1 = pf.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7])
        hdu = pf.PrimaryHDU()
        hdulist = pf.HDUList([hdu,hdu1])
        hdulist[1].header['EXTNAME'] = 'CONFIG'
        os.system('rm -f %s' % (config_file))
        hdulist.writeto(config_file)
        print "%s was generated" % (config_file)
## single exposure for each arm ##
## --> pfsArm-%6d-%1d%1s.fits ##
        C=0
        for arm_num in arm_list:
            output_file = './out/pfsArm-%06d-%1d%1s.fits' % (visit,1,arm_name(arm_num=arm_num,res_mode=res_mode))
            wav_a = np.zeros((nreal,npix_arm[C]))
            obj_a = np.zeros((nreal,npix_arm[C]))
            flx_a = np.zeros((nreal,npix_arm[C]))
            err_a = np.zeros((nreal,npix_arm[C]))
            var_a = np.zeros((nreal,npix_arm[C]))
            msk_a = np.zeros((nreal,npix_arm[C]))
            sky_a = np.zeros((nreal,npix_arm[C]))
            cov_a = np.zeros((nreal,3,npix_arm[C]))
            for Creal in range(nreal):
                j = 0
                for i in range(len(wav)):
                    if arm[i] == arm_num:
                        wav_a[Creal,j] = wav[i]
                        obj_a[Creal,j] = obj_o[Creal,i]
                        flx_a[Creal,j] = flx_o[Creal,i]
                        err_a[Creal,j] = err_o[Creal,i]
                        var_a[Creal,j] = var_o[Creal,i]
                        msk_a[Creal,j] = msk_o[Creal,i]
                        sky_a[Creal,j] = sky_o[Creal,i]
                        cov_a[Creal,0,j] = var_o[Creal,i]
                        j += 1
            Npix = npix_arm[C]
            col1 = np.array(flx_a, dtype='f4')
            col2 = np.array(cov_a, dtype='f4')
            col3 = np.array(msk_a, dtype='i4')
            col4 = np.array(wav_a, dtype='f4')
            col5 = np.array(sky_a, dtype='f4')
            col61 = pf.Column(name = 'pfsConfigId', format = 'K', array = np.array([pfsConfigId]))
            col62 = pf.Column(name = 'visit', format = 'J', array = np.array([visit]))
            hdu1 = pf.ImageHDU(data = col1, header=None)
            hdu2 = pf.ImageHDU(data = col2, header=None)
            hdu3 = pf.ImageHDU(data = col3, header=None)
            hdu4 = pf.ImageHDU(data = col4, header=None)
            hdu5 = pf.ImageHDU(data = col5, header=None)
            hdu6 = pf.BinTableHDU.from_columns([col61,col62])
            hdu = pf.PrimaryHDU()
            hdulist = pf.HDUList([hdu])
            hdulist.append(hdu1)
            hdulist.append(hdu2)
            hdulist.append(hdu3)
            hdulist.append(hdu4)
            hdulist.append(hdu5)
            hdulist.append(hdu6)
            hdulist[1].header['EXTNAME'] = 'FLUX'
            hdulist[2].header['EXTNAME'] = 'COVAR'
            hdulist[3].header['EXTNAME'] = 'MASK'
            hdulist[4].header['EXTNAME'] = 'WAVELENGTH'
            hdulist[5].header['EXTNAME'] = 'SKY'
            hdulist[6].header['EXTNAME'] = 'CONFIG'
            os.system('rm -f %s' % (output_file))
            hdulist.writeto(output_file)
            print "%s was generated" % (output_file)
            C += 1
## pfsObject (generated by Robert Lupton's code ##            
        from pfs.datamodel.pfsConfig import PfsConfig
        from pfs.datamodel.pfsArm import PfsArmSet
        from pfs.datamodel.pfsObject import PfsObject, makePfsObject
        for i in range(len(fiberid)):
            driver(pfsConfigId = pfsConfigId, tract = tract, patch = patch, fiberId = fiberid[i], objId = objid[i])
            print "pfsObject-%04d-%s-%03d-%08x-%02d-0x%08x.fits was generated" % (tract, patch, catid, objid[i], visit % 100, pfsVisitHash)
    else:
        print "FITS file was NOT generated"
        print "Please install pyfits and numpy to generate FITS format"    

