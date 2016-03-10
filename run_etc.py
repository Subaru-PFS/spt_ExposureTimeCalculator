#!/usr/bin/env python
import math,os,sys
import re
import argparse
import time
import subprocess, shlex

if __name__ == '__main__':
    ### Some parameters ###########################
    ### CAUTION: ##################################
    ### CHANGE BELOW ON YOUR OWN RESPONSIBILITY ###
    ###############################################
    ETC            = 'gsetc'
    INSTR_SETUP    = './config/PFS.dat'
    INSTR_SETUP_MR = './config/PFS.redMR.dat'
    SKYMODELS      = '11005'
    OFFSET_FIB     = 0.00
    SKY_SUB_FLOOR  = 0.01
    DIFFUSE_STRAY  = 0.02
    ###############################################
    start = time.time()
    if os.path.exists('bin') == False:
        os.mkdir('bin')
    try:
        p_gcc = subprocess.call(shlex.split("gcc ./src/%s.c -lm -O3 -DHGCDTE_SUTR -DMOONLIGHT_ -o ./bin/%s.x" % (ETC, ETC)))
        if p_gcc != 0:
            exit("Compile error: %d" % p_gcc)
    except OSError, e:
        exit("Compile error: %s (%s)" % ETC, e)
    print "ETC ready..."
    param_name=['SKYMODELS','SKY_SUB_FLOOR','DIFFUSE_STRAY','OFFSET_FIB']
    param_value={'SKYMODELS':SKYMODELS,'SKY_SUB_FLOOR':str(SKY_SUB_FLOOR),'DIFFUSE_STRAY':str(DIFFUSE_STRAY),'OFFSET_FIB':str(OFFSET_FIB)}
    parser = argparse.ArgumentParser(description='PFS ETC developed by Chris Hirata, modified by Kiyoto Yabe and Yuki Moritani')
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
#        try:
#            os.system('tr \\\\r \\\\n \<default.param \>default2.param')
#        except:
#            print "warning!"
#        dat=open(args.params,'r')
#        print dir(dat)
#        print dat.read()
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
 
    ## Medium Resolution Mode ? ##
    if param_value['MR_MODE'].lower() == 'yes' or param_value['MR_MODE'].lower() == 'y':
        param_value['INSTR_SETUP'] = INSTR_SETUP_MR
    else:
        param_value['INSTR_SETUP'] = INSTR_SETUP
    ## make continuum magnitude file ##
    try:
        mag = float(param_value['MAG_FILE'])
        if os.path.exists('tmp') == False:
            os.mkdir('tmp')
        file = open('tmp/mag.dat','w')
        file.write('300.0 %.2f\n 1300. %.2f\n'%(mag,mag))
        file.close()
        mag_file = 'tmp/mag.dat'
    except:
        mag_file = param_value['MAG_FILE']
    ## reuse noise data ? ##
    flag = '0'
    if param_value['NOISE_REUSED'].lower()=='y':
        flag='1'
    elif param_value['NOISE_REUSED'].lower()=='n':
        flag='0'
    param_value['NOISE_REUSED']=flag
    ## check file overwritten ##
    C = 0
    if param_value['OVERWRITE'].lower() == 'no' or param_value['OVERWRITE'].lower() == 'n':
        if os.path.exists(param_value['OUTFILE_NOISE']):
            print "Error: %s already exists... "%(param_value['OUTFILE_NOISE'])
            C += 1
        if os.path.exists(param_value['OUTFILE_SNC']):
            print "Error: %s already exists... "%(param_value['OUTFILE_SNC'])
            C += 1
        if os.path.exists(param_value['OUTFILE_SNL']):
            print "Error: %s already exists... "%(param_value['OUTFILE_SNL'])
            C += 1
    if param_value['OVERWRITE'].lower() == 'yes' or param_value['OVERWRITE'].lower() == 'y':
            C = 0
    ## run ETC ##
    print param_value['INSTR_SETUP']
    if C != 0:
        exit('No execution of ETC')
    try:
        p_etc = subprocess.Popen(['./bin/%s.x' % ETC], stdin = subprocess.PIPE)
        p_etc.communicate("\n".join([
                    param_value['INSTR_SETUP'],
                    param_value['SKYMODELS'],
                    param_value['SEEING'],
                    param_value['ZENITH_ANG'],
                    param_value['GALACTIC_EXT'],
                    param_value['FIELD_ANG'],
                    param_value['OFFSET_FIB'],
                    param_value['MOON_ZENITH_ANG'],
                    param_value['MOON_TARGET_ANG'],
                    param_value['MOON_PHASE'],
                    param_value['EXP_TIME'],
                    param_value['EXP_NUM'],
                    param_value['SKY_SUB_FLOOR'],
                    param_value['DIFFUSE_STRAY'],
                    param_value['NOISE_REUSED'],
                    param_value['OUTFILE_NOISE'],
                    param_value['OUTFILE_OII'],
                    param_value['OUTFILE_SNL'],
                    param_value['LINE_FLUX'],
                    param_value['LINE_WIDTH'],
                    param_value['OUTFILE_SNC'],
                    '-',
                    mag_file,
                    param_value['REFF']
                ]))
    except OSError, e:
        exit('Execution error of "%s" (%s)' % ETC, e)
## end of the script ##
    elapsed_time = time.time() - start
    print "elapsed_time:{0}".format(elapsed_time) + "[sec]"
