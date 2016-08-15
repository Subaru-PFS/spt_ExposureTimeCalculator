#!/usr/bin/env python
import math,os,sys
import re
import argparse
import time
import subprocess, shlex

#######################
offset = 0.01
HOME_DIR = "path-to-etc"
#######################

def main():
    def convert_arg_line_to_args(arg_line):
        """Make argparse handle Hirata-style parameter files"""
        for i, arg in enumerate(arg_line.split()):
            if not arg.strip():
                continue
            if arg[0] == '#':
                return
            if i == 0:
                arg = "--%s" % arg
            yield arg
    ### Some parameters ###########################
    ### CAUTION: ##################################
    ### CHANGE BELOW ON YOUR OWN RESPONSIBILITY ###
    ###############################################
    SKYMODELS      = '11005'
    OFFSET_FIB     = '0.00'
    SKY_SUB_FLOOR  = '0.01'
    DIFFUSE_STRAY  = '0.02'
    ###############################################
    start = time.time()
    parser = argparse.ArgumentParser(description='PFS ETC developed by Chris Hirata, modified by Kiyoto Yabe, Yuki Moritani, and Atsushi Shimono', fromfile_prefix_chars='@')
    parser.convert_arg_line_to_args = convert_arg_line_to_args
    parser.add_argument("--SEEING", type=str, help="seeing", default="0.80")
    parser.add_argument("--ZENITH_ANG", type=str, help="zenith angle", default="45.00")
    parser.add_argument("--GALACTIC_EXT", type=str, help="galactic extinction", default="0.00")
    parser.add_argument("--MOON_ZENITH_ANG", type=str, help="moon-zenith angle", default="30.0")
    parser.add_argument("--MOON_TARGET_ANG", type=str, help="moon-target angle", default="60.0")
    parser.add_argument("--MOON_PHASE", type=str, help="moon phase", default="0")
    parser.add_argument("--EXP_TIME", type=str, help="exposure time", default="450")
    parser.add_argument("--EXP_NUM", type=str, help="exposure number", default="8")
    parser.add_argument("--FIELD_ANG", type=str, help="field angle", default="0.675")
    parser.add_argument("--MAG_FILE", type=str, help="magnitude input file", default="22.5")
    parser.add_argument("--REFF", type=str, help="effective radius", default="0.3")
    parser.add_argument("--LINE_FLUX", type=str, help="emission line flux", default="1.0e-17")
    parser.add_argument("--LINE_WIDTH", type=str, help="emission line width (sigma)", default="70")
    parser.add_argument("--NOISE_REUSED", type=str, help="noise vector reused flag", default="N")
    parser.add_argument("--OUTFILE_NOISE", type=str, help="noise vector output file", default="out/ref.noise.dat")
    parser.add_argument("--OUTFILE_SNC", type=str, help="continuum results output file", default="out/ref.snc.dat")
    parser.add_argument("--OUTFILE_SNL", type=str, help="emission line results output file", default="out/ref.snl.dat")
    parser.add_argument("--OUTFILE_OII", type=str, help="[OII] emission line results output file", default="-")
    parser.add_argument("--MR_MODE", type=str, help="medium resolution mode on-off", default="N")
    parser.add_argument("--OVERWRITE", type=str, help="overwrite on-off", default="Y")

    args = parser.parse_args()

    if not os.path.exists(HOME_DIR):
        exit("Unable to find path; please run make and try again")

    ETC            = HOME_DIR + 'bin/gsetc.x'
    INSTR_SETUP    = HOME_DIR + 'config/PFS.dat'
    INSTR_SETUP_MR = HOME_DIR + 'config/PFS.redMR.dat'

    if not os.path.exists(HOME_DIR + 'bin'):
        os.mkdir(HOME_DIR + 'bin')
    if not os.path.exists('out'):
        os.mkdir('out')
    if not os.path.exists(ETC):
        exit("Unable to find %s; please run make and try again" % ETC)

    ## Medium Resolution Mode ? ##
    if args.MR_MODE.lower() == 'yes' or args.MR_MODE.lower() == 'y':
        INSTR_SETUP = INSTR_SETUP_MR
    else:
        INSTR_SETUP = INSTR_SETUP
    ## make continuum magnitude file ##
    try:
        mag = float(args.MAG_FILE)
        if os.path.exists(HOME_DIR + 'tmp') == False:
            os.mkdir(HOME_DIR + 'tmp')
        file = open(HOME_DIR + 'tmp/mag.dat','w')
        file.write('300.0 %.2f\n 1300. %.2f\n'%(mag,mag))
        file.close()
        mag_file = HOME_DIR + 'tmp/mag.dat'
    except:
        mag_file = args.MAG_FILE
    ## reuse noise data ? ##
    flag = '0'
    if args.NOISE_REUSED.lower()=='y':
        flag='1'
    elif args.NOISE_REUSED.lower()=='n':
        flag='0'
    args.NOISE_REUSED=flag
    ## check file overwritten ##
    C = 0
    if args.OVERWRITE.lower() == 'no' or args.OVERWRITE.lower() == 'n':
        if os.path.exists(args.OUTFILE_NOISE):
            print "Error: %s already exists... "%(args.OUTFILE_NOISE)
            C += 1
        if os.path.exists(args.OUTFILE_SNC):
            print "Error: %s already exists... "%(args.OUTFILE_SNC)
            C += 1
        if os.path.exists(args.OUTFILE_SNL):
            print "Error: %s already exists... "%(args.OUTFILE_SNL)
            C += 1
    if args.OVERWRITE.lower() == 'yes' or args.OVERWRITE.lower() == 'y':
            C = 0
    ## run ETC ##
    print INSTR_SETUP
    if C != 0:
        exit('No execution of ETC')
    try:
        p_etc = subprocess.Popen([ETC], stdin = subprocess.PIPE)
        p_etc.communicate("\n".join([
                    INSTR_SETUP,
                    SKYMODELS,
                    args.SEEING,
                    args.ZENITH_ANG,
                    args.GALACTIC_EXT,
                    args.FIELD_ANG,
                    OFFSET_FIB,
                    args.MOON_ZENITH_ANG,
                    args.MOON_TARGET_ANG,
                    args.MOON_PHASE,
                    args.EXP_TIME,
                    args.EXP_NUM,
                    SKY_SUB_FLOOR,
                    DIFFUSE_STRAY,
                    args.NOISE_REUSED,
                    args.OUTFILE_NOISE,
                    args.OUTFILE_OII,
                    args.OUTFILE_SNL,
                    args.LINE_FLUX,
                    args.LINE_WIDTH,
                    args.OUTFILE_SNC,
                    '-',
                    mag_file,
                    args.REFF
                ]))
    except OSError, e:
        exit('Execution error of "%s" (%s)' % ETC, e)
## end of the script ##
    elapsed_time = time.time() - start
    print "elapsed_time: %.1f[sec]" % (elapsed_time)
    return 0

if __name__ == '__main__':
    main()