#!/usr/bin/env python

import os
import sys
import argparse
import time
import subprocess
import numpy as np

from pfsspecsim import pfsetc


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


def main():
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
    parser.add_argument("--INFILE_OIICat", type=str, help="input catalogue for [OII] emitters", default="-")
    parser.add_argument("--OUTFILE_OIICat", type=str, help="output catalogue for [OII] emitters", default="-")
    parser.add_argument("--minSNR", type=str, help="minimum SNR for [OII] emission", default="9.0")
    args = parser.parse_args()

    etc = pfsetc.Etc()
    for arg in args._get_kwargs():
        etc.set_param(arg[0], arg[1])
    etc.run()

    return 0

if __name__ == '__main__':
    main()
