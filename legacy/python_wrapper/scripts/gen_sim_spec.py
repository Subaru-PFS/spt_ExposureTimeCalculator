#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np

from pfsspecsim import pfsspec


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
    parser = argparse.ArgumentParser(description='PFS Spectral Simulator developed by Kiyoto Yabe and Robert Lupton',
                                     fromfile_prefix_chars='@')
    parser.convert_arg_line_to_args = convert_arg_line_to_args
    parser.add_argument("--EXP_NUM", type=int, help="Number of exposures", default=4)
    parser.add_argument("--MAG_FILE", type=str, help="magnitude input file", default="22.5")
    parser.add_argument("--etcFile", type=str, help="continuum results input file", default="out/ref.snc.dat")
    parser.add_argument("--nrealize", type=int, help="the number of realization", default=1)
    parser.add_argument("--outDir", type=str, help="Directory for outputs", default="out")
    parser.add_argument("--tract", type=int, help="tract", default=0)
    parser.add_argument("--patch", type=str, help="patch", default='0,0')
    parser.add_argument("--visit0", type=int, help="the first visit number", default=1)
    parser.add_argument("--objId", type=int, help="object id of first realisation ", default=1)
    parser.add_argument("--catId", type=int, help="catalogue id", default=0)
    parser.add_argument("--spectrograph", type=int, help="spectrograph number", default=1)
    parser.add_argument("--countsMin", type=float, help="Minimum counts per pixel for noise", default=0.1)
    parser.add_argument("--asciiTable", type=str, help="simulated spectrum output ASCII file", default='None')
    parser.add_argument("--writeFits", type=str, help="Write FITS files", default="True")
    parser.add_argument("--writePfsArm", type=str, help="Write pfsArm files (writeFits must be set)",
                        default="True")
    parser.add_argument("--plotArmSet", action='store_true', help="Plot the pfsArmSet data")
    parser.add_argument("--plotObject", action='store_true', help="Plot the pfsObject data")
    args = parser.parse_args()

    sim = pfsspec.Pfsspec()
    for arg in args._get_kwargs():
        sim.set_param(arg[0], arg[1])
    sim.make_sim_spec()

    return 0


if __name__ == '__main__':
    main()
