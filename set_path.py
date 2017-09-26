#!/usr/bin/env python
import os
import subprocess

file1 = 'run_etc.py'
file2 = 'gen_sim_spec.py'

if __name__ == '__main__':
    path = os.path.abspath(os.path.dirname(__file__)) + '/'
    file = open('tmp.txt', 'w')
    for line in open(file1, 'r'):
        if "path-to-etc" in line:
            file.write(line.replace('path-to-etc', path))
        else:
            file.write(line)
    file.close()
    os.rename('tmp.txt', file1)
    file = open('tmp.txt', 'w')
    for line in open(file2, 'r'):
        if "path-to-etc" in line:
            file.write(line.replace('path-to-etc', path))
        else:
            file.write(line)
    file.close()
    os.rename('tmp.txt', file2)
    print "Set HOME_DIR to %s" % (path)
