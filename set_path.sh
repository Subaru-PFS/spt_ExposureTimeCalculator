#!/bin/sh
grep -v "HOME_DIR" run_etc.defaults > tmp.txt
echo "HOME_DIR         $PWD   # FULL PATH to the package" >> tmp.txt
mv tmp.txt run_etc.defaults
grep -v "HOME_DIR" gen_sim_spec.defaults > tmp.txt
echo "HOME_DIR         $PWD   # FULL PATH to the package" >> tmp.txt
mv tmp.txt gen_sim_spec.defaults
echo "Done!"
