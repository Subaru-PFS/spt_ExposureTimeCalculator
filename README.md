PFS Exposure Time Calculator and Spectrum Simulator 
============================================================

This package is developed by the following people
---------------------------------
The original ETC was developed and written by Christopher Hirata (Ohio State University), which is based on the package developed for WFIRST (C. Hirata; arXiv:1204.5151) and altered for use in PFS project.

The code modification, the python wrapping, and the development of the spectral simulator were done by Kiyoto Yabe, Yuki Moritani, Atsushi Shimono (Kavli IPMU) and Robert Lupton (Princeton University)


Release Note
------------
* Version 1.0  Feb. 26, 2016
* Version 1.1  Apr. 27, 2016
* Version 1.2  Feb. 05, 2021

Requirements
------------
* Standard C compiler (e.g., GCC)
* Python3 is recommended (Python2 is NOT fully supported)
* scipy   (0.14 and higher)  see http://www.scipy.org
* pyfits  (3.3 and higher)  see http://www.stsci.edu/institute/software_hardware/pyfits
* matplotlib (if you use the plotting options)
* Sufficient computing power
 (Note1) pyfits is required only for using gen_sim_spec.py and the PFS datamodel package. If you don't have these modules, please install them from the above website. The version of the module is the minimum one that we confirmed so far. If you have any updates, let me know please.
 (Note2) Standard unix system including Linux and Mac OSX is recommended. There has been reported that this code does not work properly on a Linux system mounted on a Windows drive. This package is tested under Mac OSX 10.9.5 on 2.8GHz Quad-Core Intel Xeon machine and Fedora Core 20 on Intel Core i5-4690 3.50GHz machine. Depending on the machine power, it takes **several minutes** if you run all the standard process. We recommend sufficient computing power at least similar to that we have tested. With our testing machine above, it takes about ~90 sec. (~40 sec. for noise calculation, ~25 sec. for emission line S/N calculation, and ~20 sec. for continuum S/N calculation).

Installation
------------
To install the package, get the git repository by typing the following command on the terminal (if you have git installed):
  
    git clone --recursive https://github.com/Subaru-PFS/spt_ExposureTimeCalculator.git
    cd spt_ExposureTimeCalculator
    make
    python setup.py install

Once you clone the repository, you can pull updates from the next time on the directory like this:

    git pull
    git submodule update --init
    make
    python setup.py install

You also can get the zip or tar ball from the following page:

http://sumire.pbworks.com/w/page/107534730/PFS%20ETC

Before you use the package, please reed `README.md` carefully.

Description
-----------
This package includes two parts: one is the exposure time calculator (scripts/run_etc.py) and the other one is the spectral simulator (scripts/gen_sim_spec.py). You can get S/N information of an object in a given exposure time and various conditions by using the ETC, which is based on the "Chris Hirata's simulator". By using the results from the ETC, you can get the simulated spectra in the format of the current PFS datamodel with the spectral simulator.

Exposure Time Calculator (ETC)
------------------------------
The ETC can be run as follows:

    python run_etc.py <input parameter file> [--param1=value1] [--param2=value2] ...

In the default setting, the code can be run by typing the following command at your terminal (the @ means that this is an input file):

    python run_etc.py @run_etc.defaults

Or using some arguments as follows:

    python run_etc.py @run_etc.defaults --MAG_FILE=23.0 --LINE_FLUX=5.0e-17 --LINE_WIDTH=100

All parameters can be specified in a parameter file or passed as arguments. Here are a list and the description of parameters that users can manage.

### Setup of the observation condition 
In this section, users can manage the observational conditions including seeing, zenith angle, galactic extinction, lunar condition, exposure time, and distance from FoV center.

For instance, if you want observational condition under 0.80 arcsec seeing, the target zenith angle of 45 deg., galactic extinction of 0.05 mag., new moon with a zenith angle of 30 deg. and target - moon separation of 60 deg., 8 exposures with a single exposure of 450 sec. (3600 sec in total), and observation at the edge of the FoV, the parameters in .param file should be as follows: 

| Parameter       | Default value | Description             | Unit                          |
|:---             |:---           |:---                     |:---                           |
| SEEING          | 0.80          | Seeing FWHM size        | [arcsec]                      |
| ZENITH_ANG      | 45.00         | Zenith angle            | [deg.]                        |
| GALACTIC_EXT    | 0.00          | Galactic extinction     | [ABmag.]                      |
| MOON_ZENITH_ANG | 30.0          | Moon zenith angle       | [deg.]                        |
| MOON_TARGET_ANG | 60.0          | Moon-target separation  | [deg.]                        |
| MOON_PHASE      | 0             | Moon phase              | [0=New,0.25=quarter,0.5=full] |
| EXP_TIME        | 450           | Single exposure time    | [sec.]                        |
| EXP_NUM         | 8             | The number of exposures |                               |
| FIELD_ANG       | 0.675         | Field angle             | [deg.; center=0, edge=0.675]  |
| degrade         | 1.0           | Throughput degradation factor |                          |

### Setup of the target information
In this section, users can describe the target information including target magnitude, effective radius, emission line flux and width. If you want to calculate the S/N for a flat continuum with 22.5 ABmag and an emission line with line flux of 1.0x10^-17 erg/s/cm^2, line width of sigma=70 km/s, and effective radius of 0.3 arcsec, the parameters are as follows:

| Parameter       | Default value | Description                 | Unit                  |
|:---             |:---           |:---                         |:---                   |
| MAG_FILE        | 22.5          | Magnitude or input spectrum | [ABmag] or _filename_ |        
| REFF            | 0.3           | Effective radius            | [arcsec]              |
| LINE_FLUX       | 1.0e-17       | Emission line flux          | [erg/s/cm^2]          |
| LINE_WIDTH      | 70            | Emission line width sigma   | [km/s]                |

Note: You can input your own target spectra into the ETC. Just specify the file name for `MAG_FILE`. The file should include wavelength [nm] and magnitude [ABmag] in the first and second column, respectively. Please note that ETC does not convolve the instrument resolution to the input spectrum internally, so the resolution should be considered in the input spectrum in advance. Also note that the wavelength is resampled with the sampling of 0.5A which is slightly larger than the pixel sampling of the PFS detector (~0.7A, ~0.8A, and ~0.8A for blue, red, and NIR arm, respectively). We do not recommend to include "NaN" or other non-numeric values in the input file.

### Setup of the name of the output file
The outputs of the calculation are saved to the files defined by users. `OUTFILE_NOISE` defines the file name for the outputs of the noise calculation, `OUTFILE_SNC` defines the file name for the continuum S/N, and `OUTFILE_SNL` is for the emission line S/N. If you want to measure the S/N for [OII] doublet, kindly use `OUTFILE_OII` for the output file. If you set `NOISE_REUSED` to Y, the ETC will skip the process of generating noise vector which is a time-consuming process. The process time of the ETC can be reduced to roughly half by this mode. This mode is useful if you want to calculate the S/N of objects with various magnitudes and line fluxes with the same noise assumption (zenith angle, field angle, lunar condition, and exposure time). If you use this mode, you need to specified the noise vector file that you want to use in `OUTFILE_NOISE`. Also, You can skip outputting the results of continuum S/N and emission line S/N by replace the name to '-'. If you want to use medium resolution mode in the red arm, set `MR_MODE` to Y.

| Parameter       | Default value     | Description                       | Unit       |
|:---             |:---               |:---                               |:---        |
| NOISE_REUSED    | N                 | Noise Vector Reused?              | [Y/N] (If Y, please specify file name used in OUTFILE_NOISE) |
| OUTFILE_NOISE   | out/ref.noise.dat | Output file for noise vector      | _filename_ |
| OUTFILE_SNC     | out/ref.snc.dat   | Output file for continuum S/N     | _filename_ |
| OUTFILE_SNL     | out/ref.snl.dat   | Output file for emission line S/N | _filename_ | 
| OUTFILE_OII     | -                 | Output file for [OII] doublet S/N | _filename_ |
| MR_MODE         | N                 | Medium resolution mode switch     | [Y/N]      |
| OVERWRITE       | Y                 | Overwrite switch                  | [Y/N]      |

### Descriptions on the output file
Each output file contains the following information:

a. `OUTFILE_NOISE`: noise variance per pixel

* (1) the spectrograph arm number
* (2) the pixel number
* (3) the vacuum wavelength in [nm]
* (4) the noise variance in [e^2/pix] (per one exposure without contribution from object continuum)
* (5) the sky background in [e]

Note: (4) the noise variance in this table does not include the contribution from the object continuum.

b. `OUTFILE_SNC`: S/N of continuum per pixel

* (1)  the spectrograph arm number
* (2)  the pixel number
* (3)  the vacuum wavelength in [nm]
* (4)  the continuum S/N per pixel
* (5)  the signal in one exposure in [e]
* (6)  the noise variance in [e^2/pix] (per one exposure without contribution from object continuum)
* (7)  the noise variance in [e^2/pix] (per one exposure with contribution from object continuum)
* (8)  the input spectra in [ABmag]
* (9)  the conversion factor from flux density [erg/s/cm2/Hz] to [e]
* (10) the sampling factor
* (11) the sky background in [e]

Note: (6) the noise variance in this table includes the contribution from the object continuum.

c. `OUTFILE_SNL`: S/N of emission line

* (1) the vacuum wavelength in [nm]
* (2) the fiber aperture factor
* (3) the effective collecting area [m^2]
* (4) the emission line S/N in blue arm
* (5) the emission line S/N in  red arm
* (6) the emission line S/N in  NIR arm
* (7) the emission line S/N in total

d. `OUTFILE_OII`: S/N of [OII] emission lines (here we assume line ratio is 1:1)

* (1) the redshift
* (2) the vacuum wavelength of [OII]3726 in [nm]
* (3) the vacuum wavelength of [OII]3729 in [nm]
* (4) the fiber aperture factor
* (5) the effective collecting area [m^2]
* (6) the emission line S/N in blue arm
* (7) the emission line S/N in  red arm
* (8) the emission line S/N in  NIR arm
* (9) the emission line S/N in total

Some caution again:

1. `NOISE_REUSED` option is only valid for the same condition (zenith angle, field angle, lunar condition, and exposure time). Please be careful when you use this option and DO NOT use the different values from those used when the noise file was created.
2. You can input your own spectra by using `MAG_FILE` option, but please note that ETC does not convolve the instrument resolution to the input spectrum internally, so the resolution should be considered in the input spectrum in advance. Also note that the wavelength sampling should be larger than the pixel sampling of the PFS detector (~0.7A, ~0.8A, and ~0.8A for blue, red, and NIR arm, respectively).

Spectral Simulator 
------------------
(N.b. you will need the python module numpy/scipy to run the simulator.  To write FITS files you'll also need pyfits, and to use the plotting options you'll need matplotlib)

We can generate simulated outputs for use by the 1-D pipeline using a subset of the outputs from ETC above and then running gen_sim_spec.py.  The only output that we use is `OUTFILE_SNC`, and we only use a subset of the information in that file.  In particular, the settings
    `--EXP_NUM`  (the number of exposures)
    `--MAG_FILE` (the input spectrum)
have *no* effect upon the simulator!

Parameters that do matter are
    `--MR_MODE`  (use low/medium resolution grating in the R arm)
and of course things to do with observing such as the exposure time, seeing, moon, ...
    `--EXP_TIME`
    `--SEEING`
    `--MOON_PHASE` ...
The settings used when running run_etc.py are the ones that will be used to generate an output spectrum. 
There is a default value (N, i.e. the low resolution grating), but we'll be explicit when preparing files for the simulator.

So, before running gen_sim_spec.py, you need to prepare input files with commands like:

    python run_etc.py @run_etc.defaults --MR_MODE=no  --EXP_TIME=450 --OUTFILE_SNC=out/etc-t450-lr.dat --OUTFILE_SNL=-

    python run_etc.py @run_etc.defaults --MR_MODE=yes --EXP_TIME=450 --OUTFILE_SNC=out/etc-t450-mr.dat --OUTFILE_SNL=-

You are now ready to generate simulated spectra.  The input spectrum is given by the `MAG_FILE` and may be either a file (with columns of wavelength and AB magnitude) or a floating point number (the AB magnitude of a flat-spectrum source). The simulator will add appropriate noise, taking into account the number of integrations requested (see cautions in the previous section).  Run the gen_sim_spec.py like this:

    python gen_sim_spec.py --etcFile=out/etc-t450-lr.dat

You can either specify options on the command line (use --help to list them), or by providing a file of overrides by adding e.g. "@gen_sim_spec.defaults" to the command line (the @ means that this is an input file). If you want to plot the spectra, try `--plotObject` or `--plotArm`.

So a more complete example is:

    python gen_sim_spec.py @gen_sim_spec.defaults --outDir=out --etcFile=out/etc-t450-lr.dat --asciiTable=test.sim --writeFits=False --MAG_FILE=20.0

to write the ASCII file (out/test.sim.dat) instead of the fits files, simulating a 20th magnitude flat spectrum source.  The noise is taken from the output of run_etc.py, used without specifying `OUTFILE_SNC` (which we do not recommend!). Note that, as usual, arguments such as `--asciiTable` may be abbreviated.

If you specify a value of `--nrealize` > 1 then multiple realizations of the noise are generated, each with its own objId, and the corresponding number of pfsObject files are written.

1. The parameter file (e.g. gen_sim_spec.defaults) has the same format as the input to the ETC.

| Parameter  | Default value   | Description                                                              |
|:---        |:---             |:---                                                                      |
| etcFile    | out/ref.snc.dat | Input noise file for the simulator (from ETC)                            |
| nrealize   | 1               | The number of realization                                                |
| asciiTable | None            | Output ASCII table name without extension ("None" to omit writing table) |

Here are some descriptions:

a. etcFile: Information about the background+instrumental noise per pixel; the ETC's --OUTFILE_SNC file

b. nrealize: The number of realizations of the noise, incrementing the objId for each

c. asciiTable: The file name of the output ASCII table of the simulated spectrum
The table includes the following items:

* (1)  WAVELENGTH  [nm]
* (2)  FLUX        [10^-17 erg/s/cm^2/A]
* (3)  ERROR       [10^-17 erg/s/cm^2/A]
* (4)  MASK        [1=masked]
* (5)  SKY         [10^-17 erg/s/cm^2/A]
* (6)  ARM         [0=blue,1=red,2=NIR,3=redMR]

The FITS format output will be automatically generated unless you specify "--writeFits False". The filenames are defined by the datamodel, using values:

| Parameter | Default value  | Description           |
|:---       |:---            |:---                   |
| tract     | 0              | Tract                 |
| patch     | 0,0            | Patch                 |
| visit0    | 1              | The fist visit number |
| objId     | 1              | Object ID             |
| catId     | 0              | Catalogue ID          |

The format is defined by the datamodel, which may be found at
     https://github.com/Subaru-PFS/datamodel/blob/master/datamodel.txt

Note 1: The wavelength information on HDU #1,3,5,6 is described in the header of #1 by using `CRPIX1`, `CRVAL1`, and `CD1_1`. The wavelength (lambda [nm]) corresponding to the pixel (X [pix]) is as follows:

      lambda = CRVAL1 + CD1_1*(X - CRPIX1)

Values in these HDU are resampled from the original data so that the pixel scale is 1.0 A (TBD). The flux in the original sampling can be found in `FLUXTBL` HDU or the outputs in the ASCII table.

Note 2: Currently, the format of data table described above is still discussed in the DRP team and may be changed in the future release.

Note 3: PFS configuration information, including catalogue ID, object ID, coordinate (dummy value is used), fiber flux (flux density of the input magnitude), MPS centroid parameter (dummy value is used), will be saved in "pfsConfig" file (e.g., "pfsConfig-0x00000001.fits").

#### Realization of multiple spectra
If you have many spectra you want to realize, you can do that in a single run using an input magnitude file (`MAG_FILE`) containing each spectral information (with columns like this: wavelength magnitude1 magnitude2 ... magnitude1000). Then you can get the output file of each spectrum. Please note that `--nrealize=1` when you use this mode.

Usage in your Python code or Jupyter notebook (under development)
---------------------------
These functionality can be used by importing pfsspecsim module in your own Python codes and on Jupyter notebooks like this: 

For calculation of S/N curves:

```python
from pfsspecsim import pfsetc

etc = pfsetc.Etc()
etc.set_param('EXP_TIME', 1200)
etc.set_param('EXP_NUM', 3)
etc.set_param('OUTFILE_NOISE','out/ref.noise.dat')
etc.set_param('OUTFILE_SNC','out/ref.snc.dat')
etc.set_param('OUTFILE_SNL','out/ref.snl.dat')
etc.set_param('OUTFILE_OII','out/ref.snoii.dat')
etc.run()
```

For making simulated spectra,

```
from pfsspecsim import pfsspec

sim = pfsspec.Pfsspec()
sim.set_param('ra', 150.0)
sim.set_param('dec', 2.0)
sim.set_param('etcFile', 'out/ref.snc.dat')
sim.set_param('MAG_FILE', 19.0)
sim.set_param('EXP_NUM',16)
sim.set_param('asciiTable','test')
sim.set_param('nrealize',1)
sim.set_param('plotObject','t')
sim.set_param('plotArmSet','f')
sim.make_sim_spec()
```

See example/notebooks/ETC Example.ipynb for details.
Some examples (under development)
---------------------------
#### Input spectra

There are some examples of input spectra for the ETC. If you make additional example in your own calculation, let the contact person described below know please.

##### An observation of a star-forming galaxy with a continuum and emission lines:

This example is a SDSS galaxy classified as `GALAXY`, whose spectra is redshifted to z=0.8 and scaled to mag=21 ABmag at 1000 nm. For this object with 1 hour exposure time under the seeing of 0.7 arcsec and gray lunar-phase condition, assuming that the object in a fiber on the FoV center and the elevation angle is 60 deg., type the following command:

    python run_etc.py @run_etc.defaults --MAG_FILE=./example/spec/ex_gal_sf.dat --EXP_TIME=900 --EXP_NUM=4 --REFF=0.30 --OUTFILE_NOISE=./out/ex_gal_sf.noise.dat --OUTFILE_SNC=./out/ex_gal_sf.snc.dat --OUTFILE_SNL=- --NOISE_REUSED=N --MR_MODE=N --OVERWRITE=Y --SEEING=0.70 --ZENITH_ANG=30.0 --FIELD_ANG=0.00 --MOON_PHASE=0.25 --MOON_ZENITH_ANG=30.0 --MOON_TARGET_ANG=60.0

    python gen_sim_spec.py @gen_sim_spec.defaults --etcFile=./out/ex_gal_sf.snc.dat --asciiTable=sim.ex_gal_sf --MAG_FILE=./example/spec/ex_gal_sf.dat --EXP_NUM=4 --outDir=out

##### An observation of a star-burst galaxy with (mostly) only emission lines:

This example is a SDSS galaxy classified as `GALAXY STARBUSRT`, whose spectra is redshifted to z=0.8 and scaled to mag=26 ABmag at 1000 nm. For this object with 0.5 hour exposure time under the seeing of 0.5 arcsec and bright lunar-phase condition, assuming that the object in a fiber on the FoV center and the elevation angle is 50 deg., type the following command:

    python run_etc.py @run_etc.defaults --MAG_FILE=./example/spec/ex_gal_sb.dat --EXP_TIME=900 --EXP_NUM=2 --REFF=0.30 --OUTFILE_NOISE=./out/ex_gal_sb.noise.dat --OUTFILE_SNC=./out/ex_gal_sb.snc.dat --OUTFILE_SNL=- --NOISE_REUSED=N --MR_MODE=N --OVERWRITE=Y --SEEING=0.50 --ZENITH_ANG=40.0 --FIELD_ANG=0.00 --MOON_PHASE=0.5 --MOON_ZENITH_ANG=30.0 --MOON_TARGET_ANG=60.0

    python gen_sim_spec.py @gen_sim_spec.defaults --etcFile=./out/ex_gal_sb.snc.dat --asciiTable=sim.ex_gal_sb --MAG_FILE=./example/spec/ex_gal_sb.dat --EXP_NUM=2 --outDir=out

##### An observation of a passive galaxy with a continuum and absorption lines:

This example is generated by using a CB07 stellar population synthesis model with Chabrier IMF, solar abundance, and single burst. The stellar age of the the galaxy is ~5 Gyr and with no dust extinction. The intrinsic spectra is redshifted to z=1.2 and scaled to mag=21 ABmag at 1000 nm. For this object with 5 hour exposure time under the seeing of 0.5 arcsec and bright lunar-phase condition, assuming that the object in a fiber on the FoV center and the elevation angle is 60 deg., type the following command:

    python run_etc.py @run_etc.defaults --MAG_FILE=./example/spec/ex_gal_pv.dat --EXP_TIME=900 --EXP_NUM=20 --REFF=0.30 --OUTFILE_NOISE=./out/ex_gal_pv.noise.dat --OUTFILE_SNC=./out/ex_gal_pv.snc.dat --OUTFILE_SNL=- --NOISE_REUSED=N --MR_MODE=N --OVERWRITE=Y --SEEING=0.50 --ZENITH_ANG=30.0 --FIELD_ANG=0.00 --MOON_PHASE=0.5 --MOON_ZENITH_ANG=30.0 --MOON_TARGET_ANG=60.0

    python gen_sim_spec.py @gen_sim_spec.defaults --etcFile ./out/ex_gal_pv.snc.dat --asciiTable sim.ex_gal_pv.dat --MAG_FILE=./example/spec/ex_gal_pv.dat --EXP_NUM=20 --outDir=out

#### Noise model

Because the calculation of noise models is somewhat time consuming, we prepared some examples of noise model for some cases. In ./example/noise/ directory, there are noise models varying zenith angle from 0 to 70 deg. and field angle from 0.0 to 0.7 deg. for the dark moon phase and 8 exposures of 450 sec.  (e.g.,ref.noise.drk.z30.f04.t450.dat).

Other notes
-----------
1. The other parameters that are implemented in the Chris Hirata's ETC are indicated in the run_etc.py including instrument setups and sky subtraction floor. Users can change these parameters on your own responsibility.

2. The fraction of light covered by fiber aperture may be overestimated by up to 10% depending on the field angle compared to the ray-tracing calculation by using the PFS optical model including the telescope and the wide-field corrector.

3. We assumed no fiber central position offset, additional 1% systematic sky subtraction error (on each 1D pixel), and additional 2% diffuse stray light (the entire 2D detector surface) for the noise model calculation.

4. As a fiducial sky continuum model, we use the recent sky model provided by Jim Gunn, which is basically in
consistent with the observations in the SDSS/BOSS in optical and the result obtained the MOSFIRE engineering runs. We use the sky emission line model taken from UVES visible line atlas and theoretical model (Rousselot et al. 2000, A&A, 354, 1134) in NIR. As a fiducial atmospheric transmission model, we use Kitt Peak model for <900nm and a simulated ATRAN model with 3 mm PWV at longer wavelengths.

5. The original code (gsetc.c) provided by Chris Hirata, which is partially modified by people in Kavli IPMU, can be found in the src directory. The manual for the code (Manual_v5.pdf) may be useful for understanding some assumptions in the noise calculation.

6. The current throughput model includes several actual measurements such as Telescope reflectivity, optical elements in WFC and PFI, SpS optics, dichroic mirrors, VPH gratings, CCD QE. The details about the throughput models are described elsewhere.

Have fun! and your feedback is highly appreciated!
--------------------------------------------------

Contact
-------
Kiyoto Yabe (Kavli IPMU)
e-mail: kiyoto.yabe@ipmu.jp
