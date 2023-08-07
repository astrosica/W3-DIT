# W3 DIT

This code performs data reduction and calibration for DIT Photometry of W3 and is written by Jessica Campbell.

David A. Dunlap Department of Astronomy & Astrophysics, Canadian Institute for Theoretical Astrophysics, University of Toronto

## Context
The Dunlap Institute Telescope (DIT) is a 1m telescope that observed the W3 giant molecular cloud from Mexico. Observing was led by Suresh Sivanandam and Nick Law. The DIT data was obtained using a 3x3 grid of pointings across the W3 field with a 10th pointing centered on KR 140. The data is organized into folders named YYYYMMDD and include the riz passbands. Bias measurements, dark fields, and flat fields were obtained for the use of removing instrumental artifacts through a process called data reduction.

## Data Reduction

The data reduction steps used here are as follows:

### 1. Quality Check
Flag problematic images that are not to be used (e.g., star trails in dark images) and create a mask for bad pixels (e.g., very negative pixels).

### 2. Astrometric Solution
Use [Astrometry.net](astrometry.net) to obtain an astrometric solution for each science image. The plate solution was originally done with PinPoint (via header) but this was found to be inaccurate.

### 3. Bias Subtraction
Create a master bias frame by averaging all bias frames across all nights; they all have the same exposure (0s) and do not change with wavelength. Bias measurements create an offset (~100s of counts) measured by the CCD before the start of the exposure that needs to be subtracted from all images. 

### 4. Dark subtraction
Make master dark frames for each integration time by averaging across all nights for a given exposure time and subtract from science images with corresponding exposure time (either 5s or 120s). The dark frame subtracted from a science image should have the same exposure time. Dark measurements were made in sets of 5 and 10 frames over a wide range in exposure time (10s, 30s, 60s, 120s, 300s). Check that dark counts are linear with exposure time. The dark-subtraction accounts for the bias (i.e., bias is included in the dark frame). 

### 5. Flat Fielding
Create a master flat frame for each passband by average dark-subtracted flat field images across all nights for each passband; normalize once more at the end by dividing by the median value. The flat field images should be dark subtracted using a dark frame with the same exposure time. The flat field measures the CCD response/sensitivity; this is not uniform and varies across the detector and as a function of wavelength. A flat field measurement is usually taken at dusk and/or dawn as it requires a high photon count (~60-70% saturation). 

### 6. Calibration
Use Source Extractor to obtain tabulated photometric measurements of stars in calibrated science images and cross-match to Pan-STARRS catalog. Fit a straight line to (Source Extractor - Pan-STARRS) vs. Pan-STARRS magnitudes to measure the zero point and slope in the data. Only use a science image if the slope is equal to zero within 3-sigma uncertainties. Verify that the zero-point is zero post-calibration.

### 6. Stack Images
Create final science images by stacking (e.g., median averaging) reduced science images. Science images are dithered so they are projected onto a common grid first using a simple WCS grid. This common grid is the same for all passbands for a given pointing and helps with things like cosmic rays. When stacking science images, measures the standard deviation to obtain a measure of noise across pixels. 

Note: The "dupe" in several file names does not necessarily mean that they are in fact duplicates because the observation time and header information changes.

## Usage

Requires that images are organized in the following way:
raw science images: <parent directory>/<date>/<images>/
bias images: <parent directory>/<date>/<calibration folder>/<images>/
dark images: <parent directory>/<date>/<calibration folder>/<images>/
flat images: <parent directory>/<date>/<flat folder>/<images>/

Requires that filenames are named in the following way:
raw science images: objectname_[]-[]-[].fts
bias images: Bias-[]-[]-[]-.fts
dark images: Dark-[]-[]-[].fts
flat images: AutoFlat-[]-[]-[].fts

The following is done with the raw science images:
1. dark subtracted and flat fielded
2. astronometry.net is used to obtain accurate world coordinate system
3. Source Extractor is used for photometry

Outputs the following directories
- :

Requires the following header keys in the science images' FITS files:
- OBJECT
- FILTER
- EXPTIME

Requires the following software:
    - sextractor
    - astrometry.net
Requires the following packages: os, re, yaml, docopt, numpy, astropy
Requires the following files:
    - calibration.py
    - callastrometry.py

Usage:
main [options]

Options:
  -h, --help           Show this screen.
                       [default: False]
  -l, --logoutput      Generate log outputs.
                       [default: True]
  -v, --verbose        Generate verbose output.
                       [default: True]
  -g, --generate       Generate outputs.
                       [default: True]
  -w, --wcs            Create WCS solution via Astronomy.net.
                       [default: False]
  -r, --reduce         Perform data reduction.
                       [default: False]
  -c, --calibrate      Perform zero-point calibration.
                       [default: False]
  -p, --panstarrs      Use Pan-STARRS for calibration.
                       [default: False]
  -a, --atlas          Use ATLAS for calibration.
                       [default: True]
  -s, --sextractor     Run Source Extractor for calibration step.
                       [default: False]
  -k, --checkcal       Check the zero-point calibration.
                       [default: False]
  -x, --stack          Stack calibrated images.
                       [default: False]
  -f, --finalcalcheck  Check the zero-point calibration on stacked images.
                       [default: False]
  -n, --finalsexrun    Final SExtractor run on stacked images.
                       [default: False]
  -m, --matchcats      Match source catalogues.
                       [default: False]
