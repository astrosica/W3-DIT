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
