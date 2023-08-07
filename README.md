# W3 DIT

This Python pipeline performs data reduction and calibration for DIT Photometry of W3.

David A. Dunlap Department of Astronomy & Astrophysics, Canadian Institute for Theoretical Astrophysics, University of Toronto

## About DIT
The Dunlap Institute Telescope (DIT) is a 1 m telescope that observed W3 from Mexico. Observing was led by Suresh Sivanandam and Nick Law. The DIT data was obtained using a 3x3 grid of pointings across the W3 field with a 10th pointing centered on KR 140. The data is organized into folders named YYYYMMDD. The plate solution was originally done with PinPoint (via header) but this is not accurate and was done again. The "dupe" in several file names does not necessarily mean that they are in fact duplicates because the observation time and header information changes (they might have different exposure times).

## Data Reduction

The data reductions steps are as follows:

### 1. Quality Check
Visually inspect images and flag those that will not be used (e.g., lots of dark images have star trails). Create mask for bad pixels (e.g., very negative pixels).

### 2. Astrometric Solution
Use [Astrometry.net](astrometry.net) to obtain an astrometric solution for each science image.

### 3. Bias Subtraction
Bias creates an offset (~100s of counts) measured by the CCD before the start of th expoure that needs to be subtracted from all images. Bias measurements were made using either 10 or 20 x 0s exposures. Can (median) average all bias frames across all nights since they all have the same exposure and do not change with wavelength.

### 4. Dark subtraction
Dark measurements were made using either 5 or 10 x frames over a wide range in exposure time (10s, 30s, 60s, 120s, 300s). Check that dark counts are linear with exposure time. A dark-subtraction should account for the bias as well (i.e., bias is included in the dark frame). The dark frame subtracted from a science image should have the same exposure time, so should make separate master dark frames for each integration time used in science and flat field images by (median) averaging across all nights for a given exposure time. Most science images have either 5s or 120s exposure times, other dark frame exposure times might have been for the flat frames.

### 5. Flat Fielding
The flat field measures the CCD response/sensitivity; this is not uniform and varies across the detector and as a function of wavelength. A flat field measurement is usually taken at dusk and/or dawn as it requires a high photon count (~60-70% saturation). Flat field images should be dark subtracted using a dark frame with the same exposure time. Will need a master flat frame for each passband used in science images, so create these by (median) average dark-subtracted flat field images over all nights for each passband (may need to normalize once more at the end).

Normalize each flat by dividing by the median value.

### 6. Calibration

### 7. Stack Images
Create final science images by stacking (e.g., median averaging) reduced and calibrated science images. Science images are dithered so they will need to be projected onto a common grid first (use a simple WCS grid). This common grid should be the same for all passbands for a given pointing. This will help with things like cosmic rays. When stacking science images, can measure standard deviation to obtain a measure of noise across pixels.

## Calibration
We use Pan-STARRS data to calibrate the DIT photometry.
