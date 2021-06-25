# DIT Photometry Data Reduction and Calibration

This repository contains my code for reducing and calibrating the DIT photometric observations of W3.

## About DIT
The Dunlap Institute Telescope (DIT) is a 1 m telescope that observed W3 from Mexico. Observing was led by Suresh Sivanandam and Nick Law. Natalie Price-Jones and Jielai Zhang reduced and calibrated the Dragonfly photometric data which some of this is based on.

The DIT data is hosted on the CITA server at: /mnt/raid-project/hp/campbell/DIT/DITdata

The DIT data was obtained using a 3x3 grid of pointings across the W3 field with a 10th pointing centered on KR140. The data is organized into folders named YYYYMMDD. The plate solution was originally done with PinPoint (via header) but this is not accurate and will need to be done again. The "dupe" in several file names does not necessarily mean that they are in fact duplicates because the observation time and header information changes (they might have different exposure times).

## Data Reduction

A useful resource for this is Handbook of CCD Astronomy by Steve B. Howell. 

### 1. Quality Check
Visually inspect images and flag those that will not be used (e.g., some dark images have star trails).

### 2. Astrometric Solution
Use [Astrometry.net](astrometry.net) to obtain an astrometric solution for each science image.

### 3. Bias Subtraction
Bias creates an offset (~100s of counts) measured by the CCD before the start of th expoure that needs to be subtracted from all images. Bias measurements were made using either 10 or 20 x 0s exposures.

### 4. Dark subtraction
Dark measurements were made using either 5 or 10 x frames over a wide range in exposure time (10s, 30s, 60s, 120s, 300s). Check that dark counts are linear with exposure time.

### 5. Flat Fielding
The flat field measures the CCD response/sensitivity; this is not uniform and varies across the detector and as a function of wavelength. A flat field measurement is usually taken at dusk and/or dawn as it requires a high photon count (~60-70% saturation).

Normalize each flat by dividing by the median value.

### 6. Calibration

### 7. Stack (e.g., median average) reduced and calibrated science images.
This will help with things like cosmic rays.

## Calibration
We use Pan-STARRS data to calibrate the DIT photometry.

![alt text](figures/zp_plot.png)
