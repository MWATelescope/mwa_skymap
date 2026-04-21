# mwa_skymap
All-sky maps and movies showing the MWA telescope primary beam

This package uses the  mwa_hyperbeam package (or, optionally, the old mwa_pb
package) to generate all-sky maps of one or more MWA observations, showing the 
primary beam on the sky, with the HASLAM radio image as a background.

It consists of a library (mwaplot.py) implementing the core functionality, and a
command-line tool (skymap, implemented in mwa_skymap.py) to generate all-sky 
maps and movies.

Features include:
  * Generate still frames showing the primary beam and phase centre at a single instant.
  * Generate movies showing the primary beam and phase centre moving across the sky,
    repointing as the observations progress between the given start and stop times.
  * Single frame output formats include GIF, JPEG, PNG, etc, and movie formats can be as an
    animated PNG (with variable length frames for more precise repointing times),
    or an MPEG (.mpg) file. Output formats are determined by the file extension.
  * If an MWA observation includes voltage beams, they will be plotted as well 
    as the phase centre.
  * If an observation includes primary beam subarrays (multiple sets of tiles, 
    each with a different pointing centre), they will be plotted as well, with
    a different contour color for each subarray.
  * You can optionally add the GLEAM catalog of radio sources to the map.
  * You can plot the HASLAM background as white-on-black (the default), or 
    in inverse mode, as black-on-white. Source markers and labels are automatically
    color corrected to stand out in either mode.
  * You can specify the primary beam calculation to use - one of 
    'HBA' (mwa_hyperbeam analytic), 'HBFEE' (mwa_hyperbeam FEE - slow without a GPU),
    or 'MWA_PB' (mwa_pb analytic).
  * You can specify a specific coarse channel to use for the beam model (defaults to 
    the 13th channel in the observation).

# Getting started:

* Create and activate a Python environment.
* Install the package:
  ```
  git clone https://github.com/MWATelescope/mwa_skymap.git
  cd mwa_skymap
  pip install .
  ```

# (Optional) Download the MWA primary beam data (to use the HBFEE model)

* Download the MWA primary beam data:
  ```
  skymap getbeamfile
  ```

# Command line usage

The skymap command has three subcommands:

* `skymap single` - Generate a single frame of the primary beam
* `skymap movie` - Generate a movie of the primary beam
* `skymap getbeamfile` - Download the MWA primary beam data
* `skymap beamtypes` - print the list of primary beam models available

You can get help on the command line options with:

```
skymap --help
```

and, for example:

```
skymap single --help
```

# Usage examples:

```
skymap single 1441333936  
skymap single --inverse 1441333936   # THe same observation as black-on-white
skymap single --gleamsources --cchan=57 1441333936  # Use channel 57, and show GLEAM sources

# Or an entire 4.4 hour block of EoR observations and calibrators, at 10 frames per second (each
# showing one minute of actual observing time)

skymap movie --outfile=20250908.mpg --startgps=1441386096 --stopgps=1441401784 --inverse
```
