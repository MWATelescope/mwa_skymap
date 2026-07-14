# mwa_skymap
All-sky maps and movies showing the MWA telescope primary beam, based on the 
'plot_MWAconstellations.py' script by David Kaplan.

This package uses the  mwa_hyperbeam package (or, optionally, the old mwa_pb
package) to generate all-sky maps of one or more MWA observations, showing the 
primary beam on the sky, with the HASLAM radio image as a background.

It consists of a library (mwaplot.py) implementing the core functionality, and a
command-line wrapper script (skymap, implemented in skymap.py) to generate all-sky 
maps and movies. You can use the 'mwaplot.py' library from your own code if you want to 
customise the plots (provide your own source lists, for example).

Features include:
  * Generate still frames showing the primary beam and phase centre at a single instant.
  * Generate movies showing the primary beam and phase centre moving across the sky,
    repointing as the observations progress between the given start and stop times.
  * Single frame output formats include GIF, JPEG, PNG, etc, and movie formats can be as an
    animated PNG (with variable length frames for more precise repointing times),
    or an MPEG (.mpg) file in H.264. Output formats are determined by the file extension.
  * If an MWA observation includes real-time beamformer voltage beams, they will be plotted 
    as well as the phase centre.
  * If an observation includes primary beam subarrays (multiple sets of tiles, 
    each with a different pointing centre), they will be plotted with
    a different contour colour for each subarray.
  * You can optionally add the GLEAM catalogue of radio sources to the map.
  * You can plot the HASLAM background as white-on-black (the default), or 
    in inverse mode, as black-on-white. Source markers and labels are automatically
    colour corrected to stand out in either mode.
  * You can specify the primary beam calculation to use - one of 
    'HBA' (mwa_hyperbeam analytic), 'HBFEE' (mwa_hyperbeam FEE - slow without a GPU),
    or 'MWA_PB' (mwa_pb analytic).
  * You can specify a specific coarse channel to use for the beam model (defaults to 
    the 13th channel in the observation).
  * Choice of white, black or transparent background.
  * Customise (or remove) the text shown in the plot by providing a template with 
    optional fields: '%(obsid)d', '%(viewgps_gps)d', '%(viewgps_utc)s', '%(obsname)s', 
    '%(freq_mhz)s MHz', and '%(constellation)s'. These fields will be filled in 
    with the observation ID, GPS time in GPS seconds, GPS time in UTC, observation name, 
    frequency, and constellation, respectively. When generating movies, the template fields
    will be filled in for each frame, changing for each new observation displayed.

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
```
<img src="https://ws.mwatelescope.org/plots/1441333936.png" width="512" height="512">

The same observation as black-on-white JPEG, using the hyperbeam FEE beam model:
```
skymap single --inverse 1441333936 --beam_type=HBFEE --background=white --outfile=1441333936-inverse.jpg
```
<img src="https://ws.mwatelescope.org/plots/1441333936-inverse.jpg" width="512" height="512">

Use channel 57, and show GLEAM sources:
```
skymap single --gleamsources --cchan=57 1441333936 --outfile=1441333936-57-gleam.png
```
<img src="https://ws.mwatelescope.org/plots/1441333936-57-gleam.png" width="512" height="512">

Showing an observation with four subarray primary beams in different colours:
```
skymap single 1377994280
```
<img src="https://ws.mwatelescope.org/plots/1377994280.png" width="512" height="512">

A 4.4-hour block of EoR observations and calibrators, at 10 frames per second (each
fram howing one minute of actual observing time):
```
skymap movie --outfile=20250908.mpg --startgps=1441386096 --stopgps=1441401784 --inverse
```
[View here](https://ws.mwatelescope.org/plots/20250908.mpg)

# Local observation files

By default, skymap will use the MWA web services to look for observations in the 
MWA observation database. You can also use local observation files by setting the 
`--ldir` option to a directory containing dummy observation files in JSON format.

These files have the form `<obsid>.json` and are created by calling the `local_obs` 
command, a (very simplified) version of the same command used to add real
MWA observations to the schedule database. 

Note that like real MWA observations, the phase centre will be the actual RA/Dec provided,
but the primary beam will be centered on the nearest sweet-spot grid position, and dipole
delays calcualted accordingly.

## usage

`local_obs` generates JSON files with simulated MWA telescope observations. Parameters include:

  * --ldir: Directory to write the JSON file describing the observation.
  * Either --starttime= or both --utdate= and --lst= to sepcify the observation start
  * --stoptime= to specify the observation end
  * One of --source=, or --ra= and --dec=, or --alt= and --az=
  * --freq= with a channel specifier (defaults to 121,24)
  * --obsname: Observation name.

To define the start of the observation, you can use --startime= with an argument
that's one of:
  * An integer in GPS seconds
  * A date/time in the form 'yyyy-mm-dd,hh:mm:ss'
  * A modifier, e.g. ++10s, ++1m, ++1h, relative to the current time.

Or, you can use --utdate=yyyy-mm-dd to specify the UTC date, and --lst= to specify
the Local Sidereal Time, in hours, on that date.

To define the end of an observation, you use --stoptime= with an argument that's
one of:
  * An integer in GPS seconds
  * A date/time in the form 'yyyy-mm-dd,hh:mm:ss'
  * A modifier, e.g. ++10s, ++1m, ++1h, relative to the observation start time.
    
To define the pointing, you can:
  * Specify a source name with --source= to be resolved using 
      astropy.coordinates.SkyCoord.from_name
  * Specify --ra= and --dec= using either floats (in degrees), or sexagesimal values 
    (DD:MM:SS). Sexagesimal values for RA will be interpreted as hours, not degrees.
  * Specify --alt= and --az= using floats (in degrees), or sexagesimal values (DD:MM:SS) 
    in degrees.
        
The output will be a JSON file of the form <obsid>.json in the current directory, 
or in the directory specified by the --ldir= option.

To add a (simulated) subarray pointing, call %prog once with the arguments
for the observation, then call it again with exactly the same arguments
except for a new pointing (ra/dec, alt/az or sourcename), and --rfstream=1 to
add a subarray with the new pointing. You can repeat with --rfstream=2, etc.

To add a (simulated) real-time voltage beam, call %prog once with the arguments
for the observation, then call it again with exactly the same arguments
except for a new pointing (ra/dec, alt/az or sourcename), and --voltbeam=1 to
add voltage beam with the new pointing. You can repeat with --voltbeam=2, etc.

Example: Schedule a 32-second observation of HerA starting in 8 seconds time:

`local_obs` --starttime=++8 --stoptime=++32 --source=HerA --freq=121,24

A 296s observation of CenA at a specific time, written to `dummy` directory:

`local_obs` --ldir=dummy --source=CenA --stoptime=++296 --starttime="2026-07-14,10:00:00 --obsname=cenatest --freq=145,24

A 120s observation of PicA on August 18th 2026, at an LST of 5.5 hours:

`local_obs` --ldir=dummy --source=PicA --stoptime=++120 --utdate=2026-08-18 --lst=5.5 --obsname=picatest --freq=145,24

A 296s observation of HydA starting 32 seconds from now:

`local_obs` --ra='09:18:06' --dec='-12:05:44' --stoptime=++296 --starttime=++32 --obsname=radectest --freq=145,24

Once you've generated local JSON files, you can use them with skymap:

Example: generate a skymap for the observation defined in dummy/1441333936.json

`skymap single --ldir=dummy 1441333936`

If the directory specified by --ldir contains more than one observation, you can pass a list of obsids or use
--startgps and --stopgps to define the observation time range for a 'skymap movie ...' command.