#!/usr/bin/env python

"""Script to generate fake MWA observations in JSON files, so that they can be
   plotted with skymap without actually being in the schedule database.

   New, untested, not working yet.
"""

import optparse
import os
import sys

import json

import astropy
from astropy.time import Time, TimeDelta
from astropy.coordinates import Angle, AltAz, EarthLocation, SkyCoord

from mwa_skymap import tile_geometry

usage = """insert usage help here"""
MWAPOS = EarthLocation.from_geodetic(lon="116:40:14.93485", lat="-26:42:11.94986", height=377.827)  # GDA94


######################################################################
# Command line input parsing functions

def next_mod8(gpstime=None, buf=1):
    """Takes a value 'gpstime' in GPS seconds, and calculates the next modulo-8 second time
       at least (buf + 1) seconds after that time. If gpstime is not given, use the current time.

       The value of buffer defaults to 1 second, so if 'buf' is not specified,
       the time returned will be at least 2 seconds (buf + 1) in the future.

       :param gpstime: time in GPS seconds
       :param buf: Extra buffer time, in seconds, defaults to 1
       :return: Next modulo 8 time, in GPS seconds
    """
    return int((gpstime + buf + 1) / 8.0) * 8 + 8


def convert_time(newtime, lasttime):
    """
    timeout = convert_time(newtime, lasttime)

    Converts the timestring given in newtime to a time in GPSseconds, return as timeout
    lasttime is the last time returned, used for increments
    formats for newtime:

      ++dt                 - increments lasttime by dt seconds
      yyyy-mm-dd,hh:mm:ss  - UT date/time
      t                    - GPS seconds


    :param newtime: new time string to parse
    :param lasttime: previous time passed, used for parsing increments
    :return: GPStime

    """

    timeout = newtime

    if (isinstance(newtime, str)):
        if (newtime.startswith('++')):
            try:
                dt_string = newtime.replace('++', '')
                if (dt_string.count('s')):
                    # seconds
                    dt = int(dt_string.replace('s', ''))
                elif (dt_string.count('m')):
                    # minutes
                    dt = int(60 * float(dt_string.replace('m', '')))
                elif (dt_string.count('h')):
                    # hours
                    dt = int(3600 * float(dt_string.replace('h', '')))
                else:
                    # assume seconds
                    dt = int(dt_string)
            except ValueError:
                print('ERROR: Unable to interpret time increment: %s' % newtime)
                return 0
            timeout = next_mod8(gpstime=lasttime + dt - 8)
        elif ':' in newtime:
            try:
                # Accept a seperator betwen date and time of ',' or 'T', as well as a space for true ISO format
                newtime = newtime.replace(',', ' ')
                newtime = newtime.replace('T', ' ')

                newtime = newtime.replace('/', '-')   # Allow a / between year, month, and day, as well as a '-'
                timeout = next_mod8(Time(newtime, format='iso').gps) - 8
            except ValueError:
                print('ERROR: astropy.time.Time unable to interpret ISO timestamp: %s' % newtime)
                return 0
        else:
            try:
                timeout = next_mod8(int(newtime) - 8)
            except ValueError:
                print('ERROR: Unable to interpret GPStime: %s' % newtime)
                return 0

    return timeout


def parse_coord(value='', is_hours=False):
    """
    Parse an RA or Dec option, and return a value in degrees

    :param value: argument to the option, as a string
    :param is_hours: Boolean, true if a sexagesimal value should be interpreted as hours, not degrees
    :return: Value in degrees
    """
    if ':' in value:
        try:
            if is_hours:
                x = Angle(value, unit=astropy.units.hour)
            else:
                x = Angle(value, unit=astropy.units.deg)
        except ValueError:
            raise optparse.OptionValueError("Could not parse %s" % value)
    else:
        try:
            x = Angle(value, unit=astropy.units.deg)
        except ValueError:
            raise optparse.OptionValueError("Could not parse %s" % value)

    return float(x.deg)


def nearest_gridpoint(az, el):
    """
    Given an azimuth and elevation, return a tuple of (az, el) for the nearest grid point (the point on
    the sky where the 16 dipole delay values are all integers).

    :param az: azimuth in degrees
    :param el: elevation in degrees
    :return: A tuple of (az, el) in degrees
    """
    gripoint_number, azimuth, elevation, distance_in_deg = tile_geometry.find_closest_grid_pointing(az, el)
    return azimuth, elevation


def parse_freqs(input_spec='', numchannels=24, separator=";"):
    """
    Parse the frequency specification given on the command line
    that can take the form of::

      <channel>
      <center>,<width>,<increment>
      <start>:<stop>:<increment>

    where the increments default to 1.  Multiple entries can be given separated by <separator>,
    which defaults to ;
    if using ;, entries should be enclosed in quotes

    :param input_spec: input string
    :param numchannels: number of frequency channels (int)
    :param separator: separator between the frequency channels (string)
    :return: list of frequency channels (ints)

    """
    freqs = []
    atoms = input_spec.split(separator)
    for atom in atoms:
        if not (":" in atom or "," in atom or "-" in atom):
            # just a number
            try:
                freqs.append(int(atom))
            except ValueError:
                print('ERROR: Unable to parse frequency channel: %s' % atom)
                return freqs
        elif (":" in atom):
            # assumes <channelstart>:<channelstop>:<increment>
            increment = 1
            res = atom.split(":")
            if (len(res) > 2):
                try:
                    increment = int(res[2])
                except ValueError:
                    print('ERROR: Unable to parse frequency increment: %s' % res[2])
                    return freqs
            try:
                freqstart = int(res[0])
            except ValueError:
                print('ERROR: Unable to parse frequency start: %s' % res[0])
                return freqs
            try:
                freqstop = int(res[1])
            except ValueError:
                print('ERROR: Unable to parse frequency stop: %s' % res[1])
                return freqs
            for freq in range(freqstart, freqstop + 1, increment):
                freqs.append(freq)
        elif ("-" in atom):
            # assumes <channelstart>-<channelstop>-<increment>
            increment = 1
            res = atom.split("-")
            if (len(res) > 2):
                try:
                    increment = int(res[2])
                except ValueError:
                    print('ERROR: Unable to parse frequency increment: %s' % res[2])
                    return freqs
            try:
                freqstart = int(res[0])
            except ValueError:
                print('ERROR: Unable to parse frequency start: %s' % res[0])
                return freqs
            try:
                freqstop = int(res[1])
            except ValueError:
                print('ERROR: Unable to parse frequency stop: %s' % res[1])
                return freqs
            for freq in range(freqstart, freqstop + 1, increment):
                freqs.append(freq)

        elif ("," in atom):
            # assumes <center>,<width>,<increment>
            increment = 1
            res = atom.split(",")
            if (len(res) > 2):
                try:
                    increment = int(res[2])
                except ValueError:
                    print('ERROR: Unable to parse frequency increment: %s' % res[2])
                    return freqs
            try:
                freqcenter = int(res[0])
            except ValueError:
                print('ERROR: Unable to parse frequency center: %s' % res[0])
                return freqs
            try:
                freqwidth = int(res[1])
            except ValueError:
                print('ERROR: Unable to parse frequency width: %s' % res[1])
                return freqs
            for freq in range(freqcenter - int(freqwidth / 2.0), freqcenter + int(freqwidth / 2.0 + 0.5), increment):
                freqs.append(freq)

    # remove duplicates
    origfreqs = freqs
    freqs = list(set(freqs))
    if (len(freqs) < len(origfreqs)):
        print("WARNING: Removed duplicate items from frequency list")
    # sort
    freqs.sort()
    # trim if necessary
    if len(freqs) > numchannels > 0:
        print("WARNING: Removed excess items from frequency list: %s", freqs[numchannels:])
        freqs = freqs[:numchannels]

    if (len(freqs) == 1):
        # only a single frequency
        print('WARNING: --freq=%d requested, but interpreting it as --freq=%d,24' % (freqs[0], freqs[0]))
        freqcenter = freqs[0]
        freqwidth = 24
        increment = 1
        freqs = list(range(freqcenter - int(freqwidth / 2.0), freqcenter + int(freqwidth / 2.0 + 0.5), increment))

    if len(freqs) < numchannels:  # Pad to numchannels
        if freqs[-1] <= 255 - (24 - len(freqs)):
            freqs += list(range(freqs[-1] + 1, freqs[-1] + (24 - len(freqs)) + 1))
        elif freqs[0] > (24 - len(freqs)):
            freqs = list(range(freqs[0] - (24 - len(freqs)), freqs[0], 1)) + freqs
        else:
            freqs += [x for x in range(60, 200) if x not in freqs][:(24 - len(freqs))]
            freqs.sort()

    if (min(freqs) < 0) or (max(freqs) > 255):
        print("Centre channel too close to 0 or 255, some channels would be < 0 or > 255")
        return []

    return freqs


def main():
    parser = optparse.OptionParser(usage=usage, version="%prog")

    parser.add_option("--ldir",
                      dest="dir",
                      default='.',
                      help="Local directory to write <obsid>.json dummy observation files.")
    parser.add_option("--obsname", "--name",
                      dest="obsname",
                      default='',
                      help="Base name of observation. Will have source name and centre frequency channel appended.")
    parser.add_option("--starttime",
                      dest="starttime",
                      default=None,
                      help="""Observation start time in GPSseconds or yyyy-mm-dd,hh:mm:ss, or ++<dt>[s/m/h] from now.""")
    parser.add_option("--stoptime",
                      dest="stoptime",
                      default=None,
                      help="Observation stop time as above, or ++DT[s/m/h] from starttime. If omitted, use starttime + (nsources*nfreqs*exptime)")
    parser.add_option("--utdate", "--date",
                      dest="utdate",
                      default=None,
                      help="UT Date (YYYY-MM-DD) to schedule observation - must be used with the --lst parameter, and WITHOUT specifying --starttime")
    parser.add_option("--lst", "--lsthour",
                      dest="lst",
                      type='string',
                      help="LST in hours to schedule observation - must be used with the --date parameter, and WITHOUT specifying --starttime")
    parser.add_option("--source",
                      dest="source",
                      default=None,
                      help="Name of the source to track. It will be looked up online using the CDS name resolver.")
    parser.add_option("--ra",
                      dest='ra',
                      type='string',
                      help="Source to track by RA,Dec (hh:mm:ss or decimal degrees) in ICRS reference frame.")
    parser.add_option("--dec",
                      dest='dec',
                      type='string',
                      help="Source to track by RA,Dec (dd:mm:ss or decimal degrees) in ICRS reference frame.")
    parser.add_option("--azimuth", "--az",
                      dest='az',
                      default=None,
                      type='string',
                      help="Pointing direction in Az,El (dd:mm:ss or decimal degrees).")
    parser.add_option("--elevation", "--altitude", "--el", "--alt",
                      dest='alt',
                      default=None,
                      type='string',
                      help="Pointing direction in Az,El (dd:mm:ss or decimal degrees).")
    parser.add_option("--rfstream",
                      dest='rfstream_number',
                      default=0,
                      type="int",
                      help="RF stream number [default=%default], pass value > 0 to add a subarray to a previously saved observation.")
    parser.add_option("--voltbeam",
                      dest='voltbeam_number',
                      default=0,
                      type="int",
                      help="Voltage number [default=%default], pass value > 0 to add a voltage beam to a previously saved observation.")
    parser.add_option("--freq", "--frequencies", "--frequency",
                      dest="freq",
                      default="121,24",
                      help="Frequency, some of <channel>, <channelstart>:<channelstop>:[<increment>], " +
                           "<centerchannel>,<channelwidth>,[<increment>], separated by ; " +
                           "(channels are 1.28 MHz coarse channel numbers).  [default=121,24]")

    starttimelast = int(Time.now().gps / 8.0) * 8 + 8

    (options, args) = parser.parse_args()

    lst = None
    if options.lst is not None:
        lst = parse_coord(value=options.lst, is_hours=True)

    if options.utdate is not None:
        try:
            # Convert to datetime struct, extract just the date at 0h UTC, then convert back to an astropy Time object
            utdate = Time(Time(options.utdate).datetime.date().isoformat(), location=MWAPOS)
        except:
            print('ERROR: Invalid utdate parameter - must be YYYY-MM-DD')
            sys.exit(-1)

        if (lst is None) or (options.starttime is not None):
            print('ERROR: Must provide both utdate AND lst, or starttime')
            sys.exit(-1)

        utc_starthour = (lst - utdate.sidereal_time('apparent').hour) * 23.93447 / 24.0
        if utc_starthour < 0:
            utc_starthour += 23.93447

        starttime = (utdate + TimeDelta(utc_starthour * 3600, format='sec')).gps
        starttime = 8 * int(starttime / 8)

    elif lst is not None:
        print('ERROR: Must provide both utdate AND lst, or starttime')
        sys.exit(-1)

    else:
        starttime = convert_time(options.starttime, starttimelast)

    if (not starttime):
        print("Can't determine start time in GPS seconds")
        sys.exit(-1)

    stoptime = convert_time(options.stoptime, starttime)
    obs_midpoint = Time((starttime + stoptime) / 2.0, format='gps', scale='utc')
    altaz_frame = AltAz(obstime=obs_midpoint, location=MWAPOS)

    if options.source is not None:
        source = SkyCoord.from_name(options.source)
        ra, dec = float(source.ra.value), float(source.dec.value)
        app_pos = source.transform_to(altaz_frame)
        exact_alt, exact_az = float(app_pos.alt.value), float(app_pos.az.value)
        az, alt = nearest_gridpoint(az=exact_az, el=exact_alt)
    elif options.ra is not None and options.dec is not None:
        ra = parse_coord(value=options.ra, is_hours=True)
        dec = parse_coord(value=options.dec, is_hours=False)
        body = SkyCoord(ra=ra / 15.0 * astropy.units.hour, dec=dec * astropy.units.degree)
        app_pos = body.transform_to(altaz_frame)
        exact_alt, exact_az = float(app_pos.alt.value), float(app_pos.az.value)
        az, alt = nearest_gridpoint(az=exact_az, el=exact_alt)
    elif options.alt is not None and options.az is not None:
        exact_alt = parse_coord(value=options.alt, is_hours=False)
        exact_az = parse_coord(value=options.az, is_hours=False)
        aacoord = SkyCoord(alt=exact_alt * astropy.units.degree,
                           az=exact_az * astropy.units.degree,
                           frame='altaz',
                           obstime=obs_midpoint,
                           location=MWAPOS,
                           unit=astropy.units.degree)
        coord = aacoord.transform_to('icrs')
        coord.obstime = obs_midpoint
        ra, dec = float(coord.ra.to(astropy.units.degree).value), float(coord.dec.to(astropy.units.degree).value)
        az, alt = nearest_gridpoint(az=exact_az, el=exact_alt)
    else:
        print("Must specify either source, ra/dec or alt/az")
        sys.exit(-1)

    if alt < 15:
        print("Altitude must be above 15 degrees")
        sys.exit(-1)

    freq_list = parse_freqs(input_spec=options.freq)

    filename = '%s/%d.json' % (options.dir, starttime)

    if options.rfstream_number > 0 and options.voltbeam_number == 0:
        if not os.path.exists(filename):
            print('Create obsid=%d first, then add subarray using --rfstream=%d' % (starttime, options.rfstream_number))
            sys.exit(-1)
        print('Adding subarray rfstream %d to existing obsid=%d' % (options.rfstream_number, starttime))
        f = open(filename, 'r')
        obsinfo = json.load(f)
        f.close()
        print('Adding new subarray number %d')
        obsinfo['rfstreams'][str(options.rfstream_number)] = {'az':az, 'el':alt, 'frequencies': freq_list}
    elif options.rfstream_number == 0 and options.voltbeam_number > 0:
        if not os.path.exists(filename):
            print('Create obsid=%d first, then add voltage beam using --voltbeam=%d' % (starttime, options.voltbeam_number))
            sys.exit(-1)
        print('Adding voltage beam %d to existing obsid=%d' % (options.voltbeam_number, starttime))
        f = open(filename, 'r')
        obsinfo = json.load(f)
        f.close()
        obsinfo['voltage_beams'][str(options.voltbeam_number)] = {'ra': ra, 'dec': dec, 'target_name': options.obsname}
    elif options.rfstream_number > 0 and options.voltbeam_number > 0:
        print('Can only specify one of rfstream or voltbeam')
        sys.exit(-1)
    else:
        if os.path.exists(filename):
            print('Overwriting obsid=%d' % starttime)
        obsinfo = {'starttime': starttime,
                   'stoptime': stoptime,
                   'obsname': options.obsname,
                   'rfstreams': {'0': {'az':az, 'el':alt, 'frequencies': freq_list}},
                   'ra_phase_center': ra,
                   'dec_phase_center': dec,
                   'voltage_beams': {}}

    obs_json = json.dumps(obsinfo, indent=4)

    f = open(filename, 'w')
    f.write(obs_json)
    f.close()
    print(obs_json)


if __name__ == '__main__':
    main()