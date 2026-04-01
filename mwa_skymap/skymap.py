#!/usr/bin/python

import io
import logging
import math
import os
import sys
import traceback
import warnings

import apng

import imageio.v2 as iio

warnings.filterwarnings("ignore")

logging.basicConfig()
DEFAULTLOGGER = logging.getLogger()

import numpy

import ephem   # Only used to look up what constellation a given ra/dec is in

import astropy
from astropy.coordinates import EarthLocation, get_body, Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.time import Time

import matplotlib

if 'matplotlib.backends' not in sys.modules:
    matplotlib.use('agg')
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

# This will work provided you don't install it in some wierd way that leaves the
# files in a zip-file or other compressed format. If you want to do that, you can
# work out how to do this bit, because I've fought importlib.resources long enough,
# and gave up in horror at the sheer volume of overly complicated, badly documented
# crap in you need to go through.

DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
GLEAMCAT_FILE = os.path.join(DATA_DIR, 'G4Jy_catalogue_allEGCcolumns.fits')
RADIO_IMAGE_FILE = os.path.join(DATA_DIR, 'radio408.RaDec.fits')   # Haslam image:
RADIO_IMAGE_LOW = 100     # Value in RADIO_IMAGE_FILE below which is considered black
RADIO_IMAGE_HIGH = 10000  # Value in RADIO_IMAGE_FILE above which is considered white
if os.getenv('MWA_BEAM_FILE'):
    MWA_BEAM_FILE = os.getenv('MWA_BEAM_FILE')
elif os.path.exists(os.path.join(DATA_DIR, 'mwa_full_embedded_element_pattern.h5')):
    MWA_BEAM_FILE = os.path.join(DATA_DIR, 'mwa_full_embedded_element_pattern.h5')
else:
    MWA_BEAM_FILE = None
MWAPOS = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)

# noinspection PyUnresolvedReferences
CM = plt.cm.gray
CMI = CM.reversed()
DEFAULT_PLOTSIZE = 1200
FIGSIZE = 8
DPI = 150

SKYDATA = None   # Will be an instance of SkyData() after the first use

DEFAULT_PLOT_TEXT = ("Obs ID %(obsid)d at %(viewgps_utc)s:\n" +
                     "%(obsname)s at %(freq_mhz)d MHz\n" +
                     "in the constellation %(constellation)s'")

# Used for multiple primary beams - the first beam is green contours, the second is cyan, etc
PBCOLORS = [(0.0, 1.0, 0.0),  # green
            (0.0, 1.0, 1.0),  # cyan
            (0.0, 0.0, 1.0),  # blue
            (1.0, 0.0, 1.0),  # magenta
            (1.0, 0.0, 0.0),  # red
            (1.0, 1.0, 0.0)]  # yellow

BEAMS = {}
try:
    import mwa_hyperbeam
except:
    mwa_hyperbeam = None

if mwa_hyperbeam:
    try:
        BEAMS['HBA'] = mwa_hyperbeam.AnalyticBeam(rts_behaviour=False)
    except:
        pass

    try:
        if MWA_BEAM_FILE:
            BEAMS['HBFEE'] = mwa_hyperbeam.FEEBeam(MWA_BEAM_FILE)
        else:
            BEAMS['HBFEE'] = mwa_hyperbeam.FEEBeam()
    except:
        print('Failed to load the hyperbeam FEE beam - you probably need to run "plotmwa getbeamfile" first.')

try:
    from mwa_pb import primary_beam
    BEAMS['MWA_PB'] = primary_beam
except ImportError:
    primary_beam = None


class Source():
    """
    Represents a bright source to plot on the sky map, including the dot color, font size, and label alignment code
    """
    def __init__(self, name, ra, dec, color='k', fontsize=8, align='l'):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.color = color
        self.fontsize = fontsize
        self.align = {'l':'left', 'c':'center', 'r':'right'}[align]


class SkyData(object):
    """
    Class to hold data common to all skymaps (background radio image, GLEAM sources, and a list of bright radio sources)

    An instance of this class is Cached in the SKYDATA global in this module after the first call to plot_MWA_skymap().

    To customise this data (eg, change the list of bright sources), from your code, you can preload and then alter its contents:

    ----------------
    import skymap

    ...
    skymap.SKYDATA = SkyData()
    skymap.SKYDATA.sources['SCP'] = skymap.Source(name='South CP', ra='00:00:00.0', dec='+90:00:00.0')
    ...
    ----------------
    """
    def __init__(self, logger=DEFAULTLOGGER):
        self.valid = True

        # Read the GLEAM source list
        try:
            fi = fits.open(GLEAMCAT_FILE)
            self.gleamcat = []
            for i in range(len(fi[1].data.Name)):
                name, ra, dec, flux = fi[1].data.Name[i], fi[1].data.RAJ2000[i], fi[1].data.DEJ2000[i], fi[1].data.int_flux_151[i]
                self.gleamcat.append((name, ra, dec, flux))
            self.gleamcat.sort(key=lambda x:-x[3])  # Sort from brightest to dimmest
            fi.close()
        except:
            logger.error('Could not find GLEAM data')
            self.valid = False

        try:
            self.radio_image = fits.open(RADIO_IMAGE_FILE)

            rih = self.radio_image[0].header
            # lower case skymapra and skymapdec are one-dimensional arrays with every valid RA or every valid Dec
            self.skymapra = (rih.get('CRVAL1') + (numpy.arange(1, self.radio_image[0].data[0].shape[1] + 1) - rih.get('CRPIX1')) * rih.get('CDELT1')) / 15.0
            self.skymapdec = rih.get('CRVAL2') + (numpy.arange(1, self.radio_image[0].data[0].shape[0] + 1) - rih.get('CRPIX2')) * rih.get('CDELT2')
            self.skymapdec = self.skymapdec[self.skymapdec < (90.0 - MWAPOS.lat.deg)]   # Remove all Dec values never visible from our latitude
            # Capitalised skymapRA and skymapDec are each two-dimensional arrays of either RA or Dec, that together cover every grid point
            self.skymapRA, self.skymapDec = numpy.meshgrid(self.skymapra[::2] * 15, self.skymapdec[::2])

            # noinspection PyUnresolvedReferences
            self.mapgrid = SkyCoord(ra=self.skymapRA,  # Pass in 2D grids of RA and Dec
                                    dec=self.skymapDec,
                                    equinox='J2000',
                                    unit=(astropy.units.deg, astropy.units.deg))
            self.mapgrid.location = MWAPOS
        except:
            logger.error('Cannot open Haslam image')
            traceback.print_exc()
            self.valid = False

        self.sources = {
            'HydA': Source(name='Hyd A', ra='09:18:05.65', dec='-12:05:43.9', color='b', fontsize=10, align='l'),
            'For A': Source(name='For A (double)', ra='03:22:41.52', dec='-37:12:33.5', align='l'),
            'PicA': Source(name='Pic A', ra='05:19:49.73', dec='-45:46:43.7', color='b', fontsize=10, align='l'),
            'EOR0': Source(name='EoR0', ra='00:00:00', dec='-27:00:00', color='w', fontsize=12, align='c'),
            'EOR1': Source(name='EoR1', ra='04:00:00', dec='-30:00:00', color='b', fontsize=12, align='c'),
            'EOR2': Source(name='EoR2', ra='10:20:00', dec='-10:00:00', color='b', fontsize=12, align='c'),
            'PupA': Source(name='Pup A\n(resolved)', ra='08:24:07', dec='-42:59:48'),
            '3C161': Source(name='3C 161', ra='06:27:10.09', dec='-05:53:04.7', align='r'),
            'M42': Source(name='M42/Orion', ra='05:35:17.3', dec='-05:23:28'),
            'CasA': Source(name='Cas A', ra='23:23:24', dec='+58:48:54'),
            'CygA': Source(name='Cyg A', ra='19:59:28.36', dec='+40:44:02.1'),
            '3C444': Source(name='3C 444', ra='22:14:25.75', dec='-17:01:36.3'),
            'PKS0408': Source(name='PKS 0408-65', ra='04:08:20.37884', dec='-65:45:09.0806'),
            'PKS0410': Source(name='PKS 0410-75', ra='04:08:48.4924', dec='-75:07:19.327'),
            'LMC': Source(name='LMC', ra='05:23:34.6', dec='-69:45:22'),
            'PKS2104': Source(name='PKS 2104-25', ra='21:07:25.7', dec='-25:25:46'),
            'PKS2153': Source(name='PKS 2153-69', ra='21:57:05.98061', dec='-69:41:23.6855'),
            'PKS 1932': Source(name='PKS 1932-46', ra='19:35:56.5', dec='-46:20:41', color='w'),
            'PKS1814': Source(name='PKS 1814-63', ra='18:19:35.00241', dec='-63:45:48.1926'),
            'PKS1610': Source(name='PKS 1610-60', ra='16:15:03.864', dec='-60:54:26.14', color='w'),
            'CenB': Source(name='Cen B', ra='13:46:49.0432', dec='-60:24:29.355', color='w'),
            'Cen A': Source(name='Cen A (resolved)', ra='13:25:27.61507', dec='-43:01:08.8053'),
            '3C310': Source(name='3C 310', ra='15:04:57.108', dec='+26:00:58.28'),
            '3C409': Source(name='3C 409', ra='20:14:27.74', dec='+23:34:58.4', color='w'),
            '3C433': Source(name='3C 433', ra='21:23:44.582', dec='+25:04:27.23', color='w'),
            'SgrA': Source(name='Sgr A*', ra='17:45:40.0409', dec='-29:00:28.118', color='w'),
            'HerA': Source(name='Her A', ra='16:51:08.147', dec='+04:59:33.32', align='r'),
            '3C353': Source(name='3C 353', ra='17:20:28.147', dec='-00:58:47.12'),
            '3C327': Source(name='3C 327', ra='16:02:27.39', dec='+01:57:55.7'),
            '3C317': Source(name='3C 317', ra='15:16:44.487', dec='+07:01:18.00', align='r'),
            '3C298': Source(name='3C 298', ra='14:19:08.1788', dec='+06:28:34.757', align='r'),
            'VirA': Source(name='Vir A/M87', ra='12:30:49.42338', dec='+12:23:28.0439', color='g', align='r'),
            '3C270': Source(name='3C 270', ra='12:19:23.21621', dec='+05:49:29.6948', align='r'),
            '3C273': Source(name='3C 273', ra='12:29:06.69512', dec='+02:03:08.6628', align='r'),
            'PKS2356': Source(name='PKS 2356-61', ra='23:59:04.37', dec='-60:54:59.4'),
            'M1': Source(name='M1/Crab', ra='05:34:31.93830', dec='+22:00:52.1758', color='g')
        }


def calc_delays(az=0.0, el=0.0):
    """
       Function calc_delays - used to calculate the dipole delays
       for future observations, since they can't be obtained from the
       observation info record.

       This function takes in an azimuth and zenith angle as
       inputs and creates and returns a 16-element byte array for
       delayswitches which have values corresponding to each
       dipole in the tile having a maximal coherent amplitude in the
       desired direction.

       This will return null if the inputs are
       out of physical range (if za is bigger than 90) or
       if the calculated switches for the dipoles are out
       of range of the delaylines in the beamformer.

       azimuth of 0 is north and increases clockwise
       zenith angle is the angle down from zenith
       These angles should be given in degrees

      Layout of the dipoles on the tile:

                 N

           0   1   2   3

           4   5   6   7
      W                    E
           8   9   10  11

           12  13  14  15

                 S

       :param az: Azimuth in degrees
       :param el: Elevation in degrees
       :return: 16-element integer array of delays
    """
    dip_sep = 1.10  # dipole separations in meters
    delaystep = 435.0  # Delay line increment in picoseconds
    maxdelay = 31  # Maximum number of deltastep delays
    c = 0.000299798  # C in meters/picosecond

    # Convert to radians, and zenith angle instead of elevation angle
    azr = az * math.pi / 180.0
    zar = (90 - el) * math.pi / 180.0

    # Define arrays to hold the positional offsets of the dipoles

    # Check input sanity
    if (abs(zar) > math.pi / 2.0):
        return None

    # Offsets of the dipoles are calculated relative to the
    # center of the tile, with positive values being in the north
    # and east directions

    xoffsets = [-1.5 * dip_sep, -0.5 * dip_sep, 0.5 * dip_sep, 1.5 * dip_sep] * 4
    yoffsets = [1.5 * dip_sep] * 4 + [0.5 * dip_sep] * 4 + [-0.5 * dip_sep] * 4 + [-1.5 * dip_sep] * 4

    # First, figure out the theoretical delays to the dipoles
    # relative to the center of the tile

    # calculate exact delays in picoseconds from geometry...
    delays_picoseconds = [(xoffsets[i] * math.sin(azr) + yoffsets[i] * math.cos(azr)) * math.sin(zar) / c for i in range(16)]
    mindelay_picoseconds = min(delays_picoseconds)
    # Subtract minimum delay so that all delays are positive
    delays_picoseconds = [x - mindelay_picoseconds for x in delays_picoseconds]

    # Now minimize the sum of the deviations^2 from optimal
    # due to errors introduced when rounding the delays.
    # This is done by stepping through a series of offsets to
    # see how the sum of square deviations changes
    # and then selecting the delays corresponding to the min sq dev.

    # Go through once to get baseline values to compare
    bestoffset = -0.45 * delaystep
    minsqdev = 0

    for i in range(16):
        delay_off = delays_picoseconds[i] + bestoffset
        intdel = int(round(delay_off / delaystep))

        if (intdel > maxdelay):
            intdel = maxdelay

        minsqdev += (intdel * delaystep - delay_off) ** 2

    minsqdev = minsqdev / 16

    offset = (-0.45 * delaystep) + (delaystep / 20.0)
    while offset <= (0.45 * delaystep):
        sqdev = 0
        for i in range(16):
            delay_off = delays_picoseconds[i] + offset
            intdel = int(round(delay_off / delaystep))

            if (intdel > maxdelay):
                intdel = maxdelay
            sqdev = sqdev + ((intdel * delaystep - delay_off) ** 2)

        sqdev = sqdev / 16
        if (sqdev < minsqdev):
            minsqdev = sqdev
            bestoffset = offset

        offset += delaystep / 20.0

    rdelays = [0] * 16  # The rounded delays in units of delaystep
    for i in range(16):
        rdelays[i] = int(round((delays_picoseconds[i] + bestoffset) / delaystep))
        if (rdelays[i] > maxdelay):
            if (rdelays[i] > maxdelay + 1):
                return None  # Trying to steer out of range.
            rdelays[i] = maxdelay

    return [int(rd) for rd in rdelays]


def get_beam(Alt, Az, delays, freq_Hz, beam_type=list(BEAMS.keys())[0], logger=DEFAULTLOGGER):
    """
    Get a Stokes I value for every coordinate represented by the Alt/Az array values, for the given dipole
    delays and frequency in MHz.

    Depending on the value of the beam_type argument, this function uses one of mwa_pb.primarybeam.MWA_Tile_analytic (MWA_PB),
    mwa_hyperbeam.AnalyticBeam (HBA) or mwa_hyperbeam.FEEBeam (HBFEE).

    :param Alt: numpy array of altitude values, in degrees
    :param Az: numpy array of elevation values in degrees, with the same shape as Alt
    :param delays: list of 16 dipole delay values
    :param freq_Hz: Frequency in Hz
    :param beam_type: One of 'HBA' or 'HBFEE' for Hyperbeam analytic or FEE, or 'MWA_PB' for mwa_pb.primary beam analytic
    :param logger: Optional logging.Logger object
    :return: A numpy array with the same shape as the passed Alt and Az arrays, containing the beam power at each point,
             normalised to a maximum of 1.0
    """
    za_rad = (90 - Alt) * math.pi / 180
    za_rad[za_rad >= math.pi / 2.0] = math.pi / 2.0    # Mask off below the horizon, because hyperbeam doesn't handle it
    az_rad = Az * math.pi / 180

    if beam_type == 'MWA_PB':
        # this is the response for XX and YY
        respX, respY = BEAMS[beam_type].MWA_Tile_analytic(za=za_rad,
                                                          az=az_rad,
                                                          freq=freq_Hz,
                                                          delays=numpy.array(delays),
                                                          zenithnorm=True)
    elif beam_type == 'HBA':
        # Flatten the coordinate arrays, since hyperbeam methods expect 1D arrays, btu save the original shape
        oshape = za_rad.shape
        za_rad_flat = za_rad.ravel()
        az_rad_flat = az_rad.ravel()

        # There's a bug in the hyperbeam analytic code, so we have to do az_rad = (math.pi - az_rad) before passing it in
        jones = BEAMS[beam_type].calc_jones_array(az_rad=(math.pi - az_rad_flat),   # az_rad
                                                  za_rad=za_rad_flat,   # za_rad
                                                  freq_hz=freq_Hz,      # freq_hz
                                                  delays=delays,        # delays
                                                  amps=[1.0]*16,      # amps
                                                  latitude_rad=MWAPOS.lat.radian,    # latitude_rad
                                                  norm_to_zenith=True)  # norm_to_zenith

        # Restore the output to the original coordinate array shape
        respX, respY = jones[:, 0].reshape(oshape), jones[:, 1].reshape(oshape)
    elif beam_type == 'HBFEE':
        # Flatten the coordinate arrays, since hyperbeam methods expect 1D arrays, btu save the original shape
        oshape = za_rad.shape
        za_rad_flat = za_rad.ravel()
        az_rad_flat = az_rad.ravel()

        jones = BEAMS[beam_type].calc_jones_array(az_rad=az_rad_flat,   # az_rad
                                                  za_rad=za_rad_flat,   # za_rad
                                                  freq_hz=freq_Hz,      # freq_hz
                                                  delays=delays,        # delays
                                                  amps=[1.0] * 16,      # amps
                                                  norm_to_zenith=True,  # norm_to_zenith
                                                  latitude_rad=None,    # latitude_rad
                                                  iau_order=False)      # iau_order

        # Restore the output to the original coordinate array shape
        respX, respY = jones[:, 0].reshape(oshape), jones[:, 1].reshape(oshape)
    else:
        logger.error('Unknown beam type: %s' % beam_type)
        return None

    # make a pseudo-I beam, and normalise
    r = abs(respX) ** 2 + abs(respY) ** 2
    return r / numpy.nanmax(r)


def plot_MWA_skymap(delays=None,
                    channels=None,
                    viewgps=None,
                    gleamsources=False,
                    plot_text=None,
                    ra_pc=None,
                    dec_pc=None,
                    voltage_beams=None,
                    inverse=False,
                    background='black',
                    hidenulls=False,
                    beam_type='HBA',
                    plotsize=DEFAULT_PLOTSIZE,
                    img_format='png',
                    logger=DEFAULTLOGGER):
    """
    This function plots a skymap showing a single MWA observation at a single point in time. The background image and
    overlaid sources (GLEAM sources, planets, etc) are taken from the SKYDATA global, containing an instance of
    skymap.SkyData() that can be loaded and customised. If not pre-loaded into the global, it will be loaded and
    cached the first time this function is called.

    The MWA primary beam is plotted in contours over this background.

    If the 'delays' and 'channels' parameters are given, for one or more primary MWA beams, the
    get_beam() function is used to generate one or more numpy array/s using the 'delays' and 'channels'
    parameters. The 'delays' parameter can either be a list of 16 dipole delay values, or a list of lists of 16
    dipole delay values, one per primary beam. The 'channels' parameter can either be a single MWA channel number,
    or a list of MWA channel numbers, one per primary beam.

    If the 'primary_beams' parameter is given, it must be a numpy array as returned by get_beam(), or a list of those
    numpy arrays. The 'delays' and 'channels' parameters are then only used to populate the text in the sky map plot.

    If the voltage_beams parameter is provided, it should be a dictionary of voltage beams with key=beam_number,
    each a dict with 'ra', 'dec', and 'name' items, containing beam RA and Dec in degrees, and name as a string. These
    will be plotted as white circles, and labeled with the beam number and name.

    If ra_pc and dec_pc are provided, they should be the RA and Dec of the phase center in degrees. This will be
    plotted as a cyan circle.

    The map PNG image is returned as a byte array.

    The rest of the parameters are optional, and described below.

    :param delays: Either a list of 16 dipole delays (for a single primary beam), or a list of lists of dipole delays, one per primary beam
    :param channels: Either a single MWA channel number, or a list of MWA channel numbers, one per primary beam
    :param viewgps:  The GPS time in seconds for which the plot should be generated
    :param gleamsources:  If True, show the GLEAM source list as blue dots
    :param ra_pc:  Optional Right ascension of the phase center in degrees
    :param dec_pc:  Optional Declination of the phase center in degrees
    :param voltage_beams:  Optional dictionary of voltage beams with key=number, each a dict with 'ra', 'dec', and 'name' items
    :param plot_text:  If provided, include this text at the top left of the plot
    :param inverse:  If True, plot the Haslam radio image in white-on-black
    :param background:  One of 'black', 'white', or 'transparent' - default is 'black'
    :param hidenulls:  If True, don't show the null (power = 0.001) contour lines in black
    :param beam_type:  One of 'HBA' or 'HBFEE' for Hyperbeam analytic or FEE, or 'MWA_PB' for mwa_pb.primary beam analytic
    :param plotsize:  Output plot width and height in pixels
    :param img_format: Output image format, default is 'png'
    :param logger: Optional logging.Logger() instance
    :return: An empty string (if outfile is specified) or a byte array (if outfile is not specified)
    """
    global SKYDATA
    if SKYDATA is None:
        SKYDATA = SkyData()
    if not SKYDATA.valid:
        logger.error('Unable to load star/planet data, aborting.')
        return None

    plotscale = plotsize / 1200.0

    a_viewtime = Time(viewgps, format='gps', scale='utc')
    a_viewtime.delta_ut1_utc = 0  # We don't care about IERS tables and high precision answers
    LST_hours = a_viewtime.sidereal_time(kind='apparent', longitude=MWAPOS.lon)

    fig = plt.figure(figsize=(FIGSIZE * plotscale, FIGSIZE * plotscale), dpi=DPI)
    ax1 = fig.add_subplot(1, 1, 1)

    # Create a Basemap instance in RA and Dec, where MWA's latitude and the current viewgps LST are used as the zenith
    bmap = Basemap(projection='ortho', lat_0=MWAPOS.lat.deg, lon_0=LST_hours.hour * 15 - 360, ax=ax1)

    ax1.cla()

    # Transform ra/dec grid to the map coordinate frame.
    tform_skymap = bmap.transform_scalar(SKYDATA.radio_image[0].data[0][:, ::-1],  # FITS data HDU, with last axis in reverse order
                                         SKYDATA.skymapra[::-1] * 15,  # one-dimensional RA list in reverse order
                                         SKYDATA.skymapdec,    # one dimensional DEC list
                                         len(SKYDATA.skymapra),
                                         len(SKYDATA.skymapdec),
                                         masked=True)
    if inverse:
        cmap = CMI
    else:
        cmap = CM

    # Show the radio image on the map
    bmap.imshow(numpy.ma.log10(tform_skymap[:, ::-1]), cmap=cmap, vmin=math.log10(RADIO_IMAGE_LOW), vmax=math.log10(RADIO_IMAGE_HIGH), ax=ax1)
    # This line needs to be repeated to do anything for old matplotlib versions
    bmap.imshow(numpy.ma.log10(tform_skymap[:, ::-1]), cmap=cmap, vmin=math.log10(RADIO_IMAGE_LOW), vmax=math.log10(RADIO_IMAGE_HIGH), ax=ax1)

    if delays:
        if type(delays[0] is list):   # Multiple primary beams
            if type(channels) is list:   # Multiple channel values, one per beam
                beam_list = zip(delays, channels)
            else:
                beam_list = zip(delays, [channels] * len(delays))
        else:
            beam_list = [(delays, channels)]

        SKYDATA.mapgrid.obstime = a_viewtime
        altaz = SKYDATA.mapgrid.transform_to('altaz')
        Az, Alt = altaz.az.deg, altaz.alt.deg      # Both 2D grids of Az and Alt

        bi = 0
        for beam_info in beam_list:
            this_delays, this_channel = beam_info
            base_color = PBCOLORS[bi]
            if not hidenulls:
                contours = [0.001, 0.1, 0.5, 0.90]
                beamcolor = [(0.0, 0.0, 0.0),
                             tuple((x * 0.5 for x in base_color)),
                             tuple((x * 0.75 for x in base_color)),
                             base_color]
            else:
                contours = [0.1, 0.5, 0.90]
                beamcolor = [tuple((x * 0.5 for x in base_color)),
                             tuple((x * 0.75 for x in base_color)),
                             base_color]
            # get this primary beam
            R = get_beam(Alt,
                         Az,
                         this_delays,
                         this_channel * 1.28e6,
                         beam_type=beam_type)
            # show the beam
            X, Y = bmap(SKYDATA.skymapRA, SKYDATA.skymapDec)
            CS = bmap.contour(bmap.xmax - X, Y, R, contours, linewidths=plotscale, colors=beamcolor)
            ax1.clabel(CS, inline=1, fontsize=10 * plotscale)

            bi += 1

    X0, Y0 = bmap(LST_hours.hour * 15 - 360, MWAPOS.lat.deg)

    # Plot the phase center
    if ra_pc is not None and dec_pc is not None:
        newx, newy = bmap(ra_pc, dec_pc)
        good = (newx > bmap.xmin) & (newx < bmap.xmax) & (newy > bmap.ymin) & (newy < bmap.ymax)

        if good:
            bmap.scatter(bmap.xmax - newx, newy,
                         50.0 * plotscale,
                         'turquoise',
                         edgecolor='none',
                         alpha=0.5)
            ax1.text(bmap.xmax - newx - 2e5, newy,
                     'PC',
                     horizontalalignment='right',
                     fontsize=12 * plotscale,
                     color='turquoise',
                     alpha=0.5,
                     verticalalignment='center')

    # Plot any voltage beams
    if voltage_beams:
        for number, beam in voltage_beams.items():
            ra = numpy.array([beam['ra']])
            dec = numpy.array([beam['dec']])
            name = beam['target_name']
            newx, newy = bmap(ra, dec)
            good = (newx > bmap.xmin) & (newx < bmap.xmax) & (newy > bmap.ymin) & (newy < bmap.ymax)
            if good:
                bmap.scatter(bmap.xmax - newx, newy,
                             20.0 * plotscale,
                             'white',
                             edgecolor='none',
                             alpha=0.7)
                ax1.text(bmap.xmax - newx + 2e5, newy,
                         '%d:%s' % (int(number), name),
                         horizontalalignment='left',
                         fontsize=8 * plotscale,
                         color='white',
                         verticalalignment='center')

    # plot blue dots for all the GLEAM sources
    if gleamsources:
        ra = numpy.array([x[1] for x in SKYDATA.gleamcat])
        dec = numpy.array([x[2] for x in SKYDATA.gleamcat])
        flux = numpy.array([x[3] for x in SKYDATA.gleamcat])
        newx, newy = bmap(ra, dec)
        # testx, testy = bmap(newx, newy, inverse=True)
        good = (newx > bmap.xmin) & (newx < bmap.xmax) & (newy > bmap.ymin) & (newy < bmap.ymax)
        size = flux / 1.0
        size[size <= 7] = 7
        size[size >= 60] = 60
        bmap.scatter(bmap.xmax - newx[good], newy[good],
                     size[good] * plotscale,
                     'b',
                     edgecolor='none',
                     alpha=0.7)

    # plot the Solar-System bodies
    obstime = Time(viewgps, format='gps', scale='utc')
    bodies = {'Sun':(get_body(body='Sun', time=obstime, location=MWAPOS), 120, 'yellow'),
              'Jupiter':(get_body(body='Jupiter', time=obstime, location=MWAPOS), 60, 'cyan'),
              'Moon':(get_body(body='Moon', time=obstime, location=MWAPOS), 120, 'lightgray'),
              'Mars':(get_body(body='Mars', time=obstime, location=MWAPOS), 30, 'red'),
              'Venus':(get_body(body='Venus', time=obstime, location=MWAPOS), 40, 'violet'),
              'Saturn':(get_body(body='Saturn', time=obstime, location=MWAPOS), 50, 'skyblue')}
    for bname, bdata in bodies.items():
        body, size, color = bdata
        if inverse:
            if bname == 'Moon':
                color = 'darkgoldenrod'
            elif bname == 'Jupiter':
                color = 'sienna'
            elif bname == 'Saturn':
                color = 'purple'
        ra, dec = body.ra.deg, body.dec.deg
        newx, newy = bmap(ra, dec)
        testx, testy = bmap(newx, newy, inverse=True)
        if testx < 1e30 and testy < 1e30:
            bmap.scatter(2 * X0 - newx, newy, s=size * plotscale, c=color, alpha=1.0, latlon=False, edgecolor='none')
            ax1.text(bmap.xmax - newx + 2e5, newy,
                     bname,
                     horizontalalignment='left',
                     fontsize=12 * plotscale,
                     color=color,
                     verticalalignment='center')

    # Plot and label some fixed sources
    for source in SKYDATA.sources.values():
        # noinspection PyUnresolvedReferences
        r = Longitude(angle=source.ra, unit=astropy.units.hour).hour
        # noinspection PyUnresolvedReferences
        d = Latitude(angle=source.dec, unit=astropy.units.deg).deg

        color = source.color
        if inverse and color == 'w':
            color = 'black'

        xx, yy = bmap(r * 15 - 360, d)
        try:
            if xx < 1e30 and yy < 1e30:
                ax1.text(x=bmap.xmax - xx + 2e5,
                         y=yy,
                         s=source.name,
                         horizontalalignment=source.align,
                         fontsize=source.fontsize * plotscale,
                         color=color,
                         verticalalignment='center')
        except:
            pass

    if plot_text:
        if background.lower() == 'black':
            fontcolor = 'white'
        else:
            fontcolor = 'black'

        ax1.text(x=0,
                 y=bmap.ymax - 2e5,
                 s=plot_text,
                 fontsize=10 * plotscale,
                 color=fontcolor)

    ax1.text(x=bmap.xmax, y=Y0, s='W', fontsize=12 * plotscale, horizontalalignment='left', verticalalignment='center')
    ax1.text(x=bmap.xmin, y=Y0, s='E', fontsize=12 * plotscale, horizontalalignment='right', verticalalignment='center')
    ax1.text(x=X0, y=bmap.ymax, s='N', fontsize=12 * plotscale, horizontalalignment='center', verticalalignment='bottom')
    ax1.text(x=X0,y=bmap.ymin, s='S', fontsize=12 * plotscale, horizontalalignment='center', verticalalignment='top')

    try:
        buf = io.BytesIO()
        if background.lower() == 'transparent':
            fig.savefig(buf, format=img_format, transparent=True, facecolor='none', dpi=DPI)
        else:
            fig.savefig(buf, format=img_format, transparent=False, facecolor=background, dpi=DPI)
        buf.seek(0)
        return buf.read()
    finally:
        plt.close(fig)
        del ax1
        del fig


def plot_MWA_obs_frame(obsinfo=None,
                       viewgps=None,
                       gleamsources=False,
                       plot_text_template=DEFAULT_PLOT_TEXT,
                       inverse=False,
                       background='black',
                       hidenulls=False,
                       beam_type='HBA',
                       img_format='png',
                       plotsize=DEFAULT_PLOTSIZE,
                       logger=DEFAULTLOGGER):
    """
    Given an MWA observation info structure, with one or more rf_streams (primary beams), generate a skymap of the
    observation at the single time specified by viewgps.

    Most of the work is done by plot_MWA_skymap() function, this code just extracts:
        - a list of dipole delays from the rf_streams, or if there's more than one primary beam, a list of lists of dipole delays
        - a dict of voltage beams, if any, in the observation
        - ra_pc and dec_pc, the RA and Dec of the phase center in degrees
    from the obsinfo structure and passes them to plot_MWA_skymap() along with the other parameters.

    plot_text_template is the text to be inserted at the top left of the plot. Use %(<field_name>)s, etc, to insert the value of a field.
    Valid field names are:
        obsid          Obs start time in GPS seconds
        viewgps_gps    View time in GPS seconds
        viewgps_utc    View date/time as a UTC string
        obsname        Observation name
        freq_mhz       Centre channel frequency in MHz
        constellation  Name of the constellation that the phase centre is located in.

    :param obsinfo:  An MWA observation info structure, eg from the /metadata/obs web service, tilestatus.getObservationInfo()
    :param viewgps:  The GPS time in seconds for which the plot should be generated (default is the observation midpoint time)
    :param gleamsources:  If True, show the GLEAM source list as blue dots
    :param plot_text_template:  Template for text string summarising the observation data - see field names described above.
    :param inverse:  If True, plot the Haslam radio image in white-on-black
    :param background:  One of 'black', 'white', or 'transparent' - default is 'black'
    :param hidenulls:  If True, don't show the null (power = 0.001) contour lines in black
    :param beam_type:  One of 'HBA' or 'HBFEE' for Hyperbeam analytic or FEE, or 'MWA_PB' for mwa_pb.primary beam analytic
    :param plotsize:  Output plot width and height in pixels
    :param img_format: Output image format, default is 'png'
    :param logger: Optional logging.Logger() instance
    :return: An empty string (if outfile is specified) or a byte array (if outfile is not specified)
    """
    if obsinfo:
        all_delays = []
        all_channels = []
        r_list = list(obsinfo['rfstreams'].keys())
        r_list.sort()
        for rfs_id in r_list:
            rfs = obsinfo['rfstreams'][rfs_id]
            channel = rfs['frequencies'][12]
            # If the observation is in the future, calculate what delays will be used, instead of using the recorded actual delays
            if not rfs['xdelays']:
                delays = calc_delays(az=rfs['azimuth'], el=rfs['elevation'])
                logger.debug("Calculated future delays: %s" % delays)
            else:
                delays = rfs['xdelays']
                logger.debug("Used actual delays: %s" % delays)

            all_delays.append(delays)
            all_channels.append(channel)

        voltage_beams = obsinfo['voltage_beams']

        # Find the constellation that the beam is in
        if obsinfo['ra_phase_center'] is not None:
            ra_pc = obsinfo['ra_phase_center']
            dec_pc = obsinfo['dec_phase_center']
        else:
            ra_pc = obsinfo['metadata']['ra_pointing']
            dec_pc = obsinfo['metadata']['dec_pointing']
        if (ra_pc is not None) and (dec_pc is not None):
            constellation = ephem.constellation((ra_pc * math.pi / 180.0, dec_pc * math.pi / 180.0))
        else:
            constellation = ["N/A", "N/A"]

        if not viewgps:   # Default to midpoint of observation
            viewgps = (obsinfo['starttime'] + obsinfo['stoptime']) / 2

        plot_text = plot_text_template % {'obsid':obsinfo['starttime'],
                                          'viewgps_gps':viewgps,
                                          'viewgps_utc':Time(viewgps, format='gps', scale='utc').datetime.strftime('%Y-%m-%d %H:%M UT'),
                                          'obsname':obsinfo['obsname'],
                                          'freq_mhz':obsinfo['rfstreams'][r_list[0]]['frequencies'][12] * 1.28,
                                          'constellation':constellation[1]}
    else:
        plot_text = None
        voltage_beams = None

    return plot_MWA_skymap(delays=None,
                           channels=None,
                           viewgps=viewgps,
                           gleamsources=gleamsources,
                           plot_text=plot_text,
                           ra_pc=None,
                           dec_pc=None,
                           voltage_beams=voltage_beams,
                           inverse=inverse,
                           background=background,
                           hidenulls=hidenulls,
                           beam_type=beam_type,
                           plotsize=plotsize,
                           img_format=img_format,
                           logger=logger)


def mwa_apng_adaptive(outfile=None,
                      startgps=None,
                      stopgps=None,
                      obsinfo_list=None,
                      max_frame_duration=500,  # Maximum number of milliseconds per frame
                      frame_speed=100,  # Number of milliseconds of video per minute of actual observing time
                      gleamsources=False,
                      plot_text_template=DEFAULT_PLOT_TEXT,
                      inverse=False,
                      background='black',
                      hidenulls=False,
                      beam_type='HBA',
                      plotsize=DEFAULT_PLOTSIZE,
                      logger=DEFAULTLOGGER):
    """
    Given a list of observation info structures (assumed to be consecutive observations), generate an animated PNG of the skymaps
    for each observation, showing how the sky moves and the primary beam is repointed over the course of these observations.

    Each frame in an animated PNG can have its own individual duration, so new frames are generated on every new obsid, and
    also after every five minutes of actual observing time if the observation/s are longer than five minutes.

    :param outfile: Output filename, or if not supplied, the image is returned as a byte array
    :param startgps: Start time of the video in GPS seconds (defaults to the start time of the first observation in obsinfo_list)
    :param stopgps: End time of the video in GPS seconds (defaults to the end time of the last observation in obsinfo_list)
    :param obsinfo_list:  A list of observation information structures, each describing a single observation
    :param max_frame_duration: Longest time (in milliseconds) to display a single frame - will be shorter if a repointing happens
    :param frame_speed: Number of milliseconds of video per minute of actual observing time
    :param gleamsources:  If True, show the GLEAM source list as blue dots
    :param plot_text_template:  Template for text string summarising the observation data - see field names described above.
    :param inverse:  If True, plot the Haslam radio image in white-on-black
    :param background:  Onbe of 'black', 'white', or 'transparent' - default is 'black'
    :param hidenulls:  If True, don't show the null (power = 0.001) contour lines in black
    :param beam_type:  One of 'HBA' or 'HBFEE' for Hyperbeam analytic or FEE, or 'MWA_PB' for mwa_pb.primary beam analytic
    :param plotsize:  Output plot width and height in pixels
    :param logger: Optional logging.Logger() instance
    :return: An empty string (if outfile is specified) or a byte array (if outfile is not specified)
    """
    observations = {None:{}}
    obsid_list = []
    for obs in obsinfo_list:
        observations[obs['starttime']] = obs
        obsid_list.append(obs['starttime'])

    if not startgps:   # Default to start of the first observation
        startgps = obsid_list[0]

    if not stopgps:    # Default to end of the last observation, plus 8 seconds
        if not observations[obsid_list[-1]]:
            print('No observations passed, no stopgps')
            return None
        stopgps = observations[obsid_list[-1]]['stoptime'] + 8

    obsid_list.sort()
    sec_per_frame = int(60 * max_frame_duration / frame_speed)

    frame_list = []   # Will contain tuples of (obsid, viewgps, duration_seconds)
    if startgps < obsid_list[0]:
        for tgps in range(startgps, obsid_list[0], sec_per_frame):
            frame_list.append((None, tgps, sec_per_frame))
        if frame_list[-1][1] + sec_per_frame < obsid_list[0]:
            frame_list.append((None, (frame_list[-1][1] + sec_per_frame), obsid_list[0] - (frame_list[-1][1] + sec_per_frame)))

    for i in range(len(obsid_list)):
        obsid = obsid_list[i]
        obs = observations[obsid]
        if i == len(obsid_list) - 1:  # no next observation, so go a little after the stoptime of this observation
            exptime = stopgps - obs['starttime']
        else:  # Consider this pointing to last until the start of the next pointing, in case there's a gap
            exptime = observations[obsid_list[i + 1]]['starttime'] - obs['starttime']

        if exptime < sec_per_frame:   # Add a single frame for this observation, with the appropriate duration
            frame_list.append((obsid, obs['viewgps'], exptime))
        else:   # Add multiple frames for this observation, with the appropriate duration
            t = obsid
            while t < obs['starttime'] + exptime - sec_per_frame:
                frame_list.append((obsid, t, sec_per_frame))   # Most frames have a fixed duration
                t += sec_per_frame
            frame_list.append((obsid, t, obsid + exptime - t))  # Final frame has a duration until the next pointing

    print('Framelist:\n    ' + '\n    '.join(str(f) for f in frame_list))

    o_im = apng.APNG()
    for f in frame_list:
        oid, vgps, dur_sec = f
        im_frame = plot_MWA_obs_frame(obsinfo=observations[oid],
                                      viewgps=vgps,
                                      gleamsources=gleamsources,
                                      plot_text_template=plot_text_template,
                                      inverse=inverse,
                                      background=background,
                                      hidenulls=hidenulls,
                                      beam_type=beam_type,
                                      plotsize=plotsize,
                                      img_format='png',
                                      logger=logger)
        o_im.append(apng.PNG.from_bytes(im_frame), delay=int(dur_sec * frame_speed / 60.0))

    if not outfile:
        return o_im.to_bytes()
    else:
        o_im.save(outfile)
        return ''


def mwa_mpeg(outfile=None,
             startgps=None,
             stopgps=None,
             obsinfo_list=None,
             frame_rate=2,   # frames per second
             frame_speed=100,  # Number of milliseconds of video per minute of actual observing time (default=100, or 10 minutes per second)
             gleamsources=False,
             plot_text_template=DEFAULT_PLOT_TEXT,
             inverse=False,
             background='black',
             hidenulls=False,
             beam_type='HBA',
             plotsize=DEFAULT_PLOTSIZE,
             logger=DEFAULTLOGGER):
    """
    Given a list of observation info structures (assumed to be consecutive observations), generate an animated PNG of the skymaps
    for each observation, showing how the sky moves and the primary beam is repointed over the course of these observations.

    Each frame in an animated PNG can have its own individual duration, so new frames are generated on every new obsid, and
    also after every five minutes of actual observing time if the observation/s are longer than five minutes.

    :param outfile: Output filename, or if not supplied, the image is returned as a byte array
    :param startgps: Start time of the video in GPS seconds (defaults to the start time of the first observation in obsinfo_list)
    :param stopgps: End time of the video in GPS seconds (defaults to the end time of the last observation in obsinfo_list)
    :param obsinfo_list:  A list of observation information structures, each describing a single observation
    :param frame_rate: Number of frames per second of playing time in the output video
    :param frame_speed: Number of milliseconds of video per minute of actual observing time
    :param gleamsources:  If True, show the GLEAM source list as blue dots
    :param plot_text_template:  Template for text string summarising the observation data - see field names described above.
    :param inverse:  If True, plot the Haslam radio image in white-on-black
    :param background:  Onbe of 'black', 'white', or 'transparent' - default is 'black'
    :param hidenulls:  If True, don't show the null (power = 0.001) contour lines in black
    :param beam_type:  One of 'HBA' or 'HBFEE' for Hyperbeam analytic or FEE, or 'MWA_PB' for mwa_pb.primary beam analytic
    :param plotsize:  Output plot width and height in pixels
    :param logger: Optional logging.Logger() instance
    :return: An empty string (if outfile is specified) or a byte array (if outfile is not specified)
    """
    observations = {}
    obsid_list = []
    for obs in obsinfo_list:
        observations[obs['starttime']] = obs
        obsid_list.append(obs['starttime'])

    obsid_list.sort()
    sec_per_frame = int((60000 / frame_rate) / frame_speed)   # Number of seconds of observing time per frame

    if not startgps:   # Default to start of the first observation
        startgps = obsid_list[0]

    if not stopgps:    # Default to end of the last observation, plus 8 seconds
        stopgps = observations[obsid_list[-1]]['stoptime']

    frame_list = []   # Will contain tuples of (obsid, viewgps)
    for tgps in range(startgps, stopgps, sec_per_frame):
        if tgps >= obsid_list[0]:
            oid = max([x for x in obsid_list if x <= tgps])
        else:
            oid = None
        frame_list.append((oid, tgps))

    got_vaapi = os.path.exists('/dev/dri/renderD128')
    buf = None
    output = None
    if outfile:
        if got_vaapi:
            iio.get_writer(outfile, format='FFMPEG', mode='I', fps=1,
                           codec='h264_vaapi',
                           output_params=['-vaapi_device',
                                          '/dev/dri/renderD128',
                                          '-vf',
                                          'format=gray|nv12,hwupload'],
                           pixelformat='vaapi_vld')
        else:
            output = iio.get_writer(outfile, mode='I', format='FFMPEG', fps=frame_rate)
    else:
        buf = io.BytesIO()
        if got_vaapi:
            iio.get_writer(buf, format='FFMPEG', mode='I', fps=1,
                           codec='h264_vaapi',
                           output_params=['-vaapi_device',
                                          '/dev/dri/renderD128',
                                          '-vf',
                                          'format=gray|nv12,hwupload'],
                           pixelformat='vaapi_vld')
        else:
            output = iio.get_writer(buf, mode='I', format='FFMPEG', fps=frame_rate)

    for f in frame_list:
        oid, vgps = f
        im_frame = plot_MWA_obs_frame(obsinfo=observations[oid],
                                      viewgps=vgps,
                                      gleamsources=gleamsources,
                                      plot_text_template=plot_text_template,
                                      inverse=inverse,
                                      background=background,
                                      hidenulls=hidenulls,
                                      beam_type=beam_type,
                                      plotsize=plotsize,
                                      img_format='png',
                                      logger=logger)
        output.append_data(iio.imread(im_frame))

    if outfile:
        output.close()
        return ''
    else:
        buf.seek(0)
        return buf.read()


"""
import requests
import json
import time
from mwa_skymap import skymap
# obsid = 1456747216
# obsid = 1454999592
obsid = 1456774216
obs = obs = json.loads(requests.get('https://ws.mwatelescope.org/metadata/obs?obs_id=%d' % obsid).text)
for beam_type in ['MWA_PB', 'HBA', 'HBFEE']:
    st = time.time()
    im = skymap.plot_MWA_obs_frame(beam_type=beam_type, obsinfo=obs)
    et = time.time()
    print('Beam type %s: %0.3f seconds' % (beam_type, et - st))
    f = open('/tmp/%d_%s.png' % (obsid, beam_type), 'wb')
    _ = f.write(im)
    f.close()
    
import json
import requests
from mwa_skymap import skymap
# obsid_list = [1456774216]
obsid_list = [1458492896, 1458493496, 1458494096, 1458494696, 1458495296, 1458495896, 1458496496, 1458497096, 1458497696, 
              1458498296, 1458498896, 1458499496, 1458500096, 1458500696, 1458501296, 1458501896, 1458502496, 1458503096, 
              1458503696, 1458504296, 1458504896, 1458505496, 1458506096, 1458506696]
obsinfo_list = []
for obsid in obsid_list:
    obs = json.loads(requests.get('https://ws.mwatelescope.org/metadata/obs?obs_id=%d' % obsid).text)
    obsinfo_list.append(obs)
    
skymap.mwa_apng_adaptive(obsinfo_list=obsinfo_list, outfile='/tmp/%d-skymap-animated.png' % obsid_list[0])




"""