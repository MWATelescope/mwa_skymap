#!/usr/bin/env python

"""
Plots MWA observations with primary beam contours on a radio-map sky
"""

import glob
import json
import os
import sys

import click
import requests

from mwa_skymap import mwaplot


@click.group()
def cli():
    pass


def get_observation(obsid, ldir=None):
    """
    Returns a dictionary of observation metadata for the given obsid.

    If ldir is provided, look for an <obsid>.json file in the 'ldir' directory, otherwise
    download the actual MWA observation dictionary from the MWA web service
    :param obsid: Ten-digit observation ID
    :param ldir: Local directory name to look in for <obsid>.json
    :return: An observation dictionary structure
    """
    if ldir or (type(obsid) is str and obsid.endswith('.json')):
        if obsid.endswith('.json'):
            obs_file = obsid
        elif obsid.isdigit():
            obs_file = os.path.join(ldir, '%s.json' % obsid)
        else:
            print('Invalid obsid: %s' % obsid)
            return None

        if os.path.exists(obs_file):
            with open(obs_file, 'r') as f:
                obs = json.load(f)
            return obs
        else:
            print('File %s not found' % obs_file)
            return None
    else:
        result = requests.get('https://ws.mwatelescope.org/metadata/obs?obs_id=%s' % obsid)
        if result.status_code != 200:
            print('Failed to fine observation %s in the MWA schedule' % obsid)
            return None
        obs = result.json()
    return obs


@cli.command()
def getbeamfile():
    if os.path.exists(mwaplot.MWA_BEAM_FILE):
        print('MWA beam file already exists here: %s' % mwaplot.MWA_BEAM_FILE)
        return
    print('Downloading MWA beam file...')
    f = requests.get('https://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5')
    with open(mwaplot.MWA_BEAM_FILE, 'wb') as beamf:
        beamf.write(f.content)
    print('Saved MWA beam file to: %s' % mwaplot.MWA_BEAM_FILE)


@cli.command()
def beamtypes():
    """
    Prints all the available beam type codes
    :return:
    """
    print('Available beam types: ')
    for b in mwaplot.BEAMS.keys():
        print('    %s' % b)
    print()


@cli.command()
@click.argument('obsid', type=str, default=0)
@click.option('--ldir', type=str, default=None, help='Local directory to look for <obsid>.json dummy observation files. Default to use MWA web service')
@click.option('--viewgps', type=int, default=None, help='Plot reference time in GPS seconds (defaults to observation midpoint)')
@click.option('--cchan', type=int, default=None, help='Coarse channel number (0-255, defaults to 13th channel in observation)')
@click.option('--gleamsources', is_flag=True, help='Show GLEAM sources as blue dots')
@click.option('--text', type=str, default=None, help='Text to show on plot, instead of the default')
@click.option('--inverse', is_flag=True, help='Show HASLAM map as black-on-white')
@click.option('--background', type=str, default='black', help="One of 'black', 'white', or 'transparent' - default is 'black'")
@click.option('--hidenulls', is_flag=True, help="Don't show black contours for beam nulls")
@click.option('--beam_type', type=str, default='HBA', help="One of %s" % (', '.join(mwaplot.BEAMS.keys())))
@click.option('--plotsize', type=int, default=1200, help='Plot width and height in pixels (default 1200)')
@click.option('--outfile', type=str, default=None, help='Output filename, default is <obsid>.png - extension determines file type')
def single(obsid, ldir, viewgps, cchan, gleamsources, text, inverse, background, hidenulls, beam_type, plotsize, outfile):
    """
    Plots a single observation, as a single still frame.

    If --ldir is given, this command will try to load <obsid>.json in that directory.
    Otherwise, the MWA web service is used to download the observation dictionary.
    """
    if not outfile:
        if obsid.endswith('.json'):
            outfile = '%s.png' % os.path.splitext(os.path.basename(obsid))[0]
        elif obsid.isdigit():
            outfile = '%s.png' % obsid
        else:
            print('Invalid obsid: %s' % obsid)
            sys.exit(-1)

    img_format = outfile.split('.')[-1]

    obs = get_observation(obsid, ldir=ldir)
    if obs is None:
        print('Failed to get observation %s' % obsid)
        sys.exit(-1)

    if not viewgps:
        viewgps = int((obs['starttime'] + obs['stoptime']) / 2)

    im = mwaplot.plot_MWA_obs_frame(obsinfo=obs,
                                    viewgps=viewgps,
                                    cchan=cchan,
                                    gleamsources=gleamsources,
                                    plot_text_template=text,
                                    inverse=inverse,
                                    background=background,
                                    hidenulls=hidenulls,
                                    beam_type=beam_type,
                                    img_format=img_format,
                                    plotsize=plotsize)

    with open(outfile, 'wb') as f:
        f.write(im)


@cli.command()
@click.argument('obsids', type=str, nargs=-1)
@click.option('--ldir', type=str, default=None, help='Local directory to look for <obsid>.json dummy observation files. Default to use MWA web service')
@click.option('--startgps', type=int, default=None, help='Movie start time in GPS seconds (default to start of first obsid)')
@click.option('--stopgps', type=int, default=None, help='Movie start time in GPS seconds (default to end of last obsid)')
@click.option('--cchan', type=int, default=None, help='Coarse channel number (0-255, defaults to 13th channel in observation)')
@click.option('--fps', type=int, default=10, help='Movie frames per second (default 10)')
@click.option('--mps', type=int, default=10, help='Movie speed, in minutes of observing time per second of movie (default 10)')
@click.option('--gleamsources', is_flag=True, help='Show GLEAM sources as blue dots')
@click.option('--text', type=str, default=mwaplot.DEFAULT_PLOT_TEXT, help='Text to show on plot, instead of the default')
@click.option('--inverse', is_flag=True, help='Show HASLAM map as black-on-white')
@click.option('--background', type=str, default='black', help="One of 'black', 'white', or 'transparent' - default is 'black'")
@click.option('--hidenulls', is_flag=True, help="Don't show black contours for beam nulls")
@click.option('--beam_type', type=str, default='HBA', help="One of %s" % (', '.join(mwaplot.BEAMS.keys())))
@click.option('--plotsize', type=int, default=1200, help='Plot width and height in pixels (default 1200)')
@click.option('--outfile', type=str, default=None, help='Output filename - extension determines file type')
def movie(obsids, ldir, startgps, stopgps, cchan, fps, mps, gleamsources, text, inverse, background, hidenulls, beam_type, plotsize, outfile):
    """
    Plots a movie, either animated PNG or MPEG

    If --ldir is given, this command will try to load <obsid>.json files in that directory (either from the obsid list,
    or using startgps and/or stopgps).

    Otherwise, the MWA web service is used to search for and download the observation dictionaries.
    """
    obsinfo_list = []
    if obsids:
        obsids = list(obsids)
        obsids.sort()
        for obsid in obsids:
            obs = get_observation(obsid, ldir=ldir)
            if obs is None:
                print('Failed to get observation %s' % obsid)
                sys.exit(-1)
            obsinfo_list.append(obs)
    else:
        if not (startgps and stopgps):
            print('Need --startgps and --stopgps or at least one obsid')
            return -1
        if ldir:
            flist = glob.glob(os.path.join(ldir, '*.json'))
            for f in flist:
                foid = os.path.basename(f)[:10]
                if foid.isdigit() and (startgps <= int(foid) <= stopgps):
                    obs = get_observation(obsid=foid, ldir=ldir)
                    if obs is None:
                        print('Failed to get observation %s' % foid)
                        sys.exit(-1)
                    obsinfo_list.append(obs)
        else:
            result = json.loads(requests.get('https://ws.mwatelescope.org/metadata/find?mintime=%d&maxtime=%d' % (startgps, stopgps)).text)
            for block in result:
                obs = get_observation(obsid=block[0])
                if obs is None:
                    print('Failed to get observation %s' % block[0])
                    sys.exit(-1)
                obsinfo_list.append(obs)

    if not startgps:
        if obsinfo_list:
            startgps = int(obsinfo_list[0]['starttime'])
        else:
            print('Need --startgps or at least one obsid')
            return -1

    if not stopgps:
        if obsinfo_list:
            stopgps = int(obsinfo_list[-1]['stoptime'])
        else:
            print('Need --stopgps or at least one obsid')
            return -1

    if not outfile:
        outfile = '%d.png' % startgps
    img_format = outfile.split('.')[-1].upper()

    if img_format == 'PNG':
        _m = mwaplot.mwa_apng_adaptive(outfile=outfile,
                                       obsinfo_list=obsinfo_list,
                                       startgps=startgps,
                                       stopgps=stopgps,
                                       cchan=cchan,
                                       fps=fps,
                                       mps=mps,
                                       gleamsources=gleamsources,
                                       plot_text_template=text,
                                       inverse=inverse,
                                       background=background,
                                       hidenulls=hidenulls,
                                       beam_type=beam_type,
                                       plotsize=plotsize)
    elif img_format == 'MPG':
        _m = mwaplot.mwa_mpeg(outfile=outfile,
                              obsinfo_list=obsinfo_list,
                              startgps=startgps,
                              stopgps=stopgps,
                              cchan=cchan,
                              fps=fps,
                              mps=mps,
                              gleamsources=gleamsources,
                              plot_text_template=text,
                              inverse=inverse,
                              background=background,
                              hidenulls=hidenulls,
                              beam_type=beam_type,
                              plotsize=plotsize)
    else:
        print('Unknown image format: %s' % img_format)
        return -1

    return 0


if __name__ == '__main__':
    cli()
