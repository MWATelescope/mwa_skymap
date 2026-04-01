#!/usr/bin/env python

"""
Plots MWA observations with primary beam contours on a radio-map sky
"""

import json
import os

import click
import requests

from mwa_skymap import skymap


@click.group()
def cli():
    pass


def get_observation(obsid):
    obs = requests.get('http://ws.mwatelescope.org/metadata/obs?obs_id=%d' % obsid).json()
    return obs


@cli.command()
def getbeamfile():
    if os.path.exists(skymap.MWA_BEAM_FILE):
        print('MWA beam file already exists here: %s' % skymap.MWA_BEAM_FILE)
        return
    print('Downloading MWA beam file...')
    f = requests.get('http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5')
    with open(skymap.MWA_BEAM_FILE, 'wb') as beamf:
        beamf.write(f.content)
    print('Saved MWA beam file to: %s' % skymap.MWA_BEAM_FILE)


@cli.command()
def beamtypes():
    """
    Prints all the available beam type codes
    :return:
    """
    print('Available beam types: ')
    for b in skymap.BEAMS.keys():
        print('    %s' % b)
    print()


@cli.command()
@click.argument('obsid', type=int, default=0)
@click.option('--viewgps', type=int, default=None, help='Plot reference time in GPS seconds (defaults to observation midpoint)')
@click.option('--gleamsources', is_flag=True, help='Show GLEAM sources as blue dots')
@click.option('--text', type=str, default=None, help='Text to show on plot, instead of the default')
@click.option('--inverse', is_flag=True, help='Show HASLAM map as black-on-white')
@click.option('--background', type=str, default='black', help="One of 'black', 'white', or 'transparent' - default is 'black'")
@click.option('--hidenulls', is_flag=True, help="Don't show black contours for beam nulls")
@click.option('--beam_type', type=str, default='HBA', help="One of %s" % (', '.join(skymap.BEAMS.keys())))
@click.option('--plotsize', type=int, default=1200, help='Plot width and height in pixels (default 1200)')
@click.option('--outfile', type=str, default=None, help='Output filename, default is <obsid>.png - extension determines file type')
def single(obsid, viewgps, gleamsources, text, inverse, background, hidenulls, beam_type, plotsize, outfile):
    """
    Plots a single observation, as a single still frame
    """
    if not outfile:
        outfile = '%d.png' % obsid
    img_format = outfile.split('.')[-1]

    obs = get_observation(obsid)
    if not viewgps:
        viewgps = int((obs['starttime'] + obs['stoptime']) / 2)

    im = skymap.plot_MWA_obs_frame(obsinfo=obs,
                                   viewgps=viewgps,
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
@click.argument('obsids', type=int, nargs=-1)
@click.option('--startgps', type=int, default=None, help='Movie start time in GPS seconds (default to start of first obsid)')
@click.option('--stopgps', type=int, default=None, help='Movie start time in GPS seconds (default to end of last obsid)')
@click.option('--gleamsources', is_flag=True, help='Show GLEAM sources as blue dots')
@click.option('--text', type=str, default=None, help='Text to show on plot, instead of the default')
@click.option('--inverse', is_flag=True, help='Show HASLAM map as black-on-white')
@click.option('--background', type=str, default='black', help="One of 'black', 'white', or 'transparent' - default is 'black'")
@click.option('--hidenulls', is_flag=True, help="Don't show black contours for beam nulls")
@click.option('--beam_type', type=str, default='HBA', help="One of %s" % (', '.join(skymap.BEAMS.keys())))
@click.option('--plotsize', type=int, default=1200, help='Plot width and height in pixels (default 1200)')
@click.option('--outfile', type=str, default=None, help='Output filename - extension determines file type')
def movie(obsids, startgps, stopgps, gleamsources, text, inverse, background, hidenulls, beam_type, plotsize, outfile):
    """
    Plots a movie, either animated PNG or MPEG
    """
    img_format = outfile.split('.')[-1].upper()

    obsinfo_list = []
    if obsids:
        obsids = list(obsids)
        obsids.sort()
        for obsid in obsids:
            obs = json.loads(requests.get('https://ws.mwatelescope.org/metadata/obs?obs_id=%d' % obsid).text)
            obsinfo_list.append(obs)

    if not startgps:
        if obsinfo_list:
            viewgps = int(obsinfo_list[0]['starttime'])
        else:
            print('Need --startgps or at least one obsid')
            return -1

    if not stopgps:
        if obsinfo_list:
            viewgps = int(obsinfo_list[-1]['stoptime'])
        else:
            print('Need --stopgps or at least one obsid')
            return -1

    if img_format == 'PNG':
        im = skymap.mwa_apng_adaptive(outfile=outfile,
                                      obsinfo_list=obsinfo_list,
                                      startgps=startgps,
                                      stopgps=stopgps,
                                      gleamsources=gleamsources,
                                      plot_text_template=text,
                                      inverse=inverse,
                                      background=background,
                                      hidenulls=hidenulls,
                                      beam_type=beam_type,
                                      plotsize=plotsize)
    elif img_format == 'MPG':
        im = skymap.mwa_mpeg(outfile=outfile,
                             obsinfo_list=obsinfo_list,
                             startgps=startgps,
                             stopgps=stopgps,
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
