#!/usr/bin/env python

"""
Plots MWA observations with primary beam contours on a radio-map sky
"""
from importlib.resources import files
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
    data_dir = str(files('mwa_skymap.data'))
    beam_file = os.path.join(data_dir, 'mwa_full_embedded_element_pattern.h5')
    print('Downloading MWA beam file...')
    f = requests.get('http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5')
    with open(beam_file, 'wb') as beamf:
        beamf.write(f.content)
    print('Saved MWA beam file to: %s' % beam_file)


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
                                   plot_text=text,
                                   inverse=inverse,
                                   background=background,
                                   hidenulls=hidenulls,
                                   beam_type=beam_type,
                                   img_format=img_format,
                                   plotsize=plotsize)

    with open(outfile, 'wb') as f:
        f.write(im)


if __name__ == '__main__':
    cli()
