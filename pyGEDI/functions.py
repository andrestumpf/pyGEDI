#!/usr/bin/env python

import os, sys, ast, re

try:
    import alive_progress
except ImportError:
    sys.exit("""You need following module: alive_progress """)

try:
    import pandas
except ImportError:
    sys.exit("""You need following module: pandas """)

try:
    import geopandas
except ImportError:
    sys.exit("""You need following module: geopandas """)

try:
    import shapely
except ImportError:
    sys.exit("""You need following module: shapely """)

try:
    import requests
except ImportError:
    sys.exit("""You need following module: requests """)

try:
    import h5py
except ImportError:
    sys.exit("""You need following module: h5py """)

try:
    import pandas as pd
except ImportError:
    sys.exit("""You need following module: pandas """)

try:
    import numpy as np
except ImportError:
    sys.exit("""You need following module: numpy """)

try:
    import subprocess
except ImportError:
    sys.exit("""You need following module: subprocess """)

from pyGEDI.get import *


def gediDownload(outdir, product, version, bbox, session):
    try:
        os.makedirs(outdir)
    except OSError:
        print("Creation of the subdirectory %s failed" % outdir)
    else:
        print("Created the subdirectory %s" % outdir)

    url = 'https://lpdaacsvc.cr.usgs.gov/services/gedifinder?product=' + product + '&version=' + str(
        version) + '&bbox=' + str(bbox) + '&output=json'
    content = requests.get(url)
    listh5 = content.json().get('data')

    for url in listh5:
        url_response(outdir, url, session)


def bbox2polygon(bbox):
    [ul_lat, ul_lon, lr_lat, lr_lon] = bbox
    return shapely.geometry.Polygon([[ul_lat, ul_lon],[ul_lat, lr_lon], [lr_lat, lr_lon], [lr_lat, ul_lon]])


def extract_bbox(fileh5, bbox, latlayer='geolocation/lat_lowestmode', lonlayer='geolocation/lon_lowestmode',
                 layers=None):

    """
    :param fileh5: object of class h5py._hl.files.File
    :param bbox: list defining bounding box with the format [ul_lat, ul_lon, lr_lat, lr_lon]
    :param latlayer: string specifying the latitude layer
    :param lonlayer: string specifying the latitude layer
    :param layers: list with strings matching available input layers
    :return: Pandas data frame including filename; beam, shot number and the requested layers within the bbox
    """

    bbox_polygon = bbox2polygon(bbox)
    df = pandas.DataFrame(columns=layers)

    beams = ['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
    print('Extracting requested subset of layers for each beam...')
    with alive_progress.alive_bar(len(beams)) as bar:

        for beam in beams:
            df_beam = pandas.DataFrame()
            x = fileh5[beam][latlayer][:]
            y = fileh5[beam][lonlayer][:]
            shot_number = fileh5[beam]['shot_number'][:]

            gdf = geopandas.GeoDataFrame(geometry=geopandas.points_from_xy(x, y))
            spatial_index = gdf.sindex
            matches_index = list(spatial_index.intersection(bbox_polygon.bounds))

            # check first if specified layers are present
            for layer in layers:
                if not layer in fileh5[beam].keys():
                   raise Exception(f'Layer {layer} not found for {beam} and {fileh5}!')

            df_beam['shot_number'] = shot_number[matches_index]
            df_beam['beam'] = beam
            df_beam['filename'] = fileh5.filename

            if layers:
                for layer in layers:
                    layer_values = fileh5[beam][layer][:]
                    if 'ancillary' in layer:
                        df_beam[layer] = layer_values.tolist() * len(matches_index)
                    else:
                        df_beam[layer] = layer_values[matches_index].tolist()

            df = df.append(df_beam, ignore_index=True)
            bar()

    return df.infer_objects()


def idsBox(filesh5, latlayer, lonlayer, bbox):

    bbox_polygon = bbox2polygon(bbox)
    ids = []
    for fileh5 in filesh5:

        for beam in ['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']:
            x = fileh5[beam][latlayer][:]
            y = fileh5[beam][lonlayer][:]
            shot_number = fileh5[beam]['shot_number'][:]

            gdf = geopandas.GeoDataFrame(geometry=geopandas.points_from_xy(x, y))
            spatial_index = gdf.sindex
            matches_index = list(spatial_index.intersection(bbox_polygon.bounds))
            ids += [(beam, shot) for shot in list(shot_number[matches_index])]

    return ids


def generateBoxDataFrame(filesh5, layers, idsbox):
    df = pd.DataFrame()

    # check first if specified layers are present in all files
    for layer in layers:
        for fileh5 in filesh5:
            if not layer in fileh5['BEAM0000'].keys():
                raise Exception(f'Layer {layer} not found in {fileh5}')

    for layer in layers:
        colum = []
        for ids in idsbox:
            for fileh5 in filesh5:
                for i in np.where(fileh5[ids[0]]['shot_number'][:] == ids[1])[0]:
                    if i and (layer in fileh5[ids[0]].keys()):
                        value = [fileh5[ids[0]][layer][i]]
                        colum += value
        df[layer] = colum
    return df


def generateDataFrame(filesh5, layers):
    df = pd.DataFrame()
    for layer in layers:
        colum = []
        for beam in ['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']:
            for fileh5 in filesh5:
                if layer in fileh5[beam].keys():
                    value = fileh5[beam][layer]
                    colum += value
                if layer in ['beam', 'shot_number', 'sensitivity']:
                    break
        df[layer] = colum
    return df


def url_response(outdir, url, session):
    fileh5 = url[url.rfind('/') + 1:]
    day = url[url.rfind(':') + 41:url.rfind('/') + 1]
    path = outdir + day
    try:
        os.makedirs(path)
    except OSError:
        print("Creation of the subdirectory %s failed" % path)
    else:
        print("Created the subdirectory %s" % path)
    path5 = outdir + day + fileh5

    response = session.get(url, stream=True)
    total = response.headers.get('content-length')
    file_exists = os.path.exists(path5)
    if file_exists and total == str(os.stat(path5).st_size):
        print("File " + path5 + " is fully downloaded already on disk")
    else:
        print('Start download...')
        with open(path5, 'wb') as f:
            if total is None:
                f.write(response.content)
            else:
                downloaded = 0
                total = int(total)
                for data in response.iter_content(chunk_size=max(int(total / 1000), 1024 * 1024)):
                    downloaded += len(data)
                    f.write(data)
                    done = int(100 * downloaded / total)
                    gb = float(total / 1073741824)

                    sys.stdout.write('\r' + url[url.rfind(':') + 52:] + ' | ' + str(gb)[:5] + 'GB | ' + str(
                        100 * downloaded / total) + '% [{}{}]'.format('█' * done, '.' * (100 - done)))
                    sys.stdout.flush()
    sys.stdout.write('\n')


def waveForm(shot_number, fileh5):
    beam = getBeam(shot_number, fileh5)
    shot_number_id = list(fileh5[beam]['shot_number'][:]).index(shot_number)

    elevation_bin0 = fileh5[beam]['geolocation/elevation_bin0'][()]
    elevation_lastbin = fileh5[beam]['geolocation/elevation_lastbin'][()]
    rx_sample_count = fileh5[beam]['rx_sample_count'][()]
    rx_sample_start_index = fileh5[beam]['rx_sample_start_index'][()]

    rx_sample_start_index_n = rx_sample_start_index - min(rx_sample_start_index) + 1

    rx_sample_start = int(rx_sample_start_index_n[shot_number_id])
    rx_sample_end = int(rx_sample_start_index_n[shot_number_id] + rx_sample_count[shot_number_id] - 1)

    rxwaveform = fileh5[beam]['rxwaveform'][rx_sample_start:rx_sample_end][()]

    elevation_bin0_i = elevation_bin0[shot_number_id]
    elevation_lastbin_i = elevation_lastbin[shot_number_id]

    step = (elevation_bin0_i - elevation_lastbin_i) / rx_sample_count[shot_number_id]
    elevation = np.arange(elevation_lastbin_i, elevation_bin0_i, step)[:-1]

    return rxwaveform, elevation[::-1]
