import warnings
try:
    from osgeo import gdal
except: 
    warnings.warn("Please install GDAL (see install_gdal.py). The subannagnps module will not work without it")

import pandas as pd
import numpy as np

import networkx as nx
import os


def read_annagnps_section(filename):
    df = pd.read_csv(filename, header="infer", index_col=False)
    df.reset_index(drop=True, inplace=True)
    return df


def read_esri_asc_file(filename):
    with open(filename) as file:
        dataset = gdal.Open(file)

    ncols = dataset.RasterXSize  # int
    nrows = dataset.RasterYSize  # int
    nodataval = dataset.GetRasterBand(1).GetNoDataValue()
    geoMatrix = dataset.GetGeoTransform()
    img = dataset.ReadAsArray()

    return img, geoMatrix, ncols, nrows, nodataval, dataset


def build_network(df):
    df.loc[df["Receiving_Reach"] == "OUTLET", "Receiving_Reach"] = 0
    df["Reach_ID"] = df["Reach_ID"].astype(int)
    df["Receiving_Reach"] = df["Receiving_Reach"].astype(int)

    G = nx.from_pandas_edgelist(
        df,
        source="Reach_ID",
        target="Receiving_Reach",
        edge_attr="Length",
        create_using=nx.MultiDiGraph,
    )

    return G


def get_topagnps_subgraph(G, outlet):
    # The outlet is the last "receiving reach" to consider,
    # the downstream end of it is the outlet

    subG_preorder_nodes = nx.dfs_preorder_nodes(G.reverse(), source=outlet)
    subG = nx.subgraph(G, subG_preorder_nodes)

    return subG


def relabel_topagnps_subgraph(G, outlet):
    # G must be reversed if the node direction goes from the upstream to the
    # downstream receiving reach.
    G = G.reverse()

    preorder = list(nx.dfs_preorder_nodes(G, source=outlet))

    # relabel_dict[old_key] = new_key
    relabel_dict = {}
    for old_key in preorder:
        relabel_dict[old_key] = (
            old_key - outlet + 2
        )  # +2 because the last "actual reach" with non-zero length has label 2

    # relabel_dict[1] = 'OUTLET' # This is always true and reach 1 has 0 length

    G = G.reverse()

    G = nx.relabel_nodes(G, relabel_dict, copy=True)

    return G, relabel_dict


def subset_and_rename_reach_df(df, relabel_dict):
    old_reaches = list(relabel_dict.keys())

    def rename_reach(relabel_dict, x):
        if x in relabel_dict.keys():
            return relabel_dict[x]
        else:
            return 1

    df_sub = (
        df.loc[df["Reach_ID"].isin(old_reaches), :]
        .reset_index(drop=True)
        .copy(deep=True)
    )

    df_sub["Reach_ID"] = df_sub["Reach_ID"].apply(lambda x: relabel_dict[x])
    df_sub["Receiving_Reach"] = df_sub["Receiving_Reach"].apply(
        lambda x: rename_reach(relabel_dict, x)
    )

    # df_sub.iloc[0,2:] = df_sub.iloc[0,2:] # Copy copy reach 2 (last actual reach) and insert as first line
    df_sub = pd.concat([df_sub.iloc[[0], :], df_sub], axis=0).reset_index(drop=True)

    df_sub.loc[0, "Reach_ID"] = 1
    df_sub.loc[0, "Receiving_Reach"] = "OUTLET"
    df_sub.loc[0, "Length"] = 0  # Set last reach as 0 length

    return df_sub


def subset_and_rename_cell_df(df, relabel_dict):
    old_reaches = list(relabel_dict.keys())

    def rename_reach(relabel_dict, x):
        if x in relabel_dict.keys():
            return relabel_dict[x]
        else:
            return "OUTLET"

    df_sub = (
        df.loc[df["Reach_ID"].isin(old_reaches), :]
        .reset_index(drop=True)
        .copy(deep=True)
    )

    df_sub["Reach_ID"] = df_sub["Reach_ID"].apply(
        lambda x: rename_reach(relabel_dict, x)
    )

    return df_sub


def subset_and_rename_asc_file(
    labels, path_to_in_asc, path_to_out_asc, keep_tmp_gtiff=False
):
    in_ds = gdal.Open(path_to_in_asc)  # Template
    in_band = in_ds.GetRasterBand(1)
    in_data = in_band.ReadAsArray()
    [rows, cols] = in_data.shape

    nodata = in_band.GetNoDataValue()

    if type(labels) is list:
        out_data = np.where(
            np.isin(in_data, labels), in_data, nodata
        )  # Only keep valid rasters
    elif type(labels) is dict:
        original_values = list(labels.keys())
        out_data = np.where(
            np.isin(in_data, original_values), in_data, nodata
        )  # Only keep valid rasters and relabel them
        unique_reaches, inv = np.unique(out_data, return_inverse=True)
        out_data = np.array(
            [labels.get(old_reach, nodata) for old_reach in unique_reaches]
        )[inv].reshape(out_data.shape)

    geo_transform = in_ds.GetGeoTransform()

    x_res = abs(geo_transform[1])
    y_res = abs(geo_transform[5])

    # Because GDAL doesn't allow the creation of AAIGrid ASC file, then
    # we have to write a temporary GTiff that we can then convert to .asc
    tmp_gtiff = path_to_out_asc.replace(".asc", ".gtiff")

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(tmp_gtiff, cols, rows, 1, gdal.GDT_Int32)  # 1 = one band
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(out_data)
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    out_ds.SetProjection(in_ds.GetProjection())
    out_band.FlushCache()  # This actually writes to disk

    del in_data
    del out_band
    out_ds = None
    del in_ds

    # Convert to ASC
    opts = gdal.WarpOptions(
        format="AAIGrid",
        srcNodata=nodata,
        dstNodata=nodata,
        xRes=x_res,
        yRes=y_res,
        resampleAlg=gdal.GRA_NearestNeighbour,
    )

    with open(path_to_out_asc, "w"):
        print("Converting GTiff to ASC")
        gdal.Warp(path_to_out_asc, tmp_gtiff, options=opts)
        os.remove(path_to_out_asc + ".aux.xml")

    if keep_tmp_gtiff:
        return
    else:
        print("Removing temporary GTiff file: ", tmp_gtiff)
        os.remove(tmp_gtiff)
        return


def subset_and_rename_cell_IDs_asc(
    df_cells_sub, path_to_cell_ids_asc, outpath_to_cell_ids_asc, keep_tmp_gtiff=False
):
    valid_cells = df_cells_sub["Cell_ID"].to_list()

    subset_and_rename_asc_file(
        valid_cells,
        path_to_cell_ids_asc,
        outpath_to_cell_ids_asc,
        keep_tmp_gtiff=keep_tmp_gtiff,
    )

    return


def subset_and_rename_reach_IDs_asc(
    relabel_dict, path_to_reach_ids_asc, outpath_to_reach_ids_asc, keep_tmp_gtiff=False
):
    # This function takes relabel_dict as data input because the reaches are relabeled
    # but not the cells

    # The keys of relabel_dict are the old reach IDs and the values are the new ones

    subset_and_rename_asc_file(
        relabel_dict,
        path_to_reach_ids_asc,
        outpath_to_reach_ids_asc,
        keep_tmp_gtiff=keep_tmp_gtiff,
    )

    return


def write_csv_data_section(df, out_path):
    df.to_csv(out_path, sep=",", float_format="%.5f", index=False)

    return


def draw_network(G, with_labels=True):
    df = pd.DataFrame(index=G.nodes(), columns=G.nodes())
    for row, data in nx.shortest_path_length(G):
        for col, dist in data.items():
            df.loc[row, col] = dist

    df = df.fillna(df.max().max())

    pos = nx.kamada_kawai_layout(G, dist=df.to_dict())

    nx.draw_networkx(G, pos, with_labels=with_labels)
