#### --------------------------------------------------------------------------------------------------------------------------//
#### Execution: python cluster_analysis_GW.py Control_hsc.csv Control_ct_filt.csv Control_label.csv out_dir_path
#### Overview: This is a python helper script to the "hub_pipeline.sh" script that calculates a variety of properties 
# of hubs from an input list, including interaction and enhancer/promoter counts (among other attributes).
#
#### Args: .csv list of identified hubs from cluster-tree output/hierarchical spectral clustering, .csv list of input
# PP/EE/EP interactions, .csv list of annotated enhancer/promoters (all files are for a single condition), and output directory path
#
#### Dependencies: Multiple-- see libraries import statements below
#
#### Intermediate Outputs: Multiple files containing plots of the distribution of hubs on various calculated attributes (all optional)
#
#### Final Output: .csv file ("*cluster_plot_df_reconnected.csv") that contains a datatable of all hubs annotated with various connectivity, 
# interaction count, and enhancer/promoter count properties (including several that are optional)
#### -------------------------------------------------------------------------------------------------------------------------//

#### Import all dependencies
import functools as f
import os
import sys

import altair as alt
import igraph
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

import plots

mpl.use('Agg')
plt.switch_backend('agg')

# register the custom theme under a chosen name
alt.themes.register("publishTheme", plots.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")

# Function to plot distribution of connectedness throughout clusters
# Takes in:
# arg 1: clusters .csv
# arg 2: edge table .csv
# arg 3: IDs table .csv
# arg 4: name for output directory
# Information of what vertex is in what cluster in clusters.csv
# Information on what vertices are connected (a list of edges) in edges.csv
# The symbol associated with the gene, enhancer, or promoter ID is in IDs.csv
# Make sure to put a "/" after the output file you want the results located in
# This can be run on the files with binary weights (1s) and continuous weights
try:
    clusters = pd.read_csv(
        sys.argv[1:][0], header=None, index_col=None, skiprows=1)
    tmpClustersDF = pd.read_csv(sys.argv[1:][0], index_col=None)
    clusterVerticesDF = plots.getClusterVerticesDF(
        plots.getSymbols(tmpClustersDF))
    sample_name = os.path.basename(sys.argv[1:][0]).split("_clusters.csv")[0]
    edges = pd.read_csv(
        sys.argv[1:][1], names=["from", "to", "weight"], skiprows=1)
    IDs = pd.read_csv(sys.argv[1:][2], header=None)
    ID_dict = dict(list(zip(IDs[0].tolist(), IDs[1].tolist())))
    out = sys.argv[1:][3]

    # Make output directory
    if not os.path.exists(out):
        os.makedirs(out)

    pp = PdfPages(out + sample_name + '_plots.pdf')
    ppInteractive = out + sample_name + "_plots.html"
except:
    print("Input error")
    raise

# Giant graph for all clusters
g = igraph.Graph(len(clusters))
g.vs["name"] = clusters[0]
g.vs["cluster"] = clusters[1]
g.add_edges([(i[0], i[1]) for i in edges.values])
g.es["weight"] = edges["weight"]

# Analysis for each cluster
# Fill in dataframe with regular clusters
cg = clusters.groupby(1)
metrics = [
    "cluster", "number of vertices", "number of edges", "mean degree",
    "mean strength", "clique", "edge ratio", "cohesion", "max cohesion",
    "adhesion", "max adhesion", "mean path length", "element1", "element2"
]
clusterDF = pd.DataFrame(columns=metrics)
dis_comp = []
for i in list(cg.groups.keys()):
    # create subgraph for cluster
    curr_clust = g.induced_subgraph(cg.get_group(i)[0].tolist())
    # identify disconnected components
    if not curr_clust.is_connected():
        comps = curr_clust.components()
        main_sg = comps.giant().vs["name"]
        for sg in comps.subgraphs():
            curr_sg = sg.vs["name"]
            if curr_sg != main_sg:
                dis_comp = dis_comp + [curr_sg]
    # [OPTIONAL] determine if cluster contains element1 or element2 gene
    element1 = False
    element2 = False
    # if "chr11_69453856_69454119" in curr_clust.vs["name"]:
    #     element1 = True
    #     print("CCND1 is in " + i)
    #     curr_clust.vs["label"] = [
    #         ID_dict.get(k) for k in curr_clust.vs["name"]
    #     ]
    #     igraph.plot(
    #         curr_clust,
    #         out + sample_name + "_element1.pdf",
    #         layout=curr_clust.layout("circle"),
    #         margin=100)
    # if "chr17_70117008_70117525" in curr_clust.vs["name"]:
    #     element2 = True
    #     print("element2 is in " + i)
    #     curr_clust.vs["label"] = [
    #         ID_dict.get(k) for k in curr_clust.vs["name"]
    #     ]
    #     igraph.plot(
    #         curr_clust,
    #         out + sample_name + "_YWHAZ.pdf",
    #         layout=curr_clust.layout("circle"),
    #         margin=100)
    # mean degree
    meanDeg = np.mean(curr_clust.degree())
    # mean strength
    meanStrg = np.mean(curr_clust.strength(weights="weight"))
    # mean path length
    meanPathlen = np.mean(curr_clust.shortest_paths(weights=None))
    # clique
    cliq = int(curr_clust.omega())
    # edge ratio
    num_vert = curr_clust.vcount()
    if num_vert > 1:
        edgerat = float(curr_clust.ecount()) / (num_vert * (num_vert - 1) / 2)
    else:
        edgerat = 0
    # cohesion
    coh = curr_clust.cohesion()
    # max cohesion
    maxCohesion = plots.maxCohesion(curr_clust)
    # adhesion
    adhesion = curr_clust.adhesion()
    # max adhesion
    maxAdhesion = plots.maxAdhesion(curr_clust)
    # number of edges in cluster
    num_edg = curr_clust.ecount()
    # add cluster's values to dataframe
    clusterDF.loc[len(clusterDF)] = [
        i, num_vert, num_edg, meanDeg, meanStrg, cliq, edgerat, coh,
        maxCohesion, adhesion, maxAdhesion, meanPathlen, element1, element2
    ]
# add normalized cluster size
minV = float(min(clusterDF["number of vertices"]))
maxV = max(clusterDF["number of vertices"])
clusterDF["num of vertices (norm)"] = clusterDF["number of vertices"].apply(
    lambda i: (i - minV) / (maxV - minV))
# add normalized edge ratio
clusterDF["edge ratio (norm)"] = clusterDF[
    "num of vertices (norm)"] * clusterDF["edge ratio"]

## Plotting
plotting_metrics = [
    i for i in clusterDF.columns.tolist()
    if i not in ["element1", "element2", "cluster"]
]
clusterDF[plotting_metrics] = clusterDF[plotting_metrics].apply(pd.to_numeric)
clusterDF.to_csv(out + sample_name + "_cluster_df.csv", sep=",")
#plot_order = [6, 8, 0, 1, 7, 4, 9, 3, 5, 2]
plotDF = plots.assignRelevantText(
        pd.merge(
            left=clusterDF.reset_index(),
            right=clusterVerticesDF.reset_index(),
            on="cluster"))
for col in plotting_metrics:
    currDF = clusterDF.sort_values(col).reset_index(drop=True).reset_index()
    fig = currDF.plot.scatter(
        x='index',
        y=col,
        s=0.5,
        title=sample_name + ": " + col + " of clusters",
        grid=True)
    if col != "path length":
        x = currDF[currDF['element1']]
        if len(x) > 0:
            plt.text(x.index.item(), x[col].item(), "CCND1")
        x = currDF[currDF['element2']]
        if len(x) > 0:
            plt.text(x.index.item(), x[col].item(), "element2")
    plt.xlabel("rank")
    plt.ylabel(col)
    plt.savefig(pp, format="pdf")

    plotDF[col + " index"] = clusterDF[col].argsort().argsort()

selection = alt.selection_interval()
base = plots.plotBase(plotDF, selection)
interactivePlots = [
    plots.addTextEncodingRank(plotDF, x, plots.addSortedEncoding(base, x))
    for x in plotting_metrics
]

# Join and save interactive plots
joinedPlots = plots.plotSquare(interactivePlots)
joinedPlots.save(ppInteractive)

plotDF.to_csv(out + sample_name + "_cluster_plot_df.csv", index=False)

# Second part
# Fill in dataframe with disconnected components reintegrated into appropriate clusters
for curr_dis_comp in dis_comp:
    total_neigh = []
    for curr_vert in curr_dis_comp:
        total_neigh = total_neigh + g.neighbors(curr_vert)
    if g.vs(total_neigh)["cluster"]:
        mode_clust = stats.mode(g.vs(total_neigh)["cluster"])[0][0]
        for curr_vert in curr_dis_comp:
            rep = pd.Series([mode_clust],
                            name=1,
                            index=[g.vs.find(curr_vert).index])
            clusters.update(rep)

# Perform same analysis as above on new df
# Analysis for each cluster
# Fill in dataframe with regular clusters
cg = clusters.groupby(1)
clusterDFnew = pd.DataFrame(columns=metrics)
for i in list(cg.groups.keys()):
    curr_clust = g.induced_subgraph(cg.get_group(i)[0].tolist())
    
    # Optional variables used to search which hub a given enh/prom is part of
    element1 = False
    element2 = False
    # if "chr11_69453856_69454119" in curr_clust.vs["name"]:
    #     element1 = True
    #     print("CCND1 is in " + i)
    #     curr_clust.vs["label"] = [
    #         ID_dict.get(k) for k in curr_clust.vs["name"]
    #     ]
    #     igraph.plot(
    #         curr_clust,
    #         out + sample_name + "_element1_reconnected.pdf",
    #         layout=curr_clust.layout("circle"),
    #         margin=100)
    # if "chr17_70117008_70117525" in curr_clust.vs["name"]:
    #     element2 = True
    #     print("element2 is in " + i)
    #     curr_clust.vs["label"] = [
    #         ID_dict.get(k) for k in curr_clust.vs["name"]
    #     ]
    #     igraph.plot(
    #         curr_clust,
    #         out + sample_name + "_YWHAZ_reconnected.pdf",
    #         layout=curr_clust.layout("circle"),
    #         margin=100)
    # mean degree
    meanDeg = np.mean(curr_clust.degree())
    # mean strength
    meanStrg = np.mean(curr_clust.strength(weights="weight"))
    # mean path length
    meanPathlen = np.mean(curr_clust.shortest_paths(weights=None))
    # clique
    cliq = int(curr_clust.omega())
    # edge ratio
    num_vert = curr_clust.vcount()
    if num_vert > 1:
        edgerat = float(curr_clust.ecount()) / (num_vert * (num_vert - 1) / 2)
    else:
        edgerat = 0
    # cohesion
    coh = curr_clust.cohesion()
    # max cohesion
    maxCohesion = plots.maxCohesion(curr_clust)
    # adhesion
    adhesion = curr_clust.adhesion()
    # max adhesion
    maxAdhesion = plots.maxAdhesion(curr_clust)
    # singleton stuff
    num_edg = curr_clust.ecount()
    # add cluster's values to dataframe
    clusterDFnew.loc[len(clusterDFnew)] = [
        i, num_vert, num_edg, meanDeg, meanStrg, cliq, edgerat, coh,
        maxCohesion, adhesion, maxAdhesion, meanPathlen, element1, element2
    ]
# add normalized cluster size
minV = float(min(clusterDFnew["number of vertices"]))
maxV = max(clusterDFnew["number of vertices"])
clusterDFnew["num of vertices (norm)"] = clusterDFnew[
    "number of vertices"].apply(lambda i: (i - minV) / (maxV - minV))
# add normalized edge ratio
clusterDFnew["edge ratio (norm)"] = clusterDFnew[
    "num of vertices (norm)"] * clusterDFnew["edge ratio"]

# Plotting
# Keep most "interesting"" clusters in DF
bestClusts = pd.DataFrame(columns=[i for i in plotting_metrics if i != "path"])
clusterDFnew[plotting_metrics] = clusterDFnew[plotting_metrics].apply(
    pd.to_numeric)
clusterDFnew.to_csv(out + sample_name + "_cluster_df_reconnected.csv", sep=",")
plotReconnectedDF = plots.assignRelevantText(
        pd.merge(
            left=clusterDFnew.reset_index(),
            right=clusterVerticesDF.reset_index(),
            on="cluster"))
for col in plotting_metrics:
    currDF = clusterDFnew.sort_values(col).reset_index(drop=True).reset_index()
    if col != "path":
        bestClusts[col] = currDF.tail(20)["cluster"]
    fig = currDF.plot.scatter(
        x='index',
        y=col,
        s=0.5,
        title=sample_name + ": " + col + " of clusters (reconnected)",
        grid=True)
    if col != "path length":
        x = currDF[currDF['element1']]
        if len(x) > 0:
            plt.text(x.index.item(), x[col].item(), "CCND1")
        x = currDF[currDF['element2']]
        if len(x) > 0:
            plt.text(x.index.item(), x[col].item(), "element2")
    plt.xlabel("rank")
    plt.ylabel(col)
    plt.savefig(pp, format="pdf")
    plotReconnectedDF[col + " index"] = clusterDFnew[col].argsort().argsort()

selectionReconnected = alt.selection_interval()
baseReconnected = plots.plotBase(plotReconnectedDF, selectionReconnected)
interactivePlotsReconnected = [
    plots.addTextEncodingRank(plotReconnectedDF, x,
                              plots.addSortedEncoding(baseReconnected, x))
    for x in plotting_metrics
]

# Join and save interactive plots
joinedPlotsReconnected = plots.plotSquare(interactivePlotsReconnected)
joinedPlotsReconnected.save(
    ppInteractive.replace(".html", "_reconnected.html"))

plotReconnectedDF.to_csv(
    out + sample_name + "_cluster_plot_df_reconnected.csv", index=False)

pp.close()
u, c = np.unique(bestClusts.values.flatten().tolist(), return_counts=True)
pd.DataFrame({
    "clust": u,
    "count": c
}).sort_values(
    "count", ascending=False).to_csv(
        out + sample_name + "_cluster_rankings.csv", sep=",", index=False)

