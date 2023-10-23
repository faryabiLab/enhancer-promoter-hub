# Interactive plots for HiChIP analysis
# By
# Gregory W. Schwartz

import functools as f
from math import sqrt

import altair as alt
import igraph
import mygene
import pandas as pd


def publishTheme():
    # Typography
    font = "Arial"
    labelFont = "Arial"
    sourceFont = "Arial"
    fontSize = 12  # 12 for publishing

    # Axes
    axisColor = "#000000"

    return {
        "config": {
            "view": {
                "height": 300,  # 80 for publishing
                "width": 400,  # 100 for publishing
                "strokeOpacity": 0,  # Despine
            },
            "style": {
                "guide-title": {
                    "fontSize": fontSize
                },
                "guide-label": {
                    "fontSize": fontSize
                }
            },
            "title": {
                "fontSize": fontSize,
                "font": font,
                "fontColor": "#000000",
                "fontWeight": "normal",
            },
            "axisX": {
                "grid": False,
                "domainColor": "#000000",
                "labelFont": labelFont,
                "labelFontSize": fontSize,
                "labelAngle": 0,
                "tickColor": axisColor,
                "titleFont": font,
                "titleFontSize": fontSize,
                "titleFontWeight": "normal",
            },
            "axisY": {
                "grid": False,
                "domainColor": "#000000",
                "labelFont": labelFont,
                "labelFontSize": fontSize,
                "labelAngle": 0,
                "tickColor": axisColor,
                "titleFont": font,
                "titleFontSize": fontSize,
                "titleFontWeight": "normal",
            },
        }
    }


# Plot an interactive ranked list of points.
def plotBase(df, selection):
    base = alt.Chart(df).mark_circle().encode(
        tooltip=["item"],
        opacity=alt.condition(
            selection, alt.value(1), alt.value(0.1), legend=None),
        size=alt.condition(
            selection, alt.value(20), alt.value(10), legend=None),
        color=alt.condition(
            selection, alt.value("red"), alt.value("gray"),
            legend=None)).add_selection(selection)
    return base


# Add a sorted encoding to a chart.
def addSortedEncoding(chart, col):
    return chart.encode(x=alt.X(col + " index", title="Rank"), y=col)


# Add text to a chart, specifically Element1 and Element2 here.
def addTextEncodingRank(df, col, chart):
    text = alt.Chart(df).mark_text(
        align="left", baseline="top", dx=7).encode(
            x=alt.X(col + " index", title="Rank"), y=col,
            text="text")
    #.transform_filter((alt.datum.text == "MYC") | (alt.datum.text == "SOX9"))

    return chart + text


# Chunk a list (https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks)
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


# Join plots into a square.
def plotSquare(xs):

    chunked = chunks(xs, round(sqrt(len(xs))))
    horiz = map(lambda ls: f.reduce(lambda acc, x: acc | x, ls), chunked)
    vertHoriz = f.reduce(lambda acc, ls: acc & ls, horiz)

    return vertHoriz


# Get the cluster with list of nodes within data frame.
def getClusterVerticesDF(df):
    return df.groupby("cluster").aggregate(lambda xs: f.reduce(
        lambda acc, x: acc + " " + x, xs))


# [OPTIONAL] Assign relevant text, specifically for MYC and SOX9 here.
def assignRelevantText(df):
    df["text"] = None
#     df.loc[df["myc"], "text"] = "MYC"
#     df.loc[df["sox9"], "text"] = "SOX9"

    return df


# Return True if a gene in the string is notch dependent.
def isNotchDependent(valid, x):
    return any(map(lambda g: g in valid, x.split(" ")))


# Assign notch dependent genes based on a set of valid genes.
def assignNotchDependence(valid, df):
    df["notchDependent"] = list(
        map(lambda x: isNotchDependent(valid, x), df["item"]))

    return df


# Convert Ensembl to symbols in cluster data frame. Ignore other IDs.
def getSymbols(df):
    mg = mygene.MyGeneInfo()
    symbols = [
        g["symbol"] if "symbol" in g else df["item"][i] for i, g in enumerate(
            mg.getgenes([x if "MACS" not in x else None for x in df["item"]],
                        fields="symbol"))
    ]

    df["item"] = symbols

    return df


# Get the cohesion of the max degree vertices.
def maxCohesion(g):
    degrees = g.degree()

    vertices = [
        x.index for _, x in sorted(zip(g.degree(), g.vs()), reverse=True)
    ][0:2]

    if len(degrees) < 2:
        return 0
    else:
        return g.cohesion(
            source=vertices[0], target=vertices[1], neighbors="ignore")


# Get the adhesion of the max degree vertices.
def maxAdhesion(g):
    degrees = g.degree()

    vertices = [
        x.index for _, x in sorted(zip(g.degree(), g.vs()), reverse=True)
    ][0:2]

    if len(degrees) < 2:
        return 0
    else:
        return g.adhesion(source=vertices[0], target=vertices[1])
