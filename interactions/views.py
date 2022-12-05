from django.shortcuts import render
from plotly.offline import plot
import plotly.graph_objects as go
from django.shortcuts import get_object_or_404
import networkx as nx
import json
import numpy as np
from matplotlib import cm
import matplotlib as mpl
from .models import Pathway, Receptor, Ligand, cellClas, pathwayAndCelltype, pathwayCorrelations
from django.views import generic
# Create your views here.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:45:07 2020
@author: aransil
"""

greys = cm.get_cmap('Greys', 10000)


import math
from typing import List
from itertools import chain

# n is number of points to place evenly in a circle centered around (0,0)
def circle(n):
    x = []
    y = []
    for i in range(0,n):
        angle = (i/n)*math.pi*2
        x.append(math.cos(angle))
        y.append(math.sin(angle))
    
    return x, y

# Start and end are lists defining start and end points
# Edge x and y are lists used to construct the graph
# arrowAngle and arrowLength define properties of the arrowhead
# arrowPos is None, 'middle' or 'end' based on where on the edge you want the arrow to appear
# arrowLength is the length of the arrowhead
# arrowAngle is the angle in degrees that the arrowhead makes with the edge
# dotSize is the plotly scatter dot size you are using (used to even out line spacing when you have a mix of edge lengths)
def addEdge(start, end, edge_x, edge_y, lengthFrac=1, arrowPos = None, arrowLength=0.025, arrowAngle = 30, dotSize=20):

    # Get start and end cartesian coordinates
    x0, y0 = start
    x1, y1 = end

    # Incorporate the fraction of this segment covered by a dot into total reduction
    length = math.sqrt( (x1-x0)**2 + (y1-y0)**2 )
    dotSizeConversion = .0565/20 # length units per dot size
    convertedDotDiameter = dotSize * dotSizeConversion
    lengthFracReduction = convertedDotDiameter / length
    lengthFrac = lengthFrac - lengthFracReduction

    # If the line segment should not cover the entire distance, get actual start and end coords
    skipX = (x1-x0)*(1-lengthFrac)
    skipY = (y1-y0)*(1-lengthFrac)
    x0 = x0 + skipX/2
    x1 = x1 - skipX/2
    y0 = y0 + skipY/2
    y1 = y1 - skipY/2

    # Append line corresponding to the edge
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None) # Prevents a line being drawn from end of this edge to start of next edge
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

    # Draw arrow
    if not arrowPos == None:

        # Find the point of the arrow; assume is at end unless told middle
        pointx = x1
        pointy = y1

        eta = math.degrees(math.atan((x1-x0)/(y1-y0))) if y1!=y0 else 90.0

        if arrowPos == 'middle' or arrowPos == 'mid':
            pointx = x0 + (x1-x0)/2
            pointy = y0 + (y1-y0)/2

        # Find the directions the arrows are pointing
        signx = (x1-x0)/abs(x1-x0) if x1!=x0 else +1    #verify this once
        signy = (y1-y0)/abs(y1-y0) if y1!=y0 else +1    #verified

        # Append first arrowhead
        dx = arrowLength * math.sin(math.radians(eta + arrowAngle))
        dy = arrowLength * math.cos(math.radians(eta + arrowAngle))
        edge_x.append(pointx)
        edge_x.append(pointx - signx**2 * signy * dx)
        edge_x.append(None)
        edge_y.append(pointy)
        edge_y.append(pointy - signx**2 * signy * dy)
        edge_y.append(None)

        # And second arrowhead
        dx = arrowLength * math.sin(math.radians(eta - arrowAngle))
        dy = arrowLength * math.cos(math.radians(eta - arrowAngle))
        edge_x.append(pointx)
        edge_x.append(pointx - signx**2 * signy * dx)
        edge_x.append(None)
        edge_y.append(pointy)
        edge_y.append(pointy - signx**2 * signy * dy)
        edge_y.append(None)


    return edge_x, edge_y

def add_arrows(source_x: List[float], target_x: List[float], source_y: List[float], target_y: List[float],
               arrowLength=0.025, arrowAngle=30):
    pointx = list(map(lambda x: x[0] + (x[1] - x[0]) / 2, zip(source_x, target_x)))
    pointy = list(map(lambda x: x[0] + (x[1] - x[0]) / 2, zip(source_y, target_y)))
    etas = list(map(lambda x: math.degrees(math.atan((x[1] - x[0]) / (x[3] - x[2]))),
                    zip(source_x, target_x, source_y, target_y)))

    signx = list(map(lambda x: (x[1] - x[0]) / abs(x[1] - x[0]), zip(source_x, target_x)))
    signy = list(map(lambda x: (x[1] - x[0]) / abs(x[1] - x[0]), zip(source_y, target_y)))

    dx = list(map(lambda x: arrowLength * math.sin(math.radians(x + arrowAngle)), etas))
    dy = list(map(lambda x: arrowLength * math.cos(math.radians(x + arrowAngle)), etas))
    none_spacer = [None for _ in range(len(pointx))]
    arrow_line_x = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointx, signx, signy, dx)))
    arrow_line_y = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointy, signx, signy, dy)))

    arrow_line_1x_coords = list(chain(*zip(pointx, arrow_line_x, none_spacer)))
    arrow_line_1y_coords = list(chain(*zip(pointy, arrow_line_y, none_spacer)))

    dx = list(map(lambda x: arrowLength * math.sin(math.radians(x - arrowAngle)), etas))
    dy = list(map(lambda x: arrowLength * math.cos(math.radians(x - arrowAngle)), etas))
    none_spacer = [None for _ in range(len(pointx))]
    arrow_line_x = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointx, signx, signy, dx)))
    arrow_line_y = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointy, signx, signy, dy)))

    arrow_line_2x_coords = list(chain(*zip(pointx, arrow_line_x, none_spacer)))
    arrow_line_2y_coords = list(chain(*zip(pointy, arrow_line_y, none_spacer)))

    x_arrows = arrow_line_1x_coords + arrow_line_2x_coords
    y_arrows = arrow_line_1y_coords + arrow_line_2y_coords

    return x_arrows, y_arrows

def normalize_dot_size(min_size, max_size, val, min_p, max_p):
    val = val*(max_size-min_size)
    val = val + min_size

    return val

def get_ec(lab, ct_or_p):

    cdict = {'n':'#D12626', 'i':'#51C206', 'p':'#0855CA'}
    if ct_or_p != 'cts':
        cdict = {'s':'#056E5D', 'c':'#7905D7', 'e':'#202A25'}

    return cdict[lab]

def index(request):
    """View function for home page"""

    # Generate counts of pathways and cell types
    num_pathways = Pathway.objects.filter().count()
    num_ss = Pathway.objects.filter(interaction_type="s").count()
    num_ccc = Pathway.objects.filter(interaction_type="c").count()
    num_ecmr = Pathway.objects.filter(interaction_type="e").count()


    num_cell_types = cellClas.objects.all().count()
    num_n = cellClas.objects.filter(cell_type="n").count()
    num_i = cellClas.objects.filter(cell_type="i").count()
    num_p = cellClas.objects.filter(cell_type="p").count()

    # cell type - pathway paris
    num_celltype_pathway_pairs = pathwayAndCelltype.objects.filter(hscPercent__gt=0.049).count()
    num_ligands = Ligand.objects.count()
    num_receptors = Receptor.objects.count()

    num_correlations = pathwayCorrelations.objects.count()
    pos_correlations = pathwayCorrelations.objects.filter(correlation__gt=0).count()
    neg_correlations = pathwayCorrelations.objects.filter(correlation__lt=0).count()

    for_table = pathwayAndCelltype.objects.all().order_by('-averageScore')[:20]

    context = {
        'num_pathways' : num_pathways,
        'num_cell_types' : num_cell_types,
        'num_celltype_pathway_pairs' : num_celltype_pathway_pairs,
        'num_ss': num_ss,
        'num_ccc': num_ccc,
        'num_ecmr': num_ecmr,
        'num_n': num_n,
        'num_i': num_i,
        'num_p': num_p,
        'num_ligands': num_ligands,
        'num_receptors': num_receptors,
        'num_correlations' : num_correlations,
        'pos_correlations' : pos_correlations,
        'neg_correlations' : neg_correlations,
        'for_table' : for_table
    }

    # Render the HTML template index.html with the data in the context variable
    return render(request, 'index.html', context = context)

def help(request):
    return render(request, 'about.html')

def get_evidence_list(evidences):
    
    evidenceList = evidences.split(";")
        
    evidenceList.pop()
        
    for evidence in evidenceList:
        if ", " in evidence:
            e = evidence.split(", ")
            evidenceList.remove(evidence)
            for item in e:
                evidenceList.append(item)
        
    pmids = {}
    keggs = {}
    pmcs = {}
    
    for i in evidenceList:
        i = i.replace(" ", "")
        if "PMID" in i:
           pmids[i] = "https://pubmed.ncbi.nlm.nih.gov/" + i.split(':')[1]
        elif "KEGG" in i:
            keggs[i] = "https://www.genome.jp/dbget-bin/www_bget?pathway:" + i.split(':')[-1]
        elif "PMC" in i:
            pmcs[i] = "https://www.ncbi.nlm.nih.gov/pmc/articles/" + i
    
    return pmids, keggs, pmcs


def correlations(request):

    ligand_list = Ligand.objects.all()
    genemap, paths = make_gene_map(ligand_list)
    correlations = pathwayCorrelations.objects.all()
    plot_div = make_net_graph_corr(correlations)
    # write some function that parses it like the make_net_graph_JSON function to pass it to render...
    # in the html page, write some javascript that filters the graph based on the form submission
    context = {'genemap':genemap, 'paths':paths, 'plot_div' : plot_div}
    return render(request, 'correlations.html', context)

class PathwayListView(generic.ListView):
    model = Pathway

    def get_context_data(self, **kwargs):
        context = super(PathwayListView, self).get_context_data(**kwargs)
        context['forTable'] = pathwayAndCelltype.objects.all().order_by('-averageScore')[:20]
        return context

class PathwayDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):
        
        pathway = get_object_or_404(Pathway, pk = kwargs['pk'])
        pactsS = pathwayAndCelltype.objects.filter(pathway= pathway, sorr = 's', hscPercent__range=(0.05,2))
        pactsR = pathwayAndCelltype.objects.filter(pathway= pathway, sorr = 'r', hscPercent__range=(0.05,2))
        evidenceList = pathway.evidences
        pmids, keggs, pmcs = get_evidence_list(evidenceList)

        correlatedPathways1 = pathwayCorrelations.objects.filter(pathway1 = pathway)
        correlatedPathways2 = pathwayCorrelations.objects.filter(pathway2 = pathway)
        
        ligands = pathway.ligands.all()
        receptors = pathway.receptors.all()

        if len(pactsS) > 0 or len(pactsR) > 0:
            plot_div = make_net_graph_JSON(pactsS, pactsR, 'cts')
            #plot_div = make_net_graph_spread(pactsS, pactsR, 'cts')
            context = {'pathway': pathway, 'pactsS': pactsS, 'pactsR': pactsR,
                        'plot_div': plot_div, 'ligands':ligands, 'receptors': receptors,
                        'pmids' : pmids, 'keggs' : keggs, 'pmcs' : pmcs,
                        'correlations1' : correlatedPathways1, 'correlations2' : correlatedPathways2}
            return render(request, 'interactions/pathway_detail.html', context)
        
        else:
            return render(request, 'interactions/pathway_detail.html', context = {'pathway': pathway,'ligands':ligands, 'receptors': receptors})


class CellClassListView(generic.ListView):
    model = cellClas

    def get_context_data(self, **kwargs):
        context = super(CellClassListView, self).get_context_data(**kwargs)
        context['forTable'] = pathwayAndCelltype.objects.all().order_by('-averageScore')[:20]
        return context

class CellClassDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):
        
        cellclass = get_object_or_404(cellClas, pk = kwargs['pk'])
        pactsS = pathwayAndCelltype.objects.filter(celltype= cellclass, sorr = 's', hscPercent__range=(0.05,2))
        pactsR = pathwayAndCelltype.objects.filter(celltype= cellclass, sorr = 'r', hscPercent__range=(0.05,2))
        if len(pactsS) > 0 or len(pactsR) >0:
            plot_div = make_net_graph_JSON(pactsS, pactsR, 'paths')
            #plot_div = make_net_graph_spread(pactsS, pactsR, 'paths')
            context = {'cellclass': cellclass, 'pactsS': pactsS, 'pactsR': pactsR,
                        'plot_div': plot_div}
            return render(request, 'interactions/cellClas_detail.html', context)
        
        else:
            return render(request, 'interactions/cellClas_detail.html', context = {'pathway': cellclass})

class LigandDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):

        genename = get_object_or_404(Ligand, pk=kwargs['pk'])
        pathways = genename.pathways.all()
        paired = genename.receptors.all()
        evidenceList = genename.evidences
        pmids, keggs, pmcs = get_evidence_list(evidenceList)
        genemap, paths = make_gene_map_specific(pathways)

        gene_expression_link = "https://gexc.riken.jp/models/3/genes/" + str.lower(genename.name)

        context = {"genename": genename, 'pathways': pathways, 'paired': paired, 'link': gene_expression_link,
                    'lorr':'l', 'plot_div':genemap, 'paths':paths,
                    'pmids' : pmids, 'keggs' : keggs, 'pmcs' : pmcs}

        return render(request, 'interactions/gene_detail.html', context = context)


class ReceptorDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):

        genename = get_object_or_404(Receptor, pk=kwargs['pk'])
        pathways = genename.pathways.all()
        paired = genename.ligands.all()
        evidenceList = genename.evidences
        pmids, keggs, pmcs = get_evidence_list(evidenceList)
        genemap, paths = make_gene_map_specific(pathways)

        gene_expression_link = "https://gexc.riken.jp/models/3/genes/" + str.lower(genename.name)

        context = {"genename": genename, 'pathways': pathways, 'paired': paired, 'link': gene_expression_link,
                    'lorr':'r', 'plot_div':genemap, 'paths':paths,
                    'pmids' : pmids, 'keggs' : keggs, 'pmcs' : pmcs}

        return render(request, 'interactions/gene_detail.html', context = context)



MIN_SIZE = 20
MAX_SIZE = 50

def make_net_graph(sending, receiving, ct_or_p = 'cts'):
    """ 
    View demonstrating how to display a graph object
    on a web page with Plotly. 
    """
    g = nx.Graph()

    sending_dict_avgscore_hscP = {}
    receiving_dict_avgscore_hscP = {}

    hscPs = []

    edge_weights = []
    
    title_name = ""
    center_cell = ""

    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(center_cell,cn+"_S")
        sending_dict_avgscore_hscP[cn] = [ct.averageScore]
        sending_dict_avgscore_hscP[cn].append(ct.hscPercent)
        sending_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
        

    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(cn+"_R", center_cell)
        receiving_dict_avgscore_hscP[cn] = [ct.averageScore]
        receiving_dict_avgscore_hscP[cn].append(ct.hscPercent)
        receiving_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)

    minP = min(hscPs)
    maxP = max(hscPs)
    ewMin = min(edge_weights)
    ewMax = max(edge_weights) 
    #edge_weights = [(i-ewMin)*(9)/(ewMax-ewMin) + 1 for i in edge_weights]
    x, y = circle(len(g.nodes()))
    
    pos = {}
    k = 0
    for i in g.nodes():
        if i != center_cell:
            pos[i] = [x[k], y[k]]
            k += 1
        else:
            pos[i] = [0,0]
        

    #pos = nx.planar_layout(g)

    edge_x = []
    edge_y = []
    edge_colors = []
    check = "_R" if ct_or_p == "cts" else "_S"
    for edge in g.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]


        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

        if check in edge[1]:
            edge_colors.append("orchid")
        else:
            edge_colors.append("aqua")
        '''
        if len(g.edges()) > 2:
            if check in edge[1]:
                edge_x, edge_y = addEdge([x1,y1], [x0,y0], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("orchid")
            else:
                edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("aqua")
        else:
            edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .95, None, .05, 30, 55)

    edge_trace = {}

    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*9:i*9+9],
                                                  y = edge_y[i*9:i*9+9],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')
    '''
    edge_trace = {}
    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*3:i*3+3],
                                                  y = edge_y[i*3:i*3+3],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')

    edge_trace = list(edge_trace.values())
    
    #edge_trace = go.Scatter(
    #x=edge_x, y=edge_y,
    #line=dict(width=2, color=edge_colors),
    #hoverinfo='none',
    #mode='lines')

    node_x, node_y = [], []
    for node in g.nodes():
        x,y = pos[node]
        node_x.append(x)
        node_y.append(y)

    marker_type = 'circle' if ct_or_p == "paths" else 'square'

    node_trace = go.Scatter(
        x = node_x, y = node_y,
        mode = 'markers',
        marker_symbol = marker_type,
        hoverinfo = 'text',
        marker=dict(
        showscale=True,
        colorscale='Greys',
        reversescale=False,
        color=[],
        line = dict(width = 5),
        size=20,
        colorbar=dict(
            thickness=15,
            title='Average Probability Score',
            xanchor='left',
            titleside='right'
        ),
    ))

    node_avgscore = []
    node_text = []
    node_sizes = []
    node_ecs = []
    for node in g.nodes:
        if node == center_cell:
            node_avgscore.append(max(edge_weights))
            node_sizes.append(50)
            node_text.append(center_cell)
            node_ecs.append('black')
        else:
            cell_name = node.split("_")[0]
            if "_S" in node:
                avgScore = sending_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = sending_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(sending_dict_avgscore_hscP[cell_name][2], ct_or_p))

            else:
                avgScore = receiving_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = receiving_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(receiving_dict_avgscore_hscP[cell_name][2], ct_or_p))
    
    node_trace.marker.color = node_avgscore
    node_trace.text = node_text
    node_trace.marker.size = node_sizes
    node_trace.marker.line.color = node_ecs

    #normalize_dot_size(min_size, max_size, val, min_p, max_p)

    edge_trace.append(node_trace)
    #data = [edge_trace, node_trace]
    layout={
                'title': '<b>' + title_name + ' Pathway (Hover for info)</b>',
                'titlefont_size':16,
                'showlegend':False,
                'hovermode':'closest',
                'margin':dict(b=5,l=5,r=5,t=40),
                'xaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'yaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'width':700
                
    }
    
    
    # Getting HTML needed to render the plot.
    plot_div = plot({'data': edge_trace, 'layout': layout}, 
                    output_type='div')

    return plot_div

def genes(request):
    ligand_list = Ligand.objects.all()
    receptor_list = Receptor.objects.all()
    genemap, paths = make_gene_map(ligand_list)
    context = {'ligand_list': ligand_list, 'receptor_list': receptor_list, 'plot_div':genemap, 'paths':paths}
    
    
    return render(request, 'genes.html', context= context)

def make_gene_map(ligand_list):
    
    ligands= []
    receptors = []
    pathway = []

    for l in ligand_list:
        for r in l.receptors.all():
            for p in l.pathways.all():
                ligands.append(l.name)
                receptors.append(r.name)
                pathway.append(p.name)

    toret = json.dumps([ligands, receptors, pathway, list(set(pathway))])
    upaths = list(set(pathway))
    upaths.sort()
    return toret, upaths

def make_gene_map_specific(pathway_list):
    
    ligands= []
    receptors = []
    pathway = []

    for p in pathway_list:
        for l in p.ligands.all():
            for r in l.receptors.all():
                ligands.append(l.name)
                receptors.append(r.name)
                pathway.append(p.name)

    toret = json.dumps([ligands, receptors, pathway, list(set(pathway))])
    upaths = list(set(pathway))
    upaths.sort()
    return toret, upaths

def make_net_graph_spread(sending, receiving, ct_or_p = 'cts'):
    """ 
    View demonstrating how to display a graph object
    on a web page with Plotly. 
    """
    g = nx.Graph()

    sending_dict_avgscore_hscP = {}
    receiving_dict_avgscore_hscP = {}

    hscPs = []

    edge_weights = []
    
    title_name = ""
    center_cell = ""

    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type
        g.add_node(center_cell)
        g.add_node(cn+"_S")
        sending_dict_avgscore_hscP[cn] = [ct.averageScore]
        sending_dict_avgscore_hscP[cn].append(ct.hscPercent)
        sending_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
        
    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type
        g.add_node(center_cell)
        g.add_node(cn+"_R")
        receiving_dict_avgscore_hscP[cn] = [ct.averageScore]
        receiving_dict_avgscore_hscP[cn].append(ct.hscPercent)
        receiving_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
    
    minP = min(hscPs)
    maxP = max(hscPs)

    x, y = circle(len(g.nodes())-1)
    
    pos = {}
    k = 0
    for i in g.nodes():
        if i != center_cell:
            pos[i] = [x[k], y[k]]
            k += 1
        else:
            pos[i] = [0,0]


    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(center_cell,cn+"_S")        

    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(cn+"_R", center_cell)


    edge_x = []
    edge_y = []
    edge_colors = []
    check = "_R" if ct_or_p == "cts" else "_S"
    for edge in g.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]


        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

        if check in edge[1]:
            edge_colors.append("orchid")
        else:
            edge_colors.append("aqua")
        '''
        if len(g.edges()) > 2:
            if check in edge[1]:
                edge_x, edge_y = addEdge([x1,y1], [x0,y0], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("orchid")
            else:
                edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("aqua")
        else:
            edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .95, None, .05, 30, 55)

    edge_trace = {}

    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*9:i*9+9],
                                                  y = edge_y[i*9:i*9+9],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')
    '''
    edge_trace = {}
    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*3:i*3+3],
                                                  y = edge_y[i*3:i*3+3],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')

    edge_trace = list(edge_trace.values())
    
    #edge_trace = go.Scatter(
    #x=edge_x, y=edge_y,
    #line=dict(width=2, color=edge_colors),
    #hoverinfo='none',
    #mode='lines')

    node_x, node_y = [], []
    for node in g.nodes():
        x,y = pos[node]
        node_x.append(x)
        node_y.append(y)

    marker_type = 'circle' if ct_or_p == "paths" else 'diamond'
    marker_labels = [i[:-2] for i in list(g.nodes())[1:]]
    marker_labels.insert(0, "")

    node_trace = go.Scatter(
        x = node_x, y = node_y,
        mode = 'markers+text',
        marker_symbol = marker_type,
        text = marker_labels,
        textposition="bottom center",
        hoverinfo = 'text',
        marker=dict(
        showscale=True,
        colorscale='Greys',
        reversescale=False,
        color=[],
        line = dict(width = 5),
        size=20,
        opacity = 1,
        colorbar=dict(
            thickness=15,
            title='Average Probability Score',
            xanchor='left',
            titleside='right'
        ),
    ))

    node_avgscore = []
    node_text = []
    node_sizes = []
    node_ecs = []
    for node in g.nodes:
        if node == center_cell:
            node_avgscore.append(max(edge_weights))
            node_sizes.append(50)
            node_text.append(center_cell)
            node_ecs.append('black')
        else:
            cell_name = node.split("_")[0]
            if "_S" in node:
                avgScore = sending_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = sending_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(sending_dict_avgscore_hscP[cell_name][2], ct_or_p))

            else:
                avgScore = receiving_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = receiving_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(receiving_dict_avgscore_hscP[cell_name][2], ct_or_p))
    
    node_trace.marker.color = node_avgscore
    node_trace.hovertext = node_text
    node_trace.marker.size = node_sizes
    node_trace.marker.line.color = node_ecs

    #normalize_dot_size(min_size, max_size, val, min_p, max_p)

    edge_trace.append(node_trace)
    #data = [edge_trace, node_trace]
    layout={
                'title': '<b>' + title_name + ' Pathway (Hover for info)</b>',
                'titlefont_size':16,
                'showlegend':False,
                'hovermode':'closest',
                'margin':dict(b=5,l=5,r=5,t=40),
                'xaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'yaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'width':700
                
    }
    
    
    # Getting HTML needed to render the plot.
    plot_div = plot({'data': edge_trace, 'layout': layout}, 
                    output_type='div')

    return plot_div

def make_net_graph_corr(correlations):

    p1names = []
    p1colors = []
    p1ids = []
    p1urls = []

    p2names = []
    p2colors = []
    p2ids = []
    p2urls = []

    corrs = []
    pvals = []
    edge_colors = []    

    p1names_n = []
    p1colors_n = []
    p1ids_n = []
    p1urls_n = []

    p2names_n = []
    p2colors_n = []
    p2ids_n = []
    p2urls_n = []

    corrs_n = []
    pvals_n = []
    edge_colors_n = []    

    for corr in correlations:
        pathway1_name = corr.pathway1.name
        pathway1_id = pathway1_name + "_" + corr.p1a.upper()
        pathway1_color = "purple" if corr.p1a.upper() == "S" else "green"
        pathway2_name = corr.pathway2.name
        pathway2_id = pathway2_name + "_" + corr.p2a.upper()
        pathway2_color = "purple" if corr.p2a.upper() == "S" else "green"
        corr_val = corr.correlation
        p_val = corr.pval
        edge_color = "red" if corr_val > 0 else "blue"

        if corr_val > 0:
            p1names.append(pathway1_name)
            p1colors.append(pathway1_color)
            p1ids.append(pathway1_id)
            p1urls.append(corr.pathway1.get_absolute_url)

            p2names.append(pathway2_name)
            p2colors.append(pathway2_color)
            p2ids.append(pathway2_id)
            p2urls.append(corr.pathway2.get_absolute_url)

            corrs.append(corr_val)
            pvals.append(p_val)
            edge_colors.append(edge_color)
        
        else:
            p1names_n.append(pathway1_name)
            p1colors_n.append(pathway1_color)
            p1ids_n.append(pathway1_id)
            p1urls_n.append(corr.pathway1.get_absolute_url)

            p2names_n.append(pathway2_name)
            p2colors_n.append(pathway2_color)
            p2ids_n.append(pathway2_id)
            p2urls_n.append(corr.pathway2.get_absolute_url)

            corrs_n.append(corr_val)
            pvals_n.append(p_val)
            edge_colors_n.append(edge_color)
    
    corrs, p1names, p1colors, p1ids, p1urls, p2names, p2colors, p2ids, p2urls,  pvals, edge_colors = zip(*sorted(zip(corrs, p1names, p1colors, p1ids, p1urls, p2names, p2colors, p2ids, p2urls, pvals, edge_colors), reverse=True))

    corrs_n, p1names_n, p1colors_n, p1ids_n, p1urls_n, p2names_n, p2colors_n, p2ids_n, p2urls_n, pvals_n, edge_colors_n = zip(*sorted(zip(corrs_n, p1names_n, p1colors_n, p1ids_n, p1urls_n, p2names_n, p2colors_n, p2ids_n, p2urls_n, pvals_n, edge_colors_n)))

    toret = json.dumps([list(x) for x in [corrs, p1names, p1colors, p1ids, p2names, p2colors, p2ids, pvals, edge_colors, corrs_n, p1names_n, p1colors_n, p1ids_n, p2names_n, p2colors_n, p2ids_n, pvals_n, edge_colors_n, p1urls, p2urls, p1urls_n, p2urls_n]])

    return toret

def make_net_graph_JSON(sending, receiving, ct_or_p = 'cts'):
    """ 
    View demonstrating how to display a graph object
    on a web page with Plotly. 

    list of lists

    list[0] = list of node names
    list[1] = list of "from" names
    list[2] = list of "to" names



    """
    max_score = 0
    min_score = 100
    node_names = []
    to_names = []
    from_names = []

    sending_dict_avgscore_hscP = {}
    receiving_dict_avgscore_hscP = {}

    hscPs = []

    edge_weights = []
    
    title_name = ""
    center_cell = ""

    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        if ct.averageScore > max_score:
            max_score = ct.averageScore
        if ct.averageScore < min_score:
            min_score = ct.averageScore 

        sending_dict_avgscore_hscP[cn] = [ct.averageScore]
        sending_dict_avgscore_hscP[cn].append(ct.hscPercent)
        sending_dict_avgscore_hscP[cn].append(t)
        if ct_or_p == 'cts':
            sending_dict_avgscore_hscP[cn].append(ct.celltype.get_absolute_url)
        else:
            sending_dict_avgscore_hscP[cn].append(ct.pathway.get_absolute_url)

        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
        
    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        if ct.averageScore > max_score:
            max_score = ct.averageScore
        if ct.averageScore < min_score:
            min_score = ct.averageScore

        receiving_dict_avgscore_hscP[cn] = [ct.averageScore]
        receiving_dict_avgscore_hscP[cn].append(ct.hscPercent)
        receiving_dict_avgscore_hscP[cn].append(t)
        if ct_or_p == 'cts':
            receiving_dict_avgscore_hscP[cn].append(ct.celltype.get_absolute_url)
        else:
            receiving_dict_avgscore_hscP[cn].append(ct.pathway.get_absolute_url)
        
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)


    node_names = list(set(node_names))
    
    node_names = []
    fill_color = []
    avgscore = []
    edges = []
    hscPercent = []
    height = []
    ecs = []
    urlList = []
    k = 0
    for i in sending_dict_avgscore_hscP.keys():
        node_names.append(i + "_" + str(k))
        ecs.append(get_ec(sending_dict_avgscore_hscP[i][2], ct_or_p))

        range = max_score - min_score
        a = sending_dict_avgscore_hscP[i][0]
        if range != 0:
            a = (a - min_score) / range
        range2 = 9999
        a = (a * range2)
        fill_color.append(mpl.colors.to_hex(greys(int(a))))

        height.append(sending_dict_avgscore_hscP[i][1]*60 + 10)
        
        avgscore.append(sending_dict_avgscore_hscP[i][0])
        hscPercent.append(100*sending_dict_avgscore_hscP[i][1])
        edges.append('#00ffff')
        to_names.append(i +"_" + str(k))
        from_names.append(center_cell)
        urlList.append(sending_dict_avgscore_hscP[i][3])
        k += 1
    for i in receiving_dict_avgscore_hscP.keys():
        node_names.append(i +"_" + str(k))
        ecs.append(get_ec(receiving_dict_avgscore_hscP[i][2], ct_or_p))


        range = max_score - min_score
        a = receiving_dict_avgscore_hscP[i][0]
        if range != 0:
            a = (a - min_score) / range
        range2 = 9999
        a = (a * range2)
        fill_color.append(mpl.colors.to_hex(greys(int(a))))

        height.append(receiving_dict_avgscore_hscP[i][1]*60 + 10)

        avgscore.append(receiving_dict_avgscore_hscP[i][0])
        hscPercent.append(100*receiving_dict_avgscore_hscP[i][1])
        from_names.append(i +"_" + str(k))
        edges.append('#da70d6')
        to_names.append(center_cell)
        urlList.append(receiving_dict_avgscore_hscP[i][3])
        k += 1
    node_names.append(center_cell)
    ecs.append('3 black')
    fill_color.append("black")
    avgscore.append("NaN")
    hscPercent.append("NaN")
    height.append(70)
    urlList.append("")

    toret = json.dumps([node_names, to_names, from_names, ecs, fill_color, avgscore, hscPercent, edges, height, urlList])

    return toret