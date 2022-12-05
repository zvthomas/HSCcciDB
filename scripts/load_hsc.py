from tkinter import E
from interactions.models import Pathway, Receptor, Ligand, cellClas, pathwayAndCelltype, pathwayCorrelations
import pandas as pd
import numpy as np
# import all the code to properly extract relevant data from cellchat outputs

non_barcodes = ["Chondrocytes","EC-Arteriar","EC-Arteriolar","EC-Sinusoidal",
                        "Fibroblasts","MSPC-Adipo","MSPC-Osteo","Myofibroblasts","Osteo",
                        "Osteoblasts","Pericytes","Schwann-cells","Smooth-muscle", 
                       'CLP', 'lin-', 'MEP', 'CMP', 'GMP', 'MPP-', 'HSC', 'B cell', 'Dendritic cells',
                        'Eo/Baso prog.', 'Ery prog.', 'Ery/Mk prog.', 'Erythroblasts', 'Gran/Mono prog.',
                        'LMPPs', 'Mk prog.', 'Mono prog.', 'Monocytes', 'NK cells', 'Neutro prog.', 'Neutrophils',
                        'T cells', 'large pre-B.', 'pro-B', 'small pre-B.']

file_names = ["Chondrocytes", "ECArteriar", "ECArteriolar", "ECSinusoidal","Fibroblasts",
                  "MSPCAdipo", "MSPCOsteo", "Myofibroblasts", "Osteo", "Osteoblasts",
                  "Pericytes", "Schwanncells", "Smoothmuscle","CLP","lin","MEP","CMP","GMP","MPP",
                  "HSC","B cell","Dendritic cells","EoBaso prog","Ery prog","EryMk prog",
                  "Erythroblasts","GranMono prog","LMPPs","Mk prog","Mono prog","Monocytes",
                  "NK cells","Neutro prog","Neutrophils","T cells","large preB","proB","small preB"]

def process_csv(filename):
    pp = pd.read_csv(filename, header  = 0)
    pp = pp.set_index('Unnamed: 0')
    pp.index.names = [None]
    return pp

def process_path_file(filename, paths, thresh, celltype, df, barcodes = None):
    
    samp = process_csv(filename)
    samp = samp.T
    samp = samp.drop(columns = non_barcodes)
    if barcodes:
        samp = samp[barcodes]
    samp = samp.T
    samp = samp[[p for p in paths if p in samp.columns]]
    
    # Drop Pathways with less than 15% active HSC
    for col in samp.columns:
        zc = (samp[col] == 0).sum()
        if zc/len(samp) > ( 1 - thresh ):
            samp = samp.drop([col], axis = 1)
    
    d_to_Series = {}
    
    # Calculate Non-Zero Fraction and Non-Zero Average
    for col in samp.columns:
        zc = (samp[col] == 0).sum()
        nzf = (1 - zc/len(samp))
        nza = samp[col].mean()*len(samp)/(len(samp)-zc)
        d_to_Series[col] = [nzf, nza]
    
    if celltype == "Osteo":
        celltype = "Osteoclasts"

    ser = pd.Series(data = d_to_Series, index = d_to_Series.keys(), name = celltype)
    
    df = df.append(ser)
    
    return df

def get_cell_code(rowname):
    if rowname in ["Chondrocytes","EC-Arteriar","EC-Arteriolar","EC-Sinusoidal",
                        "Fibroblasts","MSPC-Adipo","MSPC-Osteo","Myofibroblasts","Osteoclasts",
                        "Osteoblasts","Pericytes","Schwann-cells","Smooth-muscle"] :
        return 'n'
    if rowname in ['B cell', 'Dendritic cells',
                'Erythroblasts', 'Monocytes',
                'NK cells', 'Neutrophils', 'T cells', 'large pre-B.', 'pro-B', 'small pre-B.']:
        return 'i'
    return 'p'

def drop_duplicate_rows(df):
    df = (df.reset_index().drop_duplicates(subset='index',keep='last').set_index('index').sort_index())
    return df

pathAnnot = pd.read_csv("P:/zthomas/Intercellular Interactions Project/cluster_analysis/path_annotation.csv")
signal_code = {'Cell-Cell Contact':'c','Secreted Signaling':'s',"ECM-Receptor":'e'}
def get_path_annot(path):
    subdf = pathAnnot[pathAnnot['paths'] == path]
    return signal_code[list(subdf['annot'])[0]]

def run():

    # Set up pathways and cell types and their pairs

    Pathway.objects.all().delete()
    cellClas.objects.all().delete()
    pathwayAndCelltype.objects.all().delete()
    Ligand.objects.all().delete()
    Receptor.objects.all().delete()
    pathwayCorrelations.objects.all().delete()

    path_annot = pd.read_csv('P:/zthomas/Intercellular Interactions Project/cluster_analysis/path_annotation.csv', sep = ',')
    path_annot = path_annot.set_index("paths")
    path_annot.index.names = [None]

    path_path = "P:/zthomas/Intercellular Interactions Project/cluster_analysis/general_niche/all_pathways/"

    interaction_genes = process_csv("P:/zthomas/Intercellular Interactions Project/interaction_genes.csv")

    filepathS_names = [path_path + 'hsc_to_' + f + ".csv" for f in file_names]
    filepathR_names = [path_path + f + "_to_hsc.csv" for f in file_names]
    signal_code = {'Cell-Cell Contact':'c','Secreted Signaling':'s',"ECM-Receptor":'e'}
    for signal_type in ['Cell-Cell Contact','Secreted Signaling',"ECM-Receptor"]:
        
        paths = list(path_annot[path_annot['annot'] == signal_type].index)

        NZFA_S = pd.DataFrame()
        NZFA_R = pd.DataFrame()

        for i in range(len(filepathS_names)):
            NZFA_S = process_path_file(filepathS_names[i], paths, 0, non_barcodes[i], NZFA_S)
        NZFA_S = drop_duplicate_rows(NZFA_S)

        for i in range(len(filepathR_names)):
            NZFA_R = process_path_file(filepathR_names[i], paths, 0, non_barcodes[i], NZFA_R)
        NZFA_R = drop_duplicate_rows(NZFA_R)

        for row in NZFA_S.index:
            for col in NZFA_S.columns:
                nzfa = NZFA_S.loc[row, col]
                
                if nzfa[0] == 0.0:
                    nzfa[1] = 0.0
                
                pathway, _ = Pathway.objects.get_or_create(name = col,
                                                        interaction_type = signal_code[signal_type],
                                                        hscPercent = 80, evidences = "")
                
                ct, _ = cellClas.objects.get_or_create(name = row, cell_type=get_cell_code(row))
                ct.sendingPathways.add(pathway)

                pathwayAndCelltype.objects.get_or_create(pathway = pathway,
                                                celltype = ct,
                                                hscPercent = nzfa[0],
                                                averageScore = nzfa[1],
                                                sorr = 's')

        
        for row in NZFA_R.index:
            for col in NZFA_R.columns:
                nzfa = NZFA_R.loc[row, col]
                
                if nzfa[0] == 0.0:
                    nzfa[1] = 0.0

                pathway, _ = Pathway.objects.get_or_create(name = col,
                                                        interaction_type = signal_code[signal_type],
                                                        hscPercent = 80, evidences = "")
                
                ct, _ = cellClas.objects.get_or_create(name = row, cell_type=get_cell_code(row))
                ct.receivingPathways.add(pathway)

                pathwayAndCelltype.objects.get_or_create(pathway = pathway,
                                                celltype = ct,
                                                hscPercent = nzfa[0],
                                                averageScore = nzfa[1],
                                                sorr = 'r')
                
    # Set up ligands and receptors

    
    interaction_genes = process_csv("P:/zthomas/Intercellular Interactions Project/interaction_genes.csv")

    for row in interaction_genes.index:

        pairs = interaction_genes.loc[row, ]
        evidence = pairs['evidence']
        ligand, _ = Ligand.objects.get_or_create(name = pairs['ligand'])
        receptor, _ = Receptor.objects.get_or_create(name = pairs['receptor'])
        pathway, _ = Pathway.objects.get_or_create(name = pairs['pathway_name'], interaction_type = get_path_annot(pairs['pathway_name']), hscPercent = 80)
        ligand.receptors.add(receptor)
        ligand.pathways.add(pathway)

        receptor.ligands.add(ligand)
        receptor.pathways.add(pathway)

        pathway.ligands.add(ligand)
        pathway.receptors.add(receptor)

        pathway.evidences = pathway.evidences + evidence + ";"
        pathway.save()

        receptor.evidences = receptor.evidences + evidence + ";"
        ligand.evidences = ligand.evidences + evidence + ";"
        receptor.save()
        ligand.save()



    #for row in NZFA_S.index:
    #    for col in NZFA_S.columns:
    #        print(row + "-" + col)
    
    #Pathway.objects.all().delete()
    #cellClas.objects.all().delete()
    #pathwayAndCelltype.objects.all().delete()

    #pathway = Pathway.objects.get_or_create(name=,
    #                   interaction_type=,
    #                   hscPercent=,)
    #pathway.save()

    #cell_type = (n(iche), i(mmune), p(rogenitors))
    #cellClass = cellClas(name = ,
    #                       sendingPathways = ,
    #                       receivingPathways = ,
    #                       cell_type = ,)

    #cell_type_pathway_pair = pathwayAndCelltype(pathway = ,
    #                                            celltype = ,
    #                                            hscPercent = ,
    #                                            sorr = )

    pos_correlations = pd.read_csv("P:/zthomas/Intercellular Interactions Project/cluster_analysis/pos_correlations_hide.csv")
    neg_correlations = pd.read_csv("P:/zthomas/Intercellular Interactions Project/cluster_analysis/neg_correlations_hide.csv")
        
    for i in [pos_correlations, neg_correlations]:
            for row in i.index:
                cor = i.loc[row,]

                p1 = cor['node1'].split()[1:][0]
                p1a = cor['node1'].split()[0][0]
                p2 = cor['node2'].split()[1:][0]
                p2a = cor['node2'].split()[0][0]

                corr = cor['Correlation']
                pval = cor['Adj. P-val']

                pathway1, _ = Pathway.objects.get_or_create(name = p1)
                pathway2, _ = Pathway.objects.get_or_create(name = p2)

                pathwayCorrelations.objects.get_or_create(pathway1 = pathway1,
                                                        pathway2 = pathway2,
                                                        p1a = p1a,
                                                        p2a = p2a,
                                                        correlation = corr,
                                                        pval = pval)