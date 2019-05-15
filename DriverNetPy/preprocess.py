#!/usr/bin/env python

'''
build bipartite graph, three inputs needed:
1. gene-gene interactions
2. patient-mutation matrix
3. patient-outlier matrix

'''

import sys
import logging
import pandas as pd
import networkx as nx

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def read_ggi(ggi_matirx=None):
    '''
    read gene-gene interaction matrix
    '''
    ggi = pd.read_csv(ggi_matirx, index_col=0)
    mapping = {i:j for i,j in enumerate(ggi.index)}
    ggi = nx.from_numpy_matrix(ggi.values)
    ggi = nx.relabel_nodes(ggi, mapping, copy=False)
    ggi.remove_edges_from(ggi.selfloop_edges())
    return ggi


def read_sample_out(patOutMatrix=None):
    '''
    read expression outlier matrix
    '''
    outlier_dict = {}  # a dict whose keys are smaples and values are sets containing genes

    df = pd.read_csv(patOutMatrix, index_col=0, header=0)
    for ix, row in df.iterrows():
        outlier_dict[ix] = set()
        for col_name in df.columns:
            if row[col_name] == True:
                outlier_dict[ix].add(col_name)

    return outlier_dict


def read_sample_mut(patMutMatrix=None):
    '''
    read patient-gene mutation matrix
    '''
    mutation_dict = {}  # a dict whose keys are samples and values are sets containing genes

    df = pd.read_csv(patMutMatrix, index_col=0, header=0)
    for ix, row in df.iterrows():
        mutation_dict[ix] = set()
        for col_name in df.columns:
            if row[col_name] == 1:
                mutation_dict[ix].add(col_name)

    return mutation_dict


def get_bipartite(ggi=None, outlier_dict=None, mutation_dict=None):
    '''
    get bipartite graph from thee inputs: ggi, outliers, mutated genes
    '''
    edges = dict()  # edges for the DriverNet graph. keys are samples, and values are events.

    for sample in outlier_dict:

        if not outlier_dict[sample]:    # sample contains no exp outliers
            continue

        if sample not in mutation_dict: # sample not in the mutation data
            continue

        edges[sample] = list()
        for x in outlier_dict[sample]:
            if x not in ggi:    # the outlier gene not in ggi
                continue

            # x itself
            if x in mutation_dict[sample]:
                edges[sample].append((x, x))

            # x neighbors
            x_neighbors = ggi.neighbors(x)  # genes interacting with x
            for y in x_neighbors:
                if y in mutation_dict[sample]:
                    edges[sample].append((y ,x))

    return edges


def preprocess(ggi_matirx=None, patMutMatrix=None, patOutMatrix=None):
    '''
    build a bipartite graph:
    the left contains mutated genes
    the right contains exp outliers
    '''
    ggi = read_ggi(ggi_matirx)
    logger.info("Reading gene-gene interacions done.")

    outlier_dict = read_sample_out(patOutMatrix)
    logger.info("Reading sample-gene expression outliers done.")

    mutation_dict = read_sample_mut(patMutMatrix)
    logger.info("Reading sample-gene mutation done.")

    edges = get_bipartite(ggi, outlier_dict, mutation_dict)
    logger.info("Building bipartite graph done.")

    return edges, ggi, outlier_dict, mutation_dict


def test():
    from collections import Counter

    ggi_matirx = 'sampleInfluenceGraph.txt'
    sample_mutation_matrix = 'samplePatientMutationMatrix.txt'
    sample_exp_outlier = 'samplePatientOutlierMatrix.txt'
    edges, ggi, outlier_dict, mutation_dict = preprocess(ggi_matirx, sample_mutation_matrix, sample_exp_outlier)

    candidates = Counter()
    for sample in edges:
        for edge in edges[sample]:
            candidates[edge[0]] += 1

    candidates = sorted(candidates.items(), key=lambda x:x[1], reverse=True)
    print(candidates)


if __name__ == "__main__":

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

