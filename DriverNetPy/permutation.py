#!/usr/bin/env python

'''
permutation the data to get emperical cover of each 

method: shuffle the mutated genes and outlier genes of each patient (as the way implemented in DriverNet)
'''

import sys
import logging
from collections import OrderedDict
from copy import deepcopy
import networkx as nx
from multiprocessing import cpu_count, Pool
import random

from preprocess import get_bipartite, preprocess
from greedy import greedy

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def mean(data=None):
    """Return the sample arithmetic mean of data.
    http://stackoverflow.com/a/27758326/632242
    """
    n = len(data)
    return sum(data)/float(n)

def _ss(data=None):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data=None):
    """Calculates the population standard deviation."""
    n = len(data)
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5


def multiple_testing_correction(pvalues=None, correction_type="Benjamini-Hochberg"):
    """
    Copyright 2017 Francisco Pina Martins <f.pinamartins@gmail.com>

    Taken from https://stackoverflow.com/a/21739593/3091595, remove numpy dependence

    @parameter pvalues - a list of pvalues
    @parameter correction_type - pvalue correction method

    @return qvalues - a list of qvalues
    """
    n = len(pvalues)
    qvalues = [0]*n
    if correction_type == "Bonferroni":
        qvalues = n * pvalues

    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (n-rank) * pvalue

    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]

    return qvalues


def random_net(G=None):
    '''
    Random the network but keep the original node degree roughly

    @parameter G - graph

    @return GG - random graph that only contains genes
    '''
    sequence = [d for (_, d) in G.degree]
    GG = nx.configuration_model(sequence)
    GG = nx.Graph(GG)
    GG.remove_edges_from(GG.selfloop_edges())

    # re-label
    mapping = {i:j for i,j in enumerate(G.nodes)}
    GG = nx.relabel_nodes(GG, mapping, copy=False)

    for node in GG.nodes:
        GG.nodes[node]['type'] = 'gene'

    return GG


def shuffle_data(input_dict=None):
    '''
    '''
    random_dict = dict()
    random_list = list()
    input_len = list()
    for sample in input_dict:
        random_list.extend(list(input_dict[sample]))
        input_len.append(len(input_dict[sample]))

    random.shuffle(input_len)    # random size, not elegant

    for ix, sample in enumerate(list(input_dict.keys())):
        if sample not in random_dict:
            random_dict[sample] = set()
        while len(random_dict[sample]) < input_len[ix]:
            random_dict[sample].add(random.choice(random_list))

    return random_dict


def compute_p(drivers=None, random_res=None):
    '''
    compute p-values
    '''
    pvalues_dict = {key:1 for key in drivers}

    tmp_dict = {key:0 for key in drivers}
    total_random_drivers = 0
    for i in random_res:
        total_random_drivers += len(random_res[i])
        for driver in drivers:
            for random_driver in random_res[i]:
                if random_res[i][random_driver] > drivers[driver]:
                    tmp_dict[driver] += 1

    for driver in tmp_dict:
        pvalues_dict[driver] = float(tmp_dict[driver])/total_random_drivers

    qvalues = multiple_testing_correction(list(pvalues_dict.values()))

    qvalues_dict = {key:qvalues[ix] for ix, key in enumerate(list(pvalues_dict.keys()))}

    return pvalues_dict, qvalues_dict


def permutation_helper(ggi=None, outlier_dict=None, mutation_dict=None, notPurturbGraph=False, purturbData=False):
    '''
    '''
    if purturbData:
        random_outlier_dict = shuffle_data(outlier_dict)
        random_mutation_dict = shuffle_data(mutation_dict)
    else:
        random_outlier_dict = deepcopy(outlier_dict)
        random_mutation_dict = deepcopy(mutation_dict)

    if not notPurturbGraph:
        random_ggi = random_net(ggi)
    else:
        random_ggi = deepcopy(ggi)

    random_edges = get_bipartite(random_ggi, random_outlier_dict, random_mutation_dict)

    random_drivers, random_events = greedy(random_edges, random_mutation_dict)

    return random_drivers


def permutation(ggi=None, outlier_dict=None, mutation_dict=None, permutation_time=500,
                threads=None, notPurturbGraph=False, purturbData=False):
    '''
    '''
    if not threads:
        threads = cpu_count()

    logger.info('Start permutation, using %s threads.' % threads)

    pool = Pool(threads)
    results = []
    for i in range(permutation_time):
        results.append(pool.apply_async(permutation_helper, args=(ggi,outlier_dict,mutation_dict,notPurturbGraph,purturbData,)))

    pool.close()
    pool.join()

    logger.info('Permutation done.')

    random_res = {}
    for ix, res in enumerate(results):
        random_res[ix] = res.get()

    return random_res


def test():

    ggi_matirx = 'sampleInfluenceGraph.txt'
    sample_mutation_matrix = 'samplePatientMutationMatrix.txt'
    sample_exp_outlier = 'samplePatientOutlierMatrix.txt'
    edges, ggi, outlier_dict, mutation_dict = preprocess(ggi_matirx, sample_mutation_matrix, sample_exp_outlier)

    drivers, actual_events = greedy(edges, mutation_dict)

    random_res = permutation(ggi, outlier_dict, mutation_dict)

    pvalues, qvalues = compute_p(drivers, random_res)

    print(qvalues)


if __name__ == "__main__":

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

