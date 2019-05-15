#!/usr/bin/env python

'''
Use greedy algorithm to find the optimal set of genes which covers all expression-outliers

'''

import sys
import logging
from collections import Counter, OrderedDict
from copy import deepcopy

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def greedy(edges=None, mutation_dict=None, drivers=None, actual_events=None):
    '''
    select a set of genes using greedy algorithm
    '''
    edges_copy = deepcopy(edges)    # attention! use deepcopy to avoid the function changing original values

    if drivers is None:
        drivers = OrderedDict()
    if actual_events is None:
        actual_events = OrderedDict()

    candidates = Counter()
    for sample in edges_copy:
        for edge in edges_copy[sample]:
            candidates[edge[0]] += 1    # count the events covered by each gene

    if not candidates:  # end
        return drivers, actual_events

    highest = sorted(candidates.items(), key=lambda x:x[1], reverse=True)[0]
    #print(highest, len(drivers))

    drivers[highest[0]] = highest[1]
    actual_events[highest[0]] = dict()

    tmp = deepcopy(edges_copy)
    for sample in tmp:
        if highest[0] not in mutation_dict[sample]: # the hightest is not mutated
            continue
        else:
            covered = dict()    # covered by the highest
            for edge in tmp[sample]:
                if edge[0] == highest[0]:
                    actual_events[highest[0]][sample] = edge[1]
                    covered[edge[1]] = 0
            for edge in tmp[sample]:
                if edge[1] in covered:
                    edges_copy[sample].remove(edge)  # remove covered by the highest

    return greedy(edges_copy, mutation_dict, drivers, actual_events)


def test():
    from preprocess import preprocess

    ggi_matirx = 'sampleInfluenceGraph.txt'
    sample_mutation_matrix = 'samplePatientMutationMatrix.txt'
    sample_exp_outlier = 'samplePatientOutlierMatrix.txt'
    edges, ggi, outlier_dict, mutation_dict = preprocess(ggi_matirx, sample_mutation_matrix, sample_exp_outlier)

    drivers, actual_events = greedy(edges, mutation_dict)
    print('')
    print(edges)
    drivers, actual_events = greedy(edges, mutation_dict)

    #print (drivers, len(drivers))
    #for key in drivers:
    #    print (key, drivers[key])


if __name__ == "__main__":

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

