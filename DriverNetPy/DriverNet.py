#!/usr/bin/env python

'''
A Python version of DriverNet 1.22.0

DriverNet paper:
Bashashati A, Haffari G, Ding J, Ha G, Lui K, Rosner J, Huntsman DG, Caldas C, Aparicio SA, Shah SP. 2012. DriverNet: uncovering the impact of somatic driver mutations on transcriptional networks in cancer. Genome Biology. 13:R124. doi:10.1186/gb-2012-13-12-r124.

'''

import sys
import argparse
from datetime import datetime
import logging

from preprocess import preprocess
from greedy import greedy
from permutation import permutation, compute_p

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def report(drivers, actual_events, pvalues, qvalues, output_prefix):
    '''
    report results
    '''
    rank = sorted(qvalues.items(), key=lambda x:x[1])

    output = open('%s.txt' % output_prefix, 'w')
    output.write('\t'.join(['rank', 'gene', 'mut', 'covered_events', 'p-value', 'q-value']) + '\n')

    for ix, item in enumerate(rank):
        gene = item[0]
        mut = len(actual_events[gene])
        covered_events = drivers[gene]
        p = pvalues[gene]
        q = item[1]
        output.write('\t'.join([str(ix+1), gene, str(mut), str(covered_events), str(p), str(q)]) + '\n')
    output.close()


def DriverNet(args):
    '''
    main function of DriverNet
    '''
    logger.info("Start running.")
    logger.info("Parameters: %s" %(' '.join(sys.argv)))

    ggi_matirx = args.network
    sample_mutation_matrix = args.mutation
    sample_exp_outlier = args.expression
    permutation_time = args.permutation_time
    threads = args.threads
    notPurturbGraph = args.notPurturbGraph
    purturbData = args.purturbData
    output_prefix = args.output

    edges, ggi, outlier_dict, mutation_dict = preprocess(ggi_matirx, sample_mutation_matrix, sample_exp_outlier)

    drivers, actual_events = greedy(edges, mutation_dict)
    logger.info("Driver finding done. Start significance test now.")

    random_res = permutation(ggi, outlier_dict, mutation_dict, permutation_time, threads,
                            notPurturbGraph, purturbData)
    pvalues, qvalues = compute_p(drivers, random_res)
    logger.info("Significance test done.")

    report(drivers, actual_events, pvalues, qvalues, output_prefix)
    logger.info("Writing output done. Please check %s.txt." %(output_prefix))
    logger.info("All done. Cheers.")


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
                                    usage='\n\npython DriverNet.py [-h] <-n network file> <-m mutation file> <-e expression outlier file> [-p permutation_time] [-t threads] [-o output prefix]\n', \
                                    description='', epilog='A Python implementation of DriverNet.\n \
                                                            See https://github.com/YiweiNiu/DriverNetPy for details.')

    parser.add_argument('-n', '--network', required=True, type=str, help='Path to gene-gene interaction network.', metavar='', dest="network")
    parser.add_argument('-m', '--mutation', required=True, type=str, help='Path to sample-mutation file.', metavar='', dest="mutation")
    parser.add_argument('-e', '--expression', required=True, type=str, help='Path to sample-expression file.', metavar='', dest="expression")
    parser.add_argument('-p', '--permutation', default=500, type=int, help='Permutation times.', metavar='', dest="permutation_time")
    parser.add_argument('-pg', '--not-purturbGraph', action='store_true', help='Not permutation graphs.', dest="notPurturbGraph")
    parser.add_argument('-pd', '--purturbData', action='store_true', help='Permutation data.', dest="purturbData")
    parser.add_argument('-t', '--threads', type=int, help='Number of CPUs used when permutation. Using all CPUs by default.', metavar='', dest="threads")
    parser.add_argument('-o', '--output', default='DriverNetPy_res', help='Output prefix.', metavar='', dest="output")

    args = parser.parse_args()

    DriverNet(args)


if __name__ == "__main__":

    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)

