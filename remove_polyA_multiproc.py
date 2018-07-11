#!/usr/bin/env python

"""
Simple script to trim polyA stretches of a specified length from fastq.gz files
"""

import argparse
import os
import itertools
import subprocess
import multiprocessing as mp
import json

__author__ = 'Rob Moccia'
__version__ = '0.1'
__maintainer__ = 'Rob Moccia'
__email__ = 'robert.moccia@pfizer.com'
__status__ = "Development"


def remove_polyA_tail( inName, length, min_length ):
    '''
    Rewrite fastq file with polyA tail removed.
    Arguments:  input -- str, path to the fastq.gz file
                length -- int, minimum length to define a polyA tail
                min_length -- int, minimum number of reads remaining after trimming
                                otherwise read is discarded
    Output: a new fastq.gz file with the trimmed reads
    '''

    # declare accumulator for number of records processed
    num_processed = 0
    # declare accumulator for number of records deleted because too short
    num_too_short = 0
    # declare accumulator for number of reads trimmed
    num_trimmed = 0
    # build a string representing the minimum polyA tail
    polyA = 'A' * length

    # generate a filename for the output file and a log file
    directory, fname = os.path.split(inName)
    if fname.endswith('.fastq.gz'):
        basename = fname[:-len('.fastq.gz')]
    outName = os.path.join( directory, 
        basename + '.polyA_trimmed.fastq' )
    logName = os.path.join( directory, basename + '.polyA_trim_log.json' )

    with open( outName, 'w' ) as outFile:
        with os.popen( 'zcat ' + inName, 'r' ) as inFile:
            for record in itertools.zip_longest( *[ inFile ] * 4 ):
                num_processed += 1
                read = record[1]
                quality = record[3]
                try:
                    polyA_start = read.index(polyA)
                    if polyA_start <= min_length:
                        num_too_short += 1
                        next
                    else:
                        # trim polyA tail from read
                        read = read[:polyA_start] + '\n'
                        # trim the quality scores corresponding to the trimmed bases
                        quality = quality[:polyA_start] + '\n'
                        record = ( record[0], read, record[2], quality )
                        outFile.write( ''.join(record) )
                        num_trimmed += 1
                except ValueError:
                    # this will catch cases where there is no matching polyA tail in read
                    outFile.write( ''.join(record) )

    subprocess.run( [ 'gzip', '-f', outName ] )

    with open( logName, 'w' ) as log:
        json.dump( { 'file': fname,
                     'settings': { 'length': length, 'min_length': min_length },
                     'records_processed': num_processed,
                     'records_trimmed': num_trimmed,
                     'percent_trimmed': round( num_trimmed/num_processed * 100, 1 ),
                     'records_removed': num_too_short,
                     'percent_removed': round( num_too_short/num_processed * 100, 1 ) 
                    },
                    fp = log )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument( '-l', '--length', type = int, required = True, default = 10,
        help = 'minimum length of polyA stretch to define a polyA tail' )
    parser.add_argument( 'directory',
        help = 'directory containing fastq.gz files to process' )
    parser.add_argument( '-m', '--min_length', type  = int, required = True, default = 20,
        help = 'minimum length of read after trimming, otherwise read is discarded' )

    args = parser.parse_args()

    fastq_files = [ file for file in os.listdir(args.directory) if file.endswith('.fastq.gz') and 'polyA_trimmed' not in file ]
    jobs = []
    for fq in fastq_files:
        p = mp.Process( target = remove_polyA_tail, args = (fq, ),
                        kwargs = { 'length': args.length,
                                    'min_length': args.min_length } )
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()


