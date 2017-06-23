import argparse
import os
import itertools
import subprocess

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

    # generate a filename for the output file
    directory, file = os.path.split(inName)
    basename = file.strip('.fastq.gz')
    outName = os.path.join( directory, 
        basename + '.polyA_trimmed.fastq' )

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
    print( 'Reads processed: {}'.format(num_processed) )
    print( 'Reads trimmed: {}'.format(num_trimmed) )
    print( 'Percent reads trimmed: {}%'.format( round( num_trimmed/num_processed * 100, 1 ) ) )
    print( 'Reads removed (too short): {}'.format(num_too_short) )
    print( 'Percent reads removed (too short): {}%'.format( round( num_too_short/num_processed * 100, 1 ) ) )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument( '-l', '--length', type = int, required = True, default = 10,
        help = 'minimum length of polyA stretch to define a polyA tail' )
    parser.add_argument( '-i', '--input', required = True,
        help = 'path to the input fastq file' )
    parser.add_argument( '-m', '--min_length', type  = int, required = True, default = 20,
        help = 'minimum length of read after trimming, otherwise read is discarded' )

    args = parser.parse_args()
    remove_polyA_tail( inName = args.input, length = args.length, min_length = args.min_length )


