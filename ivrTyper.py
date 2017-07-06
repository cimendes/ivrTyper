#!/usr/bin/env python

import os
import sys
import time
import argparse
import modules.utils as utils
import modules.runTyper as runTyper

version = '0.1'

def ivrTyper(args):

    number_samples_successfully = 0
    samples_total_number = 0

    #creating work directory if necessary
    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    print '\n' + 'Starting...' + '\n'

    #Script path, setting environment variables and reference path
    script_path = os.path.abspath(sys.argv[0])
    utils.setPATHvariable(args.skipProvidedSoftware, script_path)
    reference = os.path.join(os.path.dirname(script_path), 'src/seq/D39_ivr_extended.fasta')
    # print reference
    # indexRun = runTyper.indexSequenceBowtie2(reference, args.threads)

    #load reads file paths
    sample_data = utils.getReadsFiles(workdir)
    samples_total_number = len(sample_data)

    #run bowtie
    runTyper.alignSamples(sample_data,reference,args.threads, workdir, script_path)

    return number_samples_successfully, samples_total_number


def main():

    parser = argparse.ArgumentParser(prog='ivrTyper.py', description='Reads mapping against pneumococcal ivr locus for typing.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_optional_general = parser.add_argument_group('General options')
    parser_optional_general.add_argument('-w', '--workdir', type=str, metavar='/path/to/workdir/directory/',
                                         help='Path to the directory containing the read files in subdirectories, one per sample (fastq.gz or fq.gz and pair-end direction coded as _R1_001 / _R2_001 or _1 / _2) or to be downloaded',
                                         required=True)
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use',
                                         required=False, default=1)
    parser_optional_general.add_argument('-u', '--skipProvidedSoftware', action='store_true', help='Do not use provided software',
                                         required=False, default=False)

    args = parser.parse_args()


    if args.workdir is None:
        parser.error('A directory containing at least one paired-end sample should be provided to --workdir.')

    start_time = time.time()

    print '\n' + "_,.-'``'-.,_,.='`` irvTyper ``'-.,_,.-'``'-.,_"
    print '\n' '-- ivr locus determination from paired-end genomic data --'

    #run ivrTyper
    number_samples_successfully, samples_total_number = ivrTyper(args)

    print '\n' + 'ivrTyper has finished!'
    print '\n' + str(number_samples_successfully) + ' samples out of ' + str(samples_total_number) + ' ran successfully'

    #TODO - Run time (source: ReMatCh)
    #time_taken = utils.runTime(start_time)
    #del time_taken

    if number_samples_successfully == 0:
        sys.exit('No samples ran successfully.')


if __name__ == "__main__":
    main()
