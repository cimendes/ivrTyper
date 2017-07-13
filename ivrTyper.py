#!/usr/bin/env python

import os
import sys
import time
import argparse
import modules.utils as utils
import modules.runTyper as runTyper


version = '0.3'


def ivrTyper(args, time):

    number_samples_successfully = 0

    #creating work directory if necessary
    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    #Script path, setting environment variables and reference path
    script_path = os.path.abspath(sys.argv[0])
    utils.setPATHvariable(args.skipProvidedSoftware, script_path)

    reference = os.path.join(os.path.dirname(script_path), 'src/seq/D39_ivr_extended.fasta')

    #load reads file paths
    sample_data = utils.getReadsFiles(workdir)
    samples_total_number = len(sample_data)

    if samples_total_number == 0:
        sys.exit('No samples found.')

    #Run typing algorithm
    for sample, files in sorted(sample_data.items()):
        print(utils.bcolors.HEADER + '\n-> ' + sample + '\n' + utils.bcolors.ENDC)
        success = runTyper.alignSamples(sample, files, reference, args.threads, workdir, script_path, args.keepFiles,
                                    args.minCoverage, args.proportionCutOff, time)
        if success:
            number_samples_successfully+=1

    return number_samples_successfully, samples_total_number


def main():

    parser = argparse.ArgumentParser(prog='ivrTyper.py', description='Reads mapping against pneumococcal ivr locus for typing.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_optional_general = parser.add_argument_group('General options')
    parser_optional_general.add_argument('-w', '--workdir', type=str, metavar='/path/to/workdir/directory/',
                                         help='Path to the directory containing the read files in subdirectories, one '
                                         'per sample (fastq.gz or fq.gz and pair-end direction coded as _R1_001 / '
                                         '_R2_001 or _1 / _2) or to be downloaded',
                                         required=True)
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use',
                                         required=False, default=1)
    parser_optional_general.add_argument('-u', '--skipProvidedSoftware', action='store_true', help='Do not use provided '
                                         'software',
                                         required=False, default=False)
    parser_optional_general.add_argument('-k', '--keepFiles', action='store_true', help='Keep alignment files',
                                         required=False, default=False)

    parser_optional_ivrTyper = parser.add_argument_group('ivrTyper module facultative options')
    parser_optional_ivrTyper.add_argument('-m','--minCoverage', type=int, metavar='N', help='minimum'
                                          'coverage depth a module to be present in the sample',
                                          required=False, default=5)
    parser_optional_ivrTyper.add_argument('-c','--proportionCutOff', type=float, metavar='N', help='Proportion cut off '
                                          'for the module 1.x to be chosen as target', required=False, default=0.8)

    args = parser.parse_args()

    if args.workdir is None:
        parser.error('A directory containing at least one paired-end sample should be provided to --workdir.')

    general_start_time = time.time()
    time_str = time.strftime("%Y%m%d-%H%M%S")

    # Start logger
    sys.stdout = utils.Logger(args.workdir, time_str)

    print(utils.bcolors.HEADER + '\n' + "_,.-'``'-.,_,.='`` irvTyper ``'-.,_,.-'``'-.,_")
    print('\n' '-- ivr locus determination from paired-end genomic data --' + utils.bcolors.ENDC)

    #run ivrTyper
    number_samples_successfully, samples_total_number = ivrTyper(args, time_str)

    print(utils.bcolors.HEADER + '\n' + 'ivrTyper has finished!' + utils.bcolors.ENDC)
    print('\n {} samples out of {} ran successfully'.format(str(number_samples_successfully),str(samples_total_number)))

    end_time = time.time()
    time_taken = end_time - general_start_time
    hours, rest = divmod(time_taken, 3600)
    minutes, seconds = divmod(rest, 60)
    print("\nRuntime: {}h:{}m:{}s\n".format(str(hours), str(minutes), str(round(seconds, 2))))

    if number_samples_successfully == 0:
        sys.exit('No samples ran successfully.')

    sys.exit(0)


if __name__ == "__main__":
    main()
