#!/usr/bin/env python2

# -*- coding: utf-8 -*-

import os
import sys
import time
import argparse
import modules.utils as utils
import modules.runTyper as runTyper
import modules.getSeqFromENA as getSeq
import modules.download as download
import shutil


version = '0.6.8'

def getListIDs(workdir, fileListIDs, taxon_name):
    searched_fastq_files = False
    listIDs = []
    if fileListIDs is None and taxon_name is None:
        listIDs = getSeq.getReadsFiles(workdir)
        searched_fastq_files = True
    elif fileListIDs is not None:
        listIDs = getSeq.getListIDs_fromFile(os.path.abspath(fileListIDs))
        listIDs=dict.fromkeys(listIDs,[])
    elif taxon_name is not None and fileListIDs is None:
        listIDs = getSeq.getTaxonRunIDs(taxon_name, os.path.join(workdir, 'IDs_list.seqFromWebTaxon.tab'))
        listIDs = dict.fromkeys(listIDs, [])

    if len(listIDs) == 0:
        sys.exit('No samples were found')

    return listIDs, searched_fastq_files


def ivrTyper(args, time):

    number_samples_successfully = 0

    # creating work directory if necessary
    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    # Start logger
    sys.stdout = utils.Logger(workdir, time)

    print '\n' + 'COMMAND: ' + ' '.join(sys.argv) + '\n'

    print '\n' + 'PRESENT DIRECTORY: ' + os.path.abspath(os.getcwd()) + '\n'

    print '\n' + 'VERSION: ' + version + '\n'

    #Script path, setting environment variables and reference path
    script_path = os.path.abspath(sys.argv[0])
    utils.setPATHvariable(args.skipProvidedSoftware, script_path)

    reference = os.path.join(os.path.dirname(script_path), 'src/seq/D39_ivr_extended.fasta')

    asperaKey = os.path.abspath(args.asperaKey.name) if args.asperaKey is not None else None

    if args.read1 is not None and args.read2 is not None:
        sample_data = {}
        sample_data[os.path.dirname(os.path.abspath(args.read1)).split('/')[-1]] = [os.path.abspath(args.read1), os.path.abspath(args.read2)]
        searched_fastq_files = True
    else:
        #load reads file paths
        sample_data, searched_fastq_files = getListIDs(workdir, args.listIDs.name if args.listIDs is not None else None,
                                               args.taxon)

    samples_total_number = len(sample_data)

    #Run download if necessary and typing algorithm for each sample
    for sample, files in sorted(sample_data.items()):

        print(utils.bcolors.HEADER + '\n-> ' + sample + '\n' + utils.bcolors.ENDC)

        if not searched_fastq_files:
            workdir_sample = os.path.join(workdir, sample)
            if not os.path.isdir(workdir_sample):
                os.makedirs(workdir_sample)
            # Download Files
            out = download.runDownload(sample, 'PAIRED', asperaKey, workdir_sample, False, args.threads,
                                       args.instrumentPlatform, args.platformModel, args.librarySource)

            if out[0]: #run successfully

                files = out[1]

                success = runTyper.alignSamples(sample, files, reference, args.threads, workdir, script_path, args.keepFiles,
                                    args.minCoverage, args.proportionCutOff, time, args.greaterThan, args.reportFile)

                if success:
                    number_samples_successfully+=1

            else:
                print "Download information doesn't match requirements. Requiring  pair-end samples from {}, {}, " \
                      "and sample is {}, {}, {}.".format(args.instrumentPlatform, args.platformModel,
                                                       out[2]['library_layout'], out[2]['instrument_platform'],
                                                       out[2]['instrument_model'])

            if (not searched_fastq_files) and (not args.keepDownloadFiles) and (not args.keepFiles):
                shutil.rmtree(os.path.join(workdir, sample), ignore_errors=True)
        else:
            success = runTyper.alignSamples(sample, files, reference, args.threads, workdir, script_path,
                                            args.keepFiles,
                                            args.minCoverage, args.proportionCutOff, time, args.greaterThan,
                                            args.reportFile)

            if success:
                number_samples_successfully += 1

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
    parser_optional_general.add_argument('-1', '--read1', metavar='sample_1.fastq.gz',
                                         help='Forward read file (fastq.gz or fq.gz and pair-end direction coded as _R1_001 / '
                                         '_R2_001 or _1 / _2)')
    parser_optional_general.add_argument('-2', '--read2', metavar='sample_2.fastq.gz',
                                         help='Forward read file (fastq.gz or fq.gz and pair-end direction coded as _R1_001 / '
                                         '_R2_001 or _1 / _2)')
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
    parser_optional_ivrTyper.add_argument('-gt', '--greaterThan', type=float, metavar='N', help='proportion cut off '
                                          'for the "greater than" column in the final report.', required=False,
                                          default=0.5)
    parser_optional_ivrTyper.add_argument('-rf', '--reportFile',type=argparse.FileType('r'),
                                          metavar='/path/to/run.*.log', help='Logfile to append run information to.',
                                          required = False)

    parser_optional_download = parser.add_argument_group('Download facultative options')
    parser_optional_download.add_argument('-a', '--asperaKey', type=argparse.FileType('r'),
                                          metavar='/path/to/asperaweb_id_dsa.openssh',
                                          help='Download fastq files from ENA using Aspera Connect. With this option, '
                                          'the path to Private-key file asperaweb_id_dsa.openssh must be provided '
                                          '(normaly found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).',
                                          required=False)
    parser_optional_download.add_argument('-kd', '--keepDownloadFiles', action='store_true',
                                          help='Keep downloaded read files', required=False, default=False)
    parser_optional_download.add_argument('-ip', '--instrumentPlatform', type=str, metavar='ILLUMINA', help='Download '
                                          'files with specific library layout (available options: %(choices)s)',
                                          choices=['ILLUMINA', 'ALL'], required=False, default='ALL')
    parser_optional_download.add_argument('-pm', '--platformModel', type=str, metavar='HiSeq', help='Filter '
                                          'download by the model of the Illumina machine.', choices=['HiSeq', 'MiSeq',
                                          'None'], required=False, default=None)
    parser_optional_download.add_argument('-ls', '--librarySource',type=str, metavar='GENOMIC', help='Filter '
                                          'download by library source (available options: %(choices)s)', choices=['GENOMIC', 'ALL'],
                                          required=False, default=None)

    parser_optional_download_exclusive = parser.add_mutually_exclusive_group()
    parser_optional_download_exclusive.add_argument('-l', '--listIDs', type=argparse.FileType('r'),
                                                    metavar='/path/to/list_IDs.txt',
                                                    help='Path to list containing the IDs to be downloaded (one per line)',
                                                    required=False, default=None)
    parser_optional_download_exclusive.add_argument('-t', '--taxon', type=str, metavar='"Streptococcus pneumoniae"',
                                                    help='Taxon name for which fastq files will be downloaded',
                                                    required=False, default=None)

    args = parser.parse_args()

    if args.workdir is None:
        if args.listIDs is None:
            parser.error('A directory containing at least one paired-end sample should be provided to --workdir.')
        else:
            parser.error('A directory to store the output should be provided to --workdir.')
    if args.platformModel is not None and args.instrumentPlatform is not 'ILLUMINA':
        parser.error("To use the --platformModel option, the --instrumentPlatform must be set as 'ILLUMINA'")
    if (args.read1 is None and args.read2 is not None) or (args.read2 is None and args.read1 is not None):
        parser.error("Both forward and reverse reads must be provided to --read1 and --read2.")

    general_start_time = time.time()
    time_str = time.strftime("%Y%m%d-%H%M%S")

    print(utils.bcolors.HEADER + '\n' + "> irvTyper ")
    print('\n' '\t- ivr locus determination from paired-end genomic data' + utils.bcolors.ENDC)

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
