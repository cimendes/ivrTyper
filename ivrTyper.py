#!/usr/bin/env python

import os
import sys
import time
import argparse
import modules.utils as utils
import modules.runTyper as typerModule

version = '0.1'


def index_fasta_samtools(fasta, region_None, region_outfile_none, print_comand_true):
    command = ['samtools', 'faidx', fasta, '', '', '']
    shell_true = False
    if region_None is not None:
        command[3] = region_None
    if region_outfile_none is not None:
        command[4] = '>'
        command[5] = region_outfile_none
        shell_true = True
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, shell_true, None, print_comand_true)
    return run_successfully, stdout

# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads):
    outFormat_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
    if sortByName_True:
        command[6] = '-n'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    if not run_successfully:
        output_file = None
    return run_successfully, output_file

# Index alignment file
def indexAlignment(alignment_file):
    command = ['samtools', 'index', alignment_file]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    return run_successfully

def runTyper(args):
    number_samples_successfully = 0
    samples_total_number = 0

    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    print '\n' + 'STARTING' + '\n'

    script_path = os.path.abspath(sys.argv[0])

    utils.setPATHvariable(False, script_path)

    sample_data = utils.getReadsFiles(workdir)
    samples_total_number = len(sample_data)

    reference = os.path.join(os.path.dirname(script_path),'src/seq/D39_ivr_extended.fasta')
    #print reference
    #indexRun = runTyper.indexSequenceBowtie2(reference, args.threads)

    #run bowtie
    for sample, files in sample_data.items():
        print sample
        runMapping, samFile = typerModule.mappingBowtie2(files,reference,args.threads,workdir,True,1,None)

        if runMapping:
            runSortAlignment, bamFile = sortAlignment(samFile, str(os.path.splitext(samFile)[0] + '.bam'), False, args.threads)

            if runSortAlignment:
                os.remove(samFile)
                # Index bam
                run_successfully = indexAlignment(bamFile)

                if run_successfully:

                    runIndexSamtools, stdout = index_fasta_samtools(reference, None, None, True)



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

    args = parser.parse_args()


    if args.workdir is None:
        parser.error('A directory containing at least one paired-end sample should be provided to --workdir.')

    start_time = time.time()

    number_samples_successfully, samples_total_number = runTyper(args)

    print '\n' + 'ivr Typer has finished'
    print '\n' + str(number_samples_successfully) + ' samples out of ' + str(samples_total_number) + ' ran successfully'
    #time_taken = utils.runTime(start_time)
    #del time_taken

    if number_samples_successfully == 0:
        sys.exit('No samples ran successfully!')


if __name__ == "__main__":
    main()
