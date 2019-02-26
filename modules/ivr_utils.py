#!/usr/bin/env python3

import shlex
import os
import sys
import subprocess
from threading import Timer
import shutil
import pickle

version = '1.0'


class bcolors:
    # The entire table of ANSI color codes - https://gist.github.com/chrisopedia/8754917
    BOLD = '\033[1m'
    HEADER = '\033[1;96m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


def getProportionsModule2(module2, allele):
    return float(module2[allele]/float(module2["2.1"]+module2["2.2"]+module2["2.3"]))


def getProportionsModule1(module1, allele):
    return float(module1[allele]/float(module1["1.1"]+module1["1.2"]))


def indexSequenceBowtie2(referenceFile, threads, write_command=True):
    # Indexing reference file using Bowtie2
    if os.path.isfile(str(referenceFile + '.1.bt2')):
        run_successfully = True
    else:
        command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
        run_successfully, stdout, stderr = runCommandPopenCommunicate(command, True, None, write_command)
    return run_successfully


def mappingBowtie2(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc, bowtieOPT,
                   write_command=True, name='alignment'):
    # Mapping with Bowtie2

    sam_file = os.path.join(outdir, str(name+'.sam'))

    # Index reference file
    run_successfully = indexSequenceBowtie2(referenceFile, threads)

    if run_successfully:
        command = ['bowtie2', '-k', str(numMapLoc), '-q', '', '--threads', str(threads), '-x', referenceFile, '',
                   '--no-unal', '', '-S', sam_file]

        if len(fastq_files) == 1:
            command[9] = '-U ' + fastq_files[0]
        elif len(fastq_files) == 2:
            command[9] = '-1 ' + fastq_files[0] + ' -2 ' + fastq_files[1]
        else:
            return False, None

        if conserved_True:
            command[4] = '--sensitive'
        else:
            command[4] = '--very-sensitive-local'

        if bowtieOPT is not None:
            command[11] = bowtieOPT

        run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, write_command)

    if not run_successfully:
        sam_file = None

    return run_successfully, sam_file


def index_fasta_samtools(fasta, region_None, region_outfile_none, print_comand_true):
    command = ['samtools', 'faidx', fasta, '', '', '']
    shell_true = False
    if region_None is not None:
        command[3] = region_None
    if region_outfile_none is not None:
        command[4] = '>'
        command[5] = region_outfile_none
        shell_true = True
    run_successfully, stdout, stderr = runCommandPopenCommunicate(command, shell_true, None, print_comand_true)
    return run_successfully, stdout


def indexAlignment(alignment_file, write_command=True):
    # Index alignment file
    command = ['samtools', 'index', alignment_file]
    run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, write_command)
    return run_successfully


def sortAlignment(alignment_file, output_file, sortByName_True, threads, write_command=True):
    # Sort alignment file
    outFormat_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
    if sortByName_True:
        command[6] = '-n'
    run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, write_command)

    if not run_successfully:
        output_file = None
    return run_successfully, output_file


def bam2fastq(bam_file, write_command=True):
    command = ['samtools', 'fastq', bam_file, '>', bam_file+'.fastq']
    run_successfully, stdout, stderr = runCommandPopenCommunicate(command, True, None, write_command)
    return run_successfully


def setPATHvariable(doNotUseProvidedSoftware, script_path):
    path_variable = os.environ['PATH']
    script_folder = os.path.dirname(script_path)
    # Set path to use provided software
    if not doNotUseProvidedSoftware:
        bowtie2 = os.path.join(script_folder, 'src', 'bowtie2-2.2.9')
        samtools = os.path.join(script_folder, 'src', 'samtools-1.3.1', 'bin')
        bcftools = os.path.join(script_folder, 'src', 'bcftools-1.3.1', 'bin')

        os.environ['PATH'] = str(':'.join([bowtie2, samtools, bcftools, path_variable]))

    # Print PATH variable
    print(bcolors.BOLD + '\n' + 'PATH variable:' + bcolors.ENDC)
    print(os.environ['PATH'])
    print('\n')


def runCommandPopenCommunicate(command, shell_True, timeout_sec_None, print_comand_True):
    run_successfully = False
    if not isinstance(command, str):
        command = ' '.join(command)
    command = shlex.split(command)

    if print_comand_True:
        print('Running: ' + ' '.join(command))

    if shell_True:
        command = ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    not_killed_by_timer = True
    if timeout_sec_None is None:
        stdout, stderr = proc.communicate()
    else:
        timer = Timer(timeout_sec_None, kill_subprocess_Popen, args=(proc, command,))
        timer.start()
        stdout, stderr = proc.communicate()
        timer.cancel()
        not_killed_by_timer = timer.isAlive()

    if proc.returncode == 0:
        run_successfully = True
    else:
        if not print_comand_True and not_killed_by_timer:
            print('Running: ' + str(command))
        if len(stdout) > 0:
            print('STDOUT')
            print(stdout.decode("utf-8"))
        if len(stderr) > 0:
            print('STDERR')
            print(stderr.decode("utf-8"))
    return run_successfully, stdout, stderr


def kill_subprocess_Popen(subprocess_Popen, command):
    print('Command run out of time: ' + str(command))
    subprocess_Popen.kill()


class Logger(object):
    def __init__(self, out_directory, time_str):

        self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
        self.terminal = sys.stdout
        self.log = open(self.logfile, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        pass

    def getLogFile(self):
        return self.logfile


def rchop(string, ending):
    if string.endswith(ending):
        string = string[:-len(ending)]
    return string


def removeDirectory(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)


def saveVariableToPickle(variableToStore, outdir, prefix):
    pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
    with open(pickleFile, 'wb') as writer:
        pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
    with open(pickleFile, 'rb') as reader:
        variable = pickle.load(reader)
    return variable


def RepresentsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False