import shlex
import os
import subprocess
from threading import Timer


def getPorpotionsModule2(module2,allele):
    return float(module2[allele]/float(module2["2.1"]+module2["2.2"]+module2["2.3"]))

def getPorpotionsModule1(module1,allele):
    return float(module1[allele]/float(module1["1.1"]+module1["1.2"]))

# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads, write_command=True):
    if os.path.isfile(str(referenceFile + '.1.bt2')):
        run_successfully = True
    else:
        command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
        run_successfully, stdout, stderr = runCommandPopenCommunicate(command, True, None, write_command)
    return run_successfully


# Mapping with Bowtie2
def mappingBowtie2(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc, bowtieOPT,
                   write_command=True, name='alignment'):
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


# Index alignment file
def indexAlignment(alignment_file, write_command=True):
    command = ['samtools', 'index', alignment_file]
    run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, write_command)
    return run_successfully


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads, write_command=True):
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
    print '\n' + 'PATH variable:\n'
    print os.environ['PATH']
    print '\n'

def getReadsFiles(workdir):

    sample_data={}

    for sample in os.listdir(workdir):
        if os.path.isdir(os.path.join(workdir,sample)):
            sample_files = os.listdir(os.path.join(workdir,sample))

            sample_file_foward = None
            sample_file_reverse = None

            for file in sample_files:
                if file.endswith('_R1_001.fastq.gz') or file.endswith('_1.fastq.gz') or file.endswith(
                   '_R1_001.fq.gz') or file.endswith('_1.fq.gz'):
                    sample_file_foward = os.path.join(workdir, sample, file)
                    # print sample_file_foward
                elif file.endswith('_R1_002.fastq.gz') or file.endswith('_2.fastq.gz') or file.endswith(
                    '_R1_002.fq.gz') or file.endswith('_2.fq.gz'):
                    sample_file_reverse = os.path.join(workdir, sample, file)
                    # print sample_file_reverse
            if sample_file_foward == None or sample_file_reverse == None:
                print 'WARNING: No files found for ' + sample
            else:
                sample_data[sample] = [sample_file_foward,sample_file_reverse]

    #for key, value in sample_data.items():
    #    print key
    #    print value

    return sample_data

def runCommandPopenCommunicate(command, shell_True, timeout_sec_None, print_comand_True):
    run_successfully = False
    if not isinstance(command, basestring):
        command = ' '.join(command)
    command = shlex.split(command)

    if print_comand_True:
        print 'Running: ' + ' '.join(command)

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
            print 'Running: ' + str(command)
        if len(stdout) > 0:
            print 'STDOUT'
            print stdout.decode("utf-8")
        if len(stderr) > 0:
            print 'STDERR'
            print stderr.decode("utf-8")
    return run_successfully, stdout, stderr


def kill_subprocess_Popen(subprocess_Popen, command):
    print 'Command run out of time: ' + str(command)
    subprocess_Popen.kill()

