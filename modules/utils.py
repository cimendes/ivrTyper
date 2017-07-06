import shlex
import os
import subprocess
from threading import Timer


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
