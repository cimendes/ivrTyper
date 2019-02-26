#!/usr/bin/env python3

import os
import argparse
import subprocess
import time

version = '2.0'


def runTyper(args):
    print('\n' + '> Restarting ivrTyper' + '\n')

    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    initialWorkdir = os.path.abspath(args.initialWorkdir)

    files_required = get_files_required(initialWorkdir)

    samples_run = get_samples_run(files_required['report']['file'])

    command, list_ids, taxon, threads, initial_present_directory = get_command(files_required['run']['file'])

    if list_ids is not None:
        total_samples = getListIDs_fromFile(list_ids)
    elif taxon:
        total_samples = getTaxonRunIDs(files_required['IDs_list.seqFromWebTaxon']['file'])
    else:
        samples_fastq = searchFastqFiles(initialWorkdir)
        total_samples = samples_fastq.keys()

    samples_to_run = list(set(total_samples) - set(samples_run))

    print(str(len(samples_to_run)) + ' samples out of ' + str(
        len(total_samples)) + ' will be analysed by ivrTyper' + '\n')

    if list_ids is not None or taxon:
        samples_to_run_file = write_samples_to_run(samples_to_run, workdir)
    else:
        setSamples_fromFolders(samples_to_run, samples_fastq, workdir)

    command.extend(['-w', workdir])
    command.extend(['-j', str(threads) if args.threads is None else str(args.threads)])
    if list_ids is not None or taxon:
        command.extend(['-l', samples_to_run_file])

    print('ivTyper will start in 5 seconds...')
    time.sleep(5)

    os.chdir(initial_present_directory)
    print(command)
    subprocess.call(command[1:])


def write_samples_to_run(samples_to_run, workdir):
    samples_to_run_file = os.path.join(workdir, 'restart_ivrTyper.samples_to_run.txt')
    with open(samples_to_run_file, 'wt') as writer:
        for sample in samples_to_run:
            writer.write(sample + '\n')
    return samples_to_run_file


def get_files_required(initialWorkdir):
    files_required = {'report': {'extension': 'csv'}, 'run': {'extension': 'log'},
                      'IDs_list.seqFromWebTaxon': {'extension': 'tab'}}

    files = sorted([f for f in os.listdir(initialWorkdir) if
                    not f.startswith('.') and os.path.isfile(os.path.join(initialWorkdir, f))])
    for file_found in files:
        file_path = os.path.join(initialWorkdir, file_found)
        file_modification = os.path.getmtime(file_path)
        for prefix, values in files_required.items():
            if file_found.startswith(prefix) and file_found.endswith('.' + values['extension']):
                if 'file' not in values:
                    files_required[prefix]['file'] = file_path
                    files_required[prefix]['modification'] = file_modification
                else:
                    if file_modification > files_required[prefix]['modification']:
                        files_required[prefix]['file'] = file_path
                        files_required[prefix]['modification'] = file_modification
    return files_required


def get_samples_run(sample_report_file):
    samples_run = []
    with open(sample_report_file, 'rtU') as reader:
        for line in reader:
            sample_info = line.split(',')[0]
            if len(sample_info) > 0: # line is not empty
                if not line.startswith('Sample') and sample_info not in samples_run:
                    samples_run.append(sample_info)

    return samples_run


def get_command(log_file):
    variables = {'command': False, 'directory': False}

    with open(log_file, 'rtU') as reader:
        for line in reader:
            if any([isinstance(value, bool) for value in variables.values()]):
                if len(line) > 0:
                    if 'COMMAND:' in line:
                        variables['command'] = line.split(' ')
                    elif 'PRESENT DIRECTORY:' in line:
                        variables['directory'] = line.split(' ')[-1].strip()
            else:
                break

    command = {'command': [], 'listIDs': None, 'taxon': False, 'threads': None}
    if all([not isinstance(value, bool) for value in variables.values()]):
        counter = 0
        while counter < len(variables['command']):
            if variables['command'][counter].startswith('-'):
                if variables['command'][counter] not in ('-t', '--taxon'):
                    if variables['command'][counter] in ('-l', '--listIDs'):
                        command['listIDs'] = variables['command'][counter + 1]
                        counter += 1
                    elif variables['command'][counter] in ('-w', '--workdir'):
                        counter += 1
                    elif variables['command'][counter] in ('-j', '--threads'):
                        command['threads'] = int(variables['command'][counter + 1])
                        counter += 1
                else:
                    command['taxon'] = True
                    for i in range(counter, len(variables['command'])):
                        if i + 1 < len(variables['command']):
                            if variables['command'][i + 1].startswith('-'):
                                counter = i
                                break
                        else:
                            counter = i
            else:
                command['command'].append(variables['command'][counter])
            counter += 1

    return command['command'], command['listIDs'], command['taxon'], command['threads'], variables['directory']


def getTaxonRunIDs(IDs_list_seqFromWebTaxon_file):
    list_ids = []
    with open(IDs_list_seqFromWebTaxon_file, 'rtU') as reader:
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')
                    list_ids.append(line[0])
    return list_ids


def getListIDs_fromFile(listIDs_file):
    list_ids = []
    with open(listIDs_file, 'rtU') as lines:
        for line in lines:
            line = line.splitlines()[0]
            if len(line) > 0:
                list_ids.append(line)
    return list_ids


def searchFastqFiles(initialWorkdir):

    filesExtensions = ['.fastq.gz', '.fq.gz']
    pairEnd_filesSeparation = [['_R1_001.f', '_R2_001.f'], ['_1.f', '_2.f']]

    list_ids = {}
    directories = [d for d in os.listdir(initialWorkdir) if
                   not d.startswith('.') and os.path.isdir(os.path.join(initialWorkdir, d, ''))]
    for directory_found in directories:
        directory_path = os.path.join(initialWorkdir, directory_found, '')

        fastqFound = []
        files = [f for f in os.listdir(directory_path) if
                 not f.startswith('.') and os.path.isfile(os.path.join(directory_path, f))]
        for file_found in files:
            if file_found.endswith(tuple(filesExtensions)):
                fastqFound.append(file_found)

        if len(fastqFound) == 1:
            list_ids[directory_found] = [os.path.join(directory_path, f) for f in fastqFound]
        elif len(fastqFound) >= 2:
            file_pair = []

            # Search pairs
            for PE_separation in pairEnd_filesSeparation:
                for fastq in fastqFound:
                    if PE_separation[0] in fastq or PE_separation[1] in fastq:
                        file_pair.append(fastq)

                if len(file_pair) == 2:
                    break
                else:
                    file_pair = []

            # Search single
            if len(file_pair) == 0:
                for PE_separation in pairEnd_filesSeparation:
                    for fastq in fastqFound:
                        if PE_separation[0] not in fastq or PE_separation[1] not in fastq:
                            file_pair.append(fastq)

                if len(file_pair) >= 1:
                    file_pair = file_pair[0]

            if len(file_pair) >= 1:
                list_ids[directory_found] = [os.path.join(directory_path, f) for f in file_pair]

    return list_ids


def setSamples_fromFolders(samples_to_run, samples_fastq, workdir):
    for sample in samples_to_run:
        sample_dir = os.path.join(workdir, sample, '')
        if not os.path.isdir(sample_dir):
            os.mkdir(sample_dir)
        for file_found in samples_fastq[sample]:
            link_path = os.path.join(sample_dir, os.path.basename(file_found))
            if os.path.islink(link_path):
                os.remove(link_path)
            if not os.path.isfile(link_path):
                os.symlink(file_found, link_path)


def main():
    parser = argparse.ArgumentParser(prog='restartTyper.py', description='Restart a ivr Typer run abruptly terminated',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-i', '--initialWorkdir', type=str, metavar='/path/to/initial/workdir/directory/',
                                 help='Path to the directory where ivr Typer was running', required=True)

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-w', '--workdir', type=str, metavar='/path/to/workdir/directory/',
                                         help='Path to the directory where irv Typer will run again', required=False,
                                         default='.')
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N',
                                         help='New number of threads to use instead', required=False)

    args = parser.parse_args()

    runTyper(args)


if __name__ == "__main__":
    main()
