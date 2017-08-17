#!/usr/bin/env python

import os
import sys
sys.path.append("..")
import time
import argparse
import shutil
import glob
import modules.getSeqFromENA as getSeq
import modules.download as download
import modules.utils as utils

version = '0.1'

def runIVRTyper(workdir_sample,threads,ivrReport,asperaKey):

    command = ['ivrTyper.py', '--workdir', workdir_sample, '--threads',str(threads), '--reportFile', ivrReport]

    if asperaKey is not None:
        command += ['--asperaKey', asperaKey]

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully

def runReMatCH_ivr():
    command = ['ivrTyper.py', '--workdir', workdir_sample, '--threads', str(threads), '--reportFile', ivrReport]

    if asperaKey is not None:
        command += ['--asperaKey', asperaKey]

    print command

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully

def runReMatCH_mlst(workdir_sample,threads,asperaKey):
    command = ['rematch.py', '--w', workdir_sample, '-j', str(threads), '--mlst', '"Streptococcus pneumoniae"',
               '--mlstConsensus', 'noMatter', '--mlstReference']

    if asperaKey is not None:
        command += ['--asperaKey', asperaKey]

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully

def main():

    parser = argparse.ArgumentParser(prog='runTyperReMatCh.py', description='Reads mapping against pneumococcal ivr '
                                     'locus for typing and locus evaluation - ReMatCh and ivrTyper integration. ',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_optional_general = parser.add_argument_group('General options')
    parser_optional_general.add_argument('-w', '--workdir', type=str, metavar='/path/to/workdir/directory/',
                                         help='Path to the output directory',
                                         required=True)
    parser_optional_general.add_argument('-l','--listIDs', type=argparse.FileType('r'),
                                                    metavar='/path/to/list_IDs.txt',
                                                    help='Path to list containing the IDs to be downloaded (one per line)',
                                                    required=False, default=None)
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use',
                                         required=False, default=1)

    parser_optional_download = parser.add_argument_group('Download facultative options')
    parser_optional_download.add_argument('-a', '--asperaKey', type=argparse.FileType('r'),
                                          metavar='/path/to/asperaweb_id_dsa.openssh',
                                          help='Download fastq files from ENA using Aspera Connect. With this option, '
                                          'the path to Private-key file asperaweb_id_dsa.openssh must be provided '
                                          '(normaly found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).',
                                          required=False)
    parser_optional_download.add_argument('-kd', '--keepDownloadFiles', action='store_true',
                                          help='Keep downloaded read files', required=False, default=False)


    args = parser.parse_args()

    if args.workdir is None:
        parser.error('An output directory should be provided to --workdir.')

    general_start_time = time.time()
    time_str = time.strftime("%Y%m%d-%H%M%S")

    workdir = os.path.abspath(args.workdir)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    asperaKey = os.path.abspath(args.asperaKey.name) if args.asperaKey is not None else None

    if args.listIDs is not None:
        IDs = getSeq.getListIDs_fromFile(os.path.abspath(args.listIDs.name))
    else:
        parser.error('The path to list containing the IDs to be downloaded (one per line) should be provided.')

    samples_total_number = len(IDs)

    #create report file -ivrTyper
    ivrReport = os.path.join(workdir, "report_" + time_str + ".csv")
    with open(ivrReport, "w") as report:
        report.write("Sample,1.1,1.2,2.1,2.2,2.3,pA,pB,pE,pD,pC,pF,Most Prevalent, gt" + str(0.5) + "\n")

    for sample in IDs:

        workdir_sample = os.path.join(workdir, sample)
        if not os.path.isdir(workdir_sample):
            os.makedirs(workdir_sample)

        download_workdir = os.path.join(workdir_sample, sample,)
        if not os.path.isdir(download_workdir):
            os.makedirs(download_workdir)

        run_successfully, downloaded_files, sequencingInformation = download.runDownload(sample, 'PAIRED', asperaKey,
                                                                     download_workdir, False, args.threads, 'ALL', None,
                                                                     None)

        #success_ivrTyper = runIVRTyper(workdir_sample,args.threads,ivrReport,asperaKey)

        #TODO
        #runReMatCH_ivr()

        success_rematch_mlst =runReMatCH_mlst(workdir_sample,args.threads,asperaKey)

        if success_rematch_mlst:
            print os.path.join(workdir_sample + '/' + "mlst_report*")
            lala=glob.glob(os.path.join(workdir_sample + '/' + "mlst_report*"))
            print lala

        #TODO
        #runPneumoCaT()



        if not args.keepDownloadFiles:
            shutil.rmtree(download_workdir, ignore_errors=True)


    end_time = time.time()
    time_taken = end_time - general_start_time
    hours, rest = divmod(time_taken, 3600)
    minutes, seconds = divmod(rest, 60)
    print("\nRuntime: {}h:{}m:{}s\n".format(str(hours), str(minutes), str(round(seconds, 2))))

    sys.exit(0)


if __name__ == "__main__":
    main()
