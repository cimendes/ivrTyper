#!/usr/bin/env python

import os
import sys
sys.path.insert(0,'..')
import time
import argparse
import shutil
import glob
import xml.etree.ElementTree as ET
import modules.getSeqFromENA as getSeq
import modules.download as download
import modules.utils as utils

version = '1.0'

def runIVRTyper(workdir_sample,threads,ivrReport,asperaKey):

    command = ['ivrTyper.py', '--workdir', workdir_sample, '--threads',str(threads), '--reportFile', ivrReport]

    if asperaKey is not None:
        command += ['--asperaKey', asperaKey]

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully

def runReMatCH_ivr(workdir_sample,threads,asperaKey, script_path):

    reference = os.path.join(os.path.dirname(script_path),'src', 'seq', 'SpnD39III_plus_flanking.fasta')
    command = ['rematch.py', '--w', workdir_sample, '-j', str(threads), '-r', reference, '--reportSequenceCoverage',
               '--minGeneCoverage', '30', '--minGeneIdentity', '80']

    if asperaKey is not None:
        command += ['--asperaKey', asperaKey]

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully

def runReMatCH_mlst(workdir_sample,threads,asperaKey):
    command = ['rematch.py', '--w', workdir_sample, '-j', str(threads), '--mlst', '"Streptococcus pneumoniae"',
               '--mlstConsensus', 'noMatter', '--mlstReference']

    if asperaKey is not None:
        command += ['--asperaKey', asperaKey]

    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully

def runPneumoCaT(reads_sample, workdir_sample):

    files=[None,None]

    files[0] = glob.glob(os.path.join(reads_sample + '/' + "*_1.fq.gz"))
    if len(files[0]) != 1:
        files[0] = glob.glob(os.path.join(reads_sample + '/' + "*_1.fastq.gz"))
        if len(files[0]) != 1:
            print "Error finding read file 1"
            return False
    else:
        newname = files[0][0].replace('fq','fastq')
        os.rename(files[0][0], newname)
        files[0][0] = newname

    files[1] = glob.glob(os.path.join(reads_sample + '/' + "*_2.fq.gz"))
    if len(files[1]) != 1:
        files[1] = glob.glob(os.path.join(reads_sample + '/' + "*_2.fastq.gz"))
        if len(files[1]) != 1:
            print "Error finding read file 2"
            return False
    else:
        newname = files[1][0].replace('fq', 'fastq')
        os.rename(files[1][0], newname)
        files[1][0] = newname

    command = ['PneumoCaT.py', '-1', files[0][0], '-2', files[1][0], '-o', workdir_sample, '--cleanup']

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

    script_path = os.path.abspath(sys.argv[0])

    if args.listIDs is not None:
        IDs = getSeq.getListIDs_fromFile(os.path.abspath(args.listIDs.name))
    else:
        parser.error('The path to list containing the IDs to be downloaded (one per line) should be provided.')

    samples_total_number = len(IDs)

    #create report file -ivrTyper
    ivrReport = os.path.join(workdir, "report_" + time_str + ".csv")
    with open(ivrReport, "w") as report:
        report.write("Sample,1.1,1.2,2.1,2.2,2.3,pA,pB,pE,pD,pC,pF,Most Prevalent, gt" + str(0.5) + "\n")

    mlstReport = os.path.join(workdir, "mlst_report_" + time_str + ".csv")
    with open(mlstReport, "w") as mlstReport_file:
        mlstReport_file.write("#sample\tReMatCh_run\tconsensus_type\tST\taroE\tgdh\tgki\trecP\tspi\txpt\tddl\n")

    locusReport = os.path.join(workdir, "locus_report_" + time_str + ".csv")
    with open(locusReport, "w") as locustReport_file:
        locustReport_file.write("#sample\t1_hrcA\t2_SPD_0457\t3_SPD_0456\t4_hsdR\t5_hsdM\t6_hsdS\t7_SPD_0452_creX\t"
                              "8_SPD_0450_hsdS_line\t9_SPD_0451_hsdS_line_line\t10_SPD_0449\t11_glnA\n")

    pneumocatReport = os.path.join(workdir, "serotype_report_" + time_str + ".csv")
    with open(pneumocatReport, "w") as pneumoCaT_file:
        pneumoCaT_file.write("#sample\tSerotype\tQC_coverage\n")

    errorReport = os.path.join(workdir, "status_" + time_str + ".csv")
    with open(errorReport, "w") as error:
        error.write("Sample,Download,ivrTyper,ReMatCh_locus,ReMatCh_mlst,PneumoCaT\n")

    for sample in IDs:

        print "Running {} out of {}".format(IDs.index(sample)+1, samples_total_number)

        workdir_sample = os.path.join(workdir, sample)
        if not os.path.isdir(workdir_sample):
            os.makedirs(workdir_sample)

        download_workdir = os.path.join(workdir_sample, sample,)
        if not os.path.isdir(download_workdir):
            os.makedirs(download_workdir)

        download_successfully, downloaded_files, sequencingInformation = download.runDownload(sample, 'PAIRED',
                                                                                             asperaKey,
                                                                     download_workdir, False, args.threads, 'ALL', None,
                                                                     None)

        if download_successfully:
            status=['PASS','FAIL','FAIL','FAIL','FAIL']

            try:
                success_ivrTyper = runIVRTyper(workdir_sample,args.threads,ivrReport,asperaKey)

                if success_ivrTyper:
                    status[1]='PASS'
                else:
                    print "Error running ivrTyper"
                    status[1] = 'FAIL'

            except:
                print "Error running ivrTyper"
                status[1] = 'FAIL'

            try:
                success_rematch_locus = runReMatCH_ivr(workdir_sample,args.threads,asperaKey, script_path)

                if success_rematch_locus:
                    locusFile = glob.glob(os.path.join(workdir_sample + '/' + "combined_report.data_by_gene.first_run.sequence_coverage.*"))
                    if locusFile >= 1:
                        with open(locusFile[-1],'r') as report1:
                            data = report1.readlines()
                            with open(locusReport, 'a') as report2:
                                report2.write(data[1])
                        status[2]='PASS'
                    else:
                        print "No locus report for {}".format(sample)
                        status[2]='ERROR'

                else:
                    print "Error running ReMatCh for ivr locus"
                    status[2] = 'FAIL'

            except:
                print "Error running ReMatCh for ivr locus"
                status[2] = 'FAIL'

            try:
                success_rematch_mlst =runReMatCH_mlst(workdir_sample,args.threads,asperaKey)

                if success_rematch_mlst:
                    mlstFile=glob.glob(os.path.join(workdir_sample + '/' + "mlst_report*"))
                    if mlstFile >= 1:
                        with open(mlstFile[-1],'r') as report1:
                            data = report1.readlines()
                            with open(mlstReport, 'a') as report2:
                                report2.write(data[1])
                        status[3]='PASS'
                    else:
                        print "No mlst report for {}".format(sample)
                        status[3]='ERROR'
                else:
                    print "Error running ReMatCh for mlst"
                    status[3] = 'FAIL'
            except:
                print "Error running ReMatCh for mlst"
                status[3] = 'FAIL'

            try:
                success_pneumoCaT = runPneumoCaT(download_workdir, workdir_sample)

                if success_pneumoCaT:
                    pneumocatXML = glob.glob(os.path.join(workdir_sample + '/' + "*.results.xml"))

                    if pneumocatXML >= 1:
                        with open(pneumocatXML[-1], 'r') as xmlRepot:
                            Serotype = 'NA'
                            QC_coverage = 'NA'
                            for line in xmlRepot:
                                if '"Serotype"' in line:
                                    Serotype = line.split('value=')[1].replace('>','').strip()
                                elif '"QC_coverage"' in line:
                                    QC_coverage = line.split('value=')[1].replace('/','').replace('>','').strip()

                            with open(pneumocatReport,'a') as pneumoCaT_file:
                                pneumoCaT_file.write(sample + '\t' + Serotype + '\t' + QC_coverage + '\n')

                            status[4] = 'PASS'
                    else:
                        print "No PneumoCaT report for {}".format(sample)
                        status[4]='ERROR'
                else:
                    print "Error running PneumoCaT"
                    status[4] = 'FAIL'
            except:
                print "Error running PneumoCaT"
                status[4] = 'FAIL'

            with open(errorReport, 'a') as statusReport:
                statusReport.write('\t'.join(status)+'\n')

        else:
            print "Error downlading {} read files.".format(sample)
            with open(errorReport, 'a') as error:
                error.write(sample+',FAIL,NA,NA,NA,NA\n')


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
