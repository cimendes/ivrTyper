#!/usr/bin/env python3

# Adapted from https://github.com/B-UMMI/ReMatCh/blob/master/rematch.py

from urllib.request import urlopen
import sys
import urllib.parse
import xml.etree.ElementTree as ET
import time
import modules.ivr_utils as utils
import os


version = '1.0'


def getReadsFiles(workdir):

    sample_data={}

    for sample in os.listdir(workdir):
        if os.path.isdir(os.path.join(workdir,sample)):
            sample_files = os.listdir(os.path.join(workdir,sample))

            sample_file_foward = None
            sample_file_reverse = None

            for file in sample_files:
                if file.endswith('_R1_001.fastq.gz') or file.endswith('_1.fastq.gz') or \
                        file.endswith('R1_001.fq.gz') or file.endswith('_1.fq.gz'):
                    sample_file_foward = os.path.join(workdir, sample, file)

                elif file.endswith('_R1_002.fastq.gz') or file.endswith('_2.fastq.gz') or \
                        file.endswith('_R1_002.fq.gz') or file.endswith('_2.fq.gz'):
                    sample_file_reverse = os.path.join(workdir, sample, file)

            if sample_file_foward is None or sample_file_reverse is None:
                print(utils.bcolors.WARNING + 'WARNING: No files found for ' + sample + utils.bcolors.ENDC)
            else:
                sample_data[sample] = [sample_file_foward, sample_file_reverse]

    return sample_data


def getListIDs_fromFile(fileListIDs):
    list_ids = []

    with open(fileListIDs, 'rtU') as lines:
        for line in lines:
            line = line.splitlines()[0]
            if len(line) > 0:
                list_ids.append(line)

    if len(list_ids) == 0:
        sys.exit('No runIDs were found in ' + fileListIDs)

    return list_ids


def getTaxonRunIDs(taxon_name, outputfile):
    runSeqFromWebTaxon(taxon_name, outputfile, True, True, True, False)

    runIDs = []
    with open(outputfile, 'rtU') as reader:
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')
                    runIDs.append(line[0])

    return runIDs


def runSeqFromWebTaxon(taxonname, outputfile, getmachine, getOmicsDataType, getLibraryType, print_True):
    print('\n' + 'Searching RunIDs for ' + taxonname)

    taxonname = urllib.parse.quote(taxonname)
    url = "http://www.ebi.ac.uk/ena/data/view/Taxon%3A" + taxonname + "&display=xml"
    try:
        content = urlopen(url)
        xml = content.read()
        tree = ET.fromstring(xml)
        taxonid = ''
    except:
        print("Ooops!There might be a problem with the ena service, try later or check if the xml is well formated at " + url)
        raise
    for child in tree:
        taxonid = child.get('taxId')
    if (taxonid):
        print("\n" + "Taxon ID found: " + taxonid)
        url = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree%28" + taxonid + "%29%22&result=read_run&display=xml"

        content = urlopen(url)
        xml = content.read()
        tree = ET.fromstring(xml)

        runid = ''
        n = 0
        with open(outputfile, "w") as f:
            f.write('#' + str(time.strftime("%d/%m/%Y")) + "\n")
            model = ''
            prjid = ''
            length_line = 0
            omics = ''
            libraryType = ''
            for child in tree:
                runid = child.get('accession')

                n += 1

                if getmachine is True or getOmicsDataType is True or getLibraryType is True:
                    for child2 in child:
                        if child2.tag == 'EXPERIMENT_REF':
                            expid = child2.get('accession')
                            url2 = "http://www.ebi.ac.uk/ena/data/view/" + expid + "&display=xml"
                            content = urlopen(url2)
                            xml = content.read()
                            tree2 = ET.fromstring(xml)
                            try:
                                for child3 in tree2:
                                    for child4 in child3:
                                        if child4.tag == 'PLATFORM':
                                            for child5 in child4:
                                                for child6 in child5:
                                                    if child6.tag == 'INSTRUMENT_MODEL':
                                                        model = child6.text
                                        elif child4.tag == 'STUDY_REF':
                                            prjid = child4.get('accession')
                                        elif child4.tag == 'DESIGN':
                                            if getOmicsDataType is True or getLibraryType is True:
                                                for child5 in child4:
                                                    if child5.tag == 'LIBRARY_DESCRIPTOR':
                                                        for child6 in child5:
                                                            if child6.tag == 'LIBRARY_SOURCE' and \
                                                                    getOmicsDataType is True:
                                                                omics = child6.text
                                                            elif child6.tag == 'LIBRARY_LAYOUT' and \
                                                                    getLibraryType is True:
                                                                libraryType = child6[0].tag
                            except:
                                model = 'not found'
                                omics = 'not found'
                                libraryType = 'not found'
                    f.write(str(runid) + "\t" + model + "\t" + prjid + "\t" + omics + "\t" + libraryType + "\n")
                    if print_True:
                        line = "run acession %s sequenced on %s from project %s for %s %s end data" % \
                               (runid, model, prjid, omics, libraryType)

                        if length_line < len(line):
                            length_line = len(line)
                        sys.stderr.write("\r" + line + str(' ' * (length_line - len(line))))
                        sys.stderr.flush()
                else:
                    f.write(str(runid) + '\t' * 4 + "\n")
                    if print_True:
                        line = "run acession %s %s" % (runid, prjid)
                        if length_line < len(line):
                            length_line = len(line)
                        sys.stderr.write("\r" + line + str(' ' * (length_line - len(line))))
                        sys.stderr.flush()
        print("\n")
        print("\nfound %s run id's" % n)

    else:
        print("taxon name does not exist")