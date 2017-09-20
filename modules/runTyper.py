import os
import os.path
import utils
import pysam
import shutil


def writeReport(workdir, sample, time, module1_reads, module2_reads, proportions, threshold, minCoverage,
                proportionCutOff, reportFile):

    fname = os.path.join(workdir, "report_" + time + ".csv")

    if reportFile is not None:
        fname = reportFile.name

    alleles=['A','B','E','D','C','F']

    maxP=0
    gt_t=0
    for p in proportions:
        if utils.RepresentsFloat(p):
            if float(p)>maxP:
                maxP = float(p)
    if maxP >  threshold:
        gt_t = maxP

    if maxP == 0:
        classifier = 'NT'
        if all(int(i) < minCoverage for i in module1_reads): # and any(i > minCoverage for i in module1_reads):
            classifier = 'NT - module 1 absent'
        else:
            if 'NA' not in module2_reads and all(int(i) < minCoverage for i in module2_reads):
                classifier = 'NT - module 2 absent'
            elif all(int(i) > minCoverage for i in module1_reads) and all(i == 'NA' for i in module2_reads):
                classifier = 'NT - inconclusive module 1 '

    elif proportions.count(str(maxP)) > 1:
        mix_alleles = []
        #print proportions.count(str(maxP))
        start = 0
        for i in range(proportions.count(str(maxP))):
            mix_alleles.append(alleles[proportions.index(str(maxP),start)])
            start = proportions.index(str(maxP),start) + 1

        classifier =' '.join(mix_alleles)
    else:
        classifier = alleles[proportions.index(str(maxP))]

    gt_allele='NT'
    if gt_t != 0:
        alleles_gt = []
        start = 0
        for i in range(proportions.count(str(gt_t))):
            alleles_gt.append(alleles[proportions.index(str(gt_t),start)])
            start = proportions.index(str(gt_t),start) + 1
        gt_allele = ' '.join(alleles_gt)


    toWrite = [sample] + module1_reads + module2_reads + proportions + [classifier] + [gt_allele]

    if os.path.isfile(fname):
        with open(fname, "a") as report:
            report.write(','.join(toWrite)+'\n')
    else:
        with open(fname, "w") as report:
            report.write("Sample,1.1,1.2,2.1,2.2,2.3,pA,pB,pE,pD,pC,pF,Most Prevalent, gt" + str(threshold) + "\n")

            report.write(','.join(toWrite)+'\n')


def getType(workdir,sample, module1, module2, minCoverage):

    toWrite=['NA']*6

    print(utils.bcolors.BOLD + "\n---> Typing:\n" + utils.bcolors.ENDC)

    if module1.keys()[0] == '1.1':
        if module2['2.1'] >= minCoverage:
            pA = utils.getProportionsModule2(module2, '2.1')
            print(utils.bcolors.OKGREEN + "Allele A (1.1 2.1): {} \n".format(format(pA, '.2f')) + utils.bcolors.ENDC)
            toWrite[0]=str(pA)
        else:
            print(utils.bcolors.WARNING + "Coverage on the 2.1 module below the minCoverage {} read "
                                          "threshold.\n".format(minCoverage) + utils.bcolors.ENDC)
            toWrite[0] = '0'

        if module2['2.2'] >= minCoverage:
            pB = utils.getProportionsModule2(module2, '2.2')
            print(utils.bcolors.OKGREEN +  "Allele B (1.1 2.2): {} \n".format(format(pB, '.2f')) + utils.bcolors.ENDC)
            toWrite[1] = str(pB)
        else:
            print(utils.bcolors.WARNING +  "Coverage on the 2.2 module below the minCoverage {} read "
                                           "threshold.\n".format(minCoverage) + utils.bcolors.ENDC)
            toWrite[1] = '0'

        if module2['2.3'] >= minCoverage:
            pE = utils.getProportionsModule2(module2, '2.3')
            print(utils.bcolors.OKGREEN +  "Allele E (1.1 2.3): {} \n".format(format(pE, '.2f')) + utils.bcolors.ENDC)
            toWrite[2] = str(pE)
        else:
            print(utils.bcolors.WARNING +  "Coverage on the 2.3 module below the minCoverage {} read "
                                           "threshold.\n".format(minCoverage) + utils.bcolors.ENDC)
            toWrite[2] = '0'

    elif module1.keys()[0] == '1.2':
        if module2['2.1'] >= minCoverage:
            pD = utils.getProportionsModule2(module2, '2.1')
            print(utils.bcolors.OKGREEN + "Allele D (1.2 2.1): {} \n".format(format(pD, '.2f')) + utils.bcolors.ENDC)
            toWrite[3] = str(pD)
        else:
            print(utils.bcolors.WARNING + "Coverage on the 2.1 module below the minCoverage {} read "
                                          "threshold.\n".format(minCoverage) + utils.bcolors.ENDC)
            toWrite[3] = '0'

        if module2['2.2'] >= minCoverage:
            pC = utils.getProportionsModule2(module2, '2.2')
            print(utils.bcolors.OKGREEN +  "Allele C (1.2 2.2): {} \n".format(format(pC, '.2f')) + utils.bcolors.ENDC)
            toWrite[4] = str(pC)
        else:
            print(utils.bcolors.WARNING + "Coverage on the 2.2 module below the minCoverage {} read "
                                          "threshold.\n".format(minCoverage) + utils.bcolors.ENDC)
            toWrite[4] = '0'

        if module2['2.3'] >= minCoverage:
            pF = utils.getProportionsModule2(module2, '2.3')
            print(utils.bcolors.OKGREEN +  "Allele F (1.2 2.3): {} \n".format(format(pF, '.2f')) + utils.bcolors.ENDC)
            toWrite[5] = str(pF)
        else:
            print(utils.bcolors.WARNING + "Coverage on the 2.3 module below the minCoverage {} read "
                                          "threshold.\n".format(minCoverage) + utils.bcolors.ENDC)
            toWrite[5] = '0'
    else:
        print "something went wrong with he module1 dictionary keys"
        return False, None

    return True, toWrite


def typeSeq_moduleTwo(files, threads, workdir, script_path):
    readCount_2_1 = None
    readCount_2_2 = None
    readCount_2_3 = None

    # 2.1
    runMapping, samFile1 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '2.1.fasta'), threads, workdir, False, 1, '--seed 42', True, "2.1")
    if runMapping:
        runSortAlignment, bamFile1 = utils.sortAlignment(samFile1, str(os.path.splitext(samFile1)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile1, False)
            if runIndex:
                samfile1 = pysam.AlignmentFile(bamFile1, "rb")
                readCount_2_1 = samfile1.mapped
                print(utils.bcolors.BOLD + "\t -> Module 2.1 - {} reads".format(readCount_2_1) + utils.bcolors.ENDC)
    else:
        print(utils.bcolors.FAIL + 'Failed 2.1 Bowtie mapping' + utils.bcolors.ENDC)
        return False, None, None

    # 2.2
    runMapping, samFile2 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '2.2.fasta'), threads, workdir, False, 1, '--seed 42', True, "2.2")
    if runMapping:
        runSortAlignment, bamFile2 = utils.sortAlignment(samFile2, str(os.path.splitext(samFile2)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile2, False)
            if runIndex:
                samfile2 = pysam.AlignmentFile(bamFile2, "rb")
                readCount_2_2 = samfile2.mapped
                print(utils.bcolors.BOLD + "\t -> Module 2.2 - {} reads".format(readCount_2_2) + utils.bcolors.ENDC)
    else:
        print(utils.bcolors.FAIL + 'Failed 2.2 Bowtie mapping' + utils.bcolors.ENDC)
        return False, None, None

    # 2.3
    runMapping, samFile3 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '2.3.fasta'), threads, workdir, False, 1, '--seed 42', True, "2.3")
    if runMapping:
        runSortAlignment, bamFile3 = utils.sortAlignment(samFile3, str(os.path.splitext(samFile3)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile3, False)
            if runIndex:
                samfile3 = pysam.AlignmentFile(bamFile3, "rb")
                readCount_2_3 = samfile3.mapped
                print(utils.bcolors.BOLD + "\t -> Module 2.3 - {} reads".format(readCount_2_3) + utils.bcolors.ENDC)
    else:
        print(utils.bcolors.FAIL +  'Failed 2.3 Bowtie mapping' + utils.bcolors.ENDC)
        return False, None, None

    module2 = {"2.1": int(readCount_2_1), "2.2": int(readCount_2_2), "2.3": int(readCount_2_3)}
    toWrite = [str(readCount_2_1),str(readCount_2_2),str(readCount_2_3)]
    return True, module2, toWrite


def getSeq_moduleTwo(first_type, bamfile, threads, workdir, script_path, sample):

    pysam_fullRef = pysam.AlignmentFile(bamfile, "rb")

    if first_type == '1.1':
        iter = pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus', 8029 - 1, 8446 - 1)
        matePairs_conserved = pysam.AlignmentFile(os.path.join(workdir, sample + "_matepairs_conserved_module2.bam"),
                                                               "wb", template=pysam_fullRef)

        #TODO - remove verification
        for x in iter:

            try:
                # read needs to be paired and mapped to the reverse strand of the reference
                if x.is_paired and x.is_reverse:
                    mate = pysam_fullRef.mate(x)
                    matePairs_conserved.write(mate)
            except:
                pass

        matePairs_conserved.close()

        # TODO - this looks ugly..
        runSortAlignment, bamFile_1_1 = utils.sortAlignment(os.path.join(workdir,
                                                            sample + "_matepairs_conserved_module2.bam"),
                                                            os.path.join(workdir,
                                                            sample + "_matepairs_conserved_module2.bam"),
                                                            False, threads, False)
        utils.indexAlignment(bamFile_1_1, False)
        utils.bam2fastq(bamFile_1_1, False)

        run, module2, report2 =typeSeq_moduleTwo([bamFile_1_1 + '.fastq'], threads, workdir, script_path)

        return run, module2, report2

    elif first_type == '1.2':
        iter = pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus', 4462 - 1, 4861 - 1)
        matePairs_conserved = pysam.AlignmentFile(os.path.join(workdir, sample + "_matepairs_conserved_module2.bam"),
                                                               "wb", template=pysam_fullRef)

        #TODO - remove verification...
        for x in iter:
            #TODO
            try:
                # read needs to be mapped to the forward strand of the reference and both mates need to be pair
                if x.is_paired and not x.is_reverse:
                    mate = pysam_fullRef.mate(x)
                    matePairs_conserved.write(mate)
            except:
                pass
        matePairs_conserved.close()

        # TODO - this looks ugly..
        runSortAlignment, bamFile_1_2 = utils.sortAlignment(os.path.join(workdir,
                                                                         sample + "_matepairs_conserved_module2.bam"),
                                                            os.path.join(workdir,
                                                                         sample + "_matepairs_conserved_module2.bam"),
                                                            False, threads, False)
        utils.indexAlignment(bamFile_1_2, False)
        utils.bam2fastq(bamFile_1_2, False)

        run, module2, report2 = typeSeq_moduleTwo([bamFile_1_2 + '.fastq'], threads, workdir, script_path)

        return run, module2, report2

    else:
        print(utils.bcolors.FAIL + "No module 2.x found. This sample is untypable." + utils.bcolors.ENDC)
        return False, None


def typeSeq_moduleOne(files, threads, workdir, script_path, minCoverage, proportionCutOff):

    readCount_1_1=None
    readCount_1_2=None

    #1.1
    runMapping, samFile1 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '1.1.fasta'), threads, workdir, False, 1, '--seed 42', True, "1.1")
    if runMapping:
        runSortAlignment, bamFile1 = utils.sortAlignment(samFile1, str(os.path.splitext(samFile1)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile1, False)
            if runIndex:
                samfile1 = pysam.AlignmentFile(bamFile1, "rb")
                readCount_1_1 = samfile1.mapped
                print(utils.bcolors.BOLD + "\t -> Module 1.1 - {} reads".format(readCount_1_1) + utils.bcolors.ENDC)
    else:
        print(utils.bcolors.FAIL + 'Failed 1.1 Bowtie mapping' + utils.bcolors.ENDC)
        return False, readCount_1_2, None

    #1.2
    runMapping, samFile2 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '1.2.fasta'), threads, workdir, True, 1, '--seed 42', True, "1.2")
    if runMapping:
        runSortAlignment, bamFile2 = utils.sortAlignment(samFile2, str(os.path.splitext(samFile2)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile2, False)
            if runIndex:
                samfile2 = pysam.AlignmentFile(bamFile2, "rb")
                readCount_1_2 = samfile2.mapped
                print(utils.bcolors.BOLD + "\t-> Module 1.2 - {} reads".format(readCount_1_2) + utils.bcolors.ENDC)
    else:
        print(utils.bcolors.FAIL + "Failed 1.2 Bowtie mapping" + utils.bcolors.ENDC)
        return False, readCount_1_2, ['Fail']*2

    if readCount_1_1 < minCoverage and readCount_1_2 < minCoverage:
        print(utils.bcolors.FAIL + "\n\t- Coverage on the 1.1 and 1.2 modules below the minCoverage {} read threshold. "
                                   "This sample cannot be typed.".format(minCoverage) + utils.bcolors.ENDC)
        return False, None, [str(readCount_1_1), str(readCount_1_2)]

    else:
        if readCount_1_1 < minCoverage:
            print(utils.bcolors.OKGREEN + "\n\t- 1.2 module found with {} reads. Using it as new target.".format(
                readCount_1_2) + utils.bcolors.ENDC)
            module1 = {'1.2': int(readCount_1_2)}
            toWrite=[str(readCount_1_1),str(readCount_1_2)]
            return True, module1, toWrite

        elif readCount_1_2 < minCoverage:
            print(utils.bcolors.OKGREEN +"\n\t- 1.1 module found with {} reads. Using it as new target.".format(
                readCount_1_1) + utils.bcolors.ENDC)
            module1 = {'1.1': int(readCount_1_1)}
            toWrite=[str(readCount_1_1), str(readCount_1_2)]
            return True, module1, toWrite
        else:
            toWrite=[str(readCount_1_1),str(readCount_1_2)]
            module1 = {'1.1': int(readCount_1_1), '1.2': int(readCount_1_2)}
            p_1_1 = utils.getProportionsModule1(module1, '1.1')
            p_1_2 = utils.getProportionsModule1(module1, '1.2')
            if p_1_1 >= proportionCutOff:
                print(utils.bcolors.OKGREEN + "\n\t- Choosing 1.1 module as new target, with {} reads (proportion = "
                      "{})".format(readCount_1_1, p_1_1) + utils.bcolors.ENDC)
                del module1['1.2']
                return True, module1, toWrite
            elif p_1_2 >= proportionCutOff:
                print(utils.bcolors.OKGREEN + "\n\t- Choosing 1.2 module as new target, with {} reads (proportion = "
                      "{}".format(readCount_1_2, p_1_2) + utils.bcolors.ENDC)
                del module1['1.1']
                return True, module1, toWrite
            else:
                print(utils.bcolors.FAIL + "\n\t- Inconclusive results on the 1.x module! Sample cannot be typed."
                      + utils.bcolors.ENDC)
                return False, None, toWrite


def alignSamples(sampleName, sampleFiles, reference, threads, workdir, script_path, keepFiles, minCoverage,
                 proportionCutOff, time, greaterThan, reportFile):

    success = False

    newWorkdir=os.path.join(workdir, sampleName, "tmp")
    if not os.path.isdir(newWorkdir):
        os.makedirs(newWorkdir)

    runMapping, samFile_fullRef = utils.mappingBowtie2(sampleFiles, reference, threads, newWorkdir, False, 1, '--seed 42',
                                                       True, sampleName)
    if runMapping:
        runSortAlignment, bamFile_fullRef = utils.sortAlignment(samFile_fullRef, str(os.path.splitext(samFile_fullRef)[0]
                                                                               + '.bam'), False, threads, True)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile_fullRef, False)
            if runIndex:

                pysam_fullRef = pysam.AlignmentFile(bamFile_fullRef, "rb")

                #fetching reads mapped to the conserved target region in hsdS
                iter=pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus',8532-1,8722-1, until_eof=True)
                matePairs_conserved = pysam.AlignmentFile(os.path.join(newWorkdir, sampleName +
                                                          "_pysam_matepairs_conserved_module1.bam"), "wb",
                                                          template=pysam_fullRef)

                for x in iter:

                    #TODO - remove exception - Throws error is mate is unmapped
                    #read needs to be mapped to the reverse strand of the reference and both mates need to be pair
                    try:
                        if x.is_paired and x.is_reverse:
                            mate=pysam_fullRef.mate(x)
                            matePairs_conserved.write(mate)
                    except:
                        pass

                matePairs_conserved.close()


                #TODO - this looks ugly..
                runSortAlignment, bam_matepairs = utils.sortAlignment(os.path.join(newWorkdir,
                                                       sampleName + "_pysam_matepairs_conserved_module1.bam"),
                                                       os.path.join(newWorkdir, sampleName +
                                                       "_pysam_matepairs_conserved_module1.bam"),
                                                       False, threads, False)
                utils.indexAlignment(bam_matepairs, False)
                utils.bam2fastq(bam_matepairs, False)

                print(utils.bcolors.BOLD + "\n--> Mapping Module 1" + utils.bcolors.ENDC)

                success1, module1, report1 = typeSeq_moduleOne([bam_matepairs+'.fastq'], threads, newWorkdir,
                                                          script_path, minCoverage, proportionCutOff)
                if success1:

                    print(utils.bcolors.BOLD + "\n--> Mapping Module 2" + utils.bcolors.ENDC)
                    success2, module2, report2 = getSeq_moduleTwo(module1.keys()[0], bamFile_fullRef, threads,
                                                                 newWorkdir, script_path, sampleName)
                    if success2:
                        successTyping, report3 = getType(workdir, sampleName, module1, module2, minCoverage)

                        if successTyping:
                            success = True
                            writeReport(workdir, sampleName, time, report1, report2, report3, greaterThan,
                                        minCoverage, proportionCutOff, reportFile)
                        else:
                            writeReport(workdir, sampleName, time, report1, report2, ['NA'] * 6, greaterThan,
                                        minCoverage, proportionCutOff, reportFile)
                    else:
                        writeReport(workdir, sampleName, time, report1, report2, ['NA'] * 6, greaterThan, minCoverage,
                                    proportionCutOff, reportFile)
                else:
                    writeReport(workdir, sampleName, time, report1, ['NA'] * 3, ['NA'] * 6, greaterThan, minCoverage,
                                proportionCutOff, reportFile)
    else:
        writeReport(workdir, sampleName, time, ['NA']*2, ['NA']*3, ['NA']*6, greaterThan, minCoverage,
                    proportionCutOff, reportFile)

    if not keepFiles:
        #TODO - this is STILL not removing /tmp/ folder in windows
        shutil.rmtree(newWorkdir+'/', ignore_errors=True)
        if os.path.exists(newWorkdir) and not os.listdir(newWorkdir):
            os.remove(newWorkdir)

    return success
