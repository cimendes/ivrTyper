import os
import os.path
import utils
import pysam
import shutil


def writeReport(workdir, sample, module1, module2, type):
    fname=os.path.join(workdir, "report.csv")
    if not os.path.isfile(fname):
        with open(fname, "w") as report:
            report.write("Sample,Module 1,Module 2,hsdS Type\n")
            report.write("%s,%s,%s,%s\n" % (sample, module1, module2, type))
    else:
        with open(fname, "a") as report:
            report.write("%s,%s,%s,%s\n" % (sample, module1, module2, type))


def typeSeq_moduleOne(files, threads, workdir, script_path, minCoverage):

    readCount_1_1=None
    readCount_1_2=None

    #1.1
    runMapping, samFile1 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '1.1.fasta'), threads, workdir, False, 1, None, True, "1.1")
    if runMapping:
        runSortAlignment, bamFile1 = utils.sortAlignment(samFile1, str(os.path.splitext(samFile1)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile1, False)
            if runIndex:
                samfile1 = pysam.AlignmentFile(bamFile1, "rb")
                readCount_1_1 = samfile1.mapped
    else:
        print 'Failed 1.1 Bowtie mapping'
        return False, readCount_1_1

    #1.2
    runMapping, samFile2 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '1.2.fasta'), threads, workdir, True, 1, None, True, "1.2")
    if runMapping:
        runSortAlignment, bamFile2 = utils.sortAlignment(samFile2, str(os.path.splitext(samFile2)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile2, False)
            if runIndex:
                samfile2 = pysam.AlignmentFile(bamFile2, "rb")
                readCount_1_2 = samfile2.mapped
    else:
        print "Failed 1.2 Bowtie mapping"
        return False, readCount_1_2

    print "\n1.1 - " + str(readCount_1_1)
    print "1.2 - " + str(readCount_1_2) + '\n'

    new_target = None

    if readCount_1_1 > readCount_1_2:
        if readCount_1_1 > minCoverage:
            new_target = '1.1'
        else:
            print "coverage on the 1.1 module below the minCoverage %s read threshold " % (minCoverage)
            print "results inconclusive"
            return False, None, None
    elif readCount_1_1 < readCount_1_2:
        if readCount_1_1 > minCoverage:
            new_target = '1.2'
        else:
            print "coverage on the 1.2 module below the minCoverage %s read threshold " % (minCoverage)
            print "results inconclusive"
            return False, None, None
    else:
        print "results inconclusive"
        return False, None, None

    print "\t- New target sequence: module %s \n" % (new_target)
    module1={'1.1': int(readCount_1_1), '1.2': int(readCount_1_2)}
    return True, new_target, module1

def typeSeq_moduleTwo(files, threads, workdir, script_path):
    readCount_2_1 = None
    readCount_2_2 = None
    readCount_2_3 = None

    # 2.1
    runMapping, samFile1 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '2.1.fasta'), threads, workdir, False, 1, None, True, "2.1")
    if runMapping:
        runSortAlignment, bamFile1 = utils.sortAlignment(samFile1, str(os.path.splitext(samFile1)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile1, False)
            if runIndex:
                samfile1 = pysam.AlignmentFile(bamFile1, "rb")
                readCount_2_1 = samfile1.mapped
    else:
        print 'Failed 2.1 Bowtie mapping'
        return False, None

    # 2.2
    runMapping, samFile2 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '2.2.fasta'), threads, workdir, False, 1, None, True, "2.2")
    if runMapping:
        runSortAlignment, bamFile2 = utils.sortAlignment(samFile2, str(os.path.splitext(samFile2)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile2, False)
            if runIndex:
                samfile2 = pysam.AlignmentFile(bamFile2, "rb")
                readCount_2_2 = samfile2.mapped
    else:
        print 'Failed 2.2 Bowtie mapping'
        return False, None

    # 2.3
    runMapping, samFile3 = utils.mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq',
                                                '2.3.fasta'), threads, workdir, False, 1, None, True, "2.3")
    if runMapping:
        runSortAlignment, bamFile3 = utils.sortAlignment(samFile3, str(os.path.splitext(samFile3)[0] + '.bam'), False,
                                                         threads, False)
        if runSortAlignment:
            runIndex = utils.indexAlignment(bamFile3, False)
            if runIndex:
                samfile3 = pysam.AlignmentFile(bamFile3, "rb")
                readCount_2_3 = samfile3.mapped
    else:
        print 'Failed 2.3 Bowtie mapping'
        return False, None

    print "\n2.1 - " + str(readCount_2_1)
    print "2.2 - " + str(readCount_2_2)
    print "2.3 - " + str(readCount_2_3) + '\n'

    module2 = {"2.1": int(readCount_2_1), "2.2": int(readCount_2_2), "2.3": int(readCount_2_3)}
    return True, module2


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

        run, module2 =typeSeq_moduleTwo([bamFile_1_1 + '.fastq'], threads, workdir, script_path)

        return run, module2

    elif first_type == '1.2':
        iter = pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus', 4462 - 1, 4861 - 1)
        matePairs_conserved = pysam.AlignmentFile(sample + "_matepairs_conserved_module2.bam", "wb",
                                                  template=pysam_fullRef)

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
        runSortAlignment, bamFile3 = utils.sortAlignment(os.path.join(workdir,
                                                         sample + "_matepairs_conserved_module2.bam"),
                                                         os.path.join(workdir, sample +
                                                         "_matepairs_conserved_module2.bam"),
                                                         False, threads, False)
        utils.indexAlignment(bamFile3, False)
        utils.bam2fastq(bamFile3, False)

        run, module2 = typeSeq_moduleTwo([bamFile3 + '.fastq'], threads, workdir, script_path)

        return run, module2

    else:
        print "this sample is untypable"
        return False, None

def getType(module1, module2, minCoverage):
    #TODO - implement percentage!! - confidence calling

    totalReadCount = sum(module1.values() + module2.values())

    if  module2["2.1"] < minCoverage:
        print "coverage on the 2.1 module below the minCoverage %s read threshold " % (minCoverage)
        module2["2.1"] = 0
    elif module2["2.2"] < minCoverage:
        print "coverage on the 2.2 module below the minCoverage %s read threshold " % (minCoverage)
        module2["2.2"] = 0
    elif module2["2.3"] < minCoverage:
        print "coverage on the 2.3 module below the minCoverage %s read threshold " % (minCoverage)
        module2["2.3"] = 0

    print "\nTyping support:\n"

    #readPercentageA = float((module1["1.1"] + module2["2.1"])) / float(totalReadCount) * 100
    readPercentageA = float(float(module1["1.1"]/float(module1["1.1"]+module1["1.2"]))+float(module2["2.1"]/float(
        module2["2.1"]+module2["2.2"]+module2["2.3"])))/2
    print "allele A (1.1 2.1): {} \n".format(format(readPercentageA, '.2f'))

    #readPercentageB = float((module1["1.1"] + module2["2.2"])) / float(totalReadCount) * 100
    readPercentageB = float(float(module1["1.1"] / module1["1.1"] + module1["1.2"]) + float(module2["2.2"] / module2[
        "2.1"] + module2["2.2"] + module2["2.3"])) / 2
    print "allele B (1.1 2.2): {} \n".format(format(readPercentageB, '.2f'))

    readPercentageC = float(float(module1["1.2"] / module1["1.1"] + module1["1.2"]) + float(module2["2.2"] / module2[
        "2.1"] + module2["2.2"] + module2["2.3"])) / 2
    print "allele C (1.2 2.2): {} \n".format(format(readPercentageC, '.2f'))

    readPercentageD = float(float(module1["1.2"] / module1["1.1"] + module1["1.2"]) + float(module2["2.1"] / module2[
        "2.1"] + module2["2.2"] + module2["2.3"])) / 2
    print "allele D (1.2 2.1): {}\n".format(format(readPercentageD, '.2f'))

    readPercentageE = float(float(module1["1.1"] / module1["1.1"] + module1["1.2"]) + float(module2["2.3"] / module2[
        "2.1"] + module2["2.2"] + module2["2.3"])) / 2
    print "allele E (1.1 2.3): {}\n".format(format(readPercentageE, '.2f'))

    readPercentageF = float(float(module1["1.2"] / module1["1.1"] + module1["1.2"]) + float(module2["2.3"] / module2[
        "2.1"] + module2["2.2"] + module2["2.3"])) / 2
    print "allele F (1.2 2.3): {}\n".format(format(readPercentageF, '.2f'))



def alignSamples(sampleFiles, reference, threads, workdir, script_path, keepFiles, minCoverage):

    for sample, files in sorted(sampleFiles.items()):
        print '\n-> ' + sample + '\n'

        newWorkdir=os.path.join(workdir, sample, "tmp")
        if not os.path.isdir(newWorkdir):
            os.makedirs(newWorkdir)

        runMapping, samFile_fullRef = utils.mappingBowtie2(files, reference, threads, newWorkdir, False, 1, None,
                                                           True, sample)

        #TODO - is sorting needed?
        if runMapping:
            runSortAlignment, bamFile_fullRef = utils.sortAlignment(samFile_fullRef, str(os.path.splitext(samFile_fullRef)[0]
                                                                                   + '.bam'), False, threads, False)
            if runSortAlignment:
                runIndex = utils.indexAlignment(bamFile_fullRef, False)
                if runIndex:
                    pysam_fullRef = pysam.AlignmentFile(bamFile_fullRef, "rb")

                    #fetching reads mapped to the conserved target region in hsdS
                    iter=pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus',8532-1,8722-1)
                    matePairs_conserved = pysam.AlignmentFile(os.path.join(newWorkdir, sample +
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
                                                           sample + "_pysam_matepairs_conserved_module1.bam"),
                                                           os.path.join(newWorkdir,sample +
                                                           "_pysam_matepairs_conserved_module1.bam"),
                                                           False, threads, False)
                    utils.indexAlignment(bam_matepairs, False)
                    utils.bam2fastq(bam_matepairs, False)
                    success1, newTarget, module1 = typeSeq_moduleOne([bam_matepairs+'.fastq'], threads, newWorkdir,
                                                              script_path, minCoverage)

                    success2, module2 = getSeq_moduleTwo(newTarget, bamFile_fullRef, threads, newWorkdir,
                                                             script_path, sample)

                    if success1 and success2:
                         getType(module1,module2, minCoverage)

                    #    print "--> Sample has a type %s ivr locus!" % (ivrType)
                    #    writeReport(workdir, sample, first_unit, second_unit, ivrType)
                    #else:
                    #    print "No ivr locus found for this sample"


        if not keepFiles:
            #TODO - this is STILL not removing /tmp/ folder
            shutil.rmtree(newWorkdir+'/', ignore_errors=True)
            if os.path.exists(newWorkdir) and not os.listdir(newWorkdir):
                os.remove(newWorkdir)
