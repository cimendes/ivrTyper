import os
import os.path
import utils
import pysam
import glob
import shutil



# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads, write_command=True):
    if os.path.isfile(str(referenceFile + '.1.bt2')):
        run_successfully = True
    else:
        command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, write_command)
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

        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, write_command)

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
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, shell_true, None, print_comand_true)
    return run_successfully, stdout


# Index alignment file
def indexAlignment(alignment_file, write_command=True):
    command = ['samtools', 'index', alignment_file]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, write_command)
    return run_successfully


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads, write_command=True):
    outFormat_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
    if sortByName_True:
        command[6] = '-n'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, write_command)
    if not run_successfully:
        output_file = None
    return run_successfully, output_file


def bam2fastq(bam_file, write_command=True):
    command = ['samtools', 'fastq', bam_file, '>', bam_file+'.fastq']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, write_command)
    return run_successfully


def typeSeq_moduleOne(files, threads, workdir, script_path):

    readCount_1_1=None
    readCount_1_2=None
    #1.1
    runMapping, samFile1 = mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq', '1.1.fasta'),
                                          threads, workdir, False, 1, None, True, "1.1")
    if runMapping:
        runSortAlignment, bamFile1 = sortAlignment(samFile1, str(os.path.splitext(samFile1)[0] + '.bam'), False,
                                                   threads, False)
        if runSortAlignment:
            runIndex = indexAlignment(bamFile1, False)
            if runIndex:
                samfile1 = pysam.AlignmentFile(bamFile1, "rb")
                readCount_1_1 = samfile1.mapped
    else:
        print 'Failed 1.1 Bowtie mapping'
        return False, readCount_1_1

    #1.2
    runMapping, samFile2 = mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq', '1.2.fasta'),
                                          threads, workdir, True, 1, None, True, "1.2")
    if runMapping:
        runSortAlignment, bamFile2 = sortAlignment(samFile2, str(os.path.splitext(samFile2)[0] + '.bam'), False,
                                                   threads, False)
        if runSortAlignment:
            runIndex = indexAlignment(bamFile2, False)
            if runIndex:
                samfile2 = pysam.AlignmentFile(bamFile2, "rb")
                readCount_1_2 = samfile2.mapped
    else:
        print "Failed 1.2 Bowtie mapping"
        return False, readCount_1_2

    print "1.1 - " + str(readCount_1_1)
    print "1.2 - " + str(readCount_1_2)

    if readCount_1_1 > readCount_1_2:
        return True, '1.1'
    elif readCount_1_1 < readCount_1_2:
        return True, '1.2'
    else:
        print "results inconclusive"
        return False, None

def typeSeq_moduleTwo(files, threads, workdir, script_path):
    readCount_2_1 = None
    readCount_2_2 = None
    readCount_2_3 = None

    # 2.1
    runMapping, samFile1 = mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq', '2.1.fasta'),
                                          threads, workdir, False, 1, None, True, "2.1")
    if runMapping:
        runSortAlignment, bamFile1 = sortAlignment(samFile1, str(os.path.splitext(samFile1)[0] + '.bam'), False,
                                                   threads, False)
        if runSortAlignment:
            runIndex = indexAlignment(bamFile1, False)
            if runIndex:
                samfile1 = pysam.AlignmentFile(bamFile1, "rb")
                readCount_2_1 = samfile1.mapped
    else:
        print 'Failed 2.1 Bowtie mapping'
        return False, readCount_2_1

    # 2.2
    runMapping, samFile2 = mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq', '2.2.fasta'),
                                          threads, workdir, False, 1, None, True, "2.2")
    if runMapping:
        runSortAlignment, bamFile2 = sortAlignment(samFile2, str(os.path.splitext(samFile2)[0] + '.bam'), False,
                                                   threads, False)
        if runSortAlignment:
            runIndex = indexAlignment(bamFile2, False)
            if runIndex:
                samfile2 = pysam.AlignmentFile(bamFile2, "rb")
                readCount_2_2 = samfile2.mapped
    else:
        print 'Failed 2.2 Bowtie mapping'
        return False, readCount_2_2

    # 2.3
    runMapping, samFile3 = mappingBowtie2(files, os.path.join(os.path.dirname(script_path), 'src', 'seq', '2.3.fasta'),
                                          threads, workdir, False, 1, None, True, "2.3")
    if runMapping:
        runSortAlignment, bamFile3 = sortAlignment(samFile3, str(os.path.splitext(samFile3)[0] + '.bam'), False,
                                                   threads, False)
        if runSortAlignment:
            runIndex = indexAlignment(bamFile3, False)
            if runIndex:
                samfile3 = pysam.AlignmentFile(bamFile3, "rb")
                readCount_2_3 = samfile3.mapped
    else:
        print 'Failed 2.3 Bowtie mapping'
        return False, readCount_2_3

    print "2.1 - " + str(readCount_2_1)
    print "2.2 - " + str(readCount_2_2)
    print "2.3 - " + str(readCount_2_3)

    if readCount_2_1 == readCount_2_2 == readCount_2_3:
        print "Results inconclusive"
        return False, None
    elif readCount_2_1 > readCount_2_2 and readCount_2_1 > readCount_2_3:
        return True, "2.1"
    elif readCount_2_2 > readCount_2_1 and readCount_2_2 > readCount_2_3:
        return True, "2.2"
    elif readCount_2_3 > readCount_2_1 and readCount_2_3 > readCount_2_2:
        return True, "2.3"


def getSeq_moduleTwo(first_type, bamfile, threads, workdir, script_path):

    pysam_fullRef = pysam.AlignmentFile(bamfile, "rb")

    if first_type == '1.1':
        iter = pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus', 8029 - 1, 8446 - 1)
        matePairs_conserved = pysam.AlignmentFile("matepairs_conserved_module2.bam", "wb", template=pysam_fullRef)

        #TODO - remove verification
        for x in iter:
            # read needs to be mapped to the reverse strand of the reference and both mates need to be pair
            if x.is_proper_pair and x.is_reverse:
                mate = pysam_fullRef.mate(x)
                matePairs_conserved.write(mate)
        matePairs_conserved.close()

        # TODO - this looks ugly..
        runSortAlignment, bamFile_1_1 = sortAlignment("matepairs_conserved_module2.bam", "matepairs_conserved_module2.bam",
                                                   False, threads, False)
        indexAlignment(bamFile_1_1, False)
        bam2fastq(bamFile_1_1, False)
        run, type =typeSeq_moduleTwo([bamFile_1_1 + '.fastq'], threads, workdir, script_path)
        return run, type

    elif first_type == '1.2':
        iter = pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus', 4462 - 1, 4861 - 1)
        matePairs_conserved = pysam.AlignmentFile("matepairs_conserved_module2.bam", "wb", template=pysam_fullRef)

        #TODO - remove verification...
        for x in iter:
            # read needs to be mapped to the forward strand of the reference and both mates need to be pair
            if x.is_proper_pair and not x.is_reverse:
                mate = pysam_fullRef.mate(x)
                matePairs_conserved.write(mate)
        matePairs_conserved.close()

        # TODO - this looks ugly..
        runSortAlignment, bamFile3 = sortAlignment("matepairs_conserved_module2.bam", "matepairs_conserved_module2.bam",
                                                   False, threads, False)
        indexAlignment(bamFile3, False)
        bam2fastq(bamFile3, False)
        run, type = typeSeq_moduleTwo([bamFile3 + '.fastq'], threads, workdir, script_path)
        return run, type

    else:
        print "this sample is untypable"
        return False, None

def getType(module1, module2):

    if module1 == '1.1':
        if module2 == "2.1":
            return "A"
        elif module2 == "2.2":
            return "B"
        elif module2 == "2.3":
            return "E"
    elif module1 == "1.2":
        if module2 == "2.1":
            return "D"
        elif module2 == "2.2":
            return "C"
        elif module2 == "2.3":
            return "F"
    else:
        return None


def alignSamples(sampleFiles, reference, threads, workdir, script_path, keepFiles):

    for sample, files in sampleFiles.items():
        print '\n-> ' + sample

        newWorkdir=os.path.join(workdir, sample, "tmp")
        if not os.path.isdir(newWorkdir):
            os.makedirs(newWorkdir)

        runMapping, samFile_fullRef = mappingBowtie2(files, reference, threads, newWorkdir, False, 1, None, True, sample)

        #TODO - is sorting needed?
        if runMapping:
            runSortAlignment, bamFile_fullRef = sortAlignment(samFile_fullRef, str(os.path.splitext(samFile_fullRef)[
                                                                                      0] + '.bam'), False,
                                                      threads, False)
            if runSortAlignment:
                runIndex = indexAlignment(bamFile_fullRef, False)
                if runIndex:
                    pysam_fullRef = pysam.AlignmentFile(bamFile_fullRef, "rb")

                    #fetching reads mapped to the conserved target region in hsdS
                    iter=pysam_fullRef.fetch('CP000410_extraction_-_Type_I_RM_system_locus',8532-1,8722-1)
                    matePairs_conserved = pysam.AlignmentFile(sample + "_pysam_matepairs_conserved_module1.bam", "wb",
                                                              template=pysam_fullRef)

                    for x in iter:
                        try:
                            mate = pysam_fullRef.mate(x)
                            matePairs_conserved.write(mate)
                        except:
                            pass

                        #read needs to be mapped to the reverse strand of the reference and both mates need to be pair
                        '''if x.is_proper_pair and x.is_reverse: #TODO - remove verification
                            mate=pysam_fullRef.mate(x)
                            matePairs_conserved.write(mate)'''
                    matePairs_conserved.close()

                    #TODO - this looks ugly..
                    runSortAlignment, bam_matepairs = sortAlignment(sample + "_pysam_matepairs_conserved_module1.bam",
                                                           sample + "_pysam_matepairs_conserved_module1.bam",
                                                               False, threads, False)
                    indexAlignment(bam_matepairs, False)
                    bam2fastq(bam_matepairs, False)
                    success1, first_unit = typeSeq_moduleOne([bam_matepairs+'.fastq'], threads, newWorkdir, script_path)

                    success2, second_unit = getSeq_moduleTwo(first_unit, bamFile_fullRef, threads, newWorkdir, script_path)

                    if success1 and success2:
                        ivrType = getType(first_unit,second_unit)
                        print "--> Sample has a type %s ivr locus!" % (ivrType)
                    else:
                        print "No ivr locus found for this sample"

        if not keepFiles:
            #TODO - this is not removing /tmp/ folder
            shutil.rmtree(newWorkdir+'/', ignore_errors=True)
            #TODO - move this to /tmp/ in the sample folder
            print script_path + '/' + sample + "*"
            for filename in glob.glob(os.path.join(script_path, sample + "*")):
                os.remove(filename)



