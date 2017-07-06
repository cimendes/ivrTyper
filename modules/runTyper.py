import os
import os.path
import utils
import pysam



# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
    if os.path.isfile(str(referenceFile + '.1.bt2')):
        run_successfully = True
    else:
        command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, True)
    return run_successfully


# Mapping with Bowtie2
def mappingBowtie2(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc, bowtieOPT):
    sam_file = os.path.join(outdir, str('alignment.sam'))

    # Index reference file
    run_successfully = indexSequenceBowtie2(referenceFile, threads)

    if run_successfully:
        command = ['bowtie2', '-k', str(numMapLoc), '-q', '', '--threads', str(threads), '-x', referenceFile, '', '--no-unal', '', '-S', sam_file]

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

        run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

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
def indexAlignment(alignment_file):
    command = ['samtools', 'index', alignment_file]
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    return run_successfully


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads):
    outFormat_string = os.path.splitext(output_file)[1][1:].lower()
    command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
    if sortByName_True:
        command[6] = '-n'
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
    if not run_successfully:
        output_file = None
    return run_successfully, output_file

def typeSeq(seq, script_path):
    from Bio import pairwise2

    #print script_path

    status = 'NA'
    #1.1
    alignments = pairwise2.align.globalxx(seq, os.path.join(os.path.dirname(script_path),'1.1.fasta'),score_only=True, one_alignment_only=True)
    print alignments





    return status
def alignSamples(sampleFiles, reference, threads, workdir, script_path):

    for sample, files in sampleFiles.items():
        print '-> ' + sample

        runMapping, samFile = mappingBowtie2(files, reference, threads, workdir, True, 1, None)

        #TODO - is sorting needed?
        if runMapping:
            runSortAlignment, bamFile = sortAlignment(samFile, str(os.path.splitext(samFile)[0] + '.bam'), False, threads)
            if runSortAlignment:
                runIndex = indexAlignment(bamFile)
                if runIndex:
                    samfile = pysam.AlignmentFile(bamFile, "rb")
                    pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)
                    for read in samfile.fetch():
                        #TODO - check if both reads are mapped
                        if read.is_paired:
                            pairedreads.write(read)
                    pairedreads.close()

                    runSortAlignment, bamFile2 = sortAlignment("allpaired.bam", "allpaired.bam", False, threads)
                    runIndex = indexAlignment(bamFile2)
                    newPaired=pysam.AlignmentFile("allpaired.bam","rb")

                    #fetching reads mapped to the conserved target
                    iter=newPaired.fetch('CP000410_extraction_-_Type_I_RM_system_locus',8532-1,8722-1)
                    for x in iter:
                        #TODO - add verification that they are reverse
                        if x.is_paired:
                            try:
                                mate=newPaired.mate(x)
                                mateseq = str(mate).split()[9]
                                #print mateseq
                                classification = typeSeq(mateseq, script_path)

                                #startpos = mate.get_reference_positions()[0]
                                #endpos = mate.get_reference_positions()[-1]
                                #TODO - not needed?

                                #if startpos>=8029-1 and endpos<=8446-1:
                                    #mate= str(mate)
                                    #print mate.split()[9]
                                    #print mate.query_sequence()



                            except:
                                print "ah!"



        '''
        if runMapping:
            runSortAlignment, bamFile = sortAlignment(samFile, str(os.path.splitext(samFile)[0] + '.bam'), False,args.threads)
            if runSortAlignment:
                os.remove(samFile)
                # Index bam
                run_successfully = indexAlignment(bamFile)
                if run_successfully:
                    runIndexSamtools, stdout = index_fasta_samtools(reference, None, None, True)
        '''