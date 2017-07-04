import multiprocessing
import utils
import functools
import sys
import pickle
import os
import traceback
import shlex
import subprocess
from threading import Timer
import shutil
import time
import os.path

# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
	if os.path.isfile(str(referenceFile + '.1.bt2')):
		run_successfully = True
	else:
		command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
		run_successfully, stdout, stderr = runCommandPopenCommunicate(command, True, None, True)
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

		run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None, True)

	if not run_successfully:
		sam_file = None

	return run_successfully, sam_file

def kill_subprocess_Popen(subprocess_Popen, command):
	print 'Command run out of time: ' + str(command)
	subprocess_Popen.kill()

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

