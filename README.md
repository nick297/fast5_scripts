# fast5_scripts

Simple python scripts for working with nanopore data. Other software is available and recommended including poretools and poRe. Some code shamelessly borrowed from poretools.

fast5Watcher
============

fast5Watcher was written for extracting fastq and timing data from MinKNOW basecalled fast5 files. Other software was available and will be better supported than this, however after a MinKNOW update they were not updated for the new internal HDF file structure and I am impatient. The file was subsequently updated to use linux/bash find command to account for the large number of files produced. This also allowed me to create two lists, one of all the files in the directory, and one of the written fastq file of all the previously seen fast5 files. The lists are then subtracted so only the newest files are processed. Linux find is very fast and is subject to cacheing so first use can be slow but subsequent uses are faster.

usage
-----

	fast5Watcher.py /path/to/fast5/dirctory/ out.fq out.timestats.csv

nanoStats
=========

Simple read length and quality statistics from fastq file(s). Can also output a histogram. One or more fastq files can be used.

usage
-----
	nanoStats.py out.fq [out1.fq out2.fq...]


