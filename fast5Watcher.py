#!/usr/bin/env python
import h5py
import sys
from Bio import SeqIO
from StringIO import StringIO
from subprocess import check_output
import glob
import os
import time
import datetime

def file2fq(fi,fo,sfo):
	f = h5py.File(fi,'r')
	try:
		key="/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
		fq = f[key][()]
	except:
		key="/Analyses/Basecall_RNN_1D_000/BaseCalled_template/Fastq"
		fq = f[key][()]
	rec = SeqIO.read(StringIO(fq), "fastq")
	rec.description=str("{0} {1}".format(rec.id,fi))
	SeqIO.write([rec], fo, "fastq")

	expstart =f["UniqueGlobalKey"]["tracking_id"].attrs.get("exp_start_time")
	# stats
	read=f["Raw"]["Reads"].keys()[0]
	#try:
	#	rstart = f["Analyses"]["Basecall_1D_000"]["Summary"]["basecall_1d_template"].attrs.get("start_time")
	#	start = int(expstart) + int(rstart)
	#except:
	#	#rstart = None #f["Analyses"]["Basecall_RNN_1D_000"]["Summary"]["basecall_1d_template"].attrs.get("start_time")
	#	rstart = None #f["Raw"]["Reads"][read].attrs.get("start_time")
	#if rstart == None:
	tstamp = f["Analyses"]["Segmentation_000"].attrs.get("time_stamp")
	start = int(time.mktime(datetime.datetime.strptime(tstamp, "%Y-%m-%dT%H:%M:%SZ").timetuple()))
	rstart=start-int(expstart)
	l=len(rec.seq)
	avQual=sum(rec.letter_annotations["phred_quality"])/float(l)
	accu = (1-(10**(-avQual/10)))*100
	sfo.write("{0},{1},{2},{3},{4},{5},{6}\n".format(rec.id,start,l,rstart,expstart,avQual,accu))


def watcher(folder):
	ol=[]
	stf=open(sys.argv[3],'a') # stat file
	print "getting old name"
	#ol=getOldNames(sys.argv[3])
	try:
		ol=getFastqFileName(folder,sys.argv[2])
	except:
		print "First time ey"
		ol=[]
	print "getting all names"
	#nl = glob.glob("{0}/*fast5".format(folder))
	#nl = os.listdir("{0}".format(folder))
	out = check_output(["find", "{0}".format(sys.argv[1]),"-name","*.fast5"])
	out=out.replace(folder,'')
	nl=out.split("\n")
	del nl[-1]
	#print nl
	#print ol
	print "subtracting names"
	#newfiles=[item for item in nl if item not in ol]
	newfiles=set(nl)-set(ol)
	outf=open(sys.argv[2],'a')
	print "extracting reads"
	n=1
	x=len(newfiles)
	for nf in newfiles:
		#file2fq("{0}/{1}".format(folder,nf),outf)
		p=int((float(n)/x)*100)
		sys.stdout.write("\r{0} reads of {1}, {2}%".format(n,x,p))
		n+=1
		try:
			file2fq("{0}/{1}".format(folder,nf),outf,stf)
		except:
			print "Failed {0}".format(nf)
	sys.stdout.write("\ndone\n")
		#ln.write(nf+'\n')
	#infile=sys.argv[1]	
	#outf=open(sys.argv[2],'a')

def getFastqFileName(fol,fq):
	files=[]
	for sr in SeqIO.parse(open(fq),'fastq'):
		#files.append("{0}/{1}".format(fol,sr.description.split(' ')[-1].split('/')[-1]))
		files.append("{0}".format(sr.description.split(' ')[-1].split('/')[-1]))
	return files


def getOldNames(f):
	olist=[]
	for line in open(f,'r'):
		olist.append(line.replace('\n',''))
	return olist

watcher(sys.argv[1])
