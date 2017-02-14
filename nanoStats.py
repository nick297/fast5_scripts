#!/usr/bin/env python
import sys
from Bio import SeqIO
import gzip
import numpy as np
from collections import Counter
from matplotlib import pyplot

class nanoStatistics:
	def __init__(self, fq):
		self.lens=[]
		self.qs=[]
		self.fq=fq
		self.getLists()
		self.lenStats()
		self.qualStats()

	def getLists(self):
		for sr in SeqIO.parse(open(self.fq),'fastq'):
			self.lens.append(len(sr.seq))
			self.qs.append(sr.letter_annotations["phred_quality"])

	def lenStats(self):
		self.bases=sum(self.lens)
		self.reads=len(self.lens)
		self.av=self.bases/float(self.reads)
		self.longest=max(self.lens)
		self.shortest=min(self.lens)
		self.median=np.median(np.array(self.lens))
		data = Counter(self.lens)
		self.mode=data.most_common(1)[0] 

	def qualStats(self):
		self.allqs=[item for sublist in self.qs for item in sublist]
		self.avQual=sum(self.allqs)/float(self.bases)
		#self.accu = 10**(-self.avQual/10)
		self.accu = (1-(10**(-self.avQual/10)))*100

def histo(lst):
	fig=pyplot.figure()
	bins=np.linspace(0,5000,100)
	for l in lst:
		pyplot.hist(l.lens,bins,alpha=0.5)#, label='')
	#pyplot.legend(loc='upper right')
	#pyplot.show()
	fig.savefig("ReadLengthsHist.pdf")

cl=[]
print "Bases,Reads,AvLen,Longest,Shortest,Median,Modal,AvPHRED,AvAccuracy%"
for i in sys.argv[1:]:
	nano=nanoStatistics(i)
	print "{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(nano.bases,nano.reads,nano.av,nano.longest,nano.shortest,nano.median,nano.mode,nano.avQual,nano.accu)
	cl.append(nano)


histo(cl)




