#!/usr/bin/python
# mk_PEAR_report.py: Looking for the success of the PEAR process.  The percent of assembled reads.  If the percent of assembled reads is low then the PEAR parameters might need to be tweaked.
# usage :
#  python mk_PEAR_report.py /home/mbonteam/MBARI/reiko/processed/Analysis_20160108_0100
#
# Input argument is the directory path to the processing logfile.txt,ie /home/mbonteam/MBARI/reiko/processed/Analysis_20160108_0100
# Parses the file looking for the R1 fastq file, writes that and the percentage of merged reads
# Output looks like:
# Forward reads file.................: /home/mbonteam/MBARI/reiko/raw/BOG_data/BOG18S12_S12_L001_R1_001.fastq
# Reverse reads file.................: /home/mbonteam/MBARI/reiko/raw/BOG_data/BOG18S12_S12_L001_R2_001.fastq
# Assembled reads ...................: 150,403 / 153,994 (97.668%)
# Discarded reads ...................: 3,591 / 153,994 (2.332%)
# Not assembled reads ...............: 0 / 153,994 (0.000%)

import sys, getopt, os

def main(argv):
	ddir = sys.argv[1]
#	print(ddir)
	os.chdir(ddir)
	print(os.getcwd())
#	searchquery = 'Forward reads file'
	strings = ["Forward reads file", "Reverse reads file", "Assembled reads ..", "Discarded reads ..", "Not assembled reads .."]
	cnt=1
	cnts=0
	asum=0
	nsum=0
	
	with open('logfile.txt') as f1:
		with open('PEAR_read_report.txt', 'a') as f2:
			lines = f1.readlines()
			for i, line in enumerate(lines):
				if line.startswith("Analysis started at"):
					f2.write(line)					
				if line.startswith(tuple(strings)):
					if cnt:
						cnt=0
						f2.write(lines[i + 2]) #PEAR parameters
						f2.write(lines[i + 3])
						f2.write(lines[i + 4])
						f2.write(lines[i + 5])
						f2.write(lines[i + 6])
						f2.write(lines[i + 7])
						f2.write(lines[i + 8])
						f2.write(lines[i + 9])
						f2.write(lines[i + 10])
						f2.write(lines[i + 11])
						f2.write(lines[i + 12])
						f2.write(lines[i + 13])
						f2.write(lines[i + 14]) # should be blank line
						f2.write(line)
						bline=lines[i + 14]
					else:
						f2.write(line)      
				if line.startswith("Discarded reads .."):
					f2.write(bline)
					nread=line.split("(")
					nread2=nread[1]
					nread3=nread2.split("%")
					nsum=nsum+float(nread3[0])
				if line.startswith("Assembled reads .."):
					aread=line.split("(")
					aread2=aread[1]
					aread3=aread2.split("%")
					asum=asum+float(aread3[0])
					cnts=cnts+1
			nmean=nsum/cnts
			amean=asum/cnts

			f2.write("%s = %f%c\n" % ("Assembled read mean", amean,"%"))
			f2.write("%s = %f%c\n" % ("Discarded read mean", nmean,"%"))
			f2.write("%s = %i\n" % ("Count", cnts))
#			f2.write(amean)
#			f2.write(nmean)
			print("PEAR read result summary")
			print(" Average assembled reads: ",amean,"%")
			print(" Average discarded  reads: ",nmean,"%")
			print(" Total samples: ",cnts)
			
if __name__ == "__main__":
   main(sys.argv[1:])