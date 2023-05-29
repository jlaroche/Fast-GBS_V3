#!/prg/python/3.5/bin/python3.5

"""
Permet de scinder les séquences de chromosomes dont la taille excède 500Mb.

Usage:

/project/fbelzile/prg/SplitLargeChromosome.py <fichier_genome.fasta>

Le fichier de sortie s'appellera: genome_split.fasta
Le renommer a votre convenance ;-)
"""

import os, sys, fileinput, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC


try:
	file = sys.argv[1]
except:
	print(__doc__)
	sys.exit(1)


dot = file.index(".")
ext = file[dot+1:]
output = file[:dot]


#print("Seq","\t","Longueur".rjust(12))
entries=[]
for record in SeqIO.parse(file, "fasta"):

#	entry=SeqRecord(record.seq,id=record.name,description=record.description)		
#	entries.append(entry)

	long=len(record.seq)
	half=round(long/2)
	print(record.name,str(round(long/2)),str(long))

	entry=SeqRecord(record.seq[0:half],id=record.name+'_part1',description=record.description)		
	entries.append(entry)
	entry=SeqRecord(record.seq[half:long],id=record.name+'_part2',description=record.description)		
	entries.append(entry)
	
		
SeqIO.write(entries,output+"_split.fasta","fasta")	
