#!/prg/python/3.5/bin/python3.5

"""
Permet de remplacer les positions des variants dans le fichier source
pour celles d'un génome scindé.

Usage:

modify_source_positions.py <fichier genome.fasta> <fichier source>

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

#./modify_source_positions.py Barley_MorexV3_pseudomolecules.fasta source_file

try:
	file = sys.argv[1]
except:
	print(__doc__)
	sys.exit(1)

try:
	source = sys.argv[2]
except:
	print(__doc__)
	sys.exit(1)

dot = source.index(".")
ext = source[dot+1:]
output = source[:dot]

# A partir du fichier du genome en format fasta, construire un dictionnaire
# qui contiendra les limites (moitiée,longueur) des séquences scindées
# C'est surtour la valeur de la moitiée de la séquence qui sera utilisée
dict={}
for record in SeqIO.parse(file, "fasta"):
	long=len(record.seq)
	half=round(long/2)
#	print(record.name,str(round(long/2)),str(long))
	dict[record.name]=[]
	dict[record.name].append(half)
	dict[record.name].append(long)

out=open(output+"_new_coordinates.vcf","w")

print("\n")
print(dict)
print("\n")

with open(source,'r') as s:
	for line in s:
		if '#' in line:
			out.write(line)
		else:
			lys=line.split('\t')
			print("\nLigne du fichier source:\n")
			print(lys)
			#print(lys[0],lys[1],dict[lys[0]][0])
			if int(lys[1])<=int(dict[lys[0]][0]):
				#print("part1")
				out.write("{0}_part1\t{1}".format(lys[0],'\t'.join(lys[1:])))

			if int(lys[1])>int(dict[lys[0]][0]):
				#print("part2")
				out.write("{0}_part2\t{1}\t{2}".format(lys[0],str(int(lys[1])-int(dict[lys[0]][0])),'\t'.join(lys[2:])))

out.close()
