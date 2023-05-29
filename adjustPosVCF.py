#!/prg/python/3.5/bin/python3.5
# -*- coding: iso-8859-1 -*-

"""
Permet d'ajuster les positions des variants obtenus sur un genome de reference scindé
pour les ramener sur le genome de reference original.

Utilisation:

adjustPosVCF.py <fichier.vcf> <genome>

<fichier.vcf>	: Le fichier vcf a corriger
<genome>	: Le fichier du genome de reference avec les chromosomes scindés en format fasta

exemple:
./adjustPosVCF.py CropSNP_QOSM_PASS_indel_biSNP_Q30.vcf Hordeum_split.fasta

Le resultat est envoyé dans un fichier avec la terminaison de nom <_new_coordinates.vcf>
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

try:
        vcf = sys.argv[1]
except:
        print(__doc__)
        sys.exit(1)

## Le fichier contenant les chromosomes scindés
try:
        genome = sys.argv[2]
except:
        print(__doc__)
        sys.exit(1)

dot = vcf.index(".")
ext = vcf[dot+1:]
output = vcf[:dot]
out=open(output+"_original_coordinates.vcf","w")

dict={}
for record in SeqIO.parse(genome, "fasta"):
	dict[record.name]=len(record.seq)

#print(dict)
"""
{'Pt': 115974, 'chr6H_part2': 291690273, 'chr4H_part2': 323530078, 'chr2H_part1': 384037500,
 'chr3H_part1': 349855560, 'chr6H_part1': 291690240, 'chr2H_part2': 384037524, 
 'chr5H_part1': 335015100, 'chr1H_part2': 279267712, 'chr7H_part2': 328611980, 
 'chr4H_part1': 323530080, 'chrUn': 249774706, 'chr7H_part1': 328612020, 
 'chr1H_part1': 279267720, 'chr3H_part2': 349855554, 'chr5H_part2': 335015060}
"""

with open(vcf,'r') as f:
	for line in f:
		if '#' in line:
			out.write(line)
		else:
			listb = line.split('\t')
#			print(listb)
			if "part1" in listb[0]:
			# conserve la meme position
				chr=listb[0].split('_part1')[0]
				pos=listb[1]
				listb[0]=chr
			# pour changer le champs ID
				listb[2]=chr+':'+pos
				out.write('\t'.join(listb))
			elif "part2" in listb[0]:
			# on ajoute la longeur de part1 a la position des SNPs situes dans part2
				chr=listb[0].split('_part2')[0]
				part=listb[0].split('_part2')[1]
				pos=listb[1]
				listb[0]=chr
				listb[1]=str(dict[chr+'_part1']+int(pos))
			# pour changer le champs ID
				listb[2]=chr+':'+str(dict[chr+'_part1']+int(pos))
				out.write('\t'.join(listb))
			else:
			# pour les chromosomes qui n'ont pas ete scindes
				chr=listb[0]
				pos=listb[1]
				listb[0]=chr
			# pour changer le champs ID
				listb[2]=chr+':'+pos
				out.write('\t'.join(listb))
out.close()
