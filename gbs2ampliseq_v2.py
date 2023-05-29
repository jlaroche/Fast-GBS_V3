#!/prg/python/3.5/bin/python3.5
# -*- coding: iso-8859-1 -*-

"""
Pour identifier les variants contenu dans les sequences ampliseq

La commande suivante:
./gbs2ampliseq_v2.py FastGBS_platypus.vcf 27marker_GTA_92.txt

Va produire le fichier suivant:
FastGBS_platypus_table.txt

"""

#./gbs2ampliseq_v2.py FastGBS_platypus.vcf Source_AlleleName.txt

import os, sys, fileinput, vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

#Fichier vcf
try:
	fichier = sys.argv[1]
except:
	print(__doc__)
	sys.exit(1)

#Fichier ampliseq
try:
	ampseq = sys.argv[2]
except:
	print(__doc__)
	sys.exit(1)

# Lecture du fichier vcf avec la librairie vcf
myvcf=vcf.Reader(open(fichier,'r'))

# Fichier final de resultats
nom=os.path.splitext(os.path.basename(fichier))
name=nom[0]
ext=nom[1]
outxv=open(name+'_allele_names.txt','w')

outsi=open(name+'_snp_interet.vcf','w')

# Creation de la liste des noms de noms de colonne
list_par=['CHROM','POS','REF','ALT','REF-S','ALT-S']
#list_par=['CHROM','POS','REF','ALT']

print("1. Capture des noms d'échantillons\n")
with open(fichier,'r') as i:
	for line in i:
		if '#CHROM' in line:
			lis=line.split()
			list_par+=lis[9:]
			break
		else:
			pass

#print(list_par)

# Ecriture de l'entete des noms de colonne dans le fichier de sortie
outxv.write('\t'.join(list_par)+'\n')

print("2. Creation du dictionnaire a partir du fichier de marqueurs ampliseq a traiter\n")
# Ce sont ces marqueurs qu'il faut rechercher la presence dans le fichier VCF
dic2={}
with open(ampseq,'r') as f:
	for line in f:
		if "#" not in line:
			lys=line.split()
#			print(lys)

# j'ai ajouté les alleles dans les clés du dictionnaire afin de
# distinguer les marqueurs ont même positions mais qui ont des génotypes différents
			if lys[0]+"_"+lys[1] not in dic2:
				dic2[lys[0]+"_"+lys[1]]={}

				dic2[lys[0]+"_"+lys[1]][lys[2]+"_"+lys[3]]=[]
				
				dic2[lys[0]+"_"+lys[1]][lys[2]+"_"+lys[3]].append(lys[4])
				dic2[lys[0]+"_"+lys[1]][lys[2]+"_"+lys[3]].append(lys[5])

			else:
				dic2[lys[0]+"_"+lys[1]][lys[2]+"_"+lys[3]]=[]
	
				dic2[lys[0]+"_"+lys[1]][lys[2]+"_"+lys[3]].append(lys[4])
				dic2[lys[0]+"_"+lys[1]][lys[2]+"_"+lys[3]].append(lys[5])

		else:
			pass
#for x in dic2:
#	print(x,dic2[x])


print("3. Lecture du fichier VCF et remplissage du fichier avec les noms d'allèles\n")
for x in myvcf:
	liste1=[]
	if x.CHROM+"_"+str(x.POS) in dic2:
		print(x.CHROM,x.POS)
		liste1.append(x.CHROM)
		liste1.append(str(x.POS))
		liste1.append(x.REF)
		liste1.append(str(x.ALT).replace('[','').replace(']',''))

		if x.REF +"_"+ str(x.ALT).replace('[','').replace(']','') in dic2[x.CHROM+"_"+str(x.POS)]:

			liste1.append('=')
			liste1.append('=')

			for sample in x.samples:
				if sample['GT'] == './.':
					liste1.append('MIS')
				elif sample['GT'] == '0/0':
					liste1.append(dic2[x.CHROM+"_"+str(x.POS)][x.REF +"_"+ str(x.ALT).replace('[','').replace(']','')][0])
					#liste1.append('A')	
				elif sample['GT'] == '1/1':
					liste1.append(dic2[x.CHROM+"_"+str(x.POS)][x.REF +"_"+ str(x.ALT).replace('[','').replace(']','')][1])
					#liste1.append('B')
				elif sample['GT'] == '0/1':
					liste1.append('HET')
				elif sample['GT'] == '1/2':
					liste1.append('HET')
				elif sample['GT'] == '2/1':
					liste1.append('HET')
				elif sample['GT'] == '1/0':
					liste1.append('HET')			
				else:
					pass

			liste1.append('///')
			for x in dic2[x.CHROM+"_"+str(x.POS)]:
				liste1.append(x.split('_')[0])
				liste1.append(x.split('_')[1])

		else:

			liste1.append('<>')
			liste1.append('<>')

			for sample in x.samples:
				if sample['GT'] == './.':
					liste1.append('MIS')
				elif sample['GT'] == '0/0':
					for y in dic2[x.CHROM+"_"+str(x.POS)]:
						liste1.append(y.split('_')[0])
				elif sample['GT'] == '1/1':
					for y in dic2[x.CHROM+"_"+str(x.POS)]:
						liste1.append(y.split('_')[1])
				elif sample['GT'] == '0/1':
					liste1.append('HET')
				elif sample['GT'] == '1/2':
					liste1.append('HET')
				elif sample['GT'] == '2/1':
					liste1.append('HET')
				elif sample['GT'] == '1/0':
					liste1.append('HET')			
				else:
					pass
					
			liste1.append('///')
			for x in dic2[x.CHROM+"_"+str(x.POS)]:
				liste1.append(x.split('_')[0])
				liste1.append(x.split('_')[1])

		outxv.write('\t'.join(liste1)+'\n')
	else:
		pass

print("4. Création du fichier VCF avec seulement les SNP d'intérêt\n")		
with open(fichier,'r') as i:
	for line in i:
		if '#' in line:
			outsi.write(line)
		else:
			lys=line.split()
			if lys[0]+"_"+lys[1] in dic2:
				outsi.write(line)
			else:
				pass

outxv.close()
outsi.close()
