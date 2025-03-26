#!/prg/python/3.5/bin/python3.5
# -*- coding: iso-8859-1 -*-

import os, sys
import sys
import vcf

'''Version1.0'''
'''Fonctions permettant d'analyser des résultats de GBS dans des fichiers vcf '''


def convert_vcf(results):
	"""Permet de convertir le fichier vcf en tableau texte.
	N'applique aucun filtre sur les entrees (eg, il y a plusieurs alleles alternatifs).
	Si c'est necessaire, utiliser la fonction readvcf_plus.
	"""
	myvcf=vcf.Reader(open(results,'r'))
	nom=os.path.splitext(os.path.basename(results))
	name=nom[0]
	ext=nom[1]
	out=open(name+'_table.txt','w')
	
	list_par=['CHROM','POS','REF','ALT']
	
	with open('vcf_header','r') as k:
		for line in k:
			lis=line.split()
	list_par+=lis[9:]

	out.write('\t'.join(list_par)+'\n')
	
	for x in myvcf:
		liste1=[]
		liste1.append(x.CHROM),liste1.append(str(x.POS)),liste1.append(x.REF),liste1.append(str(x.ALT))
		for sample in x.samples:
			if sample['GT'] == None:
				liste1.append('-')
			else:
				liste1.append(sample['GT'])
#If a SNP has no genotype, null for all samples, don't print that SNP to output file
#minus 4 is for the first four values of the list that are for CHROM, POS, REF and ALT
		if liste1.count('-') ==len(liste1)-4:
			pass
		else:
			out.write('\t'.join(liste1)+'\n')
			
	out.close()

if __name__ == '__main__':
	#en sys.argv[1] le nom du fichier des positions de genes
	convert_vcf(sys.argv[1])
