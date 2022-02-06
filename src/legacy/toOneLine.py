#!/usr/bin/python
import argparse
import os

def args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--fa_file',default=None,help='path to fa file to generate chromosome separte oneline file')
	parser.add_argument('--species',default='hg38',help='path to fa file to generate chromosome separte oneline file')
	argv = parser.parse_args()
	fa_file = argv.fa_file
	species = argv.species
	return fa_file,species

def main(fa_file,species):
	fa = open(fa_file,'r')
	if not os.path.exists('oneLine'):
		os.makedirs('oneLine')
	skip = False
	for line in fa.readlines():
		line = line.rstrip('\n')
		if('>' in line):
			chro = line.split(' ')[0]
			chro = chro[1:]
			if('chr' not in chro):
				chro = 'chr'+chro

			ww = open('oneLine/%s.%s.fa'%(species,chro),'w')
			if(len(chro)>6 or 'M' in chro):
				skip = True
				os.system('rm oneLine/%s.%s.fa'%(species,chro))
				continue
			else:
				#ww.write('%s\n'%line)
				skip = False
		else:
			if(skip):
				continue
			ww.write(line)

	ww.close()

if __name__ == "__main__":
	main(*args())
