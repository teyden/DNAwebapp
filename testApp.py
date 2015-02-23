# testApp.py

from DnaAnalyzer import *

def main():

	dna = DnaAnalyzer()
	dna.readFasta('cftr-partial.fasta')
	dna.reportOutput()
	
	# make test suite 


if __name__ == "__main__":
	main()

