# DnaAnalyzer.py

class DnaAnalyzer:

	""" 
	sequenceAnalyzer class analyzes # of nucleotides, codons, amino acids
	and GC content of our DNA

	self.codon
	self.aa

	"""

	def __init__(self):
		
		self.codonTable = 
		{
		"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
	    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
	    "UAU":"Y", "UAC":"Y", "UAA":"-", "UAG":"-",
	    "UGU":"C", "UGC":"C", "UGA":"-", "UGG":"W",

	    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
	    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",

	    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
	    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
	    
	    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
	    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
		}

		# initalize codon and AA dicts
		self.codon = {}
		self.aa = {}

		# dict contains nucleotide count
		self.nuc = { "A":0, "T":0, "G":0, "C":0, "U":0, "N":0 }  

		# add keys to codon and aa, and initialize to zero
		for key in self.codonTable:
			self.codon[key] = 0
			self.aa[self.codonTable[key]] = 0

		# list holding headers
		self.header = []
		self.ps = "%"
		self._totalBase = 0
		self._totalMB = 0.0


	def readFasta(self, fp):
		"""
		readFasta reads in fasta file, calls parseFasta
		then calls analyzeSequence to analyze DNA
		"""

		# reads in fasta
		for head, seq in self.parseFasta(fp):
			# analyzes sequence
			self.analyzeSequence(seq)
			if head == "":
				continue 
			else:
				# saves header
				self.header.append(head)


	def parseFasta(self, fp):
		"""
		parseFasta reads fasta from file, takes filepath as argument;
		outputs header and sequence
		"""

		fh = open(fp, 'r') 

		header = ""
		sequence = ""

		for line in fh:

			# hit header
			if line.startswith(">"):

				# return the header, seq
				yield header, sequence

				# grab new heder
				header = line.replace(">", "")

				# reset our sequence string
				sequence = ""

				# don't add the dna this round
				continue 

			sequence += line.strip()

			# return the final header and sequence
			yield header, sequence


	def convertDNAtoRNA(self, seq):
		""" 
		convertDNAroRNA... also parses the DNA
		"""

		temp = ""

		# take out all non dna / rna characters
		for base in seq:
			if base in self.nuc:
				temp += base

		seq = temp

		# convert dna to rna
		seq = seq.replace('T', 'U')

		return(seq)



	def analyzeSequence(self, seq):
		"""
		analyzeSequence increments the respective dictionaries 
		"""

		# convert seq to uppercase letters
		seq = seq.upper()

		# increment the bases to our nucleotide dictionary
		for base in seq:
			if base in self.nuc:
				self.nuc[base] += 1
			print("%s: %d " % (base, self.nuc[base]))

		# convert to RNA
		rna = self.convertDNAtoRNA(seq)

		# add codons to the codon dictionary
		for i in range(0, len(rna), 3):
			if rna[i:i+3] in self.codonTable:
				self.codon[rna[i:i+3]] += 1
				self.aa[self.codonTable[rna[i:i+3]]] += 1    # key = rna[i:i+3]


	def gcContent(self):

		gc = (float(self.nuc["G"] + self.nuc["C"]) / float(self.nuc["G"] + self.nuc["C"] + self.nuc["T"] + self.nuc["A"]))
		gc *= 100
		return gc


	def getSeqSize(self):
		"""
		getSeqSize calculates total # bases 
		and total # MB. call getSeqSize before reportOutput
		to get up-to-date size 

		self._totalBase is 0 upon call to getSeqSize()
		"""
		self._totalBase = 0
		for base in self.nuc:
			self._totalBase += self.nuc[base]

		self._totalMB = float(self._totalBase) / 1000000.00 


	def totalAA(self):
		 
		totalAA = 0

		for aa in self.aa:
			totalAA += self.aa[aa]

		return totalAA


	def reportAA(self):

		keys = list(self.aa.keys())
		keys.sort()

		print("Amino Acid Composition")

		for key in keys: 
			print("%s : %3.2f%s %d " % (key, 
			(float(self.aa[key])/self.totalAA()*100), 
			self.ps, self.aa[key]))

		print()


	def reportCodon(self):

		keys = list(self.codon.keys())
		keys.sort()

		print("Codon Usage")

		for key in keys: 
			print("%s : %s %3.2f%s %d " % (key, self.codonTable[key],
			(float(self.codon[key])/self.totalAA()*100), self.ps, 
			self.codon[key] ))
		
		print()


	def reportOutput(self):

		print (self.header[1]), 
		self.getSeqSize()
		print("Total Bases = %d " % self._totalBase)
		for base in self.nuc:
			per_yield = float(self.nuc[base] / float(self._totalBase))*100
			print("	%s: %d nts | %.2f%s " % (base, self.nuc[base], per_yield, self.ps))
		print("MB = %.6f " % self._totalMB)
		print("GC content %.2f%s " % (self.gcContent(), self.ps))
		print("*********************")

		self.reportAA()
		self.reportCodon()




