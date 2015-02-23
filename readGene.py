def readGene(myFile):
	g=0
	c=0
	a=0
	t=0

	gene = open(myFile, "r")
	gene.readline()

	for line in gene:
		line = line.lower()
		for char in line:
			if char == "g":
				g+=1
			if char == "a":
				a+=1
			if char == "c":
				c+=1
			if char == "t":
				t+=1

	print "Nucleotide count: \n"
	
	count = [g,c,a,t]
	nts = ["G", "C", "A", "T"]
	for i in range(4):
		print "%s: " %(nts[i]) + str(count[i])


	gc = (g+c+0.) / (a+t+c+g+0.)
	at = (a+t+0.) / (a+t+c+g+0.)
	print "GC content: " + str(gc*100) + "%"
	print "AT content: " + str(at*100) + "%"  




