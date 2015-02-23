if (typeof String.prototype.startsWith != 'function') {
  // see below for better implementation!
  String.prototype.startsWith = function (str){
    return this.indexOf(str) == 0;
  };
}

var addictions = {"rs1800497": {"description": "Addictive behaviours",
								"chromosome": 11,
								"gene": "ANKK1",
								"alleles": "AA, AG, GG",
								"tested-by": "A"}};


function execute() {
	var gene = document.getElementById('gene-input').value;
	var lines = gene.split(/\r?\n/);
	var collect = {};

	lines.forEach(function (line) {   // lines is an array, line is string element
		if (!line.startsWith("#")) {
			var lineArr = line.split(/\s+/);
			var rsid = lineArr[0],
				chromosome = lineArr[1],
				position = lineArr[2],
				genotype = lineArr[3];
			if (addictions[rsid] != null) {
				alert('MATCH! \n\nSNP ID: ' + addictions[rsid].description);
			} else {
				alert('Not matched \n\nSNP ID: ' + rsid);
				 
			}
		}	
	});
}
