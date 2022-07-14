import sys
import gzip
import intron_class
import re


def read_fasta(filename):

	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	name = None
	seqs = []
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith(">"):
			if len(seqs) > 0:
				seq = "".join(seqs)
				yield (name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield (name, "".join(seqs))
	fp.close()


def read_gff(filename, feat_type, check_source, source=""):
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	seqid = None
	features = []
	for line in fp.readlines():
		if line[0] == "#":
			continue

		field = line.split()
		if not field[0] == seqid:
			if len(features) > 0:
				yield (seqid, features)
				seqid = field[0]
				features = []
			else:
				seqid = field[0]

		correct_source = True

		if check_source and source != field[1]:
			correct_source = False

		if correct_source and field[2] == feat_type:
			feat = intron_class.Intron(
				seqid,
				int(field[3]) - 1,  # start, with offset from GFF file
				int(field[4]),  # end
				field[6] == "+",  # strand direction, fwd = True
				field[8],  # info
				None,
				None,
				None,
				None,
				[],
			)
			features.append(feat)

	yield (seqid, features)
	fp.close()


def read_ortholog(filename, introns):
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	wbgenes = set()
	for intr in introns:
		wbgenes.add(intr.wbgene)

	name = None
	new_section = False
	for line in fp:
		if line[0] == "#":
			continue

		field = line.split()
		if new_section:
			name = field[0]
			new_section = False
		elif field[0] == "=":
			new_section = True

		else:
			genus, species = field[0], field[1]
			orth, public = field[2], field[3]
			if orth[:6] == "WBGene":  # restricting to only WormBase (?)
				orth_info = [genus + "_" + species, orth]
				if name in wbgenes:
					for intr in introns:
						if name == intr.wbgene:
							intr.orth_genes.append(orth_info)

	return introns
	fp.close()


def seq_to_wbgene(filename, introns):

	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	# this function passes the introns list in so that the wbgene stays associated with the the intron
	# another way to to do this would be like with the wb_to_seqgene function
	# which only exchanges the gene names through a dictionary
	intr_genes = set()
	for intr in introns:
		gene = re.search("([\w]+\.[\d]+)", intr.gene).group(1)
		intr_genes.add(gene)
		intr.clean_gene = gene

	for line in fp.readlines():
		field = line.split(",")
		# sequence name is the identifier from the sequence
		wb, public, seq_name = field[1], field[2], field[3]
		if seq_name in intr_genes:
			for intr in introns:
				if seq_name == intr.clean_gene:
					intr.wbgene = wb

	return introns
	fp.close()

	# wb_genes returns genes in shuffled order, not corresponding to seq_gene order
	return wb_list


def wb_to_seqgene(filename, wbgenes):
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	genes = {}
	for line in fp.readlines():
		field = line.split(",")

		wb, public, seq_name = field[1], field[2], field[3]
		if wb in wbgenes:
			genes[wb] = [seq_name, []]  # array so that other attributes can be added later

	fp.close()

	# wb_genes returns genes in shuffled order, not corresponding to seq_gene order
	return genes
