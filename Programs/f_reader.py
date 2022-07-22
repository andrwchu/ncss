import gzip
import re
import sys

import bio_lib
import intron_class


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


def read_gff(filename, feat_type, species, check_source, source=""):
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
				species,
				seqid,
				int(field[3]) - 1,  # start, with offset from GFF file
				int(field[4]),  # end
				field[6] == "+",  # strand direction, fwd = True
				field[8],  # info
				None,
				None,
				None,
				None,
				None,
				[],
				[],
			)
			features.append(feat)

	yield (seqid, features)
	fp.close()

def read_exons(filename, introns):
	# introns --> genes --> exons/introns locations
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	gene_list = set()
	for intr in introns:
		gene_list.add(intr.gene)

	seqid = None
	features = {}
	for line in fp.readlines():
		if line[0] == "#":
			continue

		field = line.split()
		if not field[0] == seqid:
			seqid = field[0]

		if field[1] == "WormBase" and field[2] in ["exon", "intron"]:
			gene = field[8].replace("Parent=Transcript:", "").split(";")[0]
			if gene in gene_list:
				feat = (
					int(field[3]) - 1,  # start, with offset from GFF file
					int(field[4]),  # end
					field[6] == "+",  # strand direction, fwd = True,
					field[2] == "exon",  # exons = True, intron = False
				)
				if seqid not in features:
					features[seqid] = {}
				if gene not in features[seqid]:
					features[seqid][gene] = []
				features[seqid][gene].append(feat)
	fp.close()
	return features


def find_mask(nc_intr, FA, GFF):
	features = read_exons(GFF, nc_intr)

	masked = {}
	for header, seq in read_fasta(FA):
		seqid = header.split()[0]
		if seqid in features:
			for gene in features[seqid]:
				gene_seq = []
				for feat in features[seqid][gene]:
					beg, end, fwd, is_exon = feat
					feat_seq = seq[beg:end]
					if not is_exon:
						feat_seq = feat_seq[:2] + "nnnnn" + feat_seq[-2:]
						feat_seq = feat_seq.lower()

					if not fwd:
						feat_seq = bio_lib.rev_comp(feat_seq)
						gene_seq.insert(0, feat_seq)
					else:
						gene_seq.append(feat_seq)

				if gene not in masked:
					masked[gene] = ""
				masked[gene] = "".join(gene_seq)
				# print(gene, ''.join(gene_seq))

	return masked  # {gene : mask_seq}

def read_exons_orth(filename, gene_list):
	# introns --> genes --> exons/introns locations
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	seqid = None
	features = {}
	for line in fp.readlines():
		if line[0] == "#":
			continue

		field = line.split()
		if not field[0] == seqid:
			seqid = field[0]

		if field[1] == "WormBase" and field[2] in ["exon", "intron"]:
			gene = field[8].replace("Parent=Transcript:", "").split(";")[0]
			gene = gene.replace("Parent=Pseudogene:","")[:8]
			print(gene)
			if gene in gene_list:
				feat = (
					int(field[3]) - 1,  # start, with offset from GFF file
					int(field[4]),  # end
					field[6] == "+",  # strand direction, fwd = True,
					field[2] == "exon",  # exons = True, intron = False
				)
				if seqid not in features:
					features[seqid] = {}
				if gene not in features[seqid]:
					features[seqid][gene] = []
				features[seqid][gene].append(feat)
	fp.close()
	return features


def find_mask_orth(gene_list, FA, GFF):
	features = read_exons_orth(GFF, gene_list)

	masked = {}
	for header, seq in read_fasta(FA):
		seqid = header.split()[0]
		if seqid in features:
			for gene in features[seqid]:
				gene_seq = []
				for feat in features[seqid][gene]:
					beg, end, fwd, is_exon = feat
					feat_seq = seq[beg:end]
					if not is_exon:
						feat_seq = feat_seq[:2] + "nnnnn" + feat_seq[-2:]
						feat_seq = feat_seq.lower()

					if not fwd:
						feat_seq = bio_lib.rev_comp(feat_seq)
						gene_seq.insert(0, feat_seq)
					else:
						gene_seq.append(feat_seq)

				if gene not in masked:
					masked[gene] = ""
				masked[gene] = "".join(gene_seq)
				# print(gene, ''.join(gene_seq))

	return masked  # {gene : mask_seq}


def combine_orth(intr_orths):
	# combine genes from the same species (list --> dict)
	orth_dict = {}
	for orth in intr_orths:
		species = orth[0]
		gene = orth[1]
		if species not in orth_dict:
			orth_dict[species] = []
		orth_dict[species].append(gene)

	intr_orths = orth_dict
	if len(orth_dict.keys()) == 0:
		intr_orths = None

	return intr_orths


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
	fp.close()

	for intr in introns:
		intr.orth_genes = combine_orth(intr.orth_genes)
	return introns


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


def wb_to_seqgene(filename, introns, target):
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	# gather a set of all seqgenes under target species
	wbgenes = set()
	for intr in introns:
		if intr.orth_genes is not None:
			if target in intr.orth_genes:
				wbgenes.update(intr.orth_genes[target])  # add list to set

	genes = {}
	for line in fp.readlines():
		field = line.split(",")

		wb, public, seq_name = field[1], field[2], field[3]
		if wb in wbgenes:
			genes[wb] = seq_name

	fp.close()

	# wb_genes returns genes in shuffled order, not corresponding to seq_gene order
	return genes
