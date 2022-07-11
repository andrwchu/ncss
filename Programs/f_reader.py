import sys
import gzip
from extract_intron import Feature

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
			feat_info = [
				int(field[3]) - 1,  # start, with offset from GFF file
				int(field[4]),  # end
				field[6] == "+",  # strand direction, fwd = True
				field[8], # info
				None,
				None,
				None
			)
			features.append(feat)

	yield (seqid, features)
	fp.close()


def read_ortholog(filename):
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	name = None
	orthologs = []
	new_section = False
	for line in fp:
		if line[0] == "#":
			continue

		field = line.split()
		if new_section:
			name = field[0]
			orthologs = []
			new_section = False
		elif field[0] == "=":
			new_section = True
			if len(orthologs) > 0:
				yield (name, orthologs)
		else:
			genus, species = field[0], field[1]
			orth, public = field[2], field[3]
			if orth[:6] == "WBGene":  # restricting to only WormBase (?)
				orth_info = [genus + "_" + species, orth]
				orthologs.append(orth_info)

	yield (name, orthologs)
	fp.close()


def seq_to_wbgene(filename, seqs):

	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	wb_list = []
	for line in fp.readlines():
		field = line.split(",")

		# sequence name is the identifier from the sequence
		wb, public, seq_name = field[1], field[2], field[3]
		if seq_name in seqs:
			wb_list.append((seq_name, wb))

	fp.close()

	# wb_genes returns genes in shuffled order, not corresponding to seq_gene order
	return wb_list
