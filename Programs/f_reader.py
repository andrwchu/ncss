import sys
import gzip


def read_fasta(filename):
	name = None
	seqs = []

	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

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


def read_gff(filename, feat_type):
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
		field = line.split()
		if not field[0] == seqid:
			if len(features) > 0:
				yield (seqid, features)
				seqid = field[0]
				features = []
			else:
				seqid = field[0]

		if field[2] == feat_type:
			feat_info = [
				int(field[3]) - 1,  # start, with offset from GFF file
				int(field[4]),  # end
					field[6] == "+",  # strand direction, fwd = True
					field[5],
					field[7],
			]
			features.append(feat_info)

	yield (seqid, features)
	fp.close()
