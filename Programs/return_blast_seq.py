import f_reader


def read_blast(filename, seqid):
	fp = None
	if filename == "-":
		fp = sys.stdin
	elif filename.endswith(".gz"):
		fp = gzip.open(filename, "rt")
	else:
		fp = open(filename)

	matches = []
	for line in fp.readlines():
		if line[0] == "#":
			continue

		field = line.split()
		# Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

		qID, sID, beg_s, end_s, e = (
			field[0],
			field[1],
			int(field[8]),
			int(field[9]),
			field[10],
		)
		if sID == seqid:
			matches.append((qID, sID, beg_s, end_s, e))
	fp.close()
	return matches


fa = "../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.genomic.fa.gz"
blast = "../WormBase/briggsae/001briggsae"

for name, seq in f_reader.read_fasta(fa):
	seqid = name.split()[0]
	print(seqid, len(seq))
	for match in read_blast(blast, seqid):
		print(match)
		qID, sID, beg_s, end_s, e = match
		print(seq[beg_s:end_s])
