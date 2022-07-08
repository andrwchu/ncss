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
		wb_list, public, seq_name = field[1], field[2], field[3]
		if seq_name in seqs:
			wb_list.append(wb)
	fp.close()

	return wb_list
