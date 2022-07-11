import f_reader

gff = '../WormBase/remanei/c_remanei.PRJNA577507.WS284.annotations.gff3.gz'

features = {}
for gff_seqid, feature_i in f_reader.read_gff(gff, "gene", False, "WormBase_imported"):
	for i in range(len(feature_i)):
		feat = feature_i[i]
		start, end = feat[0], feat[1]
		fwd = feat[2]

		# cleans up info field, assuming WormBase source
		gene = feat[3].replace("ID=gene:", "").split(";")[0]
		print(gene)
