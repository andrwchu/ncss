import argparse
import f_reader

# p3 feature_masker.py --fa ../../datacore/genome_celegans/1pct_elegans.fa --gff ../../datacore/genome_celegans/1pct_elegans.gff3 --feat exon -hm


def mask_feature(fa, gff, feat_type, hard):
	features = {}
	for gff_seqid, feature_i in f_reader.read_gff(gff, feat_type, False):
		assert gff_seqid not in features.keys()
		features[gff_seqid] = feature_i

	for name, seq in f_reader.read_fasta(fa):
		fa_seqid = name.split()[0]
		new_seq = seq
		if fa_seqid in features.keys():
			feat_list = features[fa_seqid]

			for i in range(len(feat_list)):
				feat = feat_list[i]
				start, end = feat[0], feat[1]
				fwd = feat[2]
				# assumes that features may overlap; if not, should use a linear replacement
				if hard:
					new_seq = new_seq[:start] + ("N" * (end - start)) + new_seq[end:]
				else:
					new_seq = (
						new_seq[:start] + new_seq[start:end].lower() + new_seq[end:]
					)

		header = ">" + name
		yield (header, new_seq)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="mask features")
	parser.add_argument(
		"--fa",
		required=True,
		type=str,
		metavar="<str>",
		help="required string argument",
	)
	parser.add_argument(
		"--gff",
		required=True,
		type=str,
		metavar="<str>",
		help="required string argument",
	)
	parser.add_argument(
		"--feat",
		required=True,
		type=str,
		metavar="<str>",
		help="required string argument",
	)
	parser.add_argument(
		"--hm", default=False, action="store_true", help="optional boolean flag"
	)
	arg = parser.parse_args()

	for header, seq in mask_feature(arg.fa, arg.gff, arg.feat, arg.hm):
		print(header)
		print(seq)
