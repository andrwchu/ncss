import argparse
import f_reader
import bio_lib

# p3 extract_intron.py --fa ../../datacore/genome_celegans/1pct_elegans.fa --gff ../../datacore/genome_celegans/1pct_elegans.gff3

def extract_intron(fa, gff, field_start, field_end):
	features = {}
	for gff_seqid, feature_i in f_reader.read_gff(gff, "intron", True, "WormBase"):
		assert gff_seqid not in features.keys()
		features[gff_seqid] = feature_i

	for name, seq in f_reader.read_fasta(fa):
		fa_seqid = name.split()[0]

		if fa_seqid in features.keys():

			feat_list = features[fa_seqid]  # retrieve features list from dictionary

			for i in range(len(feat_list)):
				feat = feat_list[i]
				start, end = feat[0], feat[1]
				fwd = feat[2]

				# cleans up info field, assuming WormBase source
				gene = feat[3].replace("Parent=Transcript:", "").split(";")[0]

				feat_seq = seq[start:end]
				if not fwd:
					# shortcut: skips finding rev_comp of the whole seq
					feat_seq = bio_lib.rev_comp(feat_seq)

				fprime = feat_seq[:2]
				tprime = feat_seq[-2:]

				ss = fprime + tprime

				if not ss == "GTAG":
					yield (fa_seqid, start, end, ss, gene, fwd, feat_seq)[
						field_start:field_end
					]


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="find donor/acceptor splice sites")
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
	arg = parser.parse_args()

	for fields in extract_intron(arg.fa, arg.gff, 0, 7):
		print(fields)
