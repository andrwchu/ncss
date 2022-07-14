import argparse
import f_reader
import bio_lib

# p3 gen_intron_exon.py --fa ../../datacore/genome_celegans/1pct_elegans.fa --gff ../../datacore/genome_celegans/1pct_elegans.gff3

def find_exons(exons, seqid, beg_in, end_in, gene_in):
	exon_list = exons[seqid]  # retrieve features list from dictionary

	beg_i, end_i = None, None
	for feat in exon_list:
		gene_ex = feat[3].replace("Parent=Transcript:", "").split(";")[0]
		if gene_in != gene_ex:
			continue

		beg_ex, end_ex = feat[0], feat[1]

		if end_ex == beg_in:
			beg_i = feat[0]
		if end_in == beg_ex:
			end_i = feat[1]
		if beg_i is not None and end_i is not None:
			return (beg_i, end_i)
	return (beg_i, end_i)

def ex_in_ex(fa, introns, exons):
	for name, seq in f_reader.read_fasta(fa):
		fa_seqid = name.split()[0]

		if fa_seqid in introns.keys():

			intr_list = introns[fa_seqid]  # retrieve features list from dictionary

			for feat in intr_list:
				start, end = feat[0], feat[1]
				fwd = feat[2]

				# cleans up info field, assuming WormBase source
				gene = feat[3].replace("Parent=Transcript:", "").split(";")[0]

				feat_seq = seq[start:end]
				if not fwd:
					# shortcut: skips finding rev_comp of the whole seq
					feat_seq = bio_lib.rev_comp(feat_seq)

				ss = feat_seq[:2] + feat_seq[-2:]

				if not ss == "GTAG":
					beg_i, end_i = find_exons(exons, fa_seqid, start, end, gene)

					trio = seq[beg_i: end_i]
					if not fwd:
						trio = bio_lib.rev_comp(trio)

					header = (
						f">{fa_seqid}_{beg_i}..{start}..{end}..{end_i}_{fwd}_{ss}_{gene}"
					)

					yield (header, trio)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="generate ex in ex file")
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
	fa = arg.fa
	gff = arg.gff
	introns = {}
	for gff_seqid, feature_i in f_reader.read_gff(gff, "intron", True, "WormBase"):
		introns[gff_seqid] = feature_i

	exons = {}
	for gff_seqid, feature_i in f_reader.read_gff(gff, "exon", True, "WormBase"):
		exons[gff_seqid] = feature_i

	for name, seq in ex_in_ex(fa, introns, exons):
		print(name)
		print(seq)
