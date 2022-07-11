import argparse
import f_reader
import bio_lib
import extract_intron

# p3 find_proteins.py --fa ../../datacore/genome_celegans/1pct_elegans.fa --gff ../../datacore/genome_celegans/1pct_elegans.gff3


def extract_cds(gff, genes_list):
	cds_dict = {}
	for seqid, feature_i in f_reader.read_gff(gff, "CDS", True, "WormBase"):
		for feat in feature_i:

			genes = feat[3]
			genes = genes.replace("Parent=", "").split(";")[1]
			genes = genes.replace("Transcript:", "").split(",")

			for gene in genes:
				if gene in genes_list:
					start, end = feat[0], feat[1]
					fwd = feat[2]
					if seqid not in cds_dict:
						cds_dict[seqid] = {}
					if gene not in cds_dict[seqid]:
						cds_dict[seqid][gene] = []
					cds_dict[seqid][gene].append((start, end, fwd))
	return cds_dict  # {I : {gene : (1, 10, True), (11, 20, True)}}


def cds_to_protein(fa, cds_dict):
	for name, seq in f_reader.read_fasta(fa):
		fa_seqid = name.split()[0]

		if fa_seqid in cds_dict.keys():

			cds_list = cds_dict[fa_seqid]

			for gene in cds_list.keys():
				indices = cds_list[gene]
				cds_seq = ""

				for field in indices:
					start, end, fwd = field[0], field[1], field[2]

					feat_seq = seq[start:end]
					if not fwd:
						feat_seq = bio_lib.rev_comp(feat_seq)
						cds_seq = feat_seq + cds_seq
					else:
						cds_seq += feat_seq

				yield (fa_seqid, gene, bio_lib.translate(cds_seq))


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="find proteins")
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

	genes_list = []
	for fields in extract_intron.extract_intron(arg.fa, arg.gff, 4, 5):
		genes_list.append(fields[0])

	cds = extract_cds(arg.gff, genes_list)
	for protein in cds_to_protein(arg.fa, cds):
		print(f">{protein[0]} {protein[1]}")
		print(protein[2])
