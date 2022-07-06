import argparse
import f_reader


def rev_comp(seq):
	rc = ""
	for c in seq:
		if c == "A":
			rc += "T"
		elif c == "T":
			rc += "A"
		elif c == "C":
			rc += "G"
		else:
			rc += "C"
	return rc[::-1]


parser = argparse.ArgumentParser(description="find donor/acceptor splice sites")
parser.add_argument(
	"--fa", required=True, type=str, metavar="<str>", help="required string argument"
)
parser.add_argument(
	"--gff", required=True, type=str, metavar="<str>", help="required string argument"
)
arg = parser.parse_args()

features = {}
for gff_seqid, feature_i in f_reader.read_gff(arg.gff, "intron"):
	assert gff_seqid not in features.keys()
	features[gff_seqid] = feature_i

full_ss_dict = {}
for name, seq in f_reader.read_fasta(arg.fa):
	fa_seqid = name.split()[0]

	if fa_seqid in features.keys():

		feat_list = features[fa_seqid]  # retrieve features list from dictionary
		ss_dict = {}

		for i in range(len(feat_list)):
			feat = feat_list[i]
			start = feat[0]
			end = feat[1]
			fwd = feat[2]

			feat_seq = seq[start:end]

			if not fwd:
				# shortcut: skips finding rev_comp of the whole seq
				feat_seq = rev_comp(feat_seq)

			fprime = feat_seq[:2]
			tprime = feat_seq[-2:]

			ss = fprime + tprime

			if ss not in ss_dict:
				ss_dict[ss] = 0
			if ss not in full_ss_dict:
				full_ss_dict[ss] = 0

			ss_dict[ss] += 1
			full_ss_dict[ss] += 1

		print(fa_seqid, ss_dict, "\n")

full_ss_dict = {
	k: v for k, v in sorted(full_ss_dict.items(), key=lambda x: x[1], reverse=True)
}

total_ss = sum(full_ss_dict.values())
for seq, count in full_ss_dict.items():
	print(f"{seq}\t{count}\t{(count / total_ss):.4f}")
