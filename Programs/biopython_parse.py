from Bio.Blast.Applications import NcbiblastnCommandline
import os


def gen_file_dir(directory, start, end):
	i = start
	while True:
		num_i = str(i).rjust(4, "0")
		celegans = f"{directory}{num_i}celegans.fa"
		if os.path.exists(celegans):
			with open(celegans) as f:
				elegans_id = f.readline().rstrip()

			j = 0
			while True:
				num_j = str(j).rjust(3, "0")
				orth_id = f"{num_i}orth_{num_j}"
				orth = f"{directory}{num_i}orth_{num_j}.fa"
				if os.path.exists(orth):
					yield (celegans, elegans_id, orth, orth_id)
				else:
					break
				j += 1
		else:
			break
		i += 1
		if i >= end:
			break

directory = "../out/bl2seq4/"
beg = 0
end = 1506
matches = []
for q, elegans_id, s, orth_id in gen_file_dir(directory, beg, end):
	print(orth_id)
	xml = directory + "xml/" + orth_id + ".xml"

	cline = NcbiblastnCommandline(
		query=q,
		subject=s,
		reward=1,
		penalty=-1,
		gapopen=3,
		gapextend=1,
		perc_identity=50,
		word_size=7,
		out=xml,
		outfmt=5,
	)()
	# calling cline() returns (stdout, stderr)
