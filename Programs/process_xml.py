from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os
import re


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


def find_match(qry, sub, e_ss):
	indices = set()
	regex = r"([\w]{2})(NNNNN)([\w]{2})"
	matches_q = re.finditer(regex, qry)
	for match in matches_q:
		indices.add(match.span())

	pairs = set()
	intron_skip = 0
	for i, j in indices:
		if (
			qry[i + 2 : j - 2] == sub[i + 2 : j - 2] and qry[i + 2 : j - 2] == "NNNNN"
		):  # check for intron alignment
			ss1 = qry[i : i + 2] + qry[j - 2 : j]
			ss2 = sub[i : i + 2] + sub[j - 2 : j]
			if e_ss == ss1:
				pairs.add((ss1, ss2))
		elif qry[i + 2 : j - 2] == "NNNNN" and "-" in sub[i + 2 : j - 2]:
			intron_skip += 1
	return pairs, intron_skip, len(indices)


def fuzzy(eleg, orth):
	counter = 0
	for i in range(len(eleg)):
		if eleg[i] == orth[i]:
			counter += 1
	return counter


def stats(matches):
	pair_dict = {}
	for match in matches:
		pair = match[0]

		if pair not in pair_dict:
			pair_dict[pair] = 0
		pair_dict[pair] += 1
	return pair_dict


directory = "../out/bl2seq3/"
beg = 0
end = 1491
matches = []
orth_spec = {
	"Caenorhabditis_briggsae": 0,
	"Caenorhabditis_brenneri": 0,
	"Caenorhabditis_remanei": 0,
}
skips = {
	"Caenorhabditis_briggsae": [0, 0],
	"Caenorhabditis_brenneri": [0, 0],
	"Caenorhabditis_remanei": [0, 0],
}
orth_matches = {}
for q, elegans_id, s, orth_id in gen_file_dir(directory, beg, end):
	# print(orth_id)
	xml = directory + "xml/" + orth_id + ".xml"

	result_handle = open(xml)
	blast_records = NCBIXML.parse(result_handle)

	E_VALUE_THRESH = 0.01
	SEQ_DISPLAY = 75
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					(
						e_num,
						e_spec,
						e_seqid,
						e_beg,
						e_end,
						e_gene,
						e_ss,
					) = elegans_id.split(",")

					(o_num, o_spec, o_gene,) = alignment.hit_id.split(",")
					orth_spec[
						o_spec
					] += 1  # count up total successful aligns for each species

					pairs, intron_skip, tot_intron = find_match(
						hsp.query, hsp.sbjct, e_ss
					)
					skips[o_spec][0] += intron_skip
					skips[o_spec][1] += tot_intron

					if len(pairs) > 0:
						# print(pairs)
						for pair in pairs:
							score = fuzzy(pair[0], pair[1])
							matches.append(((pair[0], pair[1]), score, o_spec))



						if pair[0] == pair[1] and pair[0] in (
							"GCAG",
							"GTGG",
							"ATAA",
							"GGAG",
							"GTAT",
							"ATAC",
						):
							if o_num[:4] not in orth_matches:
								orth_matches[o_num[:4]] = []
							orth_matches[o_num[:4]].append((o_num, o_spec, pair[0]))
	result_handle.close()

for i in orth_matches:
	matches = orth_matches[i]
	if len(matches) > 2:
		print(matches)
"""
print("overall")
pairs = stats(matches)
pairs = {k: v for k, v in sorted(pairs.items(), key=lambda x: x[1], reverse=True)}

print(f"total alignment attempts: ", sum(orth_spec.values()))

for key in pairs:
	print(key, pairs[key])

for i in [
	"Caenorhabditis_briggsae",
	"Caenorhabditis_brenneri",
	"Caenorhabditis_remanei",
]:
	print(i)
	print("total alignment attempts: ", orth_spec[i])
	print(f"intron skips: {(skips[i][0] / skips[i][1]):.4f}")
	spec = []
	for j in matches:
		if j[-1] == i:
			spec.append(j)
	pairs = stats(spec)
	pairs = {k: v for k, v in sorted(pairs.items(), key=lambda x: x[1], reverse=True)}

	for key in pairs:
		print(f"{key}\t{pairs[key]}")
"""
