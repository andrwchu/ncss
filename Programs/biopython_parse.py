from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os
import re


def gen_file_dir(directory):
	i = 0
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
		if i > 5:
			break

		i += 1


def find_match(qry, sub):
	indices = set()
	regex = r"([\w]{2})(NNNNN)([\w]{2})"
	matches_q = re.finditer(regex, qry)
	for match in matches_q:
		indices.add(match.span())

	matches_s = re.finditer(regex, sub)
	for match in matches_s:
		indices.add(match.span())

	pairs = set()
	for i, j in indices:
		ss1 = qry[i : i + 2] + qry[j - 2 : j]
		ss2 = sub[i : i + 2] + sub[j - 2 : j]
		pairs.add((ss1, ss2))
	print(pairs)


directory = "../out/bl2seq2/"
for q, elegans_id, s, orth_id in gen_file_dir(directory):
	xml = directory + "xml/" + orth_id + ".xml"

	cline = NcbiblastnCommandline(
		query=q,
		subject=s,
		reward=1,
		penalty=-1,
		gapopen=3,
		gapextend=1,
		out=xml,
		outfmt=5,
	)()

	# calling cline() returns (stdout, stderr)
	# print(cline)

	result_handle = open(xml)
	blast_records = NCBIXML.parse(result_handle)

	E_VALUE_THRESH = 1
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
					(
						o_num,
						o_spec,
						o_seqid,
						o_beg,
						o_end,
						o_gene,
						o_ss,
					) = alignment.hit_id.split(",")

					em_beg = int(e_beg) + hsp.query_start - 1
					em_end = int(e_beg) + hsp.query_end
					om_beg = int(o_beg) + hsp.sbjct_start - 1
					om_end = int(o_beg) + hsp.sbjct_end

					print("q: ", elegans_id)
					print("s:  ", alignment.hit_id)
					print("length:", hsp.align_length)
					print("e value:", hsp.expect)
					print("bits:", hsp.bits)

					print(f"{em_beg}..{em_end}")
					print(f"{om_beg}..{om_end}")

					print(hsp.query[0:SEQ_DISPLAY] + "...")
					print(hsp.match[0:SEQ_DISPLAY] + "...")
					print(hsp.sbjct[0:SEQ_DISPLAY] + "...")

					find_match(hsp.query, hsp.sbjct)
					print()

	result_handle.close()
