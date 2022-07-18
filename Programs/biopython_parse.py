from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os




def gen_file_dir(directory):
	i = 0
	while True:
		num_i = str(i).rjust(3, "0")
		celegans = f"{directory}{num_i}celegans.fa"
		if os.path.exists(celegans):

			j = 0
			while True:
				num_j = str(j).rjust(2, "0")
				orth_id = f"{num_i}orth_{num_j}"
				orth = f"{directory}{num_i}orth_{num_j}.fa"
				if os.path.exists(orth):
					yield (celegans, orth, orth_id)
				else:
					break
				j += 1
		else:
			break
		i += 1


directory = "../out/bl2seq/"
for q, s, orth_id in gen_file_dir(directory):
	print(orth_id)

	xml = directory + "xml/" + orth_id + ".xml"

	# calling cline() returns (stdout, stderr)
	cline = NcbiblastnCommandline(query=q, subject=s, out=xml, outfmt=5)()

	result_handle = open(xml)
	blast_records = NCBIXML.parse(result_handle)

	E_VALUE_THRESH = 0.05
	SEQ_DISPLAY = 75
	for blast_record in blast_records:
		print(blast_record.descriptions)
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					print("****Alignment****")
					print("subject:", alignment.hit_id)
					print("length:", alignment.length)
					print("e value:", hsp.expect)
					print("bits:", hsp.bits)
					print(hsp.query[0:SEQ_DISPLAY] + "...")
					print(hsp.match[0:SEQ_DISPLAY] + "...")
					print(hsp.sbjct[0:SEQ_DISPLAY] + "...")
