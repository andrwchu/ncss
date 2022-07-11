import re

import f_reader
import extract_intron


def wbgene_introns(fa, gff, geneID):
	# from fasta, return intron info in wbgene form
	seq_genes = []
	intron_info = []
	introns = extract_intron.extract_intron(fa, gff, 0, 5)
	for field in introns:
		# there are some extracted introns that are located on the same gene
		# since seq_genes will have duplicates, this means len(seq_genes) > len(wb_genes)

		gene = field[4]
		pat = "([\w]+\.[\d]+)"
		gene = re.search(pat, gene).group(1)

		seq_genes.append(gene)
		intron_info.append([field[0], field[1], field[2], field[3], gene])

	wb_genes = {}
	for seq, wb in f_reader.seq_to_wbgene(geneID, seq_genes):
		wb_genes[seq] = wb

	# match up wb with intron info
	# this method allows seq_to_wbgene to look for all wb in one pass
	# alternative would be to have seq_to_wbgene look for a specific wb every time it is called
	for i in range(len(intron_info)):
		info = intron_info[i]
		seq = info[4]
		if seq in wb_genes:
			intron_info[i].append(wb_genes[seq])
	return intron_info  # [seqid, start, end, ss, seq_gene, wb_gene]



file1 = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.orthologs.txt.gz"
file2 = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.geneIDs.txt.gz"
fa1 = "../../datacore/genome_celegans/1pct_elegans.fa"
gff1 = "../../datacore/genome_celegans/1pct_elegans.gff3"

intron_info = wbgene_introns(fa1, gff1, file2)
for name, orthologs in f_reader.read_ortholog(file1):
	for intron in intron_info:
		gene = intron[5]
		if name == gene:
			print(intron)
			for orth_i in orthologs:
				print(f'{orth_i[0]}\t{orth_i[1]}')
			print()
# misses intron at I; 94929, WBGene00189949
print(count)

