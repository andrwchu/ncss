import f_reader2
import intron_class
import bio_lib
import re

fa = "../../datacore/genome_celegans/1pct_elegans.fa"
gff = "../../datacore/genome_celegans/1pct_elegans.gff3"
elegans_id = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.geneIDs.txt.gz"
orth = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.orthologs.txt.gz"
"""
fa = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.genomic.fa.gz"
gff = "../WormBase/elegans/annotations.gff3"
"""
# read gff to find intron locations
chromosome = {}
for seqid, feature_i in f_reader2.read_gff(gff, "intron", True, "WormBase"):
	chromosome[seqid] = feature_i

# read fasta to find intron sequences
nc_intr = []
for header, seq in f_reader2.read_fasta(fa):
	seqid = header.split()[0]
	if seqid in chromosome.keys():

		intr_list = chromosome[seqid]  # retrieve features list from dictionary
		for intr in intr_list:
			intr.gene = intr.gene.replace("Parent=Transcript:", "").split(";")[0]

			intr_seq = seq[intr.beg : intr.end]
			if len(intr_seq) < 4:
				continue  # should small introns be skipped?

			if not intr.fwd:
				intr_seq = bio_lib.rev_comp(intr_seq)
			intr.seq = intr_seq

			intr.ss = intr_seq[:2] + intr_seq[-2:]
			if intr.ss != "GTAG":
				nc_intr.append(intr)

# destroy chromosomes? holds canonical introns, and it is not destroyed

# find wbgene version of names, then find orthologs
nc_intr = f_reader2.seq_to_wbgene(elegans_id, nc_intr)
nc_intr = f_reader2.read_ortholog(orth, nc_intr)


def combine_orth(intr_orths):
	# combine genes from the same species (list --> dict)
	orth_dict = {}
	for orth in intr_orths:
		species = orth[0]
		gene = orth[1]
		if species not in orth_dict:
			orth_dict[species] = []
		orth_dict[species].append(gene)

	intr_orths = orth_dict
	if len(orth_dict.keys()) == 0:
		intr_orths = None

	return intr_orths


for intr in nc_intr:
	intr.orth_genes = combine_orth(intr.orth_genes)


target = "Caenorhabditis_briggsae"
target_id = "../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.geneIDs.txt.gz"
target_gff = "../WormBase/briggsae/annotations.gff3"

# gather a list of all seqgenes under target species
wbgenes = set()
for intr in nc_intr:
	if intr.orth_genes is not None:
		if target in intr.orth_genes:
			wbgenes.update(intr.orth_genes[target])
orth_info = f_reader2.wb_to_seqgene(target_id, wbgenes)
'''
orth_chromosome = {}
for gff_seqid, feature_i in f_reader2.read_gff(target_gff, "intron", True, "WormBase"):
	for intr in feature_i:
		gene = intr.gene.replace("Parent=Transcript:", "").split(";")[0].split('.')[0]
		if gene in orth_info.values(): # TODO: create a set for this
			for k, v in orth_info:
				if v[0] == gene:
					print('match')
'''
print(orth_info)

# links orth_info back to introns
"""
for intr in nc_intr:
	if intr.orth_genes is not None:
		if target in intr.orth_genes:
			for i in range(len(intr.orth_genes[target])):
				wb = intr.orth_genes[target][i]
				if wb in seqgenes:
					intr.orth_genes[target][i] = seqgenes[wb]
"""
