import f_reader
import intron_class
import bio_lib
import re

FA = "../../datacore/genome_celegans/1pct_elegans.fa"
GFF = "../../datacore/genome_celegans/1pct_elegans.gff3"
ELEGANS_ID = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.geneIDs.txt.gz"
ORTH = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.orthologs.txt.gz"

FA = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.genomic.fa.gz"
GFF = "../WormBase/elegans/annotations.gff3"


# read gff to find intron locations
chromosome = {}
for seqid, feature_i in f_reader.read_gff(
	GFF, "intron", "Caenorhabditis_elegans", True, "WormBase"
):
	chromosome[seqid] = feature_i

# read fasta to find intron sequences
nc_intr = []
for header, seq in f_reader.read_fasta(FA):
	seqid = header.split()[0]
	if seqid in chromosome:

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
nc_intr = f_reader.seq_to_wbgene(ELEGANS_ID, nc_intr)
nc_intr = f_reader.read_ortholog(ORTH, nc_intr)


target1 = "Caenorhabditis_briggsae"
target_id1 = "../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.geneIDs.txt.gz"
target_fa1 = "../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.genomic.fa.gz"
target_gff1 = "../WormBase/briggsae/annotations.gff3"

target2 = "Caenorhabditis_brenneri"
target_id2 = "../WormBase/brenneri/c_brenneri.PRJNA20035.WS284.geneIDs.txt.gz"
target_fa2 = "../WormBase/brenneri/c_brenneri.PRJNA20035.WS284.genomic.fa.gz"
target_gff2 = "../WormBase/brenneri/annotations.gff3"

target3 = "Caenorhabditis_japonica"
target_id3 = "../WormBase/japonica/c_japonica.PRJNA12591.WS284.geneIDs.txt.gz"
target_fa3 = "../WormBase/japonica/c_japonica.PRJNA12591.WS284.genomic.fa.gz"
target_gff3 = "../WormBase/japonica/annotations.gff3"

TARGETS = [
	(target1, target_id1, target_fa1, target_gff1),
	(target2, target_id2, target_fa2, target_gff2),
	(target3, target_id3, target_fa3, target_gff3),
]


def target_orth(nc_intr, target, target_id, target_fa, target_gff):
	# dict of {wb : seqgene}
	seqgenes = f_reader.wb_to_seqgene(target_id, nc_intr, target)

	# swap out wbgenes for seqgenes, create target_genes set
	target_genes = set()
	for intr in nc_intr:
		if intr.orth_genes is not None:
			if target in intr.orth_genes:
				for i in range(len(intr.orth_genes[target])):
					wb = intr.orth_genes[target][i]
					if wb in seqgenes:
						intr.orth_genes[target][i] = seqgenes[wb]
						target_genes.add(seqgenes[wb])

	# create introns for target from gff
	target_introns = {}  # {seqid : intr1, intr2}
	for gff_seqid, feature_i in f_reader.read_gff(
		target_gff, "intron", target, True, "WormBase"
	):
		for intr in feature_i:
			gene = intr.gene.replace("Parent=Transcript:", "").split(";")[0]
			intr.gene = gene
			intr.clean_gene = gene.split(".")[0]
			if intr.clean_gene in target_genes:
				if intr.seqid not in target_introns:
					target_introns[intr.seqid] = []
				target_introns[intr.seqid].append(intr)

	# filter introns for ncss, create {gene : [intr1, intr2]}
	target_nc = {}
	for header, seq in f_reader.read_fasta(target_fa):
		seqid = header.split()[0]
		if seqid in target_introns:
			for intr in target_introns[seqid]:
				intr_seq = seq[intr.beg : intr.end]
				if len(intr_seq) < 4:
					continue

				if not intr.fwd:
					intr_seq = bio_lib.rev_comp(intr_seq)
				intr.seq = intr_seq

				intr.ss = intr_seq[:2] + intr_seq[-2:]

				if intr.ss not in ["GTAG", "CTAC", "NNNN"]:
					# need to check because some introns have genes both directions
					# e.g. briggsae, III:3772625..3772815
					if intr.clean_gene not in target_nc:
						target_nc[intr.clean_gene] = []
					target_nc[intr.clean_gene].append(intr)

	# match up target introns with elegans introns using orth_genes
	for intr in nc_intr:
		if intr.orth_genes is not None and target in intr.orth_genes:
			for orth_gene in intr.orth_genes[target]:
				if orth_gene in target_nc:
					for orth in target_nc[orth_gene]:
						intr.orth_intrs.append(orth)
	return nc_intr


for i in range(len(TARGETS)):
	species, gene_id, fa, gff = TARGETS[i]
	nc_intr = target_orth(nc_intr, species, gene_id, fa, gff)


for intr in nc_intr:
	if len(intr.orth_intrs) > 0:
		print(
			f">{intr.species}\t{intr.seqid}\t\t{intr.beg}..{intr.end}\t{intr.gene}\t{intr.ss}"
		)
		for orth in intr.orth_intrs:
			print(
				f"{orth.species}\t{orth.seqid}\t\t{orth.beg}..{orth.end}\t{orth.gene}\t{orth.ss}",
				end="",
			)
			if orth.ss == intr.ss:
				print(" ***")
			else:
				print()
		print()
