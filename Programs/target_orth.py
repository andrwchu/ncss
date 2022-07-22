import f_reader
import bio_lib


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


	# {gene : mask_seq}
	masked = f_reader.find_mask_orth(target_genes, target_fa, target_gff)
	# match up target introns with elegans introns using orth_genes
	for intr in nc_intr:
		if intr.orth_genes is not None and target in intr.orth_genes:
			for orth_gene in intr.orth_genes[target]:
				print(orth_gene)
				if orth_gene in masked:
					intr.orth_intrs.append((target, orth_gene, masked[orth_gene]))
	return nc_intr


def orth_intrs(nc_intr):
	target1 = "Caenorhabditis_briggsae"
	target_id1 = "../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.geneIDs.txt.gz"
	target_fa1 = "../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.genomic.fa.gz"
	target_gff1 = "../WormBase/briggsae/annotations.gff3.gz"

	target2 = "Caenorhabditis_brenneri"
	target_id2 = "../WormBase/brenneri/c_brenneri.PRJNA20035.WS284.geneIDs.txt.gz"
	target_fa2 = "../WormBase/brenneri/c_brenneri.PRJNA20035.WS284.genomic.fa.gz"
	target_gff2 = "../WormBase/brenneri/annotations.gff3.gz"

	target3 = "Caenorhabditis_japonica"
	target_id3 = "../WormBase/japonica/c_japonica.PRJNA12591.WS284.geneIDs.txt.gz"
	target_fa3 = "../WormBase/japonica/c_japonica.PRJNA12591.WS284.genomic.fa.gz"
	target_gff3 = "../WormBase/japonica/annotations.gff3.gz"

	target4 = "Caenorhabditis_remanei"
	target_id4 = "../WormBase/remanei/c_remanei.PRJNA53967.WS284.geneIDs.txt.gz"
	target_fa4 = "../WormBase/remanei/c_remanei.PRJNA53967.WS284.genomic.fa.gz"
	target_gff4 = "../WormBase/remanei/annotations.gff3.gz"

	target5 = "Pristionchus_pacificus"
	target_id5 = "../WormBase/pacificus/p_pacificus.PRJNA12644.WS284.geneIDs.txt.gz"
	target_fa5 = "../WormBase/pacificus/p_pacificus.PRJNA12644.WS284.genomic.fa.gz"
	target_gff5 = "../WormBase/pacificus/annotations.gff3.gz"

	# orth counts {'Caenorhabditis_remanei': 1459, 'Caenorhabditis_briggsae': 1427, 'Caenorhabditis_brenneri': 1371, 'Caenorhabditis_japonica': 1251, 'Brugia_malayi': 791, 'Onchocerca_volvulus': 778, 'Pristionchus_pacificus': 1030, 'Strongyloides_ratti': 737, 'Trichuris_muris': 689}

	TARGETS = [
		(target1, target_id1, target_fa1, target_gff1),
		(target2, target_id2, target_fa2, target_gff2),
		(target3, target_id3, target_fa3, target_gff3),
		(target4, target_id4, target_fa4, target_gff4),
		(target5, target_id5, target_fa5, target_gff5),
	]

	for i in range(len(TARGETS)):
		species, gene_id, fa, gff = TARGETS[i]
		nc_intr = target_orth(nc_intr, species, gene_id, fa, gff)

	return nc_intr
