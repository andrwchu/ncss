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
	target_intr_list = []  # used to find masking
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

				# if intr.ss not in ["GTAG", "CTAC", "NNNN"]:
				if intr.ss not in ["NNNN"]:
					# need to check because some introns have genes both directions
					# e.g. briggsae, III:3772625..3772815
					if intr.clean_gene not in target_nc:
						target_nc[intr.clean_gene] = []
					target_nc[intr.clean_gene].append(intr)
					target_intr_list.append(intr)

	masked = f_reader.find_mask(target_intr_list, target_fa, target_gff)

	# match up target introns with elegans introns using orth_genes
	for intr in nc_intr:
		if intr.orth_genes is not None and target in intr.orth_genes:
			for orth_gene in intr.orth_genes[target]:
				if orth_gene in target_nc:
					for orth in target_nc[orth_gene]:
						orth.masked = masked[orth.gene]
						intr.orth_intrs.append(orth)
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
		# (target4, target_id4, target_fa4, target_gff4),
		# (target5, target_id5, target_fa5, target_gff5),
	]

	for i in range(len(TARGETS)):
		species, gene_id, fa, gff = TARGETS[i]
		nc_intr = target_orth(nc_intr, species, gene_id, fa, gff)

	return nc_intr
