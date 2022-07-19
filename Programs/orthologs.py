import bio_lib
import f_reader
import intron_class

from target_orth import orth_intrs
from gen_bl2seq import bl2seq_fasta

FA = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.genomic.fa.gz"
GFF = "../WormBase/elegans/annotations.gff3.gz"
ELEGANS_ID = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.geneIDs.txt.gz"
ORTH = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.orthologs.txt.gz"


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

# create masked versions of gene sequences
masked = f_reader.find_mask(nc_intr, FA, GFF)
for intr in nc_intr:
	if intr.gene in masked:
		intr.masked = masked[intr.gene]

# find wbgene version of names, then find orthologs
nc_intr = f_reader.seq_to_wbgene(ELEGANS_ID, nc_intr)
nc_intr = f_reader.read_ortholog(ORTH, nc_intr)

nc_intr = orth_intrs(nc_intr)


if __name__ == "__main__":
	bl2seq_fasta(nc_intr, "../out/bl2seq2/")

	'''
	orths = {}
	for intr in nc_intr:
		if intr.orth_genes is not None:
			for worm in intr.orth_genes:
				if worm not in orths: orths[worm] = 0
				orths[worm] += 1
	print(orths)
	'''
