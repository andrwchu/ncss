import f_reader
import bio_lib
fa = '../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.genomic.fa.gz'
gff = '../WormBase/briggsae/annotations.gff3'
briggsae = {

'WBGene00026674',
'WBGene00026668',
'WBGene00086870',
'WBGene00087470',
'WBGene00086870',
'WBGene00087470',
'WBGene00086870',
'WBGene00087470',
'WBGene00029166',
'WBGene00028714',
'WBGene00035874',
'WBGene00035874',
'WBGene00025898',
'WBGene00029156',
'WBGene00037226',
'WBGene00029160',
'WBGene00029159',
'WBGene00024846',
'WBGene00024847',
'WBGene00024846',
'WBGene00024847',
'WBGene00026678'
}

diff = []
for wb, seqgene in f_reader.wb_to_seqgene('../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.geneIDs.txt.gz', briggsae):
	diff.append(seqgene)
print(diff)
features = {}
for gff_seqid, feature_i in f_reader.read_gff(gff, "intron", True, "WormBase"):
	for intron in feature_i:
		gene = intron[3].replace("Parent=Transcript:", "").split(";")[0].split('.')[0]
		if gene in diff:
			if gff_seqid not in features:
				features[gff_seqid] = {}
			if gene not in features[gff_seqid]:
				features[gff_seqid][gene] = []
			features[gff_seqid][gene].append(intron[:3])
print(features)

for name, seq in f_reader.read_fasta(fa):
	for seqid in features:
		if name == seqid:
			for gene in features[seqid]:
				for feat in features[seqid][gene]:
					start, end = feat[0], feat[1]
					fwd = feat[2]

					feat_seq = seq[start:end]
					if not fwd:
						# shortcut: skips finding rev_comp of the whole seq
						feat_seq = bio_lib.rev_comp(feat_seq)

					fprime = feat_seq[:2]
					tprime = feat_seq[-2:]

					ss = fprime + tprime
					print(name, start, end)
					print(ss)

