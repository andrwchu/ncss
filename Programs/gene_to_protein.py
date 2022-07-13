import f_reader

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
diff = f_reader.wb_to_seqgene('../WormBase/briggsae/c_briggsae.PRJNA10731.WS284.geneIDs.txt.gz', briggsae)
print(diff)
features = {}
for gff_seqid, feature_i in f_reader.read_gff(gff, "intron", True, "WormBase"):
	print(gff_seqid)
	for intron in feature_i:
		gene = intron[3].replace("Parent=Transcript:", "").split(";")[0]
		if gene in diff:
			print(intron)
