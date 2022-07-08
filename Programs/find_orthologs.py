import f_reader

genes = ["WBGene00022276", "WBGene00000812", "WBGene00189949"]
file1 = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.orthologs.txt.gz"

for name, ortholog in f_reader.read_ortholog(file1):
	print(name)
	print(ortholog)
