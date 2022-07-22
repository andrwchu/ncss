from dataclasses import dataclass


@dataclass
class Gene:
	__slots__ = (
		"species",
		"seqid",
		"beg",
		"end",
		"name",
		"wbgene",
		"masked",
		"introns",
		"exons",
		"orth_genes",
	)
	species: str
	seqid: str
	beg: int
	end: int
	name: str
	wbgene: str
	masked: str
	introns: []
	exons: []
	orth_genes: []  # list of Genes
