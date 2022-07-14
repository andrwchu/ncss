from dataclasses import dataclass
import re


@dataclass
class Intron:
	__slots__ = (
		"species",
		"seqid",
		"beg",
		"end",
		"fwd",
		"gene",
		"clean_gene",
		"wbgene",
		"ss",
		"seq",
		"orth_genes",
		"orth_intrs",
	)
	species: str
	seqid: str
	beg: int
	end: int
	fwd: bool
	gene: str
	clean_gene: str
	wbgene: str
	ss: str
	seq: str
	orth_genes: []  # list of lists, then dict of {worm : [wbgene1, wbgene2]}
	orth_intrs: []
