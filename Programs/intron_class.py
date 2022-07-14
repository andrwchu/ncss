from dataclasses import dataclass
import re

@dataclass
class Intron:
	__slots__ = ("seqid", "beg", "end", "fwd", "gene", "clean_gene", "wbgene", "ss", "seq", "orth_genes")
	seqid: str
	beg: int
	end: int
	fwd: bool
	gene: str
	clean_gene: str
	wbgene: str
	ss: str
	seq: str
	orth_genes: [] # list of [worm, wbgene]
