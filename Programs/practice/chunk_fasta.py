import argparse
import f_reader

# p3 chunk_fasta.py --fa ../../datacore/genome_celegans/1pct_elegans.fa --c 1100 --o 100

def gen_lens(start, end, chunk, overlap):
	indices = []
	for i in range(start, end, chunk - overlap):
		if i + chunk > end:
			indices.append((i, end))
			break
		else:
			indices.append((i, i + chunk))
	return indices

parser = argparse.ArgumentParser(description="split fasta files")
parser.add_argument(
	"--fa",
	required=True,
	type=str,
	metavar="<str>",
	help="required string argument",
)
parser.add_argument(
	"--c",
	required=True,
	type=int,
	metavar="<int>",
	help="required integer argument",
)
parser.add_argument(
	"--o",
	required=True,
	type=int,
	metavar="<int>",
	help="required integer argument",
)
arg = parser.parse_args()

for name, seq in f_reader.read_fasta(arg.fa):
	indices = gen_lens(0, len(seq), arg.c, arg.o)
	for i in indices:
		beg = i[0]
		end = i[1]
		print(f'>{name.split()[0]} {beg + 1}..{end}')
		print(seq[beg:end])
