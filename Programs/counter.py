import bio_lib
import f_reader
import intron_class

FA = "../WormBase/elegans/c_elegans.PRJNA13758.WS284.genomic.fa.gz"
GFF = "../WormBase/elegans/annotations.gff3.gz"

# read gff to find intron locations
chromosome = {}
for seqid, feature_i in f_reader.read_gff(
	GFF, "intron", "Caenorhabditis_elegans", True, "WormBase"
):
	chromosome[seqid] = feature_i

# read fasta to find intron sequences
ss = {}
full = {}
donor = {}
acceptor = {}
counter = 0
for header, seq in f_reader.read_fasta(FA):
	seqid = header.split()[0]
	if seqid in chromosome:

		intr_list = chromosome[seqid]
		for intr in intr_list:
			intr_seq = seq[intr.beg : intr.end]
			if len(intr_seq) < 4:
				continue

			if not intr.fwd:
				intr_seq = bio_lib.rev_comp(intr_seq)
			a = intr_seq[:2]
			b = intr_seq[-2:]
			intr.ss = a + b
			if intr.ss != 'GTAG':
				print(f'{intr.seqid}:{intr.beg}..{intr.end} {intr.ss}')
			'''
			if seqid not in ss:
				ss[seqid] = {}
			if intr.ss not in ss[seqid]:
				ss[seqid][intr.ss] = 0
			ss[seqid][intr.ss] += 1

			if intr.ss not in full:
				full[intr.ss] = 0
			full[intr.ss] += 1

			if a not in donor:
				donor[a] = 0
			donor[a] += 1

			if b not in acceptor:
				acceptor[b] = 0
			acceptor[b] += 1

			counter += 1
			'''
"""
for seqid in ss:
	print(seqid)
	ss[seqid] = {k: v for k, v in sorted(ss[seqid].items(), key=lambda x: x[1])}
	for site in ss[seqid]:
		print(site, ss[seqid][site])

print("Caenorhabiditis elegans")
full = {k: v for k, v in sorted(full.items(), key=lambda x: x[1])}
for site in full:
	print(site, full[site])
"""
'''
donor = {k: v for k, v in sorted(donor.items(), key=lambda x: x[1])}
acceptor = {k: v for k, v in sorted(acceptor.items(), key=lambda x: x[1])}
print(sum(donor.values()))
print("donor")
count = 0
for site in donor:
	print(site, donor[site])
	count += donor[site]
print(count)
print("acceptor")
for site in acceptor:
	print(site, acceptor[site])
print(counter)
'''
