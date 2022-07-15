def bl2seq_fasta(nc_intr, filepath):
	intr_counter = 0
	for intr in nc_intr:
		if len(intr.orth_intrs) > 0:
			counter_s = f"{intr_counter:03d}"
			f_name_celegans = f"{counter_s}celegans"
			f_celegans = open(filepath + f_name_celegans + ".fa", "w")

			f_celegans.write(
				f">{counter_s}_{intr.species}_{intr.seqid}_{intr.beg}..{intr.end}_{intr.gene}_{intr.ss}\n"
			)

			for i in range(0, len(intr.masked), 80):
				f_celegans.write(f"{intr.masked[i : i + 80]}\n")
			f_celegans.close()

			orth_counter = 0
			for orth in intr.orth_intrs:
				orth_counter_s = f"{orth_counter:02d}"
				f_name_orth = f"{counter_s}orth_{orth_counter_s}"
				f_orth = open(filepath + f_name_orth + ".fa", "w")

				f_orth.write(
					f">{counter_s}orth_{orth_counter_s}_{orth.species}_{orth.seqid}_{orth.beg}..{orth.end}_{orth.gene}_{orth.ss}"
				)
				if orth.ss == intr.ss:
					f_orth.write("_***\n")
				else:
					f_orth.write("\n")

				for i in range(0, len(orth.masked), 80):
					f_orth.write(f"{orth.masked[i : i + 80]}\n")
				f_orth.close()
				orth_counter += 1
			intr_counter += 1
