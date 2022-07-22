def bl2seq_fasta(nc_intr, filepath):
	intr_counter = 0
	for intr in nc_intr:
		if len(intr.orth_intrs) > 0:
			counter_s = f"{intr_counter:04d}"
			f_name_celegans = f"{counter_s}celegans"
			f_celegans = open(filepath + f_name_celegans + ".fa", "w")
			f_celegans.write(
				f">{counter_s},{intr.species},{intr.seqid},{intr.beg},{intr.end},{intr.gene},{intr.ss}\n"
			)

			for i in range(0, len(intr.masked), 80):
				f_celegans.write(f"{intr.masked[i : i + 80]}\n")
			f_celegans.close()

			orth_counter = 0
			for spec, orth_gene, masked in intr.orth_intrs:
				orth_counter_s = f"{orth_counter:03d}"
				f_name_orth = f"{counter_s}orth_{orth_counter_s}"
				f_orth = open(filepath + f_name_orth + ".fa", "w")

				f_orth.write(
					f">{counter_s}orth_{orth_counter_s},{spec},{orth_gene}\n"
				)

				for i in range(0, len(masked), 80):
					f_orth.write(f"{masked[i : i + 80]}\n")
				f_orth.close()
				orth_counter += 1
			intr_counter += 1
