file = open('TreeToReads/example_out/var_site_matrix', 'r')
sample_dict = {}
for line in file:
	k = ""
	ln = line.split()
	if ln[0] not in sample_dict:
		sample_dict[ln[0]] = ""
	sample_dict[ln[0]] = sample_dict[ln[0]] + (ln[1])
#	print sample_dict
unique_seqs = {}

for key in sample_dict:
	seq = sample_dict[key]
	if seq not in unique_seqs:
		unique_seqs[seq] = 1
	else:
		unique_seqs[seq] = unique_seqs[seq] + 1

for key in unique_seqs:
	print key + " ---- " + str(unique_seqs[key])
	for sample_key in sample_dict:
		if sample_dict[sample_key] == key:
			print sample_key
				
	
