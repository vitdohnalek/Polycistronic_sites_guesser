from Bio import SeqIO


def predict_possible_polycistronic_genes(gene_dictionary, chromosome):

	chromosome_sequence = ""
	for seq_rec in SeqIO.parse("/home/vit/Desktop/Elimu/giardia_genome/GCF_000002435.2_UU_WB_2.1_genomic.fna", "fasta"):
		if f"chromosome {chromosome}" in seq_rec.description:
			chromosome_sequence = str(seq_rec.seq).upper()
	#Sequence for minus orientation
	translation_table = str.maketrans("ATCG", "TAGC")
	reverse_chromosome_sequence = chromosome_sequence[::-1].translate(translation_table)

	#Sort the dictionary by the start positions
	sorted_genes = sorted(gene_dictionary.items(), key=lambda x: x[1][0])

	close_genes = []

	for i in range(len(sorted_genes) - 1):
		
		distance_check = False
		CAAT_check = False
		AT_check = True
	    
		current_gene = sorted_genes[i]
		next_gene = sorted_genes[i + 1]
	    
	    #Distance check
		current_end = current_gene[1][1]
		next_start = next_gene[1][0]
	    
		distance = next_start - current_end
		if 0 < distance <= 25:
			distance_check = True

	    #CAAT check
		if current_gene[1][2] == "plus":
			upstream_sequence = chromosome_sequence[next_start-300:next_start-1]
			if not "CAAT" in upstream_sequence:
				CAAT_check = True
		else:
			upstream_sequence = reverse_chromosome_sequence[next_start-300:next_start-1]
			if not "CAAT" in upstream_sequence:
				CAAT_check = True

	    #AT check
		if current_gene[1][2] == "plus":
			at_required_site = chromosome_sequence[next_start-300:next_start-1]
		else:
			at_required_site = reverse_chromosome_sequence[next_start-300:next_start-1]

		def check_sites(sequence):
			AT_content = True
			window_sizes = [11] # Possible window sizes
			seq_length =len(sequence)

			def check_subsequence(subseq):
				a_count = subseq.count("A")
				t_count = subseq.count("T")
				if (a_count + t_count) >= 10:
					return True

			for size in window_sizes:
				if size > seq_length:
					continue
				for start in range(seq_length - size + 1):
					subseq = sequence[start:start + size]
					if check_subsequence(subseq):
						AT_content = False
			return(AT_content)

		AT_check = check_sites(at_required_site)

		#Conditions check
		if distance_check == True and CAAT_check == True and AT_check == True:
		    close_genes.append((current_gene[0], next_gene[0]))

	return close_genes


def get_info(chromosome_n):

	genes_info_plus = {}
	genes_info_minus = {}
	gl_IDs = {}

	with open("giardia_genes_info.txt", "r") as f:
		for l in f:
			info = l.strip().split("\t")
			if l.startswith("5741") and "GL50803_" in l and info[-4].isdigit():
				gene_name = info[2]
				orientation = info[-2]
				start_position = int(info[-4]) - 1 #-1 is neccessary offset for seq_rec.seq to get the correct first nucleotide
				end_position = int(info[-3])
				chromosome = int(info[-6])
				if chromosome == chromosome_n:
					if orientation == "plus":
						genes_info_plus[gene_name] = [start_position, end_position, "plus"]
					elif orientation == "minus":
						genes_info_minus[gene_name] = [start_position, end_position, "minus"]
	return genes_info_plus, genes_info_minus

gl_IDs = {}
with open("giardia_genes_info.txt", "r") as f:
	for l in f:
		if l.startswith("5741") and "GL50803_" in l:
			gl_IDs[l.split("\t")[2]] = l.split("\t")[5].split()[0]

plus_results = ""
minus_results = ""

for i in range(5):
	plus_results += f"Chromosome {i + 1}:\n"
	minus_results += f"Chromosome {i + 1}:\n"
	genes_info_plus, genes_info_minus = get_info(i + 1)
	plus = predict_possible_polycistronic_genes(genes_info_plus, i + 1)
	minus = predict_possible_polycistronic_genes(genes_info_minus, i + 1)

	for gene_pair in plus:
		plus_results += f"{gl_IDs[gene_pair[0]]}-{gl_IDs[gene_pair[1]]}\n"
	for gene_pair in minus:
		minus_results += f"{gl_IDs[gene_pair[0]]}-{gl_IDs[gene_pair[1]]}\n"

results = f"Plus orientation:\n{plus_results}\nMinus orientation:\n{minus_results}"
with open("Possible_sites.txt", "w") as f:
	f.write(results)
