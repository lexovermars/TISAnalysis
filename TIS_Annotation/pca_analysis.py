import sys, getopt
import os
import scipy
from operator import itemgetter
import pickle


def main(argv):
        inputfile = ''
        output_name = ''
        try:
                opts, args = getopt.getopt(argv,"hi:o:f:",["input_file=","out_name=","fasta_file="])
        except getopt.GetoptError:
                print 'usage: assess_TIS_annotation.py -i <inputpttfile> -f <fastafile> -o <outputfile>'
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print 'usage: assess_TIS_annotation.py -i <inputfile> -f <fastafile> -o <outputfile>'
                        sys.exit()
                elif opt in ("-i", "--input_file"):
                        inputfile = arg
                elif opt in ("-o", "--out_name"):
                        output_name = arg
                elif opt in ("-f", "--fasta_file"):
                        fasta_file = arg
        if len(inputfile) == 0 or len(output_name) == 0:
                print 'usage: assess_TIS_annotation.py -i <inputfile> -f <fastafile> -o <out_name>'
                sys.exit()
	if not os.path.exists(output_dir):
        	os.makedirs(output_dir)
        if not os.path.exists(r_out_dir):
                os.makedirs(r_out_dir)
        return inputfile,fasta_file,output_name


## Analysis parameters and variables #
window_up = 30
window_down = 18
codon_search_window = 198
output_dir = "output/"
first_round_postfix = "_3_max_length.txt"
all_round_postfix = "_all_starts.txt"
r_out_dir = "R_output/"
r_script_dir = "R_scripts/"
r_pca_scores_postfix = "_pca_scores.txt"
r_pca_scores_projected_postfix = "_pca_scores_projected.txt"
r_pca_loadings_postfix = "_pca_loadings.txt"
start_codons = ["ATG", "GTG", "TTG"]
stop_codons = ["TAA", "TAG", "TGA"]
sorted_results_pickle_dump = "_pc_scores_pickle_dump.txt"

## Global variables for performance analysis 
not_top_scoring_dict = {}
not_top_scoring_initial_dict = {}
not_top3_scoring_dict = {}
sorted_pca_score_dict = {}
pca_loading_dict = {}


def make_ptt_dict(ptt_file):
	ptt_dict = {}
	ptt_data = open(ptt_file,"r")
	ptt_lines = ptt_data.readlines()
	for index,line in enumerate(ptt_lines):
		if index < 3:
			continue
		line = line.split("\t")
		strand = line[1]
		location = [int(line[0].split("..")[0]),int(line[0].split("..")[1])]
		orf_length = location[1]-location[0]
		identical_strand = "no"
		if strand == "+":
			#Prev ORF
			if index-1<3:
				continue
			prev_line = ptt_lines[index-1].split("\t")
			prev_location = [int(prev_line[0].split("..")[0]),int(prev_line[0].split("..")[1])]
			prev_strand = prev_line[1]
			intergenic_dist = location[0]-prev_location[1]
			if prev_strand == "+":
				identical_strand = "yes"
		elif strand == "-":
			if index+1>=len(ptt_lines):
				continue
			#Next ORF
			next_line = ptt_lines[index+1].split("\t")
			next_location = [int(next_line[0].split("..")[0]),int(next_line[0].split("..")[1])]
			next_strand = next_line[1]
			intergenic_dist = next_location[0]-location[1]
			if next_strand == "-":
				identical_strand = "yes"
				
				
		ptt_dict[line[5]] = [location,orf_length,line[3],line[7],line[8],intergenic_dist,identical_strand]
	return ptt_dict
	

def get_genome_seq(fna_file):
	genome_seq = ""
	fna = open(fna_file,"r")
	for line in fna:
		line = line.strip()
		if line[0] != ">":
			genome_seq += line
	fna.close()
	return genome_seq
	
def reverse_sequence(sequence):
	sequence = sequence[::-1]
	sequence = sequence.replace("A","X")
	sequence = sequence.replace("T","A")
	sequence = sequence.replace("X","T")
	sequence = sequence.replace("C","X")
	sequence = sequence.replace("G","C")
	sequence = sequence.replace("X","G")
	return sequence

def read_ptt(genome):
	genome_orfs = {}
	ptt_file = open(genome,"r")
	lines = ptt_file.readlines()
	ptt_file.close()
	for line in lines[3:]:
		line = line.strip().split("\t")
		loc = line[0].split("..")
		strand = line[1]
		locus_tag = line[5]
		genome_orfs[locus_tag] = [int(loc[0]),int(loc[1]),strand]
		sorted_pca_score_dict[locus_tag] = {}
	return genome_orfs
	
def get_codon_search_seqs(genome_orfs):
	candidate_starts_per_orf = {}
	initial_pca_keys = []
	candidate_start_codons = []
	all_candidate_start_codons = []
	count = 0
	for gene in sorted(genome_orfs.keys()):
		dist_from_longest_orf = 0
		orf_candidate_start_codons = []
		start,stop,strand = genome_orfs[gene]
		start -= 1
		annotated_start_codon = ""
		## Get the codon search sequence
		if strand == "+":						
			up_seq = genome_seq[start-codon_search_window:start]
			down_seq = genome_seq[start+3:start+codon_search_window+3]
			annotated_start_codon = genome_seq[start:start+3]
		if strand == "-":
			start = stop
			up_seq = genome_seq[stop:stop+codon_search_window]
			down_seq = genome_seq[stop-codon_search_window-3:stop-3]
			up_seq = reverse_sequence(up_seq)
			down_seq = reverse_sequence(down_seq)
			annotated_start_codon = genome_seq[stop-3:stop]
			annotated_start_codon = reverse_sequence(annotated_start_codon)
		if len(up_seq)<1 or len(down_seq)<1:
			continue
			
		## Create the dictionary entry for this ORF
		candidate_starts_per_orf[gene] = {}
		
		## Search for candidate starts upstream
		candidate_starts_upstream = find_candidate_starts(up_seq,"backward",strand,start)

		## Candidate starts found, get the upstream and downstream sequence
		for candidate_start_upstream in candidate_starts_upstream:
			relative_position,absolute_position,codon = candidate_start_upstream
			#print codon
			candidate_up_seqs = get_up_down_seqs(genome_seq,absolute_position,strand)
			orf_candidate_start_codons.append([gene,absolute_position,relative_position,strand,"upstream",candidate_up_seqs[0]+codon[0],candidate_up_seqs[1]])
			candidate_starts_per_orf[gene][(gene,absolute_position)] = [gene,absolute_position,relative_position,strand,"upstream",candidate_up_seqs[0]+codon[0],candidate_up_seqs[1]]
			
			
		## Add the annotated start
		annotated_start_seqs = get_up_down_seqs(genome_seq,start,strand)
		relative_position = 0
		orf_candidate_start_codons.append([gene,start,relative_position,strand,"annotated",annotated_start_seqs[0]+annotated_start_codon[0],annotated_start_seqs[1]])
		candidate_starts_per_orf[gene][(gene,start)] = [gene,start,relative_position,strand,"annotated",annotated_start_seqs[0]+annotated_start_codon[0],annotated_start_seqs[1]]
			
		## Search for candidate starts downstream
		candidate_starts_downstream = find_candidate_starts(down_seq,"forward",strand,start)
		
		## Candidate starts found, get the upstream and downstream sequence
		for candidate_start_downstream in candidate_starts_downstream:
			relative_position,absolute_position,codon = candidate_start_downstream	
			candidate_down_seqs = get_up_down_seqs(genome_seq,absolute_position,strand)
			#print codon
			orf_candidate_start_codons.append([gene,absolute_position,relative_position,strand,"downstream",candidate_down_seqs[0]+codon[0],candidate_down_seqs[1]])
			candidate_starts_per_orf[gene][(gene,absolute_position)] = [gene,absolute_position,relative_position,strand,"downstream",candidate_down_seqs[0]+codon[0],candidate_down_seqs[1]]
		
		for index,orf_candidate in enumerate(orf_candidate_start_codons):
			candidate_start_codons.append(orf_candidate)
			initial_pca_keys.append((orf_candidate[0],orf_candidate[1]))
			#limit to the first 3
			if index > 1:
				break
		for orf_candidate in orf_candidate_start_codons:
			all_candidate_start_codons.append(orf_candidate)

			
	return candidate_starts_per_orf, initial_pca_keys
	
def get_up_down_seqs(genome_seq,start_codon_pos,strand):
	if strand == "+":
		up_seq = genome_seq[start_codon_pos-window_up:start_codon_pos]
		down_seq = genome_seq[start_codon_pos+3:start_codon_pos+window_down+3]
	if strand == "-":
		up_seq = genome_seq[start_codon_pos:start_codon_pos+window_up]
		down_seq = genome_seq[start_codon_pos-window_down-3:start_codon_pos-3]
		up_seq = reverse_sequence(up_seq)
		down_seq = reverse_sequence(down_seq)
	return up_seq,down_seq
	
def find_candidate_starts(sequence,direction,strand,start):
	relative_position = None
	alternative_starts = []
	codon = None
	if direction == "backward":
		for i in range(len(sequence)-3,0,-3):
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i-len(sequence)
				absolute_position = relative_position+start
				if strand == "-":
					absolute_position = start-relative_position
				alternative_starts.append([relative_position,absolute_position,codon])
			if sequence[i:i+3] in stop_codons:
				break
	if direction == "forward":
		for i in range(0,len(sequence),3):
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i
				absolute_position = relative_position+start+3
				relative_position += 3
				if strand == "-":
					absolute_position = start-relative_position
				alternative_starts.append([relative_position,absolute_position,codon])
			if sequence[i:i+3] in stop_codons:
				break
	return alternative_starts

			
def make_vectors(candidate_starts_per_orf,key_selection,file_name):
	output = open(file_name,"w")
	## Write output header
	output.write("Gene\tPos_start\tPos_rel\tStart_label\tUpstream_seq\tDownstream_seq")
	for i in range(0,window_up):
		sequence_pos = str(0-window_up+i)
		output.write("\tup_T"+sequence_pos+"\tup_A"+sequence_pos+"\tup_G"+sequence_pos+"\tup_C"+sequence_pos)
	#First nt of startcodon
	output.write("\tstart_NT_T_start\tstart_NT_A\tstart_NT_G\tstart_NT_C")
	for i in range(0,window_down):
		sequence_pos = str(i)
		output.write("\tdown_T"+sequence_pos+"\tdown_A"+sequence_pos+"\tdown_G"+sequence_pos+"\tdown_C"+sequence_pos)
	output.write("\n")
	
	## Write binary vectors for each included start to file
	for orf in candidate_starts_per_orf.keys():
		for orf_start in candidate_starts_per_orf[orf].keys():
			#If key_selection == None: All starts included Else: Check if orf start is in seleciton
			if key_selection != None:
				if orf_start not in key_selection :
					continue
			gene,start_pos,rel_pos,strand,label,upstream_seq,downstream_seq = candidate_starts_per_orf[orf][orf_start]
			score_vector_T_up = score_nucleotide(upstream_seq,"T")
			score_vector_A_up = score_nucleotide(upstream_seq,"A")
			score_vector_G_up = score_nucleotide(upstream_seq,"G")
			score_vector_C_up = score_nucleotide(upstream_seq,"C")
			score_vector_T_down = score_nucleotide(downstream_seq,"T")
			score_vector_A_down = score_nucleotide(downstream_seq,"A")
			score_vector_G_down = score_nucleotide(downstream_seq,"G")
			score_vector_C_down = score_nucleotide(downstream_seq,"C")
			output.write(gene+"\t"+str(start_pos)+"\t"+str(rel_pos)+"\t"+label+"\t"+upstream_seq+"\t"+downstream_seq)	
			for index,score in enumerate(score_vector_T_up):
				output.write("\t"+str(score_vector_T_up[index])+"\t"+str(score_vector_A_up[index])+"\t"+str(score_vector_G_up[index])+"\t"+str(score_vector_C_up[index]))
			for index,score in enumerate(score_vector_T_down):
				output.write("\t"+str(score_vector_T_down[index])+"\t"+str(score_vector_A_down[index])+"\t"+str(score_vector_G_down[index])+"\t"+str(score_vector_C_down[index]))			
			output.write("\n")
	output.close()

		
def score_nucleotide(sequence,nucleotide):
	score_vector = []
	nt_dict = {"A":"T","T":"A","G":"C","C":"G"}
	for nt in sequence:
		if nt.upper() == nucleotide:
			score_vector.append(1)
		else:
			score_vector.append(0)
	return score_vector
	
	
def PCA_analysis_R(file_name):
	Rcommand = "Rscript --no-save --no-restore --verbose R_scripts/pca_scores_and_loadings.R "+str(file_name)+ "> output.Rout 2>&1"
	#print Rcommand
	command_out = os.system(Rcommand)
	
	
def analyze_subset_PCA(file_name,round_count=1):
	file_name = file_name.split(".")[0]
	pca_score_dict = {}
	pca_scores_file = open(r_out_dir+file_name+r_pca_scores_postfix,"r")
	#pca_loadings_file = open(r_out_dir+file_name+r_pca_loadings_postfix,"r")
	scores_annotated = []
	scores_non_annotated = []

	data_pca_scores_file = pca_scores_file.readlines()
	for pca_score in data_pca_scores_file[1:len(data_pca_scores_file)]:
		pca_data = pca_score.split("\t")
		pca_data_clean = []
		for data in pca_data:
			try:
				pca_data_clean.append(float(data.strip()))
			except:
				pca_data_clean.append(data.strip())
		gene,relative_pos,pos_label,pc1_score,pc2_score = pca_data_clean
		#print pca_data_clean
		if pos_label == "annotated":
			scores_annotated.append(float(pc1_score))
		else:
			scores_non_annotated.append(float(pc1_score))
		if gene not in pca_score_dict.keys():
			pca_score_dict[gene] = []
		pca_score_dict[gene].append([pc1_score,relative_pos,pos_label])
	mean_annotated = scipy.mean(scores_annotated)
	mean_non_annotated = scipy.mean(scores_non_annotated)		

	if mean_annotated > mean_non_annotated:
		pca_position_annotated = "plus"
	else:
		pca_position_annotated = "min"
	correct_count = 0
	annotated_included_count = 0
	for gene in sorted(pca_score_dict.keys()):

		if pca_position_annotated == "plus":
			sorted_scores = sorted(pca_score_dict[gene], key=itemgetter(0),reverse=True)
		else:
			sorted_scores = sorted(pca_score_dict[gene], key=itemgetter(0))
		for score in sorted_scores:
			if "annotated" in score:
				annotated_included_count += 1
		if sorted_scores[0][2] == "annotated":
			"correct"
			correct_count += 1
		#occurence analysis
		else:
			if gene not in not_top_scoring_initial_dict.keys():
				not_top_scoring_initial_dict[gene] = 1
			else:
				not_top_scoring_initial_dict[gene] += 1

	# PC1 scores annotated \t PC1 scores non-annotated \t # genes \t # Annotated in PCA \t # Annotated top scores
	print "new pca\t",mean_annotated,"\t",mean_non_annotated,"\t",len(pca_score_dict.keys()),"\t",annotated_included_count,"\t",correct_count

	#Parse the PCA loadings and load them in a dict
	parse_pca_loadings(pca_position_annotated)
	
	return pca_position_annotated,mean_annotated
	
	
def analyze_projected_PCA(file_name,pca_position_annotated,round_count=1,mean_annotated_score=None):
	print "Analyzing projected scores.."
	#print "pca_position_annotated:", pca_position_annotated
	pca_score_dict = {}
	subset_keys = []
	#Number of ORFs for which the annotated start is in the top 3 max. scores.
	annotations_in_top = 0
	#Number of ORFs for which the annotated start has the max. score
	annotated_best_score = 0
	#First Parse the PCA output
	pca_scores_file = open(r_out_dir+file_name+r_pca_scores_projected_postfix,"r")
	data_pca_scores_file = pca_scores_file.readlines()
	number_of_starts = len(data_pca_scores_file)
	for pca_score in data_pca_scores_file[1:len(data_pca_scores_file)]:
		pca_data = pca_score.split("\t")
		pca_data_clean = []
		for index,data in enumerate(pca_data):
			try:
				# adjust pc score if reversed
				if index == 4 and pca_position_annotated != "plus":
					adjusted_score = float(data.strip())*-1
					pca_data_clean.append(adjusted_score)
				else:
					pca_data_clean.append(float(data.strip()))
			except:
				pca_data_clean.append(data.strip())	
		gene,absolute_pos,relative_pos,pos_label,pc1_score,pc2_score,pc3_score = pca_data_clean
		#print pca_data_clean
		if gene not in pca_score_dict.keys():
			pca_score_dict[gene] = []
		pca_score_dict[gene].append([pos_label,pc1_score,absolute_pos,relative_pos,pc2_score,pc3_score])
	#Then determine the top 3 scoring starts per ORF
	for gene in pca_score_dict.keys():
		#print gene,pca_score_dict[gene]
		annotated_in_top = False
		sorted_scores = sorted(pca_score_dict[gene], key=itemgetter(1),reverse=True)
		#Analysis loop
		for index,score in enumerate(sorted_scores):
			if index < 3:
				subset_keys.append((gene,score[2]))
			rank = index
			pc_score = score[1]
			if score[2] not in sorted_pca_score_dict[gene].keys():
				sorted_pca_score_dict[gene][score[2]] = []
			sorted_pca_score_dict[gene][score[2]].append([round_count,score[0],pc_score,rank,score[3],str(round(score[4],2)),str(round(score[5],2))])
	return subset_keys,annotated_best_score
	
def parse_pca_loadings(pca_position_annotated="plus"):
	pca_loadings_file = open(r_out_dir+output_name+r_pca_loadings_postfix,"r")
	pca_data = pca_loadings_file.readlines()
	pca_loadings_file.close()
	for data in pca_data[1:len(pca_data)]:
		data = data.strip().split("\t")
		nucleotide_id = data[0].strip()
		pca_loading1 = round(float(data[1].strip()),2)
		pca_loading2 = round(float(data[2].strip()),2)
		pca_loading3 = round(float(data[3].strip()),2)
		if pca_position_annotated != "plus":
			pca_loading1 = pca_loading1*-1
			pca_loading2 = pca_loading2*-1
			pca_loading3 = pca_loading3*-1
		if nucleotide_id not in pca_loading_dict.keys():
			pca_loading_dict[nucleotide_id] = []
		pca_loading_dict[nucleotide_id].append([pca_loading1,pca_loading2,pca_loading3])
		



if __name__ == "__main__":
        input_file,fasta_file,output_name = main(sys.argv[1:])  

	ptt_dict = make_ptt_dict(input_file)
	genome_seq = get_genome_seq(fasta_file)
	genome_orfs = read_ptt(input_file)
	candidate_starts_per_orf, initial_pca_keys = get_codon_search_seqs(genome_orfs)

	#Make initial subset binary vector file
	make_vectors(candidate_starts_per_orf,initial_pca_keys,output_dir+output_name+first_round_postfix)
	#Make complete binary vector file (all potential starts)
	make_vectors(candidate_starts_per_orf,None,output_dir+output_name+all_round_postfix)

	#Do initial PCA analysis (R script performs PCA on subset and projected PCA on complete set
	prev_annotated_best_score = 0
	print "Initial PCA"
	PCA_analysis_R(output_name)
	pca_position_annotated,mean_annotated = analyze_subset_PCA(output_name,0)
	subset_keys,annotated_best_score = analyze_projected_PCA(output_name,pca_position_annotated,0,mean_annotated)

	#Do iterative PCA on subsets
	round_count = 0
	#while annotated_best_score > prev_annotated_best_score:
	while 1:
		prev_annotated_best_score = annotated_best_score
		round_count += 1
		print "\nPCA Iteration ",round_count
		#Create new subset file
		make_vectors(candidate_starts_per_orf,subset_keys,output_dir+output_name+first_round_postfix)
		PCA_analysis_R(output_name)
		pca_position_annotated,mean_annotated = analyze_subset_PCA(output_name,round_count)
		subset_keys,annotated_best_score = analyze_projected_PCA(output_name,pca_position_annotated,round_count,mean_annotated)
		if round_count > 8:
			break

	#Delete temporary files
	try:
		os.remove(output_dir+output_name+first_round_postfix)
		os.remove(output_dir+output_name+all_round_postfix)
	except:
		pass

	#Write sorted scores to a dic dump
	dic_file = open(output_dir+output_name+sorted_results_pickle_dump,"w")
	pickle.dump(sorted_pca_score_dict,dic_file)
	dic_file.close()

	#Write pca loadings to a dic dump
	dic_file = open(output_dir+output_name+"_pca_loadings_pickle_dump.txt","w")
	pickle.dump(pca_loading_dict,dic_file)
	dic_file.close()
	
	

