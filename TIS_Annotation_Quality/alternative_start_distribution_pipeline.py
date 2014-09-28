import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import pickle
import sys, getopt

def main(argv):
	inputfile = ''
	outputfile_name = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'pca_pipeline.py -i <inputpttfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'pca_pipeline.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile_name = arg
	#print 'Input file is "', inputfile
	#print 'Output file is "', outputfile_name
	return inputfile,outputfile_name
   
if __name__ == "__main__":
	inputfile, outputfile_name = main(sys.argv[1:])

genomes_dir = "/data/genomes/"
genome = inputfile


window_up = 18
window_down = 18
codon_search_window = 198
start_codons = ["ATG", "GTG", "TTG"]
stop_codons = ["TAA", "TAG", "TGA"]
relative_scores = True

correlation_matrix_adjusted = open("final_output/correlation_matrix_adjusted_refset1.txt","a")
correlation_matrix = open("final_output/correlation_matrix_refset1.txt","a")
#correlation_matrix_adjusted = open("final_output/small_correlation_matrix_adjusted_refset.txt","a")
#correlation_matrix = open("final_output/small_correlation_matrix_refset.txt","a")

tax_dict = pickle.load(open("taxonomy_dic.txt"))

def init_genome(file_location):
	genome_seq,genome_gc = get_genome_seq(file_location)
	genome_orfs,name = read_ptt(file_location)
	if len(genome_orfs) < 500:
		sys.exit()
	adjusted_genome_orfs = []
	name = name.strip()
	uid_code = file_location.split("/")[-1].split("_")[-1].strip("uid")
	print uid_code
	print name
	#Analyze adjusted annotation
	if os.path.exists('final_output/'+name+"_adjusted_annotation.txt"):
		adjusted_genome_orfs = read_adjusted_annotation(outputfile_name,genome_orfs)	
		candidate_starts_per_orf, initial_pca_keys = get_codon_search_seqs(adjusted_genome_orfs,genome_seq,name,uid_code,genome_gc,mode="adjusted")
	#Analyze original annotation
	candidate_starts_per_orf, initial_pca_keys = get_codon_search_seqs(genome_orfs,genome_seq,name,uid_code,genome_gc,mode="annotation")
	print "# orfs in adjusted: ",len(adjusted_genome_orfs),"# orfs in refseq: ",len(genome_orfs)
	return None
        


def get_genome_seq(genome):
	genome_seq = ""
	fna_file = genome.strip(".ptt")+".fna"
	fna = open(fna_file,"r")
	for line in fna:
		line = line.strip()
		if line[0] != ">":
			genome_seq += line
	fna.close()
	genome_seq = genome_seq.upper()
	genome_length = float(len(genome_seq))
	gc = (genome_seq.count("G")+genome_seq.count("C"))/genome_length
	print "Genome sequence done.."
	return genome_seq,gc
	
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
	name = lines[0].split(",")[0]
	for line in lines[3:]:
		line = line.strip().split("\t")
		loc = line[0].split("..")
		strand = line[1]
		locus_tag = line[5]
		genome_orfs[locus_tag] = [int(loc[0]),int(loc[1]),strand]
	print "PTT file done.."
	#print genome_orfs
	return genome_orfs,name
	
def read_adjusted_annotation(outputfile_name,genome_orfs):
	file = open("final_output/"+outputfile_name+"_adjusted_annotation.txt","r")
	lines = file.readlines()
	adjusted_genome_orfs = {}
	for line in lines:
		data = line.strip().split("\t")
		locus_tag = data[0].strip()
		abs_position = int(float(data[3]))
		start,stop,strand = genome_orfs[locus_tag]
		adjusted_genome_orfs[locus_tag] = [abs_position,stop,strand]
	return adjusted_genome_orfs
		
	
	
def get_codon_search_seqs(genome_orfs,genome_seq,name,uid_code,genome_gc,mode):
	candidate_starts_per_orf = {}
	initial_pca_keys = []
	candidate_start_codons = []
	all_candidate_start_codons = []
	start_freqs_up = {}
	start_freqs_down = {}
	count = 0
	
	#coding start chance
	coding_starts = 0
	coding_triplets = 0
	upstream_longest_orf_starts = 0
	upstream_longest_orf_triplets = 0

	for i in range(codon_search_window*-1,0,3):
		start_freqs_up[i] = 0
	for i in range(0,codon_search_window+1,3):
		start_freqs_down[i] = 0
			
	for gene in sorted(genome_orfs.keys()):
		dist_from_longest_orf = 0
		#if gene != "lp_0001":
		#	continue
		orf_candidate_start_codons = []
		start,stop,strand = genome_orfs[gene]
		start -= 1
		#print gene,strand,start,stop
		end_seq_window = codon_search_window
		if (stop-start)<end_seq_window:
			end_seq_window = stop-start
		## Get the codon search sequence
		if strand == "+":						
			up_seq = genome_seq[start-codon_search_window:start]
			down_seq = genome_seq[start+3:start+codon_search_window+3]
			end_seq = genome_seq[stop-end_seq_window:stop]
		if strand == "-":
			start_non_reversed = start
			start = stop
			up_seq = genome_seq[stop:stop+codon_search_window]
			down_seq = genome_seq[stop-codon_search_window-3:stop-3]
			up_seq = reverse_sequence(up_seq)
			down_seq = reverse_sequence(down_seq)
			end_seq = reverse_sequence(genome_seq[start_non_reversed:start_non_reversed+end_seq_window])
		if len(up_seq)<1 or len(down_seq)<1:
			#print "Could not retrieve sequence for ", gene
			continue
			
		#Search for alternative starts in end of gene
		alternative_starts_in_coding = find_candidate_starts(end_seq,"forward",strand,start)
		coding_starts += len(alternative_starts_in_coding)
		coding_triplets += len(end_seq)/3
		
		#Create the dictionary entry for this ORF
		candidate_starts_per_orf[gene] = {}
		
		## Search for candidate starts upstream
		candidate_starts_upstream = find_candidate_starts(up_seq,"backward",strand,start)
		for candidate_start_upstream in candidate_starts_upstream:
			relative_position = int(candidate_start_upstream[0])
			start_freqs_up[relative_position] += 1
		#print candidate_starts_upstream
		
		#Get longest ORF (for upstream start determination
		if len(candidate_starts_upstream)>0:
			longest_orf = candidate_starts_upstream[-1][1]
			#print gene,strand,name,longest_orf,candidate_starts_upstream[-1]
		else:
			longest_orf = start
		# Don't include the first 30 nt
		if strand == "+":
			longest_orf_up_seq = genome_seq[longest_orf-codon_search_window:longest_orf]
			#longest_orf_up_seq = genome_seq[longest_orf-codon_search_window:longest_orf-30]
		if strand == "-":
			longest_orf_up_seq = genome_seq[longest_orf:longest_orf+codon_search_window]
			#longest_orf_up_seq = genome_seq[longest_orf+30:longest_orf+codon_search_window]
			longest_orf_up_seq = reverse_sequence(longest_orf_up_seq)

		orfs = find_candidate_starts_no_stop_check(longest_orf_up_seq,"backward",strand,start)
		upstream_longest_orf_starts += len(orfs)
		upstream_longest_orf_triplets += len(longest_orf_up_seq)/3
		
				
		candidate_starts_downstream = find_candidate_starts(down_seq,"forward",strand,start)
		for candidate_start_downstream in candidate_starts_downstream:
			relative_position = int(candidate_start_downstream[0])
			start_freqs_down[relative_position] += 1
	#print start_freqs_up
	#print start_freqs_down
	
	combined_dict = start_freqs_up
	for key in start_freqs_down.keys():
		combined_dict[key] = start_freqs_down[key]

	coding_alt_start_freq = coding_starts/float(coding_triplets)
	upstream_alt_start_freq = upstream_longest_orf_starts/float(upstream_longest_orf_triplets)
	print upstream_alt_start_freq
	print coding_alt_start_freq,upstream_alt_start_freq,genome_gc,mode
	plot_data(combined_dict,name,len(genome_orfs.keys()),uid_code,coding_alt_start_freq,upstream_alt_start_freq,genome_gc,mode)

	return candidate_starts_per_orf, initial_pca_keys
	
def get_up_down_seqs(genome_seq,start_codon_pos,strand):
	for gene in genome_orfs.keys():
		if strand == "+":
			up_seq = genome_seq[start_codon_pos-window_up:start_codon_pos]
			down_seq = genome_seq[start_codon_pos+3:start_codon_pos+window_up+3]
		if strand == "-":
			up_seq = genome_seq[start_codon_pos:start_codon_pos+window_down]
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
			#print sequence[i:i+3], i-len(sequence)
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i-len(sequence)
				absolute_position = relative_position+start
				if strand == "-":
					absolute_position = start-relative_position
				#print absolute_position
				#print "Identified start: ", sequence[i:i+3], "position: ", i
				alternative_starts.append([relative_position,absolute_position,codon])
				#break
			if sequence[i:i+3] in stop_codons:
				break
	if direction == "forward":
		for i in range(0,len(sequence),3):
			#print sequence[i:i+3], i
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i
				absolute_position = relative_position+start+3
				#if strand == "+":
				#	relative_position += 3
				relative_position += 3
				if strand == "-":
					absolute_position = start-relative_position
				#print absolute_position
				#print "Identified start:", sequence[i:i+3], "position: ", i
				alternative_starts.append([relative_position,absolute_position,codon])
				#break
			if sequence[i:i+3] in stop_codons:
				break
	return alternative_starts
	#return relative_position,codon

def find_candidate_starts_no_stop_check(sequence,direction,strand,start):
	relative_position = None
	alternative_starts = []
	codon = None
	if direction == "backward":
		for i in range(len(sequence)-3,0,-3):
			#print sequence[i:i+3], i-len(sequence)
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i-len(sequence)
				absolute_position = relative_position+start
				if strand == "-":
					absolute_position = start-relative_position
				#print absolute_position
				#print "Identified start: ", sequence[i:i+3], "position: ", i
				alternative_starts.append([relative_position,absolute_position,codon])
				#break
			#if sequence[i:i+3] in stop_codons:
			#	break
	if direction == "forward":
		for i in range(0,len(sequence),3):
			#print sequence[i:i+3], i
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i
				absolute_position = relative_position+start+3
				if strand == "+":
					relative_position += 3
				if strand == "-":
					absolute_position = start-relative_position-3
				#print absolute_position
				#print "Identified start:", sequence[i:i+3], "position: ", i
				alternative_starts.append([relative_position,absolute_position,codon])
				#break
			#if sequence[i:i+3] in stop_codons:
			#	break
	return alternative_starts
	#return relative_position,codon


def plot_data(combined_dict,name,number_of_orfs,uid_code,coding_alt_start_freq,upstream_alt_start_freq,genome_gc,mode):
	if 1:
		N = len(combined_dict.keys())
		values = []
		function_values = []
		keys = sorted(combined_dict.keys())
		for label in keys:
			#print type(label),up_freqs[label]
			values.append(combined_dict[label])
			if label<0:
				#function_values.append(number_of_orfs*(3/float(64))*(1-1/float(32))**(abs(label)/3))
				#with coding start freqs
				#function_values.append(number_of_orfs*(coding_alt_start_freq)*(1-1/float(32))**(abs(label)/3))
				function_values.append(number_of_orfs*(upstream_alt_start_freq)*(1-3/float(64))**(abs(label)/3))
				
			else:
				#function_values.append(number_of_orfs*(3/64.0))
				function_values.append(number_of_orfs*(coding_alt_start_freq))
		ind = np.arange(N)  # the x locations for the groups
		width = 1       # the width of the bars
		fig = plt.figure()
		fig.set_size_inches(18.5,10.5)
		ax = fig.add_subplot(111)
		ax.bar(ind, values, width, color='CornflowerBlue',edgecolor = "black") #RoyalBlue?
		#ax.set_title(name)
		ax.set_xticks(ind)
		ax.set_xticklabels(keys,rotation='vertical')

		fig.set_size_inches(18.5,10.5)
		name = name.replace("\\","")
		name = name.replace("/","")
		ax.plot(function_values,color="r",linewidth=2)
		plt.ylim([0,500])
		try:
			tax_data = tax_dict[uid_code]
		except:
			tax_data = [""]*10
		tax_string = "\t".join(tax_data)
		sum_of_square_dif_up = np.sum((np.asarray(values[0:65])-np.asarray(function_values[0:65]))**2)
		#sum_of_square_dif_up = np.sum((np.asarray(values[0:59])-np.asarray(function_values[0:59]))**2)
		prefix = ""
		try:
			correlation = stats.spearmanr(values,function_values)
			correlation_up = stats.spearmanr(values[0:65],function_values[0:65])
			#correlation_up = stats.spearmanr(values[0:59],function_values[0:59])
			if mode == "adjusted":
				correlation_matrix_adjusted.write(name+"\t"+str(genome_gc)+"\t"+uid_code+"\t"+tax_string+"\t"+str(number_of_orfs)+"\t"+str(round(correlation[0],2))+"\t"+str(round(correlation_up[0],2))+"\t"+str(sum_of_square_dif_up)+"\n")
			else:
				correlation_matrix.write(name+"\t"+str(genome_gc)+"\t"+uid_code+"\t"+tax_string+"\t"+str(number_of_orfs)+"\t"+str(round(correlation[0],2))+"\t"+str(round(correlation_up[0],2))+"\t"+str(sum_of_square_dif_up)+"\n")
			
			#uncomment to change low scoring file names
			#if correlation_up[0]<0.9:
			#	prefix = "low_"

			#corr_string = "Correlation complete & up: "+str(round(correlation[0],2))+" & "+str(round(correlation_up[0],2))
			#ax.text(3, 200, corr_string, fontsize=16)

			#title including correlation
			#ax.set_title(name+" ("+str(round(correlation_up[0],2))+")")

			#title excluding correlation
			ax.set_title(name)
		except:
			pass
		try:
			for label in ax.get_xticklabels():
				label.set_fontsize(6)
			if mode == "adjusted":
				fig.savefig('start_freqs_200nt_trendline/'+prefix+name+'_adjusted'+'.png')
			else:
				fig.savefig('start_freqs_200nt_trendline/'+prefix+name+'.png')
			plt.clf()
		except:
			print "plot for", name,"save failed..."
			plt.clf()
	else:
		print "plot for", name,"failed..."
	return None
	
init_genome(genome)
