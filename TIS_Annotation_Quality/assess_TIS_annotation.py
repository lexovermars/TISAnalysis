import os
import sys, getopt
import numpy as np
from scipy import stats
plot_result = True
try:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
except:
	print "Could not find matplotlib; will continue without distribution plot"
	plot_result = False

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
	return inputfile,fasta_file,output_name


   

## Analysis parameters ##
out_dir = "output/"
codon_search_window = 198
start_codons = ["ATG", "GTG", "TTG"]
stop_codons = ["TAA", "TAG", "TGA"]
relative_scores = True
min_number_of_orfs = 500


def init_genome(input_file,fasta_file,outname):
	genome_seq,genome_gc = get_genome_seq(fasta_file)
	genome_orfs,name = read_ptt(input_file)
	if len(genome_orfs) < min_number_of_orfs:
		print "Number of ORFs below threshold, exiting"
		sys.exit()
	candidate_starts_per_orf, initial_pca_keys = get_codon_search_seqs(genome_orfs,genome_seq,name,genome_gc)
	return None
        
def get_genome_seq(fasta_file):
	genome_seq = ""
	fna = open(fasta_file,"r")
	for line in fna:
		line = line.strip()
		if line[0] != ">":
			genome_seq += line
	fna.close()
	genome_seq = genome_seq.upper()
	genome_length = float(len(genome_seq))
	gc = (genome_seq.count("G")+genome_seq.count("C"))/genome_length
	print "Loading genome sequence done.."
	print "GC%:\t",round(gc,3)
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
	print "Loading ORF annotation done..."
	return genome_orfs,name
	
def get_codon_search_seqs(genome_orfs,genome_seq,name,genome_gc):
	candidate_starts_per_orf = {}
	initial_pca_keys = []
	candidate_start_codons = []
	all_candidate_start_codons = []
	start_freqs_up = {}
	start_freqs_down = {}
	count = 0
	
	#coding start probality
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
		orf_candidate_start_codons = []
		start,stop,strand = genome_orfs[gene]
		start -= 1
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
		
		#Get longest ORF (for upstream start determination
		if len(candidate_starts_upstream)>0:
			longest_orf = candidate_starts_upstream[-1][1]
		else:
			longest_orf = start
		# Don't include the first 30 nt
		if strand == "+":
			longest_orf_up_seq = genome_seq[longest_orf-codon_search_window:longest_orf]
		if strand == "-":
			longest_orf_up_seq = genome_seq[longest_orf:longest_orf+codon_search_window]
			longest_orf_up_seq = reverse_sequence(longest_orf_up_seq)

		orfs = find_candidate_starts_no_stop_check(longest_orf_up_seq,"backward",strand,start)
		upstream_longest_orf_starts += len(orfs)
		upstream_longest_orf_triplets += len(longest_orf_up_seq)/3
		
				
		candidate_starts_downstream = find_candidate_starts(down_seq,"forward",strand,start)
		for candidate_start_downstream in candidate_starts_downstream:
			relative_position = int(candidate_start_downstream[0])
			start_freqs_down[relative_position] += 1
	
	combined_dict = start_freqs_up
	for key in start_freqs_down.keys():
		combined_dict[key] = start_freqs_down[key]

	coding_alt_start_freq = coding_starts/float(coding_triplets)
	upstream_alt_start_freq = upstream_longest_orf_starts/float(upstream_longest_orf_triplets)
	plot_data(combined_dict,name,len(genome_orfs.keys()),coding_alt_start_freq,upstream_alt_start_freq,genome_gc)

	return candidate_starts_per_orf, initial_pca_keys
	
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

def find_candidate_starts_no_stop_check(sequence,direction,strand,start):
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
	if direction == "forward":
		for i in range(0,len(sequence),3):
			if sequence[i:i+3] in start_codons:
				codon = sequence[i:i+3]
				relative_position = i
				absolute_position = relative_position+start+3
				if strand == "+":
					relative_position += 3
				if strand == "-":
					absolute_position = start-relative_position-3
				alternative_starts.append([relative_position,absolute_position,codon])
	return alternative_starts


def plot_data(combined_dict,name,number_of_orfs,coding_alt_start_freq,upstream_alt_start_freq,genome_gc):
	N = len(combined_dict.keys())
	values = []
	function_values = []
	keys = sorted(combined_dict.keys())
	at = (1-genome_gc)/2
	gc = genome_gc/2
	stop_probability = at*at*at + at*gc*at + at*gc*at
	for label in keys:
		values.append(combined_dict[label])
		if label<0:
			#Upstream
			function_values.append(number_of_orfs*(upstream_alt_start_freq)*(1-stop_probability)**(abs(label)/3))
		else:
			#Coding
			function_values.append(number_of_orfs*(coding_alt_start_freq))
	if plot_result == True:
		ind = np.arange(N)  # the x locations for the groups
		width = 1       # the width of the bars
		fig = plt.figure()
		fig.set_size_inches(18.5,10.5)
		ax = fig.add_subplot(111)
		ax.bar(ind, values, width, color='CornflowerBlue',edgecolor = "black") #RoyalBlue?
		ax.set_xticks(ind)
		ax.set_xticklabels(keys,rotation='vertical')

		fig.set_size_inches(18.5,10.5)
		name = name.replace("\\","")
		name = name.replace("/","")
		ax.plot(function_values,color="r",linewidth=2)
		plt.ylim([0,500])
		ax.set_title(name)
	try:
		correlation = stats.spearmanr(values,function_values)
		correlation_up = stats.spearmanr(values[0:65],function_values[0:65])
		print "Quality correlation:\t",round(correlation_up[0],3)
		if not os.path.exists(out_dir):
                	os.makedirs(out_dir)
		output_file = open(out_dir+output_name+"_correlation.txt","w")
                output_file.write("Name\tGC-percentage\t#ORFs\tCorrelation Complete\tCorrelation Upstream\n")
		output_file.write(name+"\t"+str(round(genome_gc,2))+"\t"+str(number_of_orfs)+"\t"+str(round(correlation[0],2))+"\t"+str(round(correlation_up[0],2))+"\n")
		output_file.close()
		print "TIS correlation file generated.."
	except:
		print "Could not write output file.."
		pass
	if plot_result == True:
		try:
			for label in ax.get_xticklabels():
				label.set_fontsize(6)
			if not os.path.exists(out_dir):
				os.makedirs(out_dir)
			fig.savefig(out_dir+name+'_distribution.png')
			plt.clf()
			print "TIS distribution plot generated..."
		except:
			print "plot for", name,"save failed..."
			plt.clf()
	return None

if __name__ == "__main__":
        input_file,fasta_file,output_name = main(sys.argv[1:])	
	init_genome(input_file,fasta_file,output_name)
