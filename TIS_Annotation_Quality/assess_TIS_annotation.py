import os
import sys
import argparse
import numpy as np
from scipy import stats

try:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	plot_result = True
except:
	print "Could not find matplotlib; will continue without distribution plot"
	plot_result = False	


def get_args():
	"""This function parses and return arguments passed in"""
	# Assign description to the help doc
	parser = argparse.ArgumentParser(description='Script determines the TIS annotation quality of a bacterial genome')
	# Add arguments
	parser.add_argument('-i', '--pttfile', type=str, help='PTT file', required=True)
	parser.add_argument('-f', '--fastafile', type=str, help='Fasta file', required=True)
	parser.add_argument('-o', '--output', type=str, help='Output name', required=False, default="TIS_annotation_out")
	args = parser.parse_args()
	# Assign args to variables
	ptt_file = args.pttfile
	fasta_file = args.fastafile
	output_name = args.output
	# Return all variable values
	return ptt_file, fasta_file, output_name


def reverse_sequence(sequence):
	"""Return the reverse complement sequence of the genome sequence as a string."""
	sequence = sequence[::-1]
	sequence = sequence.replace("A","X")	
	sequence = sequence.replace("T","A")
	sequence = sequence.replace("X","T")
	sequence = sequence.replace("C","X")
	sequence = sequence.replace("G","C")
	sequence = sequence.replace("X","G")
	return sequence
	
	
def plot_alternative_start_distribution(Genome):
	"""Creates a png file with a plot of the distribution of the observed alternative start codons vs the expexted number of alternative start codons"""
	N = len(Genome.start_codon_freqs.keys())
	keys = sorted(Genome.start_codon_freqs.keys())
	ind = np.arange(N)  # the x locations for the groups
	width = 1       # the width of the bars	
	fig = plt.figure()
	fig.set_size_inches(18.5,10.5)
	ax = fig.add_subplot(111)
	ax.bar(ind, Genome.values, width, color='CornflowerBlue',edgecolor = "black") 
	ax.set_xticks(ind)
	ax.set_xticklabels(keys,rotation='vertical')
	fig.set_size_inches(18.5,10.5)
	name = Genome.name.replace("\\","")
	name = name.replace("/","")
	ax.plot(Genome.expected_values,color="r",linewidth=2)
	plt.ylim([0,500])
	ax.set_title(name)

	for label in ax.get_xticklabels():
		label.set_fontsize(6)

	fig.savefig(name+'_distribution.png')
	plt.clf()
	print "TIS distribution plot generated..."
	return None


class Genome:
    
	min_number_of_orfs = 500
	
	def __init__(self, ptt_file,fasta_file, out_name):
		self.fasta_file = fasta_file
		self.ptt_file = ptt_file
		self.genome_seq = ""
		self.gc = 0
		self.genome_seq,self.genome_gc = self.read_fasta_file()
		self.genome_orfs,self.name = self.read_ptt_file()
		if len(self.genome_orfs) < self.min_number_of_orfs:
			print "Number of ORFs below threshold, exiting"
			sys.exit()

		self.start_codon_freqs = ORF.start_freqs_up
		for key in ORF.start_freqs_down.keys():
			self.start_codon_freqs[key] = ORF.start_freqs_down[key]

		self.coding_alt_start_freq = ORF.coding_starts/float(ORF.coding_triplets)		
		self.upstream_alt_start_freq = ORF.upstream_longest_orf_starts/float(ORF.upstream_longest_orf_triplets)
		
		self.expected_values,self.values = self.determine_expected_values()	
		self.correlation,self.correlation_up = self.calculate_expected_vs_observed_correlation()
		

	def read_ptt_file(self):
		"""Parse the input ptt file and return a dictionary with locus tag as a key and a list with the start, stop and strand as a value."""
		genome_orfs = []
		lines = [line.rstrip('\n') for line in open(self.ptt_file)]
		name = lines[0].split(",")[0]
		for line in lines[3:]:
			line = line.split("\t")
			loc = line[0].split("..")		
			strand = line[1]
			genome_orfs.append(ORF(self,int(loc[0]),int(loc[1]),strand))
		print "Loading ORF annotation done..."
		return genome_orfs, name

		
	def read_fasta_file(self):
		"""Parse the input fasta file and return the genome sequence as a string."""
		with open(self.fasta_file) as fna:
			for line in fna:
				if line[0] != ">":
					self.genome_seq += line.strip()
		self.genome_seq = self.genome_seq.upper()
		genome_length = float(len(self.genome_seq))
		self.gc = (self.genome_seq.count("G")+self.genome_seq.count("C"))/genome_length
		print "Loading genome sequence done.."
		print "GC%:\t",round(self.gc,3)
		return self.genome_seq,self.gc
		
		
	def determine_expected_values(self):
		"""Calculates the expected number of starts per relative position based on the GC% of the genome and returns these stop_probabilities in a list"""
		N = len(self.start_codon_freqs.keys())
		values = []
		function_values = []
		keys = sorted(self.start_codon_freqs.keys())
		at = (1-self.gc)/2
		gc = self.gc/2
		stop_probability = at*at*at + at*gc*at + at*gc*at
		for label in keys:
			values.append(self.start_codon_freqs[label])
			if label<0:
				#Upstream
				function_values.append(ORF.orfCount*(self.upstream_alt_start_freq)*(1-stop_probability)**(abs(label)/3))
			else:
				#Coding
				function_values.append(ORF.orfCount*(self.coding_alt_start_freq))
		return function_values,values
	
				
	def calculate_expected_vs_observed_correlation(self):	
		correlation = stats.spearmanr(self.values,self.expected_values)
		correlation_up = stats.spearmanr(self.values[0:65],self.expected_values[0:65])
		print "Quality correlation:\t",round(correlation_up[0],3)
		return correlation,correlation_up


class ORF():

	orfCount = 0
	codon_search_window = 198
	start_codons = ["ATG", "GTG", "TTG"]
	stop_codons = ["TAA", "TAG", "TGA"]
	start_freqs_up,start_freqs_down={},{}
	coding_starts = 0
	coding_triplets = 0
	upstream_longest_orf_starts = 0
	upstream_longest_orf_triplets = 0
	
	for i in range(codon_search_window*-1,0,3):
		start_freqs_up[i] = 0
	for i in range(0,codon_search_window+1,3):
		start_freqs_down[i] = 0
	
	def __init__(self, genome,start,stop,strand):
		ORF.orfCount += 1
		self.genome = genome
		self.start = start
		self.stop = stop
		self.strand = strand
		
		# Retrieve sequences upstream, downstream and and the end of the ORF		
		self.sequence_upstream,self.sequence_downstream,self.sequence_end_coding = self.get_codon_search_seqs()

		# Identify the start codons in these sequences		
		self.candidate_starts_upstream = self.find_candidate_starts("backward",self.sequence_upstream)
		self.candidate_starts_downstream = self.find_candidate_starts("forward",self.sequence_downstream)		
		self.alternative_starts_in_coding = self.find_candidate_starts("forward",self.sequence_end_coding)	
		
		# Store the upstream and downstream positions of the start in the corresponding dictionaries
		for candidate_start_upstream in self.candidate_starts_upstream:
			relative_position = int(candidate_start_upstream[0])
			self.start_freqs_up[relative_position] += 1
				
		for candidate_start_downstream in self.candidate_starts_downstream:
			relative_position = int(candidate_start_downstream[0])
			self.start_freqs_down[relative_position] += 1
		
		# Determine number of observed start and observed triplets in coding sequencing
		ORF.coding_starts += len(self.alternative_starts_in_coding)
		ORF.coding_triplets += len(self.sequence_end_coding )/3
		
		# Based on the identified start codons, retrieve the sequence of the longest possible orf
		self.longest_orf_up_seq = self.get_longest_orf_seq()
		# Retrieve the start codons upstream of the longest ORF (regardless of any stop codons)
		self.orfs = self.find_candidate_starts("backward",self.longest_orf_up_seq,no_stop=False)
		
		# Determine number of observed start and observed triplets in upstream sequence
		ORF.upstream_longest_orf_starts += len(self.orfs)
		ORF.upstream_longest_orf_triplets += len(self.longest_orf_up_seq)/3
		
				
	def get_codon_search_seqs(self):
		up_seq,down_seq,end_seq = "","",""
		self.start -= 1
		end_seq_window = self.codon_search_window
		if (self.stop-self.start)<end_seq_window:
			end_seq_window = self.stop-self.start
		
		# Get the codon search sequence
		if self.strand == "+":					
			up_seq = self.genome.genome_seq[self.start-self.codon_search_window:self.start]
			down_seq = self.genome.genome_seq[self.start+3:self.start+self.codon_search_window+3]
			end_seq = self.genome.genome_seq[self.stop-end_seq_window:self.stop]
		if self.strand == "-":
			start_non_reversed = self.start
			self.start = self.stop
			up_seq = reverse_sequence(self.genome.genome_seq[self.stop:self.stop+self.codon_search_window])
			down_seq = reverse_sequence(self.genome.genome_seq[self.stop-self.codon_search_window-3:self.stop-3])
			end_seq = reverse_sequence(self.genome.genome_seq[start_non_reversed:start_non_reversed+end_seq_window])
		return up_seq,down_seq,end_seq

		
	def find_candidate_starts(self,direction,sequence,no_stop=True):
		relative_position = None
		alternative_starts = []
		codon = None
		if direction == "backward":
			for i in range(len(sequence)-3,0,-3):
				if sequence[i:i+3] in self.start_codons:
					codon = sequence[i:i+3]
					relative_position = i-len(sequence)
					absolute_position = relative_position+self.start
					if self.strand == "-":
						absolute_position = self.start-relative_position
					alternative_starts.append([relative_position,absolute_position,codon])
				if sequence[i:i+3] in self.stop_codons and no_stop:
					break
		if direction == "forward":
			for i in range(0,len(sequence),3):
				if sequence[i:i+3] in self.start_codons:
					codon = sequence[i:i+3]
					relative_position = i
					absolute_position = relative_position+self.start+3
					relative_position += 3
					if self.strand == "-":
						absolute_position = self.start-relative_position
					alternative_starts.append([relative_position,absolute_position,codon])
				if sequence[i:i+3] in self.stop_codons and no_stop:
					break
		return alternative_starts
	
		
	def get_longest_orf_seq(self):
		# Get the longest possible ORF
		if len(self.candidate_starts_upstream)>0:
			longest_orf = self.candidate_starts_upstream[-1][1]
		else:
			longest_orf = self.start	
		# Exclude the first 30 nt
		if self.strand == "+":
			longest_orf_up_seq = self.genome.genome_seq[longest_orf-self.codon_search_window:longest_orf]
		if self.strand == "-":
			longest_orf_up_seq = reverse_sequence(self.genome.genome_seq[longest_orf:longest_orf+self.codon_search_window])	
		return longest_orf_up_seq
		
	
def main():
	input_file,fasta_file,output_name = get_args()
	genome_TIS_data = Genome(input_file,fasta_file,output_name)
	if plot_result:
		plot_alternative_start_distribution(genome_TIS_data)
		

if __name__ == "__main__":
	main()
        
