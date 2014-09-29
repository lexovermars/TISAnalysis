import pickle
import scipy
import operator
import numpy
from operator import itemgetter
import sys, getopt

base_threshold = 0
complete_and_upstream = False


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

genome = inputfile

def make_ptt_dict(ptt_file):
	ptt_dict = {}
	ptt_data = open(ptt_file,"r")
	ptt_lines = ptt_data.readlines()
	for index,line in enumerate(ptt_lines[3:len(ptt_lines)]):
		line = line.split("\t")
		location = [int(line[0].split("..")[0]),int(line[0].split("..")[1])]
		ptt_dict[line[5]] = [location,line[1],line[7],line[8]]
	return ptt_dict


ptt_dict = make_ptt_dict(genome)

# Load the iterative pca results (scores+ranks)
result_dic_file = open("script_output/"+outputfile_name+"_pc_scores_pickle_dump.txt","r")
sorted_pca_scores = pickle.load(result_dic_file)
result_dic_file.close()

if complete_and_upstream == True:
	try:
		upstream_result_dic_file = open("script_output/"+outputfile_name+"_pc_scores_pickle_dump_upstream.txt","r")
		upstream_sorted_pca_scores = pickle.load(upstream_result_dic_file)
		upstream_result_dic_file.close()
	except:
		print "Could not find upstream data file... skipping upstream data"

# Define the output matrix
output_new_annotation = open("final_output/"+outputfile_name+"_adjusted_annotation.txt","w")
output = open("final_output/"+outputfile_name+"_5_rows_iteration_matrix.txt","w")
output_discarded = open("final_output/"+outputfile_name+"_matrix_discarded.txt","w")

# First determine the topscores for each alternative start (PC1, complete PCA)
top_scores = {}
summed_ranks = {}
for orf in sorted(sorted_pca_scores.keys()):
	for alt_start in sorted_pca_scores[orf].keys():
		score_data = sorted_pca_scores[orf][alt_start]
		#print score_data
		top_score = sorted(score_data, key=lambda score: score[2])[-1][2]
		if orf not in top_scores.keys():
			top_scores[orf] = {}
			summed_ranks[orf] = {}
		summed_rank = 0
		for round_index,score in enumerate(score_data):
			if round_index > 1:
				summed_rank += score[3]
		top_scores[orf][alt_start] = top_score
		summed_ranks[orf][alt_start] = summed_rank
		#print top_score
		
for orf in top_scores.keys():
	# sort candidates based on top scores
	#sorted_orf_candidates = sorted(top_scores[orf].iteritems(), key=operator.itemgetter(1), reverse=True)
	sorted_orf_candidates_rank = sorted(summed_ranks[orf].iteritems(), key=operator.itemgetter(1))
	
	
	#write topscore to new annotation file
	i = 0
	candidate_pos, top_score = sorted_orf_candidates_rank[i]
	alt_start_data = sorted_pca_scores[orf][candidate_pos]
	if complete_and_upstream == True:
		upstream_alt_start_data = upstream_sorted_pca_scores[orf][candidate_pos]
	#First write the start info (locus, label and position)
	output_new_annotation.write(orf+"\t"+str(alt_start_data[0][1])+"\t"+str(alt_start_data[0][4])+"\t"+str(candidate_pos)+"\t"+str(i))
	for round_index,score_data in enumerate(alt_start_data):
		output_new_annotation.write("\t"+str(score_data[0])+"\t"+str(round(score_data[2],3))+"\t"+str(score_data[5])+"\t"+str(score_data[6])+"\t"+str(score_data[3]))
		if complete_and_upstream == True:
			score_data_upstream = upstream_alt_start_data[round_index]
			output_new_annotation.write("\t"+str(round(score_data_upstream[2],3))+"\t"+str(score_data_upstream[5])+"\t"+str(score_data_upstream[6])+"\t"+str(score_data_upstream[3]))
	output_new_annotation.write("\n")

		
	#write top 5 scores to matrix
	for i in range(5):
		try:
			candidate_pos, top_score = sorted_orf_candidates_rank[i]
			alt_start_data = sorted_pca_scores[orf][candidate_pos]
			if complete_and_upstream == True:
				upstream_alt_start_data = upstream_sorted_pca_scores[orf][candidate_pos]
			#First write the start info (locus, label and position)
			output.write(orf+"\t"+str(alt_start_data[0][1])+"\t"+str(alt_start_data[0][4])+"\t"+str(candidate_pos)+"\t"+str(i))
			for round_index,score_data in enumerate(alt_start_data):
				output.write("\t"+str(score_data[0])+"\t"+str(round(score_data[2],3))+"\t"+str(score_data[5])+"\t"+str(score_data[6])+"\t"+str(score_data[3]))
				if complete_and_upstream == True:
					score_data_upstream = upstream_alt_start_data[round_index]
					output.write("\t"+str(round(score_data_upstream[2],3))+"\t"+str(score_data_upstream[5])+"\t"+str(score_data_upstream[6])+"\t"+str(score_data_upstream[3]))
			
			output.write("\n")
		except:
			output.write(orf+"\n")
			#print "Only ",len(sorted_orf_candidates_rank),"scores"
	
	
	#write all discarded scores
	if len(sorted_orf_candidates_rank)>5:
		for i in range(5,len(sorted_orf_candidates_rank)):
			try:
				candidate_pos, top_score = sorted_orf_candidates_rank[i]
				alt_start_data = sorted_pca_scores[orf][candidate_pos]
				if complete_and_upstream == True:
					upstream_alt_start_data = upstream_sorted_pca_scores[orf][candidate_pos]
				#First write the start info (locus, label and position)
				output_discarded.write(orf+"\t"+str(alt_start_data[0][1])+"\t"+str(alt_start_data[0][4])+"\t"+str(candidate_pos)+"\t"+str(i))
				for round_index,score_data in enumerate(alt_start_data):
					output_discarded.write("\t"+str(score_data[0])+"\t"+str(round(score_data[2],3))+"\t"+str(score_data[5])+"\t"+str(score_data[6])+"\t"+str(score_data[3]))
					if complete_and_upstream == True:
						score_data_upstream = upstream_alt_start_data[round_index]
						output_discarded.write("\t"+str(round(score_data_upstream[2],3))+"\t"+str(score_data_upstream[5])+"\t"+str(score_data_upstream[6])+"\t"+str(score_data_upstream[3]))
			
				output_discarded.write("\n")
			except:
				output_discarded.write(orf+"\n")
				#print "Only ",len(sorted_orf_candidates_rank),"scores"			
			

output.close()
output_discarded.close()