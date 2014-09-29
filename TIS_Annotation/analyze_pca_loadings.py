import pickle
import scipy
import operator
import numpy
from operator import itemgetter
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



# Load the iterative pca results (scores+ranks)
pca_loading_dic_file = open("script_output/"+outputfile_name+"_pca_loadings_pickle_dump.txt","r")
pca_loading_dic = pickle.load(pca_loading_dic_file)
pca_loading_dic_file.close()


# Create an output file handler
output_file = open("final_output/"+outputfile_name+"_formatted_pca_loadings.txt","w")
# Output: PC1 loading round 1 \t PC2 loading round 1 \t PC3 loading round 1 \t PC1 loading round 2......

for nucleotide_id in pca_loading_dic.keys():
	pca_rounds = pca_loading_dic[nucleotide_id]
	output_file.write(nucleotide_id)
	for pca_round in pca_rounds:
		for pca_loading in pca_round:
			output_file.write("\t"+str(pca_loading).strip())
	output_file.write("\n")

output_file.close()