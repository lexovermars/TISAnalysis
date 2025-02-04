import os
import sys
import getopt

def main(argv):
	inputfile = ''
	output_name = ''
	fasta_file = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:f:",["input_file=","out_name=","fasta_file="])
	except getopt.GetoptError:
		print 'usage: run_pca_pipeline.py -i <inputpttfile> -f <fastafile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'usage: run_pca_pipeline -i <inputfile> -f <fastafile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--input_file"):
			inputfile = arg
		elif opt in ("-o", "--out_name"):
			output_name = arg
                elif opt in ("-f", "--fasta_file"):
                        fasta_file = arg
	if len(inputfile) == 0 or len(output_name) == 0 or len(fasta_file) == 0:
		print 'usage: run_pca_pipeline.py -i <inputfile> -f <fastafile> -o <out_name>'
		sys.exit()
	return inputfile,fasta_file,output_name
				
def perform_pipeline(input_file,fasta_file,out_name):
	pca_command = "python pca_analysis.py -i "+str(input_file)+" -f "+fasta_file+" -o "+ out_name
	print pca_command
	pca_run = os.system(pca_command)
	print "PCA iterations done.."		
	score_analysis_command = "python analyze_create_new_annotation.py -i "+str(input_file)+" -o " + out_name
	score_run = os.system(score_analysis_command)	
	print "Analyzing PCA results done.."
	print "Generated output matrix with 5 best scoring potential TISs for each ORF.. " 
	print "Generated output table with best scoring TISs.."


if __name__ == "__main__":
	input_file,fasta_file,output_name = main(sys.argv[1:])	
	perform_pipeline(input_file,fasta_file,output_name)
