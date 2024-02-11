### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
	python rRNA_check.py
	--fastq <FASTQ_INPUT_FILE>
	--rRNA <rRNA_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	optional:
	--kmer <KMER_SIZE>[21]
	--cutoff <NUM_OF_HITS>[3]
	
	bug reports and feature requests: b.pucker@tu-bs.de
					"""

import sys, os, subprocess, gzip
import matplotlib.pyplot as plt

# --- end of imports --- #

def generate_kmers( seq, kmer ):
	"""! @brief generate k-mers """
	
	return [ seq[ i:i + kmer ] for i in range( 0, len( seq ), kmer ) ]


def revcomp( inputseq ):
	"""! @brief generate reverse complement of sequences """
	
	revcomp_kmers = []
	complements = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }
	newseq = []
	for nt in inputseq:
		try:
			newseq.append( complements[nt] )
		except KeyError:
			newseq.append( "N" )
	return "".join( newseq[::-1] )


def load_rRNA_kmers( rRNA_fasta_file, kmer ):
	"""! @brief load k-mers from given sequences """
	
	kmers = []
	with open( rRNA_fasta_file, "r" ) as f:
		f.readline().strip()[1:]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				kmers += generate_kmers( "".join( seq ).upper(), kmer )
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		kmers += generate_kmers( "".join( seq ).upper(), kmer )
	
	clean_kmers = []
	for each in kmers:
		if len( each ) == kmer:
			clean_kmers.append( each )
	revcomp_kmers = []
	for each in clean_kmers:
		revcomp_kmers.append( revcomp( each ) )
	
	clean_kmers = list( set( clean_kmers+revcomp_kmers ) )
	
	return clean_kmers


def generate_summary_figure( data_documentation_file, summary_figure ):
	"""! @brief load data from summary file and generate figure """
	
	# --- load data --- #
	labels, values = [], []
	with open( data_documentation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			ID = parts[0]
			if "/" in ID:
				ID = ID.split('/')[-1]
			if '.fastq.gz' in ID:
				ID = ID.replace( ".fastq.gz", "" )
			labels.append( ID )
			values.append( int( parts[1] ) )
			line = f.readline()
	
	# --- generate figure --- #
	fig, ax = plt.subplots()
	
	ax.bar( labels, values )
	
	ax.set_xticklabels( labels, rotation=90)
	
	ax.set_ylabel( "rRNA content [# reads]" )
	
	plt.tight_layout()
	
	fig.savefig( summary_figure, dpi=300 )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	rRNA_fasta_file = arguments[ arguments.index('--rRNA')+1 ]
	fastq_file = arguments[ arguments.index('--fastq')+1 ]
	if "," in fastq_file:
		fastq_files = fastq_file.split(',')
	else:
		fastq_files = [ fastq_file ]
		
	if not output_folder[-1] == "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--kmer' in arguments:
		kmer = int( arguments[ arguments.index('--kmer')+1 ] )
	else:
		kmer = 21
	
	if '--cutoff' in arguments:
		cutoff = int( arguments[ arguments.index('--cutoff')+1 ] )
	else:
		cutoff = 3

	rRNA_kmers = load_rRNA_kmers( rRNA_fasta_file, kmer )
	
	data_documentation_file = output_folder + "data_documentation.txt"
	if not os.path.isfile( data_documentation_file ):
		with open( data_documentation_file, "w" ) as out:
			for idx, fastq in enumerate( fastq_files ):
				sys.stdout.write( "processing: " + fastq +"\n" )
				sys.stdout.flush()
				extension = fastq.split('.')[-1]
				hits = []
				if extension in [ "FQ", "fq", "FASTQ", "fastq" ]:	#uncompressed file
					with open( fastq, "r" ) as f:
						line = f.readline()	#header
						while line:
							seq = f.readline()
							matches = []
							for k in rRNA_kmers:
								if k in seq:
									matches.append( k )
							unique_matches = list( set( matches ) )
							if len( unique_matches ) > cutoff:
								hits.append( seq )
							f.readline()	#useless line
							f.readline()	#quality line
							line = f.readline()	#header
				else:	#compressed input file
					with gzip.open( fastq, "rt" ) as f:
						line = f.readline()	#header
						while line:
							seq = f.readline()
							matches = []
							for k in rRNA_kmers:
								if k in seq:
									matches.append( k )
							unique_matches = list( set( matches ) )
							if len( unique_matches ) > cutoff:
								hits.append( seq )
							f.readline()	#useless line
							f.readline()	#quality line
							line = f.readline()	#header
				sys.stdout.write( fastq + ' --- hits: ' + str( len( hits ) ) + "\n" )
				sys.stdout.flush()
				out.write( fastq + "\t" + str( len( hits ) ) + "\n" )
	
	# --- generate summary figure --- #
	summary_figure = output_folder + "summary_figure.png"
	generate_summary_figure( data_documentation_file, summary_figure )


if '--fastq' in sys.argv and '--rRNA' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
