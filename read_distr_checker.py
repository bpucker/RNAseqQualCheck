### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.16 ###

#coverage file construction taken from: Pucker & Brockington, 2018: https://doi.org/10.1186/s12864-018-5360-z

__usage__ = """
	python read_distr_checker.py
	--bam <BAM_INPUT_FILE> | --cov <COVERAGE_INPUT_FILE>
	--gff <GFF_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	optional:
	--bam_is_sorted <PREVENTS_SORTING_OF_BAM>[False]
	--sample <NAME_OF_EACH_SAMPLE>
	--samtools <SAMTOOLS_PATH>[samtools]
	--bedtools <BED_TOOLS_PATH>[genomeCoverageBed]
	--chunks <NUMBER_OF_CHUNKs>[100]
	--minexpcut <MIN_CUMULATIVE_COVERAGE_TO_CONSIDER_TRANSCRIPT>[100]
	
	bug reports and feature requests: b.pucker@tu-bs.de
					"""

import sys, os, subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- end of imports --- #

def construct_cov_file( bam_file, cov_file, samtools, bedtools, bam_sorted ):
	"""! @brief construct a coverage file """
	
	if bam_sorted:
		sorted_bam_file = bam_file
	else:
		sys.stdout.write( "sorting BAM file ...\n")
		sys.stdout.flush()
		sorted_bam_file = cov_file.split('.cov')[0] + ".sorted.bam"
		cmd = samtools + " sort -m 5000000000 --threads 8 " + bam_file + " > " + sorted_bam_file
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
	
	# --- calculate read coverage depth per position --- #
	sys.stdout.write( "calculating coverage per position ....\n" )
	sys.stdout.flush()
	cmd = bedtools + " -d -split -ibam " + sorted_bam_file + " > " + cov_file
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()


def load_cov_from_file( cov_file ):
	"""! @brief load content of coverage file into dictionary """
	
	coverage = {}	#sequences are keys and coverage values are stored in lists of values
	with open( cov_file, "r" ) as f:
		first_line = f.readline().strip().split('\t')
		prev_seq = first_line[0]
		pos_counter = 1
		covs = []
		while pos_counter < int( first_line[1] ):
			covs.append( 0 )
			pos_counter += 1
		covs.append( float( first_line[1] ) )
		pos_counter += 1
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_seq:	#checks if moving from one sequence to next
				coverage.update( { prev_seq: covs } )
				prev_seq = parts[0] + ""
				covs = []
				pos_counter = 1
				while pos_counter < float( parts[1] ):
					covs.append( 0 )
					pos_counter += 1
				covs.append( float( parts[2] ) )
				pos_counter += 1
			else:
				while pos_counter < float( parts[1] ):
					covs.append( 0 )
					pos_counter += 1
				covs.append( float( parts[2] ) )
				pos_counter += 1
			line = f.readline()
		coverage.update( { prev_seq: covs } )	#adds the coverage values of last sequence
		
	return coverage


def load_transcript_structures_from_gff( gff_file, id_tag="ID" ):
	"""! @brief load exon ranges from given GFF file """
	
	transcript_structures = {}	#dictionary with transcripts as key; dictionaries as values: exon borders + orientation + chromosome name
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA", "transcript" ]:
					ID = parts[-1].split( id_tag + "=" )[1]
					if ";" in ID:
						ID = ID.split(';')[0]
					transcript_structures.update( { ID: { 'pos': [], 'chr': parts[0], 'orientation': parts[6] } } )
				elif parts[2] == "exon":
					parent = parts[-1].split( "Parent=" )[1]
					if ";" in parent:
						parent = parent.split(';')[0]
					try:
						transcript_structures[ parent ]['pos'].append( [ int( parts[3] ), int( parts[4] ) ] )
					except KeyError:
						sys.stdout.write( "EXON-ERROR: " + line )
						sys.stdout.flush()
						pass
			line = f.readline()
	return transcript_structures


def get_cov_values_per_transcript( coverage, transcript_structures ):
	"""! @brief collect all coverage values per transcript in one list """
	
	covs_per_transcript = {}
	for trans in list( transcript_structures.keys() ):
		try:
			# --- append coverage values of all exon positions to one list --- #
			cov_collection = []
			chromosome_coverage = coverage[ transcript_structures[ trans ][ 'chr' ] ]
			for exon in transcript_structures[ trans ][ 'pos' ]:
				for i in range( exon[0]-1, exon[1] ):	#generates an index list for one exon
					cov_collection.append( chromosome_coverage[ i ] )
		
			# --- check orientation and flip exon order if necessary --- #
			if transcript_structures[ trans ][ 'orientation' ] == '-':
				cov_collection = cov_collection[::-1]
			covs_per_transcript.update( { trans: cov_collection } )
		except:
			pass
	return covs_per_transcript


def write_cov_to_file( cov_input, out_file ):
	"""! @brief write coverage values per transcript into output file """
	
	with open( out_file, "w" ) as out:
		for key in list( sorted( cov_input.keys() ) ):
			values = cov_input[ key ]
			out.write( key + "\t" + ",".join( list( map( str, values ) ) ) + "\n" )


def normalize_cov_per_transcript( coverages_per_transcript, chunks, minexpcut ):
	"""! @brief normalize coverage by analyzing per chunk """
	
	normalized_cov_per_transcript = {}
	rel_normalized_cov_per_transcript = {}
	
	for key in list( coverages_per_transcript.keys() ):
		values = coverages_per_transcript[ key ]
		if len( values ) >= chunks:	#only consider trascripts with more positions than chunks
			val_per_chunk = [ ]
			for i in range( chunks ):	#generates one empty list per chunk in larger list
				val_per_chunk.append( [] )
			factor = len( values ) / float( chunks )
			for idx, val in enumerate( values ):
				cidx = int( idx/factor ) #chunk index to a append value to correct inner list
				val_per_chunk[ cidx ].append( val )
			avg_per_chunk = []	#calculate average value per chunk
			for each in val_per_chunk:
				if sum( each ) == 0:	#if there are no values, the average is 0
					avg_per_chunk.append( 0 )
				else:	#if there is at least one value, average calculation is possible
					avg_per_chunk.append( sum( each ) / float( len( each ) ) )
			normalized_cov_per_transcript.update( { key: avg_per_chunk } )
			if sum( avg_per_chunk ) > minexpcut:
				rel_avg_per_chunk = []	#generating values between 0 and 1 for all chunks
				normalizer_per_transcript = max( avg_per_chunk )
				for each in avg_per_chunk:
					rel_avg_per_chunk.append( each / normalizer_per_transcript )
				rel_normalized_cov_per_transcript.update( { key: rel_avg_per_chunk } )
	
	#Example:
	#chunks = 100
	#length = 250
	#factor = 250/100 = 2.5
	#int( idx/factor ) = cindx
		
	return normalized_cov_per_transcript, rel_normalized_cov_per_transcript


def summarize_across_transcripts( normalized_cov_per_transcript, summary_data_output_file, summary_fig_outout_file, chunks ):
	"""! @brief summarize data across all transcripts """
	
	# --- summarize all values --- #
	summary_per_position = [ ] 
	for i in range( chunks ):
		summary_per_position.append( [] )
	
	for key in list( sorted( normalized_cov_per_transcript.keys() ) ):
		values = normalized_cov_per_transcript[ key ]
		for idx, val in enumerate( values ):
			summary_per_position[ idx ].append( val )
	
	# --- write summary into file --- #
	with open( summary_data_output_file, "w" ) as out:
		for idx, vals in enumerate( summary_per_position ):
			out.write( str( idx ) + "\t" + ",".join( list( map( str, vals ) ) ) + "\n" )
	
	# --- generate summary figure --- #
	sns.set(style="whitegrid")
	fig, ax = plt.subplots()
	sns.violinplot(data=summary_per_position, ax=ax)
	ax.set(xlabel='position in transcript', ylabel='RNA-seq coverage', title='RNA-seq coverage across transcripts')
	
	fig.savefig( summary_fig_outout_file, dpi=300 )
	
	median_per_position = []
	for each in summary_per_position:
		median_per_position.append( np.median( each ) )
	return median_per_position


def generate_comparative_plot( collected_data, final_fig_file, chunks ):
	"""! @brief generate comparative plot """
	
	fig, ax = plt.subplots()
	
	x_values = list( range( 1, chunks+1 ) )
	for key in list( sorted( collected_data.keys() ) ):
		y_values = collected_data[ key ]
		ax.plot( x_values, y_values, label=key )
	
	ax.legend()
	ax.set_xlabel( "position in transcript" )
	ax.set_ylabel( "RNA-seq read coverage" )
	
	fig.savefig( final_fig_file, dpi=300 )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	
	if not output_folder[-1] == "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--samtools' in arguments:
		samtools = arguments[ arguments.index('--samtools')+1 ]
	else:
		samtools = "samtools"
	
	if '--bedtools' in arguments:
		bedtools = arguments[ arguments.index('--bedtools')+1 ]
	else:
		bedtools = "genomeCoverageBed"
	
	if '--bam_is_sorted' in arguments:
		bam_sorted = True
	else:
		bam_sorted = False
	
	if '--chunks' in arguments: #numer of chunks per transcript to analyze coverage
		chunks=int( arguments[ arguments.index('--chunks')+1 ] )
		chunks = max( 10, chunks )	#check that this is not smaller than 10
	else:
		chunks=100
	
	if '--minexpcut' in arguments: #how many bases of a transcript were sequenced in total (sum of coverage of individual positions)
		minexpcut = int( arguments[ arguments.index('--minexpcut')+1 ] )
	else:
		minexpcut=100
	
	if '--sample' in arguments:
		sample = arguments[ arguments.index('--sample')+1 ]
		if "," in sample:
			samples = sample.split(',')
	else:
		if '--bam' in arguments:
			number = arguments[ arguments.index('--bam')+1 ].count(',')+1
			samples = list( map( str, range( number ) ) )
		else:
			number = arguments[ arguments.index('--cov')+1 ].count(',')+1
			samples = list( map( str, range( number ) ) )
	
	if '--bam' in arguments:
		bam_file = arguments[ arguments.index('--bam')+1 ]	#path to BAM input file
		#generate coverage file / also allow start from coverage file
		if "," in bam_file:
			bam_files = bam_file.split(',')
		else:
			bam_files = [ bam_file ]
		cov_files = []
		for idx, bam_file in enumerate( bam_files ):
			cov_file = output_folder + samples[idx] + ".cov"
			if not os.path.isfile( cov_file ):
				sys.stdout.write( "constructing coverage file "+ samples[ idx ] +" ...\n" )
				sys.stdout.flush()
				construct_cov_file( bam_file, cov_file, samtools, bedtools, bam_sorted  )
				sys.stdout.write( "...done\n" )
				sys.stdout.flush()
			cov_files.append( cov_file )
	else:
		cov_file = arguments[ arguments.index('--cov')+1 ]
		if "," in cov_file:
			cov_files = cov_file.split(',')
		else:
			cov_files = [ cov_file ]
	

	# --- load gene structures from GFF (exon ranges of representative transcript) --- #
	sys.stdout.write( "loading transcript structures ...\n" )
	sys.stdout.flush()
	transcript_structures = load_transcript_structures_from_gff( gff_file )
	sys.stdout.write( "...done\n" )
	sys.stdout.flush()
	
	# --- iteration over all samples --- #
	collected_data = {}
	for oidx, cov_file in enumerate( cov_files ):	#iterate over all coverage files that have been created in the previous step
		# --- load coverage from file --- #
		sys.stdout.write( "loading coverage from COV file "+ samples[ oidx ] +" ...\n" )
		sys.stdout.flush()
		coverage = load_cov_from_file( cov_file )
		sys.stdout.write( "...done\n" )
		sys.stdout.flush()
		
		#collect coverage per position per transcript
		sys.stdout.write( "collecting coverage per transcript "+ samples[ oidx ] +" ...\n" )
		sys.stdout.flush()
		coverages_per_transcript = get_cov_values_per_transcript( coverage, transcript_structures )
		sys.stdout.write( "...done\n" )
		sys.stdout.flush()
		
		# --- write coverage per transcript into output file --- #
		cov_per_transcript_out_file = output_folder + samples[ oidx ] + ".cov_per_transcript.txt"
		write_cov_to_file( coverages_per_transcript, cov_per_transcript_out_file )
				
		#split each transcript into X chunks and calculate average coverage for each of them
		normalized_cov_per_transcript, rel_normalized_cov_per_transcript = normalize_cov_per_transcript( coverages_per_transcript, chunks, minexpcut )
		
		# --- write normalized coverage per transcript into output file --- #
		norm_cov_per_transcript_out_file = output_folder + samples[ oidx ] + ".norm_cov_per_transcript.txt"
		write_cov_to_file( normalized_cov_per_transcript, norm_cov_per_transcript_out_file )
				
		# --- write relative normalized coverage per transcript into output file (all value between 0 and 1) --- #
		rel_norm_cov_per_transcript_out_file = output_folder + samples[ oidx ] + ".rel_norm_cov_per_transcript.txt"
		write_cov_to_file( rel_normalized_cov_per_transcript, rel_norm_cov_per_transcript_out_file )
				
		#summarize the coverage distribution over corresponding chunks in different transcripts
		summary_data_output_file = output_folder + samples[ oidx ] + ".summary.txt"
		summary_fig_outout_file = output_folder + samples[ oidx ] + ".summary.png"
		summary = summarize_across_transcripts( rel_normalized_cov_per_transcript, summary_data_output_file, summary_fig_outout_file, chunks )
		
		# --- store summary information --- #
		collected_data.update( { samples[ oidx ]: summary } )
	
	# --- generate final comparative figure --- #
	final_fig_file = output_folder + "comparative_plot.png"
	generate_comparative_plot( collected_data, final_fig_file, chunks )


if '--bam' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
