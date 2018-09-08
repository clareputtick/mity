#!/Users/clareputtick/anaconda/bin/python3
#
# mity-merge
#
# Usage:
# merge.py -h
#
# Notes
# This code assumes that the VCFs are of the form specified in 
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# i.e. VCF v4.2. In particular it assumes that
# 
# If the nuclear VCF and the mity VCF are different versions (e.g. 4.1 and 
# 4.2), the default behaviour is to report a warning, but set the resulting 
# VCF to the nuclear VCF version. This behaviour can be changes with TODO
#
# The resulting VCF will be the nuclear VCF name, with a mity suffix added.

import sys
import gzip
import argparse
import subprocess

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='''Run mity-merge. Merges
		a nuclear VCF and a mity VCF, such that the MT variants in the nuclear
		VCF are replaces with the mity variants, and the headers are merged.''')
	parser.add_argument('--mity', action='store', help = 'mity VCF', required=True)
	parser.add_argument('--nuclear', action='store', help = 'Nuclear VCF', required=True)

	args = parser.parse_args()

	nuclear_vcf = args.nuclear
	mity_vcf = args.mity
	new_vcf_prefix = nuclear_vcf.split(".vcf.gz")[0]
	new_vcf_name = new_vcf_prefix + ".mity.vcf"
	
	nuclear_file = gzip.open(nuclear_vcf, 'rt')

	# split the header and the variants into two seperate lists
	nuclear_header = []
	nuclear_variants = []
	for line in nuclear_file:
		if line[0] == "#":
			line = line.strip()
			nuclear_header.append(line)
		else:
			line = line.strip()
			# line = line.split('\t')
			nuclear_variants.append(line)
	
	nuclear_col_names = nuclear_header[-1]
	del nuclear_header[-1]

	mity_file = gzip.open(mity_vcf, 'rt')

	# split the header and the variants into two seperate lists
	# TODO: to speed up, do we actually need to get nuclear variants?
	mity_header = []
	mity_variants = []
	for line in mity_file:
		if line[0] == "#":
			line = line.strip()
			mity_header.append(line)
		else:
			line = line.strip()
			line = line.split('\t')
			mity_variants.append(line)
	
	mity_col_names = mity_header[-1]
	del mity_header[-1]

	# print(mity_header)
	# print(nuclear_header)

	# remove file format (e.g. 4.1/4.2) for mity and nuclear add manually to the merged header
	del mity_header[0]
	del nuclear_header[0]

	# remove phasing from mity header
	mity_header = [x for x in mity_header if "phasing=" not in x]

	# remove file date from mity header
	mity_header = [x for x in mity_header if "##fileDate" not in x]

	# get mity reference
	mity_ref = [x for x in mity_header if "##reference=" in x]
	
	# should only be one line in the mity header that has reference in it
	if len(mity_ref) > 1:
		err_msg = "Mity header has more than lines that define a reference: \n"
		for x in mity_ref:
			err_msg = err_msg + x
			err_msg = err_msg + '\n'
		sys.exit(err_msg)

	mity_ref = mity_ref[0].split('##reference=')[1]
	# print(mity_ref)
	# remove mity reference from header because we will add it with nuclear reference later
	mity_header = [x for x in mity_header if "##reference" not in x]

	# get nuclear reference then remove from header
	nuclear_ref = [x for x in nuclear_header if "##reference" in x]
	# should only be one line in the nuclear header that has reference in it
	if len(nuclear_ref) > 1:
		err_msg = "HC header has more than lines that define a reference: \n"
		for x in nuclear_ref:
			err_msg = err_msg + x
			err_msg = err_msg + '\n'
		sys.exit(err_msg)

	nuclear_ref = nuclear_ref[0].split('##reference=')[1]
	# remove nuclear reference from header because we will add it with mity reference later
	nuclear_header = [x for x in nuclear_header if "##reference" not in x]

	# make new reference line and add to merged header
	new_ref_line = "##reference=If CHR=MT: " + mity_ref + ". Otherwise: " + nuclear_ref
	
	# remove all the lines in the mity header that are also in the HC header
	# this should remove contig lines if they are the same
	mity_header = [x for x in mity_header if x not in nuclear_header]

	# now if there are two lines in mity_header and nuclear_header that 
	# that have the same eg "##INFO=<ID=SRR,", then they have different definitions for the same 
	# ID. need to make nuclear and mity specific definitions.
	sep = ","
	merged_header = []
	mity_ids = [x.split(sep)[0] for x in mity_header]
	# print(mity_ids)
	# loop through the nuclear_header
	for nuclear_line in nuclear_header:
		nuclear_id = nuclear_line.split(sep, 1)[0]
		if nuclear_id in mity_ids:
			# print(nuclear_line)
			mity_line_idx = mity_ids.index(nuclear_id)
			mity_line = mity_header[mity_line_idx]
			# print(mity_line)

			# check the number is the same - if not put "."
			mity_number = mity_line.split("Number=")[1]
			mity_number = mity_number.split(",")[0]
			# print(mity_number)

			nuclear_number = nuclear_line.split("Number=")[1]
			nuclear_number = nuclear_number.split(",")[0]
			# print(nuclear_number)

			if str(nuclear_number) == str(mity_number):
				new_number = nuclear_number
			else:
				new_number = "."

			# check the type is the same
			# if not will have to stop the app with an error
			# no way of saying unknown type
			# and they really should be the same anyway
			mity_type = mity_line.split("Type=")[1]
			mity_type = mity_type.split(",")[0]
			# print(mity_type)

			nuclear_type = nuclear_line.split("Type=")[1]
			nuclear_type = nuclear_type.split(",")[0]
			# print(nuclear_type)

			if str(nuclear_type) == str(mity_type):
				new_type = nuclear_type
			else:
				new_type = nuclear_type
				# print("")
				# sys.exit("TODO")

			nuclear_description = nuclear_line.split('Description="')[1]
			# remove the last ">
			nuclear_description = nuclear_description[:-2]
			# print(nuclear_description)


			mity_description = mity_line.split('Description="')[1]
			# remove the last ">
			mity_description = mity_description[:-2]
			# print(mity_description)

			new_header_line = nuclear_id + ",Number=" + new_number + ",Type=" + new_type + ",Description=\"If CHR=MT: \'" + mity_description + "\'. Otherwise: \'" + nuclear_description + "\' \">"
			# print(new_header_line)
			merged_header.append(new_header_line)

			# remove the line from the mity_header
			del mity_header[mity_line_idx]
			# print(len(mity_header))
		else:
			merged_header.append(nuclear_line)
		
	# now add the rest of the mity_header to the merge header
	merged_header = merged_header + mity_header
	merged_header = sorted(merged_header)
	merged_header.insert(0, new_ref_line)
	merged_header.insert(0, '##fileformat=VCFv4.1')
	


	# for line in merged_header:
	# 	print(line)

	#now check that the sample names are in the right order

	# get nuclear and mity sample names
	# assumes that sample names always start from 10th column
	nuclear_col_names = nuclear_col_names.split('\t')
	nuclear_samples = nuclear_col_names[9:]
	# print(nuclear_samples)
	
	mity_col_names = mity_col_names.split('\t')
	mity_samples = mity_col_names[9:]
	# print(mity_samples)

	if nuclear_samples == mity_samples:
		sys.stderr.write("samples in same order - can just add mity variants to nuclear vcf")
				# add the nuclear variants with the mity variants - will be ordered later
		nuclear_col_names = ('\t').join(nuclear_col_names)

		new_mity_variants = []
		for line in mity_variants:
			new_line = ('\t').join(line)
			new_mity_variants.append(new_line)

		new_vcf = merged_header + [nuclear_col_names] + new_mity_variants + nuclear_variants
		for line in new_vcf:
			print(line)
	else:
		# check that they have the same sampels
		if set(nuclear_samples) != set(mity_samples):
			sys.exit("VCFs have different sample names")
		mity_variants_no_sample = [x[:9] for x in mity_variants]	
		mity_variants_sample = [x[9:] for x in mity_variants]	
		# print(mity_variants_sample)
		# sys.stderr.write("samples not in the same order - will reorder mity variants")

		# get the order of the nuclear sampels
		# print("nuclear samples")
		# print(nuclear_samples)
		# print("mity samples")
		# print(mity_samples)
		# mity_samples_idx = [nuclear_samples.index(x) for x in mity_samples]
		mity_samples_idx = [mity_samples.index(x) for x in nuclear_samples]

		# print(mity_samples_idx)
		# sys.exit()
		mity_variants_reordered_sample = []
		for line in mity_variants_sample:
			mity_variants_reordered_sample.append([line[x] for x in mity_samples_idx])
		# print(mity_variants_reordered_sample)
		# sys.exit()
		# add the newly ordered samples back onto the mity_variants_no_sample
		new_mity_variants = []
		for i in range(0, len(mity_variants_no_sample)):
			new_line = mity_variants_no_sample[i] + mity_variants_reordered_sample[i]
			# print(new_line)
			new_line = ('\t').join(new_line)
			new_mity_variants.append(new_line)

		# add the nuclear variants with the mity variants - will be ordered later
		nuclear_col_names = ('\t').join(nuclear_col_names)
		new_vcf = merged_header + [nuclear_col_names] + new_mity_variants + nuclear_variants

		new_vcf_file = open(new_vcf_name, 'w')
		for line in new_vcf:
			# print(line)
			new_vcf_file.write(line + "\n")
		new_vcf_file.close()

		bgzip_call = 'bgzip ' + new_vcf_name
		subprocess.run(bgzip_call, shell=True )

		tabix_call = 'tabix ' + new_vcf_name + ".gz"
		subprocess.run(tabix_call, shell=True )
		




