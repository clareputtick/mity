"""
Run mity-merge.
This merges a nuclear VCF and a mity VCF, such that the MT variants in the
nuclear VCF are replaces with the mity variants, and the headers are merged
"""

# Notes
#
# If the nuclear VCF and the mity VCF are different versions (e.g. 4.1 and 
# 4.2), the default behaviour is to report a warning, but set the resulting 
# VCF to the nuclear VCF version. This behaviour can be changes with TODO
#
# The resulting VCF will be the nuclear VCF name, with a mity suffix added.

# @TODO
# Instead of removing mity date, why not keep it and append with mity? 
# The contigs in the header could be made more simple if we remove the MT
# contig from the nuclear header, and only add the MT contig from the mity 
# header. Difficulty is knowing which contig line is the MT (chrM, or MT), 
# depends on the reference used

import sys
import gzip
import subprocess

def do_merge(mity_vcf, hc_vcf, prefix=None):
    
    hc_vcf
    mity_vcf
    prefix = [hc_vcf.split(".vcf.gz")[0], prefix][prefix is not None]
    new_vcf_name = prefix + ".mity.vcf"
    
    hc_file = gzip.open(hc_vcf, 'rt')
    
    # split the nuclear VCF header and variants into separate lists
    hc_header = []
    hc_variants = []
    for line in hc_file:
        if line[0] == "#":
            line = line.strip()
            hc_header.append(line)
        else:
            line = line.strip()
            hc_variants.append(line)
    
    hc_col_names = hc_header[-1]
    del hc_header[-1]
    
    mity_file = gzip.open(mity_vcf, 'rt')
    
    # split the mity VCF header and the variants into separate lists
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
    
    # get the column names from the mity VCF
    mity_col_names = mity_header[-1]
    del mity_header[-1]
    
    # remove file format (e.g. 4.1/4.2) for mity and nuclear 
    # will add manually to the merged header
    del mity_header[0]
    del hc_header[0]
    
    # remove phasing from mity header
    mity_header = [x for x in mity_header if "phasing=" not in x]
    
    # remove file date from mity header
    mity_header = [x for x in mity_header if "##fileDate" not in x]
    
    # get mity reference
    mity_ref = [x for x in mity_header if "##reference=" in x]
    
    # should only be one line in the mity header that has reference in it
    if len(mity_ref) > 1:
        err_msg = '''mity VCF header has more than lines that define a reference: \n'''
        for x in mity_ref:
            err_msg = err_msg + x
            err_msg = err_msg + '\n'
        sys.exit(err_msg)
    
    mity_ref = mity_ref[0].split('##reference=')[1]
    
    # remove mity reference from header because we will add it with nuclear 
    # reference later
    mity_header = [x for x in mity_header if "##reference" not in x]
    
    # get nuclear reference then remove from header
    hc_ref = [x for x in hc_header if "##reference" in x]
    # should only be one line in the nuclear header that has reference in it
    if len(hc_ref) > 1:
        err_msg = '''Nuclear VCF header has more than lines that define a reference: \n'''
        for x in hc_ref:
            err_msg = err_msg + x
            err_msg = err_msg + '\n'
        sys.exit(err_msg)
    
    hc_ref = hc_ref[0].split('##reference=')[1]
    # remove nuclear reference from header because we will add it with mity 
    # reference later
    hc_header = [x for x in hc_header if "##reference" not in x]
    
    # make new reference line and add to merged header
    new_ref_line = ("##reference=If CHR=MT: " + mity_ref + ". Otherwise: " +
                    hc_ref)
    
    # remove all the lines in the mity header that are also in the nuclear 
    # header. This should remove contig lines if they are the same
    mity_header = [x for x in mity_header if x not in hc_header]
    
    # If there is a line in mity_header and hc_header that 
    # has the same eg "##INFO=<ID=SRR,", then they have different 
    # definitions for the same  ID (because otherwise if the lines were 
    # exactly the same they would have already been removed from the mity 
    # header) So we need to make nuclear and mity specific definitions.
    sep = ","
    merged_header = []
    mity_ids = [x.split(sep)[0] for x in mity_header]
    # print(mity_ids)
    # loop through the hc_header
    for hc_header_line in hc_header:
        hc_id = hc_header_line.split(sep, 1)[0]
        if hc_id in mity_ids:
            # print("hc_id in mity_ids")
            # print(hc_header_line)
            # mity_line_idx = mity_ids.index(hc_id)
            # mity_line = mity_header[mity_line_idx]
            # mity_line = [line for line in mity_header if line.startswith(
            # hc_id)][0]
            mity_header_idx = ([idx for idx, s in enumerate(mity_header)
                                if s.startswith(hc_id)][0])
            mity_header_line = mity_header[mity_header_idx]
            print(mity_header_line)
            
            # check the number is the same - if not put "."
            # as the VCF format specifies for unknown number
            mity_number = mity_header_line.split("Number=")[1]
            mity_number = mity_number.split(",")[0]
            # print(mity_number)
            
            hc_number = hc_header_line.split("Number=")[1]
            hc_number = hc_number.split(",")[0]
            # print(hc_number)
            
            if str(hc_number) == str(mity_number):
                new_number = hc_number
            else:
                new_number = "."
            
            # check the type is the same, if not will set it to the nuclear 
            # type as there is no way of saying unknown type in the VCF format
            mity_type = mity_header_line.split("Type=")[1]
            mity_type = mity_type.split(",")[0]
            # print(mity_type)
            
            hc_type = hc_header_line.split("Type=")[1]
            hc_type = hc_type.split(",")[0]
            # print(hc_type)
            
            if str(hc_type) == str(mity_type):
                new_type = hc_type
            else:
                new_type = hc_type
            
            hc_description = hc_header_line.split('Description="')[1]
            # remove the last ">
            hc_description = hc_description[:-2]
            # print(hc_description)
            
            mity_description = mity_header_line.split('Description="')[1]
            # remove the last ">
            mity_description = mity_description[:-2]
            # print(mity_description)
            
            new_header_line = (hc_id + ",Number=" + new_number +
                               ",Type=" + new_type + ",Description=\"If "
                                                     "CHR=MT: \'" +
                               mity_description + "\'. Otherwise: \'" + 
                               hc_description +
                               "\' \">")
            print(new_header_line)
            merged_header.append(new_header_line)
            
            # remove the line from the mity_header
            del mity_header[mity_header_idx]
            print(mity_header)
        else:
            merged_header.append(hc_header_line)
    
    # sys.exit()
    # now add the rest of the mity_header to the merge header
    merged_header = merged_header + mity_header
    merged_header = sorted(merged_header)
    merged_header.insert(0, new_ref_line)
    merged_header.insert(0, '##fileformat=VCFv4.1')
    
    for line in merged_header:
        print(line)
    # sys.exit()
    
    # now check that the sample names are in the right order
    
    # get nuclear and mity sample names
    # assumes that sample names always start from 10th column
    hc_col_names = hc_col_names.split('\t')
    hc_samples = hc_col_names[9:]
    # print(hc_samples)
    
    mity_col_names = mity_col_names.split('\t')
    mity_samples = mity_col_names[9:]
    # print(mity_samples)
    
    if hc_samples == mity_samples:
        sys.stderr.write(
            "samples in same order - can just add mity variants to nuclear vcf")
        # add the nuclear variants with the mity variants - will be ordered 
        # later
        hc_col_names = ('\t').join(hc_col_names)
        
        new_mity_variants = []
        for line in mity_variants:
            new_line = ('\t').join(line)
            new_mity_variants.append(new_line)
        
        new_vcf = merged_header + [
            hc_col_names] + new_mity_variants + hc_variants
        for line in new_vcf:
            print(line)
    else:
        # check that they have the same sampels
        if set(hc_samples) != set(mity_samples):
            sys.exit("VCFs have different sample names")
        mity_variants_no_sample = [x[:9] for x in mity_variants]
        mity_variants_sample = [x[9:] for x in mity_variants]
        # print(mity_variants_sample)
        # sys.stderr.write("samples not in the same order - will reorder mity variants")
        
        # get the order of the nuclear sampels
        # print("nuclear samples")
        # print(hc_samples)
        # print("mity samples")
        # print(mity_samples)
        # mity_samples_idx = [hc_samples.index(x) for x in mity_samples]
        mity_samples_idx = [mity_samples.index(x) for x in hc_samples]
        
        # print(mity_samples_idx)
        # sys.exit()
        mity_variants_reordered_sample = []
        for line in mity_variants_sample:
            mity_variants_reordered_sample.append(
                    [line[x] for x in mity_samples_idx])
        # print(mity_variants_reordered_sample)
        # sys.exit()
        # add the newly ordered samples back onto the mity_variants_no_sample
        new_mity_variants = []
        for i in range(0, len(mity_variants_no_sample)):
            new_line = mity_variants_no_sample[i] + \
                       mity_variants_reordered_sample[i]
            # print(new_line)
            new_line = ('\t').join(new_line)
            new_mity_variants.append(new_line)
        
        # add the nuclear variants with the mity variants - will be ordered later
        hc_col_names = ('\t').join(hc_col_names)
        new_vcf = merged_header + [
            hc_col_names] + new_mity_variants + hc_variants
        
        new_vcf_file = open(new_vcf_name, 'w')
        for line in new_vcf:
            # print(line)
            new_vcf_file.write(line + "\n")
        new_vcf_file.close()
        
        bgzip_call = 'bgzip ' + new_vcf_name
        subprocess.run(bgzip_call, shell=True)
        
        tabix_call = 'tabix ' + new_vcf_name + ".gz"
        subprocess.run(tabix_call, shell=True)
