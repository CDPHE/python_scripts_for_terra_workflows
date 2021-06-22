#! /usr/bin/env python


import argparse
import sys
import pandas as pd
import numpy as np
from datetime import date
import re
# import os


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument('--sample_list')
    parser.add_argument("--bam_file_list",  help= "txt file with list of bam file paths")
    parser.add_argument('--percent_cvg_file_list', help = 'txt file with list of percent cvg file paths')
    
    parser.add_argument('--pangolin_lineage_csv', help = 'csv output from pangolin')
    parser.add_argument('--pangolin_version')
 
    
    parser.add_argument('--nextclade_clades_csv', help = 'csv output from nextclade parser')
    parser.add_argument('--nextclade_variants_csv')
    parser.add_argument('--nextclade_version')
    
    parser.add_argument('--seq_run', help = 'seq_run name; e.g. COVSEQ_0000 or COVMIN_0000')
    options = parser.parse_args(args)
    return options


def concat_samtools(bam_file_list):

    # get the input paths from text file into a python list
    with open(bam_file_list, 'r') as f:
        file_paths = []
        for line in f:
            file_paths.append(line.strip())

     # initiate dataframe for concatenation
    df = pd.DataFrame()
    accession_id_list = []
    samtools_numreads_list = []
    samtools_depth_list = []
    samtools_baseq_list = []
    samtools_mapq_list = []

    #loop through bam file stats files and pull data
    for file_path in file_paths:
        d = pd.read_csv(file_path, sep = '\t')
        if re.search('barcode', file_path): 
            # for nanopore runs
            sequence_name = re.findall('/([0-9a-zA-Z_\-\.]+)_barcode', file_path)[0]
        else: 
            # for illumina runs
            sequence_name = re.findall('/([0-9a-zA-Z_\-\.]+)_coverage.txt', file_path)[0]
        
        # pull data from samtools output
        num_reads = d.numreads[0]
        depth = d.meandepth[0]
        baseq = d.meanbaseq[0]
        mapq = d.meanmapq[0]

        accession_id_list.append(sequence_name)
        samtools_numreads_list.append(num_reads)
        samtools_depth_list.append(depth)
        samtools_baseq_list.append(baseq)
        samtools_mapq_list.append(mapq)

    df['accession_id'] = accession_id_list
    df['num_reads'] = samtools_numreads_list
    df['mean_depth'] = samtools_depth_list
    df['mean_base_quality'] = samtools_baseq_list
    df['mean_map_quality'] = samtools_mapq_list

    return df

def concat_percent_cvg(percent_cvg_file_list):
    
    # get the input paths from text file into a python list
    with open(percent_cvg_file_list, 'r') as f:
        file_paths = []
        for line in f:
            file_paths.append(line.strip())
    
    df_list = []
    for file in file_paths:
        d = pd.read_csv(file, dtype = {'accession_id' : object})
        df_list.append(d)
   
    df = pd.concat(df_list)
    
    return df


def get_df_spike_mutations(variants_csv):

    def get_accession_id(fasta_header):
        accession_id = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return accession_id
    
    variants = pd.read_csv(variants_csv, dtype = {'accession_id' : object})
    
    variants = variants.rename(columns = {'accession_id' : 'fasta_header'})
    
    accession_id = variants.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
    variants.insert(value = accession_id, loc = 0, column = 'accession_id')
    variants = variants.drop(columns = 'fasta_header')
    
    crit = variants.gene == 'S'
    critRBD = (variants.codon_position >= 461) & (variants.codon_position <= 509)
    critPBCS = (variants.codon_position >= 677) & (variants.codon_position <= 694)
    crit732 = variants.codon_position == 732
    crit452 = variants.codon_position == 452
    crit253 = variants.codon_position == 253
    crit13 = variants.codon_position == 13
    critdel = variants.variant_name.str.contains('del')

    variants = variants[crit & (critRBD | critPBCS | crit732 | critdel | crit452 | crit253 | crit13)]

    accession_ids = variants.accession_id.unique().tolist()

    df = pd.DataFrame()
    accession_id_list = []
    variant_name_list = []

    seperator = '; '

    for accession_id in accession_ids:
        accession_id_list.append(accession_id)

        crit = variants.accession_id == accession_id
        f = variants[crit]
        f = f.reset_index()

        mutations = []
        for row in range(f.shape[0]):
            mutations.append(f.variant_name[row])
        mutations_string = seperator.join(mutations)
        variant_name_list.append(mutations_string)

    df['accession_id'] = accession_id_list
    df['spike_mutations'] = variant_name_list
    
    return df
    

def concat_results(sample_list, samtools_df, percent_cvg_df, spike_mut_df, nextclade_clades_csv, 
                   pangolin_lineage_csv, next_version, pangolin_version, seq_run):
    
    
    # get the list of samples and create a df using the list of samples
    all_sample_accession_ids = []
    with open(sample_list, 'r') as f:
        for line in f:
            all_sample_accession_ids.append(line.strip())
    
    print(all_sample_accession_ids)
    df = pd.DataFrame(all_sample_accession_ids)
    df = df.rename(columns = {0:'accession_id'})
    df = df.set_index('accession_id')
    print(df)

    def get_accession_id(fasta_header):
        accession_id = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return accession_id
    
    def create_fasta_header(accession_id):
        return 'CO-CDPHE-%s' % accession_id
    
    # read in nextclade clade results
    pangolin = pd.read_csv(pangolin_lineage_csv, dtype = {'taxon' : object})
    pangolin = pangolin.rename(columns = {'lineage': 'pangolin_lineage', 
                              'taxon' : 'fasta_header', 
                              'note': 'pangolin_notes', 
                              'conflict': 'pangolin_probability',
                              'status': 'pangolin_status'})
    
    accession_id = pangolin.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
    pangolin.insert(value = accession_id, column = 'accession_id', loc = 0)
    pangolin = pangolin.drop(columns = 'fasta_header')
    
    pangolin['pangolin_version'] = pangolin_version
    pangolin = pangolin.set_index('accession_id')
    
    # pull out the pangoLEARN_version....
    pango_learn_version = pangolin.pangoLEARN_version[0]
    
    
    
    # read in nextclade csv
    nextclade = pd.read_csv(nextclade_clades_csv, dtype = {'accession_id' : object})
    nextclade = nextclade.rename(columns = {'accession_id' : 'fasta_header'})
    
    
    accession_id = nextclade.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
    nextclade.insert(value = accession_id, column = 'accession_id', loc = 0)
    nextclade = nextclade.drop(columns = 'fasta_header')
    
    nextclade['nextclade_version'] = next_version
    nextclade = nextclade.set_index('accession_id')
    
    
    # set index on the samtools_df and percent_cvg_df and variants_df
    samtools_df = samtools_df.set_index('accession_id')
    percent_cvg_df = percent_cvg_df.set_index('accession_id')
    spike_mut_df = spike_mut_df.set_index('accession_id')
    
    
    # join
    j = df.join(percent_cvg_df, how = 'left')
    j = j.join(samtools_df, how = 'left')
    j = j.join(nextclade, how = 'left')
    j = j.join(pangolin, how = 'left')
    j = j.join(spike_mut_df, how = 'left')
    j = j.reset_index()
    
    # add fasta header and seq run
    fasta_header = j.apply(lambda x:create_fasta_header(x.accession_id), axis=1)
    j.insert(value = fasta_header, column = 'fasta_header', loc = 0)
    j['seq_run'] = seq_run
    
    col_order = [ 'accession_id',
                 'fasta_header',
                 'spike_mutations',
                 'nextclade',
                 'pangolin_lineage',
                 
                 'percent_non_ambigous_bases',
                 'mean_depth',

                 'number_aligned_bases',
                 'number_non_ambigous_bases', 
                 
                 'number_seqs_in_fasta',
                 
                 'num_reads',
                 'mean_base_quality', 
                 'mean_map_quality',
                 'number_N_bases', 
                    
                 'nextclade_version',
                 'pangolin_probability', 
                 'pangolin_version',
                 'pangoLEARN_version', 
                 'pangolin_status',
                 'pangolin_notes', 
                 'seq_run']
    
    j = j[col_order]
                   
    # add in 'failed assembly" in missing columns
    j.spike_mutations = j.spike_mutations.fillna(value = '')
    j.nextclade = j.nextclade.fillna(value = '')
    j.nextclade_version = j.nextclade_version.fillna(value = next_version)
    j.pangolin_lineage = j.pangolin_lineage.fillna(value = 'None')
    j.percent_non_ambigous_bases = j.percent_non_ambigous_bases.fillna(value = 0)
    j.mean_depth = j.mean_depth.fillna(value = 0)
    j.number_aligned_bases = j.number_aligned_bases.fillna(value = 0)
    j.number_seqs_in_fasta = j.number_seqs_in_fasta.fillna(value = 0)
    j.num_reads = j.num_reads.fillna(value = 0)
    j.mean_base_quality = j.mean_base_quality.fillna(value = 0)
    j.mean_map_quality = j.mean_map_quality.fillna(value = 0)
    j.number_N_bases = j.number_N_bases.fillna(value = 29903)
    j.pangolin_version = j.pangolin_version.fillna(value = pangolin_version )
    j.pangoLEARN_version = j.pangoLEARN_version.fillna(value = pango_learn_version )
    
    outfile = '%s_sequencing_results.csv' % seq_run
    j.to_csv(outfile, index = False) 

    return j
                      

def make_assembly_metrics_csv(results_df, seq_run):
    
    drop_col = [ 'fasta_header', 'spike_mutations',
        'nextclade', 'pangolin_lineage',
        'nextclade_version', 'pangolin_probability',
       'pangoLEARN_version', 'pangolin_version', 'pangolin_status',
       'pangolin_notes' ]
    results_df = results_df.drop(columns = drop_col)
    
    outfile = '%s_sequence_assembly_metrics.csv' % seq_run
    results_df.to_csv(outfile, index = False)
    
    
def make_wgs_horizon_output (results_df, seq_run):
        
    col_drop = [ 'fasta_header',
                 'nextclade',
                 'mean_depth',
                 'number_aligned_bases',
                 'number_non_ambigous_bases',               
                 'num_reads',
                 'mean_base_quality', 
                 'mean_map_quality',
                 'number_N_bases', 
                 'pangolin_probability', 
                 'pangolin_status',
                 'pangolin_notes', 
                 'seq_run']

    d = results_df.drop(columns = col_drop)

    d = d.rename(columns = {'percent_non_ambigous_bases' : 'percent_coverage'})

    def report_to_epi(pangolin_lin, percent_coverage):
        VOCs = [re.compile('B.1.1.7'), re.compile('B.1.351'), re.compile('P.1'), 
                re.compile('B.1.429'), re.compile('B.1.427')]

        VUIs = [re.compile('B.1.525'), re.compile('B.1.526'), re.compile('P.2'), 
                re.compile('B.1.617'), re.compile('P.2')]
        
 
        if any(VOC.match(pangolin_lin) for VOC in VOCs) and float(percent_coverage) >= 90:
            return 'VOC assigned- %s lineage' % pangolin_lin

        elif any(VUI.match(pangolin_lin) for VUI in VUIs) and float(percent_coverage) >= 90:
            return 'VUI assigned- %s lineage' % pangolin_lin

        elif ((not any(VOC.match(pangolin_lin) for VOC in VOCs) or not any(VUI.match(pangolin_lin) for VUI in VUIs)) and 
              float(percent_coverage) >= 90):
            return 'VOC/VUI not assigned'

        elif float(percent_coverage) < 90 and any(VOC.match(pangolin_lin) for VOC in VOCs):
            return 'sequence coverage not met but assigned to VOC- %s lineage' % pangolin_lin

        elif float(percent_coverage) < 90 and any(VUI.match(pangolin_lin) for VUI in VUIs):
            return 'sequence coverage not met but assigned to VUI- %s lineage' % pangolin_lin

        elif  float(percent_coverage) < 90:
            return 'sequence coverage not met'

    VOC = d.apply(lambda x:report_to_epi(x.pangolin_lineage, x.percent_coverage), axis = 1)
    d.insert(loc =1, column = 'report_to_epi', value = VOC)
    d['Run_Date'] = str(date.today())

    col_order = ['accession_id', 'percent_coverage', 'pangolin_lineage', 'pangolin_version', 
                 'report_to_epi', 'Run_Date', 'pangoLEARN_version']
    d = d[col_order]
#     d.pangolin_lineage = d.pangolin_lineage.fillna(value = 'sample failed assembly')
#     d.pangolin_version = d.pangolin_version.fillna(value = 'sample failed assembly')
#     d.pangoLEARN_version = d.pangoLEARN_version.fillna(value = 'sample failed assembly')
    
    outfile = "%s_wgs_horizon_report.csv" % seq_run
    d.to_csv(outfile, index = False)
    
    
if __name__ == '__main__':
    
    options = getOptions()
    
    sam_df = concat_samtools(bam_file_list = options.bam_file_list)
    
    percent_cvg_df = concat_percent_cvg(percent_cvg_file_list = options.percent_cvg_file_list)
    
    spike_mut_df = get_df_spike_mutations(variants_csv = options.nextclade_variants_csv)
    
    results_df = concat_results(sample_list = options.sample_list,
                                samtools_df = sam_df, 
                               percent_cvg_df = percent_cvg_df, 
                               spike_mut_df = spike_mut_df, 
                               nextclade_clades_csv = options.nextclade_clades_csv, 
                               pangolin_lineage_csv = options.pangolin_lineage_csv, 
                               next_version = options.nextclade_version, 
                               pangolin_version = options.pangolin_version, 
                               seq_run = options.seq_run)
    
    make_assembly_metrics_csv(results_df = results_df, seq_run = options.seq_run)
    
    make_wgs_horizon_output(results_df = results_df, seq_run = options.seq_run)
    
    
    

    
    
    