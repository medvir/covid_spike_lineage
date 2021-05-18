#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Detect interesting lineages from a part of Spike Sanger sequences:
    The mutation in the range of (aa 403 to aa 681) are considered.
    It is possible to have overlapping sequences for one patient. The files 
    beloning to the same patinet should start with a patient number in the following 
    format:
        [PATINET_NUMBER]_???.ab1 or [PATINET_NUMBER]_???.fasta
    If file extension is ab1:
    For each ab1 file in current (or given directory):
        - Create fastq file and check sequence quality and trim low quality bases
    - If needed reverse complement the sequences
    - check if the sequences are coming from the same patient
        - add these sequences in one fasta file
        - create a consensus sequence allowing gaps
    - Run Nextclade docker to get the alignment, clade and other information. 
    - For nextclade quality criteria adjust the thresholds
    - generate the report in pdf format with the file name result_[DIRNAME].pdf
    in the report you get the information about:
        lineage, quality, alignment covering all the important mutations, 
        aa substitutions and nc substitutions 
    
    All output files including the pdf report will be generated in the given or 
    current folder. 
"""

import os
import logging
import logging.handlers
import subprocess
import glob
import time
import pandas as pd
import re
import shutil
import sys

CRITICAL_POS = [417,439,452,453,477, 478,484,501,518,570,614,655,677,681,701,716]

def run_cmd(cmd, exe='/bin/bash'):
    '''check_output if needed gives the error code'''
    try:
        output = subprocess.check_output(cmd, universal_newlines=True,
        shell=True, executable=exe,
        stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as ee:
        print("Execution of %s failed with returncode %d: %s" % (cmd, ee.returncode, ee.output))
        output = None
    return output

def is_aligned_fun(row):
    
    # check if alignment covers mutation 417 to 690
    # start of S-gene is at 21563. 
     
    start_aa = int((int(row['alignmentStart']) - 21563)/3)+1
    end_aa = int((int(row['alignmentEnd']) - 21563)/3)+1
    return '%s-%s' %(start_aa , end_aa)

def calc_quality(row):
    if (row['qc.mixedSites.totalMixedSites'] < 20 and
        row['qc.privateMutations.total'] < 10 and
        row['qc.snpClusters.totalSNPs'] < 5 and
        row['qc.missingData.totalMissing'] < 10):
        return 'good'
    else:
        return 'bad'

def create_md(nextclade_page1, nextclade_page2, nextclade_page3, nextclade_page4):
    logging.info('Writing report in markdown')
    md_rh = open('result.md', 'w')
    full_dir_name = os.getcwd()
    dir_name = os.path.basename(full_dir_name)    
    print('SARS-CoV-2 S-gene Sanger sequencing %s report \n' %dir_name , file=md_rh)
    date_title = time.strftime("%d-%m-%Y")
    print('date: %s' %date_title, file=md_rh)
    print('==================\n', file=md_rh)
    print(nextclade_page1.to_markdown() , file=md_rh)
    print('\n', file=md_rh)
    
    #print('\pagebreak', file=md_rh)
    print('\n', file=md_rh)
    print(nextclade_page2.to_markdown() , file=md_rh)
    print('\n', file=md_rh)
    
    print('\pagebreak', file=md_rh)
    print('\n', file=md_rh)
    print(nextclade_page3.to_markdown() , file=md_rh)
    print('\n', file=md_rh)
    
    #print('\pagebreak', file=md_rh)
    print('\n', file=md_rh)
    print(nextclade_page4.to_markdown() , file=md_rh)
    print('\n', file=md_rh)
    
    md_rh.close()
        
def convert_2_pdf(sample_id='', version='unknown'):
    """Convert markdown file to pdf with pandoc, filling sample and version info.
    :param sample_id: sample name
    :param version: MinVar version
    """

    logging.info('Converting markdown to pdf with pandoc')
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    template_full_path = "%s/aux/template.tex" % (base_dir)
    full_dir_name = os.getcwd()
    dir_name = os.path.basename(full_dir_name)
    file_name = 'result_%s.pdf' %dir_name
    pand_cml = 'pandoc --template=%s result.md -o %s' %(template_full_path, file_name)
    logging.debug(pand_cml)
    subprocess.call(pand_cml, shell=True)
    logging.debug(pand_cml)
    subprocess.call(pand_cml, shell=True)
    
def range_nc_2_aa(pos_ls):
    amb_flag = 'n'
    aa_pos_ls = []
    for nc_pos in pos_ls:
        if '-' in nc_pos:
            range_ls = nc_pos.split('-')
            aa_low = int((int(range_ls[0]) - 21563)/3)+1
            aa_high = int((int(range_ls[1]) - 21563)/3)+1
            if any(aa_pos in CRITICAL_POS for aa_pos in range(aa_low,aa_high+1)):
                amb_flag = 'y'
            aa_str = '%s-%s' %(aa_low,aa_high)
            aa_pos_ls.append(aa_str)
        else:
            aa_pos = int((int(nc_pos) - 21563)/3)+1
            aa_pos_ls.append(aa_pos)
            if aa_pos in CRITICAL_POS:
                amb_flag = 'y'
            
    return aa_pos_ls, amb_flag
    
def nc_2_aa(row):
    """" convert nucleotide positions to aa positions
    if ambiguous characters are in important positions, set the amb_flag"""
    aa_amb_pos_ls = []
    amb_flag = 'n'
    amb_flag1 = 'n'
    amb_flag2 = 'n'
    if not pd.isnull(row['nonACGTNs']):
        amb_string = row['nonACGTNs']
        amb_pos = re.sub('[A-Z]:', '', amb_string)
        nc_pos_ls = [nc_pos for nc_pos in amb_pos.split(',')]
        aa_amb_pos_ls, amb_flag1 = range_nc_2_aa(nc_pos_ls)
    
    aa_missing_pos_ls = []
    if not pd.isnull(row['missing']):
        missing_pos = row['missing']
        nc_pos_ls = [nc_pos for nc_pos in missing_pos.split(',')]
        aa_missing_pos_ls, amb_flag2 = range_nc_2_aa(nc_pos_ls)
    
    if amb_flag1 == 'y' or amb_flag2 == 'y':
        amb_flag = 'y'
    aa_pos_ls = []
    aa_pos_ls = aa_amb_pos_ls + aa_missing_pos_ls
    return aa_pos_ls, amb_flag

def detect_lineage(row):
    
    lineage=''
    
    if not pd.isnull(row['aaSubstitutions']):
        if re.search('S:[A-Z]614G', row['aaSubstitutions']):
            if re.search('S:[A-Z]484K', row['aaSubstitutions']):
                lineage = 'Brazil:P.2'
                if re.search('S:[A-Z]677H', row['aaSubstitutions']):
                    lineage = 'B.1.525'
                
                elif re.search('S:[A-Z]501Y', row['aaSubstitutions']):
                    if (re.search('S:[A-Z]681H', row['aaSubstitutions'])):
                        lineage = 'Brazil:P.3'
                    
                    elif (re.search('S:[A-Z]417N', row['aaSubstitutions']) and 
                        re.search('S:[A-Z]701V', row['aaSubstitutions']) ): 
                        lineage = 'SA:B.1.351'
                        
                    elif (re.search('S:[A-Z]417T', row['aaSubstitutions']) and 
                          re.search('S:[A-Z]655Y', row['aaSubstitutions'])):
                        lineage = 'Brazil:P.1'
                        
                elif (re.search('S:[A-Z]477N', row['aaSubstitutions']) and 
                      re.search('S:[A-Z]701V', row['aaSubstitutions'])):
                    lineage = 'B.1.526'
                    
            elif (re.search('S:[A-Z]501Y', row['aaSubstitutions']) and 
                  re.search('S:[A-Z]570D', row['aaSubstitutions']) and
                  re.search('S:[A-Z]716I', row['aaSubstitutions'])):
                lineage = 'B.1.1.7'
                
            elif re.search('S:[A-Z]452R' , row['aaSubstitutions']):
                lineage = 'Cal B.1.427/B.1.429'
                
                if re.search('S:[A-Z]681R', row['aaSubstitutions']):
                    lineage = 'B.1.617 Indian'
                    if re.search('S:[A-Z]484Q' , row['aaSubstitutions']):
                        lineage = 'B.1.617.1/B.1.617.3 Indian'
                    elif(re.search('S:[A-Z]478K', row['aaSubstitutions'])):
                        lineage = 'B.1.617.2 Indian'

        else:
            lineage = 'not concerned lineage'
        
        if not lineage:
            lineage = 'important mutations, no lineage'
        return lineage
    else:
        lineage = 'n'
        return lineage

def create_consensus(patient_multiple_sequences_dict):
    
    from io import StringIO
    from Bio import AlignIO
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import AlignInfo
    
    # add all sequences that comes from one Patient in one file
    for key in patient_multiple_sequences_dict:
        in_file_pattern = '%s*_final.fasta' %key
        for file in glob.glob(in_file_pattern):
            cmd = 'cat %s >> %s_patient.fasta' %(file, key)
            run_cmd(cmd)
            concated_files = '%s.concated' %file
            os.rename(file, concated_files )
        #perform multiple sequence alignemnet
        from Bio.Align.Applications import MafftCommandline
        in_file = '%s_patient.fasta' %key
        mafft_cline = MafftCommandline(input=in_file)
        stdout, stderr = mafft_cline()
        aligned_file = '%s_aligned.fasta' %key
        with open(aligned_file, "w") as handle:
            handle.write(stdout)
        
        alignment = AlignIO.read(StringIO(stdout), "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus(threshold=0.2, ambiguous='N')
        seq_id = '%s_cons' %key
        my_cons = SeqRecord(consensus, id=seq_id, description='')

        out_file = '%s_final.fasta' %key
        
        with open(out_file, "w") as handle:
            SeqIO.write(my_cons, out_file,'fasta')
    
    return

def main(seq_dir_path, ext):
    
    try:
        os.path.isdir(seq_dir_path)
        os.chdir(seq_dir_path)
    except:
        logging.error("The sanger directory does not exists.")
    
    logging.info('convert ab1 files to fastq, quality trim, if needed reverse complement')
    patient_dict = {}
    
    ext_string = '*.%s' %ext
    if not glob.glob(ext_string):
        logging.error('file with the extension %s does not exists.' %ext)
        sys.exit()
        
    for seq_file in glob.glob(ext_string):
        
        molis_number = seq_file.split('_')[0]
        if molis_number in patient_dict:
            patient_dict[molis_number] += 1
        else:
            patient_dict[molis_number] = 1
        
        file_name = os.path.splitext(seq_file)[0]
        
        if ext == 'ab1':
            fastq_file = '%s.fastq' %file_name
            cmd = 'seqret -sformat abi -osformat fastq  -auto -stdout -sequence %s > %s' % (seq_file, fastq_file)
            run_cmd(cmd)
            fastq_trim_file = '%s_trim.fastq' %file_name
            cmd = 'seqtk trimfq %s > %s' % (fastq_file ,fastq_trim_file)
            run_cmd(cmd)
            
            cmd = 'seqret -sequence %s_trim.fastq -outseq %s_trim.fasta' %(file_name, file_name)
            run_cmd(cmd)
        
            if '_R' in file_name:
                cmd = 'revseq -sequence %s_trim.fasta -outseq %s_final.fasta -tag N' %(file_name,file_name)
                run_cmd(cmd)
            elif '_F' in file_name:
                source_file = '%s_trim.fasta' %(file_name)
                dest_file = '%s_final.fasta' %(file_name)
                shutil.copy(source_file, dest_file)
        else:
            if '_R' in file_name:
                cmd = 'revseq -sequence %s.fasta -outseq %s_final.fasta -tag N' %(file_name,file_name)
                run_cmd(cmd)
            elif '_F' in file_name:
                source_file = '%s.fasta' %(file_name)
                dest_file = '%s_final.fasta' %(file_name)
                shutil.copy(source_file, dest_file)
    
    patient_multiple_sequences_dict = dict((k, v) for k, v in patient_dict.items() if v > 1)
    if bool(patient_multiple_sequences_dict):
        logging.info('Several sequences belong to one patient, create consensus sequence')
        create_consensus(patient_multiple_sequences_dict)
    
    cmd = 'cat *_final.fasta > all_final.fasta'
    run_cmd(cmd)
    
    timestr = time.strftime("%Y%m%d_%H%M%S")
    result_file = 'res_%s.tsv' %timestr
    
    # run nextclade
    logging.info('Run nextclade')
    cmd = "docker run -t --rm -u 1000 --volume='%s:/seq' neherlab/nextclade \
        nextclade --input-fasta '/seq/all_final.fasta' --output-tsv '/seq/%s' \
            --user '$(id -u):$(id -g)'" %(seq_dir_path, result_file)
    run_cmd(cmd)
    
    nextclade_df = pd.read_csv(result_file, sep='\t')
    columns = ['seqName', 'clade', 'aaSubstitutions', 'alignmentEnd', 
               'alignmentStart', 'qc.mixedSites.totalMixedSites', 
               'qc.privateMutations.total', 'qc.missingData.totalMissing', 
               'qc.snpClusters.totalSNPs','substitutions', 'nonACGTNs', 'missing']
    
    nextclade_res_df = nextclade_df[columns]
    
    #Check if the important positions are sequenced
    nextclade_res_df['Alignment'] = nextclade_res_df.apply(is_aligned_fun, axis=1)
    #Calculate the quality
    nextclade_res_df['quality'] = nextclade_res_df.apply(calc_quality, axis=1)
    # Identify lineages using mutations in S-gene:403 to 760
    nextclade_res_df['lineage'] = nextclade_res_df.apply(detect_lineage, axis=1)
    nextclade_res_df[['ambiguous_pos','amb_flag']] = nextclade_res_df.apply(nc_2_aa, axis=1, result_type="expand")
    
    #Modifications for the report
    nextclade_res_df['ambiguous_pos'] = nextclade_res_df['ambiguous_pos'].astype(str).str.wrap(30)
    nextclade_res_df['aaSubstitutions'] = nextclade_res_df['aaSubstitutions'].astype(str).str.wrap(30)
    nextclade_res_df['substitutions'] = nextclade_res_df['substitutions'].astype(str).str.wrap(30)
    
    col_drops = ['alignmentEnd','alignmentStart', 'qc.mixedSites.totalMixedSites', 
                 'qc.privateMutations.total', 'qc.missingData.totalMissing', 
                 'qc.snpClusters.totalSNPs', 'clade', 'nonACGTNs', 'missing']
    nextclade_res_df.drop(col_drops, axis=1, inplace=True)
    
    nextclade_res_page1_df = nextclade_res_df.drop([ 'aaSubstitutions','substitutions', 'ambiguous_pos'], axis=1)    
    col_drop_page2 = [ 'lineage','Alignment', 'quality', 'substitutions', 'ambiguous_pos','amb_flag']
    nextclade_res_page2_df = nextclade_res_df.drop(col_drop_page2, axis=1)
    col_drop_page4 = [ 'lineage','Alignment', 'quality', 'aaSubstitutions', 'substitutions','amb_flag']
    nextclade_res_page3_df = nextclade_res_df.drop(col_drop_page4, axis=1)
    col_drop_page3 = [ 'lineage','Alignment', 'quality', 'aaSubstitutions', 'ambiguous_pos','amb_flag']
    nextclade_res_page4_df = nextclade_res_df.drop(col_drop_page3, axis=1)    
    
    create_md(nextclade_res_page1_df,nextclade_res_page2_df, nextclade_res_page3_df, nextclade_res_page4_df)
    convert_2_pdf()
    
if __name__ == "__main__":
    seq_path = os.path.abspath(os.getcwd())
    main(seq_dir_path=seq_path, ext = 'ab1')