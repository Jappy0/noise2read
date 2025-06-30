# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-07 17:22:18

from Bio import SeqIO
import gzip
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
import logging
from colorlog import ColoredFormatter
import copy
import random
import editdistance
import os
import datetime
import tracemalloc
# import gc

class MemoryMonitor:
    def __init__(self, logger):
        self.current_memory = 0
        self.peak_memory = 0
        self.traceback = None
        self.logger = logger

    def start(self):
        tracemalloc.start()

    def stop(self):
        tracemalloc.stop()

    def measure(self):
        current, peak = tracemalloc.get_traced_memory()
        self.current_memory = current
        if peak > self.peak_memory:
            self.peak_memory = peak

        self.logger.info(f"Current memory usage: {self.current_memory / (1024 * 1024)} MB")
        self.logger.info(f"Recent peak memory usage: {self.peak_memory / (1024 * 1024)} MB")

    def get_memory_traceback(self, obj):
        self.traceback = tracemalloc.get_object_traceback(obj)
        self.logger.debug("Memory Allocation Traceback:")
        for filename, lineno, funcname, line in self.traceback:
            self.logger.debug(f"  File '{filename}', line {lineno}, in {funcname}")
            if line:
                self.logger.debug(f"    {line.strip()}")


def custom_logger(root_name, debug_mode) -> logging.Logger: 
    logger = logging.getLogger(root_name)    
    
    formatter = ColoredFormatter(
        "%(green)s[%(asctime)s] %(blue)s%(name)s %(log_color)s%(levelname)-8s%(reset)s %(message)s",
        datefmt=None,
        reset=True,
        log_colors={
            'DEBUG':    'cyan',
            'INFO':     'green',
            'WARNING':  'yellow',
            'ERROR':    'red',
            'CRITICAL': 'red, bg_white',
        },
        secondary_log_colors={},
        style='%'
    )
    if debug_mode:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    # Output full log
    file_handler = logging.FileHandler( datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_") + 'noise2read.log')

    file_handler.setLevel(logging.INFO)

    # formatter = logging.Formatter(log_format)
    logger.addHandler(file_handler)

    # # Output warning log
    # file_handler = logging.FileHandler('noise2read.Warning.log')
    # file_handler.setLevel(logging.WARNING)
    # # formatter = logging.Formatter(log_format)
    # file_handler.setFormatter(formatter)
    # logger.addHandler(file_handler)

    # # Output error log
    # file_handler = logging.FileHandler('noise2read.Error.log')
    # file_handler.setLevel(logging.ERROR)
    # # formatter = logging.Formatter(log_format)
    # file_handler.setFormatter(formatter)
    # logger.addHandler(file_handler)

    return logger

def parse_file_type(data_set):
    items = data_set.split(".")
    ext = items[-1]
    if ext == 'fa' or ext == 'fasta':
        f_type = 'fasta'
    elif ext == 'fq' or ext == 'fastq':
        f_type = 'fastq'
    elif ext == 'gz':
        f_type = items[-2] + '.' + ext
    return f_type

def parse_data(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        prefix = file_type.split('.')[0]
        if prefix == 'fa' or 'fasta':
            ff_type = 'fasta'
        elif prefix == 'fq' or 'fastq':
            ff_type = 'fastq'
        handle = gzip.open(data_set, 'rt')
        record_iterator = SeqIO.parse(handle, ff_type)
        return record_iterator, ff_type
    else:
        record_iterator = SeqIO.parse(data_set, file_type) 
        return record_iterator, file_type

def parse_data_index(data_set):
    file_type = parse_file_type(data_set)

    if file_type in {'fastq.gz', 'fq.gz', 'fa.gz', 'fasta.gz'}:
        prefix = file_type.split('.')[0]
        if prefix == 'fa' or 'fasta':
            ff_type = 'fasta'
        elif prefix == 'fq' or 'fastq':
            ff_type = 'fastq'        
        # Convert list to dictionary with record IDs as keys
        with gzip.open(data_set, 'rt') as handle:
            records = {rec.id: rec for rec in SeqIO.parse(handle, ff_type)}

    else:
        records = SeqIO.index(data_set, file_type)  # Dictionary-like object

    return records, file_type

def parse_data_dict(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        prefix = file_type.split('.')[0]
        if prefix == 'fa' or 'fasta':
            ff_type = 'fasta'
        elif prefix == 'fq' or 'fastq':
            ff_type = 'fastq'

        handle = gzip.open(data_set, 'rt')
        records_dict = SeqIO.to_dict(SeqIO.parse(handle, ff_type))
        return records_dict, ff_type
    else:
        records_dict = SeqIO.to_dict(SeqIO.parse(data_set, file_type))   
        return records_dict, file_type
'''
def extract_seq(shared_objects, name):
    records_dict, file_type = shared_objects
    if file_type == "fastq":
        return SeqRecord(records_dict[name].seq, id=records_dict[name].id, description=records_dict[name].description, letter_annotations=records_dict[name].letter_annotations)
    elif file_type == "fasta":
        return SeqRecord(records_dict[name].seq, id=records_dict[name].id, description=records_dict[name].description)

def extract_seqs(num_workers, name_lst, data_set, sub_dataset, logger):
    records_dict, file_type = parse_data_dict(data_set)
    sub_records = []
    shared_objects = records_dict, file_type
    with WorkerPool(num_workers, shared_objects=shared_objects, start_method='fork') as pool:
        with tqdm(total=len(name_lst), desc=logger.info("Extract Seqs")) as pbar:   
            for tmp_rec in pool.imap(extract_seq, name_lst):
                # print(tmp_rec)
                sub_records.append(tmp_rec)
                pbar.update() 
    # print(type(sub_records[0]), sub_records[0])
    with open(sub_dataset, "w") as handle:
        SeqIO.write(sub_records, handle, file_type)
    return 
'''
def extract_records(working_dir, name_lst, data_set, sub_dataset):
    f_name_lst = os.path.join(working_dir, 'tmp_name.lst')
    with open(f_name_lst, "w") as outfile:
        outfile.write("\n".join(name_lst))
    os.system("seqtk subseq %s %s > %s" % (data_set, f_name_lst, sub_dataset))
    os.system("rm %s" % f_name_lst)
    # record_iterator, file_type = parse_data(data_set)
    # records = [rec for rec in record_iterator if str(rec.id) in name_lst]
    # SeqIO.write(records, sub_dataset, file_type)
    #gc.collect()
    return
#####################################################################################
def sub_base(base):
    if base == 'A':
        return 'TCG'
    elif base == 'T':
        return 'ACG'
    elif base == 'C':
        return 'ATG'
    elif base == 'G':
        return 'ACT'
    elif base == 'N':
        return 'N'

def replace_char(seq, char, index):
    seq[index] = char
    return ''.join(seq)

def seq2substitution(read):
    """
    enumerate all the substitutions for one read 
    Args:
        read (str): a sequence
    Returns:
        set: a set contains reads
    """
    editdis1_list = []
    # substitution
    raw_seq = list(read)
    n = len(raw_seq)
    for i in range(n):
        seq = copy.deepcopy(raw_seq)
        temp_base = seq[i]
        #print(temp_base)
        temp_sub_base = list(sub_base(temp_base))
        #print(temp_sub_base)
        for b in temp_sub_base:
            if b != 'N':
                sub_seq = replace_char(seq, b, i)
                editdis1_list.append(sub_seq)
    #gc.collect()
    return set(editdis1_list)

def seq2deletion(read):
    """
    enumerate all the deletions for one read 
    Args:
        read (str): a sequence
    Returns:
        set: a set contains reads
    """
    editdis1_list = []
    seq = list(read)
    n = len(seq)
    # deletion
    for i in range(n):
        del_seq = read[:i] + read[(i+1):]
        editdis1_list.append(del_seq)
    #gc.collect()
    return set(editdis1_list)

def seq2insertion(read):
    """
    enumerate all the insertions for one read 
    Args:
        read (str): a sequence
    Returns:
        set: a set contains reads
    """
    editdis1_list = []
    seq = list(read)
    n = len(seq)
    # insertion
    bases = ['A', 'G', 'C', 'T']
    for i in range(n+1):
        for b in bases:
            raw_seq = copy.deepcopy(seq)
            raw_seq.insert(i, b)
            editdis1_list.append(''.join(raw_seq))
    #gc.collect()
    return set(editdis1_list)

def enumerate_ed1_seqs(read):
    possible_ed1 = []
    possible_ed1.extend(seq2deletion(read))
    possible_ed1.extend(seq2substitution(read))
    possible_ed1.extend(seq2insertion(read))
    #gc.collect()
    return set(possible_ed1)

def enumerate_ed2_seqs(read):
    # possible_ed1 = self.enumerate_ed1_seqs(read)
    possible_ed1 = seq2substitution(read)
    possible_ed2 = []
    for seq in possible_ed1:
        possible_ed2.extend(list(set(seq2substitution(seq)) - possible_ed1))
    #gc.collect()
    return set(possible_ed2)
    
def random_ed2_seq(read, total_reads, num):
    '''
    return num reads that each have two bases difference from given read and these generated reads not exist in total_reads
    '''
    editdis2_list = []
    # substitution
    raw_seq = list(read)
    n = len(raw_seq)
    pos_lst = []
    for i in range(num):
        pos_lst.append(random.sample(range(n), 2))
    for item in pos_lst:
        seq = copy.deepcopy(raw_seq)
        i = item[0]
        j = item[1]
        i_base = seq[i]
        j_base = seq[j]
        i_sub_base = random.sample(list(sub_base(i_base)), 1)[0]    
        j_sub_base = random.sample(list(sub_base(j_base)), 1)[0]  
        sub_seq = replace_char(seq, i_sub_base, i)
        sub_seq = replace_char(seq, j_sub_base, j)
        if sub_seq not in total_reads:
            editdis2_list.append(sub_seq)
    #gc.collect()
    return editdis2_list

def error_type_classification(read1, read2):
    f_len = len(read1)
    s_len = len(read2)  
    # position = -1  
    dis = editdistance.eval(read1, read2)
    # print(dis)
    if dis == 1:              
        if f_len == s_len:
            position = -1
            for index in range(f_len):
                if read1[index] == read2[index]:
                    continue
                else:
                    position = index
                    # first = f_seq[index]
                    # second = s_seq[index]
                    break
            first = read1[position]
            second = read2[position]
            if position == 0:
                f_kmer = read1[0:2]
                s_kmer = read2[0:2]
            elif position == f_len:
                f_kmer = read1[-2:] 
                s_kmer = read2[-2:]
            else:
                f_kmer = read1[position-1 : position+2]
                s_kmer = read2[position-1 : position+2]  

        elif f_len < s_len:
            position = -1
            num = 0
            for index in range(f_len):
                if read1[index] == read2[index]:
                    num = num + 1
                else:
                    position = index
                    break
            if num == f_len:
                position = f_len
            first = 'X'
            second = read2[position]  
            if position == 0:
                f_kmer = 'X' + read1[0]
                s_kmer = read2[0:2]
            elif position == f_len:
                f_kmer = read1[-1] + 'X'
                s_kmer = read2[-2:]
            else:
                f_kmer = read1[position-1] + 'X' + read1[position]
                s_kmer = read2[position-1 : position+2]                    
        elif f_len > s_len:
            position = -1
            num = 0
            for index in range(s_len):
                if read1[index] == read2[index]:
                    num = num + 1
                else:
                    position = index
                    break
            if num == s_len:
                position = s_len
            first = read1[position]
            second = 'X'
            if position == 0:
                f_kmer = read1[0:2]
                s_kmer = 'X' + read2[0]
            elif position == s_len:
                f_kmer = read1[-2:] 
                s_kmer = read2[-1] + 'X'
            else:
                f_kmer = read1[position-1 : position+2]
                s_kmer = read2[position-1] + 'X' + read2[position]    

        errorType = first + '-' + second
        #gc.collect()
        return [errorType, f_kmer, s_kmer]
    else:
        raise ValueError("The editdistance of two reads in the input list must equal to one!")

def usage():
    # print(" ")
    print("noise2read Usage:")
    print("   Mandatory:")
    print("     -m|--module                   module selection")
    print("   Modules: [correction, amplicon_correction, umi_correction, mimic_umi, real_umi, evaluation, simulation]")
    # print(" ")
    print("1. Using config file")
    print("     noise2read -m|--module <module_name> -c <path_configuration_file>")
    print("   Mandatory:")
    print("     -c|--config                   input configuration file")
    # print(" ")
    print("2. Using command line with the default parameters")
    print("     noise2read -m|--module <module_name> -i <path_raw_data.fastq|fasta|fa|fq>")
    print("   Mandatory:")
    print("     -i|--input                    input raw data to be corrected")
    print("   Options:")
    print("     -d|--directory                set output directory")
    print("     -a|--high_ambiguous           predict high ambiguous errors using machine learning when set true, defaut true")
    print("     -t|--true                     input ground truth data if you have")
    print("     -r|--rectification            input corrected data when using module evaluation")
    print("     -p|--parallel                 use multiple cpu cores, default total cpu cores - 2")
    print("     -g|--tree_method              use gpu for training and prediction, default auto, (options gpu_hist, hist, auto)")
    # print("     -o|--over_sampling            use over sampling or downsampling for negative samples, default True")
    print("     -h|--help                     show this help")
    # print("########################################################################################################################")

#####################################################################################
# def split_seqs(num_workers, name_set, data_set, sub_dataset_1, sub_dataset_2):
#     records_dict, file_type = parse_data_dict(data_set)
#     sub_records_1 = []
#     sub_records_2 = []
#     shared_objects = records_dict, file_type

#     with WorkerPool(num_workers, shared_objects=shared_objects, start_method='fork') as pool:
#         with tqdm(total=len(name_lst), desc="Extract Seqs") as pbar:   
#             for tmp_rec in pool.imap(extract_seq, name_lst):
#                 # print(tmp_rec)
#                 sub_records.append(tmp_rec)
#                 pbar.update() 

#     for item in record_iterator:
#         qual = {}
#         name = item.id
#         if file_type == "fastq":
#             qual['phred_quality'] = item.letter_annotations['phred_quality']
#             tmp_rec = SeqRecord(item.seq, id=name, description=item.description, letter_annotations=qual) 
#         elif file_type == "fasta":
#             tmp_rec = SeqRecord(item.seq, id=name, description=item.description)    
#         if name in name_set:             
#             sub_records_1.append(tmp_rec)
#         elif name in non_name_set:                        
#             sub_records_2.append(tmp_rec)
#     with open(sub_dataset_1, "w") as handle:
#         SeqIO.write(sub_records_1, handle, file_type)
#     with open(sub_dataset_2, "w") as handle:
#         SeqIO.write(sub_records_2, handle, file_type)
#     return 

if __name__ == "__main__":
    # input_file = '/home/pping/Data/Repo/data/noise2read_data/group1/raw/new_tcr.seq.real_SRR1543964.fastq'
    
    # records_dict, file_type = parse_data_dict(input_file)
    # sub_records = []
    # for name in records_dict:
    #     if file_type == "fastq":
    #         # qual = {}
    #         # qual['phred_quality'] = records[name].letter_annotations["phred_quality"]
    #         tmp_rec = SeqRecord(records_dict[name].seq, id=records_dict[name].id, description=records_dict[name].description, letter_annotations=records_dict[name].letter_annotations) 
    #         print(type(records_dict[name].seq))
    #         print(type(records_dict[name].id))
    #         print(type(records_dict[name].description))
    #         print(type(records_dict[name].letter_annotations))
    #     elif file_type == "fasta":
    #         tmp_rec = SeqRecord(records_dict[name].seq, id=records_dict[name].id, description=records_dict[name].description)      
    #         print(records_dict[name].seq, records_dict[name].id)
    #     break
    num = 3
    total_reads = ['AACGT', 'GACGT', 'CACGT', 'TACGT', 'AACGT', 'AGCGT','ACCGT', 'ATCGT', 'ACAGT', 'ACGGT', 'ACCGT', 'ACTGT', 'ACGAT']
    read = 'AACGT'
    seqs = random_ed2_seq(read, total_reads, num)
    print(seqs)