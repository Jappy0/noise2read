# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-01-16 15:52:44
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-02-16 11:10:55

import editdistance
import networkx as nx
import os
import csv
from tqdm import tqdm
from noise2read.utils import *
from mpire import WorkerPool
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
from collections import Counter
import sys
import pandas as pd

class DataGneration():
    """
    A class to generate genuine and ambiguous errors from 1nt/2nt-edit-distance-based graphs construted from short reads dataset
    """
    def __init__(self, logger, config):
        """
        initialize the DataGneration class

        Args:
            logger (class): customized logging
            config (class): parameters below setting using configparser
                parameters:
                config.num_workers (int): the number of subprocesses to use.
                config.input_file (str): The filename including path to be corrected.
                config.result_dir (str): Full pathname where the results will be saved.
                config.high_freq_thre (int, optional): The threshold of the high frequency. Defaults to 5.
                config.max_error_freq (int, optional): The highest frequency of error sequences. Defaults to 3.
                config.save_graph (bool, optional): If true, noise2read will save the construted graph to file. Defaults to False. 
                config.graph_visualization (bool, optional): If true, noise2read will draw the connected subgraphs of the construted  graph. Defaults to False.
                config.drawing_graph_num (int, optional): The number of subgraphs will be drawing. Defaults to 50.
                config.high_ambiguous (bool, optional): If true, noise2read will predict high ambiguous errors. Defaults to False.
                config.min_iters (int, optional): Minimum progress display update interval in iterations for the module tqdm. Defaults to 100.
                config.verbose (bool, optional): If true, noise2read will save the genuine, ambiguous errors and negative reads to csv.
        """
        self.logger = logger
        self.config = config
        if os.path.exists(self.config.result_dir):
            self.logger.warning("Directory '% s' already exists, noise2read will use it." % self.config.result_dir)
        else:
            os.makedirs(self.config.result_dir)
            self.logger.info("Directory '% s' created" % self.config.result_dir)
        
        self.file_type = parse_file_type(self.config.input_file)
        if ".gz" in self.file_type:
            self.out_file_tye = self.file_type.split(".gz")[0] 
        else:
            self.out_file_tye = self.file_type

    def graph_summary(self, graph):
        """
        print the graph info
        """
        edge_nums = graph.number_of_edges()
        node_nums = graph.number_of_nodes()
        self.logger.info("Graph Summary: Edge Number: {}, Nodes Number: {}".format(edge_nums, node_nums))
        return
        
    def err_type_classification(self, line):
        """
        calculate the error type, position, and 2- or 3-mer consisting of two bases before and after the error base.

        Args:
            line (list): [read1, read1_frequency, read1_degree, read2, read2_frequency, read2_degree]

        Raises:
            ValueError: The function raise an error when the edit distance of two reads in the input list not equal to one.

        Returns:
            list: [read1, read1_frequency, read1_degree, errorType, errorPosition, f_kmer, s_kmer, read2, read2_frequency, read2_degree]
        """
        read1 = line[0]
        read1_fre = line[1]
        read2 = line[3]   
        read2_fre = line[4] 

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
            new_line_list = [read1, read1_fre, line[2], errorType, position, f_kmer, s_kmer, read2, read2_fre, line[5]]
        else:
            raise ValueError("The edit distance of two reads in the input list must equal to one!")
        return new_line_list

    def real_ed1_seqs(self, total_seqs, read):
        """
        given a read, generate all its 1nt-edit-distance read counterparts existing in the dataset to form the edges  

        Args:
            total_seqs (list): The list consisting of all reads in the sequencing dataset.
            read (str): A DNA/RNA sequence.

        Returns:
            list: list of tuples of read pairs with only one base different
        """
        possible_ed1 = enumerate_ed1_seqs(read)
        real_seqs =  total_seqs.intersection(possible_ed1)
        read_lst = [read] * len(real_seqs)
        edges = list(zip(read_lst, real_seqs))
        return edges

    def extract_err_samples(self, graph, edit_dis):
        """
        extract genuine errors, error-free records, ambiguous errors or high ambiguous errors from read graph

        Args:
            graph (object): A graph constructed using NetworkX.
            edit_dis (int): set edit distance 1 or 2 to search edges for constructing graph

        Returns:
            DataFrame: three or four pandas dataframes saving genuine, error-free records, ambiguous errors or high ambiguous errors
        """
        subgraphs = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
        # genuine errors and ambiguous errors
        genuine_df, ambiguous_df = self.extract_genuine_ambi_errs(subgraphs, edit_dis)
        # isolated negative samples
        negative_df = self.extract_isolated_negatives(graph, edit_dis)
        # high ambiguous errors
        if edit_dis == 1 and self.config.high_ambiguous: #
            high_ambiguous_df = self.extract_high_ambiguous_errs(subgraphs)
            return genuine_df, negative_df, ambiguous_df, high_ambiguous_df
        else:
            return genuine_df, negative_df, ambiguous_df

    def extract_umi_genuine_errs(self):  
        """
        extract genuine errors from umi graph

        Returns:
            DataFrame: one pandas dataframe saving genuine errors
        """
        graph, seqs_lens_lst, seqs2id_dict, unique_seqs = self.generate_graph(self.config.input_file, edit_dis=1)
        subgraphs = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
        genuine_df = pd.DataFrame(columns= ["StartRead","StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
        for sub_graph in tqdm(subgraphs, desc=self.logger.info("Extract samples with genuine errors"), miniters=int(len(subgraphs)/self.config.min_iters)):
            edges_lst = [e for e in sub_graph.edges()]
            if len(edges_lst) > 0:
                nodes_lst = list(sub_graph.nodes)
                for node in nodes_lst:
                    node_count = sub_graph.nodes[node]['count']
                    node_degree = sub_graph.degree[node]
                    line = []
                    if node_degree >= 1 and node_degree <= 4 and node_count <= self.config.max_error_freq and not sub_graph.nodes[node]['flag']:
                        node_neis = [n for n in sub_graph.neighbors(node)]
                        # nei2count = []
                        nei_degree_count = []
                        for nei in node_neis:
                            nei_count = sub_graph.nodes[nei]['count']
                            nei_degree = sub_graph.degree[nei]
                            # nei2count.append((nei, nei_count))
                            tt = nei_degree * 0.5 + nei_count * 0.5
                            nei_degree_count.append((nei, tt, nei_count))
                        # nei2count.sort(key=lambda x:x[1], reverse=True)
                        nei_degree_count.sort(key=lambda x:x[1], reverse=True)
                        # first_nei, first_nei_count = nei2count[0]
                        first_nei, tt, first_nei_count = nei_degree_count[0]
                        first_nei_degree = sub_graph.degree[first_nei]
                        if first_nei_count >= self.config.high_freq_thre:
                            line = [first_nei, first_nei_count, first_nei_degree, node, node_count, node_degree]
                            newline = self.err_type_classification(line)
                            genuine_df.loc[len(genuine_df)] = newline
                            sub_graph.nodes[node]['flag'] = True
                    else:
                        continue   
        if self.config.verbose:
            genuine_csv = os.path.join(self.config.result_dir, "umi_gnuine.csv")
            genuine_df.to_csv(genuine_csv, index=False) 
        return genuine_df

    def extract_genuine_ambi_errs(self, subgraphs, edit_dis):
        """
        extract genuine and ambiguous errors from read graph

        Args:
            subgraphs (class): Subgraphs of graph constructed using NetworkX.
            edit_dis (int): set edit distance 1 or 2 to search edges for constructing graph

        Returns:
            DataFrame: two pandas dataframes saving genuine and ambiguous errors
        """
        if edit_dis == 1:
            genu_columns = ["StartRead","StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"]
            ambiguous_df= pd.DataFrame(columns=["idx", "StartRead","StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
            genuine_csv = self.config.result_dir + "genuine1.csv"
            ambiguous_csv = self.config.result_dir + "ambiguous1.csv"
        elif edit_dis == 2 or edit_dis == 3:
            genu_columns = ["StartRead","StartReadCount", "StartDegree", "EndRead", "EndReadCount", "EndDegree"]
            ambiguous_df= pd.DataFrame(columns=["idx", "StartRead", "StartReadCount", "StartDegree", "EndRead", "EndReadCount", "EndDegree"])
            genuine_csv = self.config.result_dir + "genuine2.csv"
            ambiguous_csv = self.config.result_dir + "ambiguous2.csv"  

        genuine_lst = []
        ambiguous_lst = []
        subgraph_num = len(subgraphs)
        shared_obs = subgraphs, edit_dis
        idx = 0
        with WorkerPool(self.config.num_workers, shared_objects=shared_obs, start_method='fork') as pool:
            for genu_ambi_lst in pool.imap(self.extract_genuine_ambi_errs_subgraph, range(subgraph_num), progress_bar=self.config.verbose):
                if genu_ambi_lst[0]:
                    genuine_lst.extend(genu_ambi_lst[0])
                    ambiguous_lst.extend(genu_ambi_lst[1])
                    
        # with WorkerPool(self.config.num_workers, shared_objects=shared_obs, start_method='fork') as pool:
        #     with tqdm(total=subgraph_num, desc=self.logger.info("Extract samples with genuine errors"), miniters=int(subgraph_num/self.config.min_iters)) as pbar:   
        #         for genu_ambi_lst in pool.imap(self.extract_genuine_ambi_errs_subgraph, range(subgraph_num)):
        #             if genu_ambi_lst[0]:
        #                 genuine_lst.extend(genu_ambi_lst[0])
        #                 ambiguous_lst.extend(genu_ambi_lst[1])
        #             pbar.update()

        genuine_df = pd.DataFrame(genuine_lst, columns=genu_columns)

        for a_item in ambiguous_lst:
            for sub_item in a_item:
                sub_item.insert(0, idx)
                ambiguous_df.loc[len(ambiguous_df)] = sub_item
            idx += 1
                
        if self.config.verbose:
            genuine_df.to_csv(genuine_csv, index=False)  
            ambiguous_df.to_csv(ambiguous_csv, index=False)    
        return genuine_df, ambiguous_df
    
    def add_genu_sample(self, shared_obs, i):
        """
        add samples to pd DataFrame

        Args:
            item_lst (list): a list containing samples

        Returns:
            DataFrame: pd dataframe saving one sample one row
        """
        item_lst, edit_dis = shared_obs
        if edit_dis == 1:
            tmp_df = pd.DataFrame(columns=["StartRead","StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
        elif edit_dis == 2:
            tmp_df = pd.DataFrame(columns=["StartRead","StartReadCount", "StartDegree", "EndRead", "EndReadCount", "EndDegree"])
        # for item in item_lst:
        tmp_df.loc[len(tmp_df)] = item_lst[i]
        return tmp_df

    def add_ambi_sample(self, idx, item_lst):
        """
        add ambiguous error samples to pd DataFrame

        Args:
            item_lst (list): a list containing samples

        Returns:
            DataFrame: pd dataframe saving one sample one row
        """
        tmp_df = pd.DataFrame()
        for a_item in item_lst:
            for sub_item in a_item:
                sub_item.insert(0, idx)
                tmp_df.loc[len(tmp_df)] = sub_item
            idx += 1
        return tmp_df

    def extract_genuine_ambi_errs_subgraph(self, shared_obs, ii):
        gen_lst = []
        ambi_lst = []
        sub_graphs, edit_dis = shared_obs 
        sub_graph = sub_graphs[ii]
        edges_lst = [e for e in sub_graph.edges()]
        if len(edges_lst) > 0:
            nodes_lst = list(sub_graph.nodes)
            for node in nodes_lst:
                node_count = sub_graph.nodes[node]['count']
                node_degree = sub_graph.degree[node]
                if node_count <= self.config.max_error_freq and not sub_graph.nodes[node]['flag']:
                    node_neis = [n for n in sub_graph.neighbors(node)]
                    if node_degree == 1:
                        nei = node_neis[0]
                        nei_count = sub_graph.nodes[nei]['count']
                        nei_degree = sub_graph.degree[nei]
                        if nei_count >= self.config.high_freq_thre:
                            line = [nei, nei_count, nei_degree, node, node_count, node_degree]
                            if edit_dis == 1:
                                newline = self.err_type_classification(line)
                                # gen_df.loc[len(gen_df)] = newline
                                gen_lst.append(newline)
                                del line, newline 
                            else:
                                # gen_df.loc[len(gen_df)] = line 
                                gen_lst.append(line)
                                del line                                   
                    elif node_degree <= self.config.ambiguous_error_node_degree:
                        high_num = 0
                        nei2count = []
                        for nei in node_neis:
                            nei_count = sub_graph.nodes[nei]['count']
                            nei_degree = sub_graph.degree[nei]
                            if nei_count >= self.config.high_freq_thre:
                                high_num += 1
                                nei2count.append((nei, nei_count))
                        if high_num == 1:             
                            nei2count.sort(key=lambda x:x[1], reverse=True)
                            first_nei, first_nei_count = nei2count[0]
                            first_nei_degree = sub_graph.degree[first_nei]
                            if first_nei_count >= self.config.high_freq_thre:
                                line = [first_nei, first_nei_count, first_nei_degree, node, node_count, node_degree]
                                if edit_dis == 1:
                                    newline = self.err_type_classification(line)
                                    # gen_df.loc[len(gen_df)] = newline
                                    gen_lst.append(newline)
                                    del line, newline  
                                else:
                                    # gen_df.loc[len(gen_df)] = line 
                                    gen_lst.append(line)
                                    del line
                        else:
                            # ambiguous errors
                            tmp_lst = []
                            for cre_nei, cur_nei_count in nei2count:
                                if cur_nei_count >= self.config.high_freq_thre:
                                    cur_nei_degree = sub_graph.degree[nei]
                                    line = [cre_nei, cur_nei_count, cur_nei_degree, node, node_count, node_degree]
                                    if edit_dis == 1:
                                        newline = self.err_type_classification(line)
                                    else:
                                        newline = line
                                    # newline.insert(0, idx)
                                    # ambi_df.loc[len(ambi_df)] = newline
                                    tmp_lst.append(newline)
                                    del newline 
                            if tmp_lst:
                                ambi_lst.append(tmp_lst)  
                            # idx += 1 
                    sub_graph.nodes[node]['flag'] = True        
                else:
                    continue
        return [gen_lst, ambi_lst]

    def data_files(self, edit_dis):
        """
        function to return the results produced by DataGneration class

        Args:
            edit_dis (int): set edit distance 1 or 2 to search edges for constructing graph

        Returns:
            MultiVariables: MultiVariables for next step error correction
        """
        # 1nt-edit-distance-based graph
        graph, seqs_lens_lst, seqs2id_dict, unique_seqs = self.generate_graph(self.config.input_file, edit_dis)
        seq_max_len = max(seqs_lens_lst)
        seq_min_len = min(seqs_lens_lst)
        self.logger.debug(seqs_lens_lst)
        self.logger.debug("Reads Max length: {}".format(seq_max_len))
        self.logger.debug("Reads Min length: {}".format(seq_min_len))
        if edit_dis == 1 and self.config.high_ambiguous:
            genuine_df, negative_df, ambiguous_df, high_ambiguous_df = self.extract_err_samples(graph, edit_dis)
        elif edit_dis == 2 or edit_dis == 1:
            genuine_df, negative_df, ambiguous_df = self.extract_err_samples(graph, edit_dis)
            
        isolates_file, non_isolates_file = self.extract_isolates(graph, unique_seqs, seqs2id_dict)
        del graph
        if self.config.high_ambiguous:
            return isolates_file, non_isolates_file, unique_seqs, seq_max_len, seq_min_len, genuine_df, negative_df, ambiguous_df, high_ambiguous_df
        else:
            return isolates_file, non_isolates_file, unique_seqs, seq_max_len, seq_min_len, genuine_df, negative_df, ambiguous_df

    def extract_isolated_negatives(self, graph, edit_dis):
        """
        extract isolated high frequency reads from read graph

        Args:
            graph (object): A graph constructed using NetworkX.
            edit_dis (int): set edit distance 1 or 2 to search edges for constructing graph

        Returns:
            DataFrame: one pandas dataframe saving isolated high frequency reads
        """
        if not nx.is_connected(graph):
            self.logger.debug("G is a connected graph: {}".format(nx.is_connected(graph)))
            isolates = set(list(nx.isolates(graph)))        
        else:
            self.logger.debug("G is a connected graph: {}".format(nx.is_connected(graph)) )

        self.logger.debug("Isolated Nodes Number: {}".format(len(isolates)))
        # save isolated nodes as negative samples
        if edit_dis == 1:
            negative_csv = self.config.result_dir + "negative1.csv"
        elif edit_dis == 2:   
            negative_csv = self.config.result_dir + "negative2.csv"  
        negative_df = pd.DataFrame(columns=["StartRead","StartReadCount", "StartDegree"])
        for k in isolates:
            k_count = graph.nodes[k]['count']
            if k_count >= self.config.high_freq_thre:
                k_degree = graph.degree[k]
                line = [k, k_count, k_degree]
                negative_df.loc[len(negative_df)] = line  
        if self.config.verbose:
            negative_df.to_csv(negative_csv, index=False) 
        return negative_df

    def extract_isolates(self, graph, unique_seqs, seqs2id_dict):
        """
        split the source file into two files of isolates and non-isolates based on constructed graph

        Args:
            graph (object): A graph constructed using NetworkX.
            unique_seqs (list): All the unique reads in the dataset
            seqs2id_dict (dict): key: read, id: sequencing id

        Returns:
            files: two files of isolates and non-isolates based on constructed graph
        """
        if not nx.is_connected(graph):
            self.logger.debug("G is a connected graph: {}".format(nx.is_connected(graph)))
            isolates = set(list(nx.isolates(graph)))
        else:
            self.logger.debug("G is a connected graph: {}".format(nx.is_connected(graph)) )

        self.logger.debug("Isolated Nodes Number: {}".format(len(isolates)))

        name_lst = []
        non_name_lst = []
        nonisolated_seqs = unique_seqs - isolates 
        # save isolated nodes as negative samples
        for k in isolates:
            name_lst.extend(seqs2id_dict[k])                  
                
        for s in nonisolated_seqs:
            non_name_lst.extend(seqs2id_dict[s])

        self.logger.debug("isolates name list: {}".format(len(name_lst)))
        self.logger.debug("non-isolates name list: {}".format(len(non_name_lst)))

        bases = self.config.input_file.split('/')[-1]
        base = bases.split('.')

        if self.file_type == 'fastq' or self.file_type == 'fq' or self.file_type == 'fastq.gz' or self.file_type == 'fq.gz':
            isolates_file = self.config.result_dir + base[0] + '_isolates.fastq'
            non_isolates_file = self.config.result_dir + base[0] + '_non_isolates.fastq'
        elif self.file_type == 'fasta' or self.file_type == 'fa' or self.file_type == 'fasta.gz' or self.file_type == 'fas.gz':
            isolates_file = self.config.result_dir + base[0] + '_isolates.fasta'
            non_isolates_file = self.config.result_dir + base[0] + '_non_isolates.fasta' 

        extract_records(self.config.result_dir, name_lst, self.config.input_file, isolates_file)
        extract_records(self.config.result_dir, non_name_lst, self.config.input_file, non_isolates_file)

        self.logger.info("Isolated nodes extraction completed.")
        return isolates_file, non_isolates_file

    def generate_graph(self, data_set, edit_dis):
        """
        construct 1nt- or 2nt-edit-distance-based read graph

        Args:
            data_set (str): The filename including path to be corrected.
            edit_dis (int): set edit distance 1 or 2 to search edges for constructing graph

        Returns:
            MultiVariables: Multi Variables after constructing edit-distance-based read graph
        """
        # if not os.path.exists(data_set):
        #     self.logger.error("No input file!")
        #     # os._exit(0)
        # else:
        self.logger.info(data_set)

        record_iterator, file_type = parse_data(data_set)
        seqs2id_dict = {}
        total_seqs = []
        seq_lens_set = set()
        for item in record_iterator:
            seq = str(item.seq)
            ll = len(seq)
            seq_lens_set.add(ll)
            total_seqs.append(seq)
            seqs2id_dict.setdefault(seq, []).append(str(item.id))
        unique_seqs = set(total_seqs)

        graph = nx.Graph()
        read_count = Counter(total_seqs)
        high_freq = []
        low_freq = []
        for read, frequency in tqdm(read_count.items(), desc=self.logger.info("Read Counts"), miniters=int(len(read_count)/self.config.min_iters)):
            if not graph.has_node(read):
                graph.add_node(read, count = frequency, flag=False)  
            if frequency >= self.config.high_freq_thre:
                high_freq.append(read)
            else:
                low_freq.append(read)
        if len(high_freq) == 0:
            self.logger.error("Error Correction Failed as no high-frequency reads detected.")
            sys.exit(1)
        self.logger.debug(len(read_count))
        ######################################################
        edges_lst = []
        if edit_dis == 1:
            shared_unique_seqs = unique_seqs
        elif edit_dis == 2:
            shared_unique_seqs = low_freq

        with WorkerPool(self.config.num_workers, shared_objects=shared_unique_seqs, start_method='fork') as pool:
            if edit_dis == 1:
                for edge_lst in pool.imap(self.real_ed1_seqs, high_freq, progress_bar=self.config.verbose):
                    edges_lst.extend(edge_lst)
            elif edit_dis == 2:
                for edge_lst in pool.imap(self.real_ed2_seqs, high_freq, progress_bar=self.config.verbose):
                    edges_lst.extend(edge_lst)

        # with WorkerPool(self.config.num_workers, shared_objects=shared_unique_seqs, start_method='fork') as pool:
        #     with tqdm(total=len(high_freq), desc=self.logger.info("Searching Edges"), miniters=int(len(high_freq)/self.config.min_iters)) as pbar:   
        #         if edit_dis == 1:
        #             for edge_lst in pool.imap(self.real_ed1_seqs, high_freq):
        #                 edges_lst.extend(edge_lst)
        #                 pbar.update()
        #         elif edit_dis == 2:
        #             for edge_lst in pool.imap(self.real_ed2_seqs, high_freq):
        #                 edges_lst.extend(edge_lst)
        #                 pbar.update()
        if len(edges_lst) > 0:
            self.logger.debug(len(edges_lst))
            self.logger.debug(edges_lst[0])
            graph.add_edges_from(edges_lst)

        ########################################################
        # save graphs
        if self.config.graph_visualization or self.config.save_graph:
            if edit_dis == 1:
                subdir = self.config.result_dir + "graph1/"                
            elif edit_dis == 2:
                subdir = self.config.result_dir + "graph2/"
            if os.path.exists(subdir):
                self.logger.info("Directory '% s' already exists" % subdir)
            else:
                os.makedirs(subdir)
                self.logger.info("Directory '% s' created" % subdir)
        if self.config.save_graph:
            graph_filename = subdir + 'graph.gexf'
            nx.write_gexf(graph, graph_filename)
            self.logger.info("Graph file *.gexf saved!")
        if self.config.graph_visualization:
            self.draw_graph(graph, subdir, self.config.drawing_graph_num)
            self.logger.info("Graph file *.png saved!")
        if edit_dis == 1:
            return graph, seq_lens_set, seqs2id_dict, unique_seqs
        else:
            return graph, unique_seqs

    def real_ed2_seqs(self, low_seqs, read):
        """
        given a read, generate all its 2nt-edit-distance substitution counterparts existing in the dataset to form the edges

        Args:
            low_seqs (list): The list consisting of all reads with low-frequency in the sequencing dataset.
            read (str): A DNA/RNA sequence.

        Returns:
            list: list of tuples of read pairs with two bases different
        """
        # possible_ed2 = set(self.ed2_seqs(read))
        possible_ed2 = enumerate_ed2_seqs(read)
        real_ed2_seqs = possible_ed2.intersection(low_seqs)
        # real_ed2_seqs = []
        # for seq in possible_ed2:
        #     if seq in low_seqs:
        #         real_ed2_seqs.append(seq)
        read_lst = [read] * len(real_ed2_seqs)
        edges = list(zip(read_lst, real_ed2_seqs))     
        return edges   

    def extract_ed2_errors(self, data_set): 
        """
        extract genuine errors, error-free records, ambiguous errors and unique reads from 2nt-edit-distance-based read graph

        Args:
            data_set (str): The filename including path for 2nt errors correction.

        Returns:
            DataFrame: three pandas dataframes of genuine_df, negative_df and ambiguous_df
            list: unique_seqs save unique reads of input dataset
        """
        self.logger.debug(data_set)
        # if self.seq_min_len > 30:   
        self.logger.info("Constructing 2nt-edit-distance based graph.")
        graph, unique_seqs = self.generate_graph(data_set, edit_dis=2)
        self.graph_summary(graph)
        genuine_df, negative_df, ambiguous_df = self.extract_err_samples(graph, edit_dis=2)
        del graph
        return genuine_df, negative_df, ambiguous_df, unique_seqs

    def extract_amplicon_err_samples(self, data_set, amplicon_low_freq=50, amplicon_high_freq=1500):
        """
        extract amplicon sequencing errors from read graph

        Args:
            data_set (str): The filename including path for further errors correction for amplicon sequencing data.
            amplicon_low_freq (int, optional): The threshold of low frquency determining a read for further rectification for amplicon sequencing ddata. Defaults to 50.
            amplicon_high_freq (int, optional): The threshold of high frquency for amplicon sequencing error correction. Defaults to 1500.

        Returns:
            DataFrame: amplicon_df, a dataframe saving the amplicon sequencing errors
            List: unique_seqs, a list saving unique reads of input dataset
        """
        self.logger.debug(data_set)
        self.logger.info("Constructing the second 1nt-edit-distance based graph for further amplicon sequencing data correction.")
        graph, seq_lens_set, seqs2id_dict, unique_seqs = self.generate_graph(data_set, edit_dis=1)
        self.graph_summary(graph)
        
        subgraphs = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
        
        amplicon_df = pd.DataFrame(columns=["idx", "StartRead","StartReadCount", "StartDegree", "ErrorTye", "ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
        for sub_graph in tqdm(subgraphs, desc=self.logger.info("Extract samples with genuine errors"), miniters=int(len(subgraphs)/self.config.min_iters)):
            edges_lst = [e for e in sub_graph.edges()]
            if len(edges_lst) > 0:
                nodes_lst = list(sub_graph.nodes)
                for node in nodes_lst:
                    # reset node visit status
                    sub_graph.nodes[node]['flag'] = False
        idx = 0            
        for sub_graph in tqdm(subgraphs, desc=self.logger.info("Extract samples with genuine errors"), miniters=int(len(subgraphs)/self.config.min_iters)):
            edges_lst = [e for e in sub_graph.edges()]
            if len(edges_lst) > 0:
                nodes_lst = list(sub_graph.nodes)
                for node in nodes_lst:
                    node_count = sub_graph.nodes[node]['count']
                    node_degree = sub_graph.degree[node]
                    line = []
                    if node_count <= amplicon_low_freq and not sub_graph.nodes[node]['flag'] and node_degree <= self.config.amplicon_error_node_degree:
                        node_neis = [n for n in sub_graph.neighbors(node)]
                        # nei2count = []
                        nei_degree_count = []
                        for nei in node_neis:
                            nei_count = sub_graph.nodes[nei]['count']
                            nei_degree = sub_graph.degree[nei]
                            if nei_count >= amplicon_high_freq:
                                line = [nei, nei_count, nei_degree, node, node_count, node_degree]
                                new_line = self.err_type_classification(line) 
                                new_line.insert(0, idx) 
                                # writer.writerow(new_line)
                                amplicon_df.loc[len(amplicon_df)] = new_line
                        idx += 1
                        sub_graph.nodes[node]['flag'] = True
                    else:
                        continue     
        if self.config.verbose:
            amplicon_csv = self.config.result_dir + "amplicon.csv"
            amplicon_df.to_csv(amplicon_csv, index=False)                                            
        return amplicon_df     

    def draw_graph(self, graph, sub_dir, drawing_graph_num = 50):
        """
        draw subgraphs from edit-distance-based read graph

        Args:
            graph (object): A graph constructed using NetworkX.
            sub_dir (str): Directory for saving subgraphs.
            drawing_graph_num (int, optional): The number of subgraphs to be drawed. Defaults to 50.
        """
        
        # layout graphs with positions using graphviz neato
        # pos = nx.nx_agraph.graphviz_layout(graph, prog="neato")
        # color nodes the same in each connected subgraph
        i = 0
        S = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
        for s in S:
            if nx.number_of_nodes(s) >= 10:
                labels = nx.get_node_attributes(s, 'count')
                # A = nx.nx_agraph.to_agraph(s)
                # print(type(A))
                # A.layout(prog="sfdp", args="-Nshape=point")
                # A.draw(graph_name + str(i) + ".png",  args="")
                nodes_num = s.number_of_nodes()
                fig_size = int(((nodes_num // 100)) * 1.5) + 6
                plt.figure(1, figsize=(fig_size, fig_size))
                nx.draw(s, pos = graphviz_layout(s, prog='sfdp', args="-Nshape=circle"), labels=labels, font_color="black", node_color="#fbb4ae", edge_color="#66c2a5")
                plt.savefig(sub_dir + str(i) + "_" + str(nodes_num) +  ".svg", format="svg")
                i += 1
                del nodes_num
                plt.close()
                if i == drawing_graph_num:
                    break
        return

    def extract_high_ambiguous_errs(self, subgraphs):
        """
        extract high ambiguous errors from read graph

        Args:
            subgraphs (class): Subgraphs of graph constructed using NetworkX.

        Returns:
            DataFrame: One pandas dataframe saving high ambiguous errors
        """
        high_ambiguous_df = pd.DataFrame(columns=["idx", "StartRead", "StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
        idx = 0
        for s in tqdm(subgraphs, desc=self.logger.info("Extract samples of high ambiguous errors from 1nt-edit-distance graph"), miniters=int(len(subgraphs)/self.config.min_iters)):
            edges_lst = [e for e in s.edges()]
            if len(edges_lst) > 0:
                for (a, b) in edges_lst:
                    a_count = s.nodes[a]['count']
                    b_count = s.nodes[b]['count']
                    a_degree = s.degree[a]
                    b_degree = s.degree[b]
                    if a_count > self.config.high_freq_thre and b_count > self.config.high_freq_thre:
                        a2b = [a, a_count, a_degree, b, b_count, b_degree]
                        new_a2b = self.err_type_classification(a2b) 
                        new_a2b.insert(0, idx)     
                        high_ambiguous_df.loc[len(high_ambiguous_df)] = new_a2b 
                        b2a = [b, b_count, b_degree, a, a_count, a_degree]
                        new_b2a = self.err_type_classification(b2a)
                        new_b2a.insert(0, idx)
                        high_ambiguous_df.loc[len(high_ambiguous_df)] = new_b2a
                        idx += 1
            del edges_lst
        if self.config.verbose:
            high_ambiguous_csv = self.config.result_dir + "high_ambiguous_1nt.csv"
            high_ambiguous_df.to_csv(high_ambiguous_csv, index=False)  
        return high_ambiguous_df

'''
    def extract_genuine_ambi_errs(self, subgraphs, edit_dis):
        """
        extract genuine and ambiguous errors from read graph

        Args:
            subgraphs (class): Subgraphs of graph constructed using NetworkX.
            edit_dis (int): set edit distance 1 or 2 to search edges for constructing graph

        Returns:
            DataFrame: two pandas dataframes saving genuine and ambiguous errors
        """
        if edit_dis == 1:
            genuine_df = pd.DataFrame(columns=["StartRead","StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
            ambiguous_df= pd.DataFrame(columns=["idx", "StartRead","StartReadCount", "StartDegree", "ErrorTye","ErrorPosition", "StartErrKmer", "EndErrKmer", "EndRead", "EndReadCount", "EndDegree"])
            genuine_csv = self.config.result_dir + "genuine1.csv"
            ambiguous_csv = self.config.result_dir + "ambiguous1.csv"
        elif edit_dis == 2 or edit_dis == 3:
            genuine_df = pd.DataFrame(columns=["StartRead","StartReadCount", "StartDegree", "EndRead", "EndReadCount", "EndDegree"])
            ambiguous_df= pd.DataFrame(columns=["idx", "StartRead", "StartReadCount", "StartDegree", "EndRead", "EndReadCount", "EndDegree"])
            genuine_csv = self.config.result_dir + "genuine2.csv"
            ambiguous_csv = self.config.result_dir + "ambiguous2.csv"  
        idx = 0
        for sub_graph in tqdm(subgraphs, desc=self.logger.info("Extract samples with genuine errors"), miniters=int(len(subgraphs)/self.config.min_iters)):
            edges_lst = [e for e in sub_graph.edges()]
            if len(edges_lst) > 0:
                nodes_lst = list(sub_graph.nodes)
                for node in nodes_lst:
                    node_count = sub_graph.nodes[node]['count']
                    node_degree = sub_graph.degree[node]
                    if node_count <= self.config.max_error_freq and not sub_graph.nodes[node]['flag']:
                        node_neis = [n for n in sub_graph.neighbors(node)]
                        if node_degree == 1:
                            nei = node_neis[0]
                            nei_count = sub_graph.nodes[nei]['count']
                            nei_degree = sub_graph.degree[nei]
                            if nei_count >= self.config.high_freq_thre:
                                line = [nei, nei_count, nei_degree, node, node_count, node_degree]
                                if edit_dis == 1:
                                    newline = self.err_type_classification(line)
                                    genuine_df.loc[len(genuine_df)] = newline
                                    del line, newline 
                                else:
                                    genuine_df.loc[len(genuine_df)] = line 
                                    del line                                   
                        else:
                            high_num = 0
                            nei2count = []
                            for nei in node_neis:
                                nei_count = sub_graph.nodes[nei]['count']
                                nei_degree = sub_graph.degree[nei]
                                if nei_count >= self.config.high_freq_thre:
                                    high_num += 1
                                    nei2count.append((nei, nei_count))
                            if high_num == 1:             
                                nei2count.sort(key=lambda x:x[1], reverse=True)
                                first_nei, first_nei_count = nei2count[0]
                                first_nei_degree = sub_graph.degree[first_nei]
                                if first_nei_count >= self.config.high_freq_thre:
                                    line = [first_nei, first_nei_count, first_nei_degree, node, node_count, node_degree]
                                    if edit_dis == 1:
                                        newline = self.err_type_classification(line)
                                        genuine_df.loc[len(genuine_df)] = newline
                                        del line, newline  
                                    else:
                                        genuine_df.loc[len(genuine_df)] = line 
                                        del line
                            else:
                                # ambiguous errors
                                for cre_nei, cur_nei_count in nei2count:
                                    if cur_nei_count >= self.config.high_freq_thre:
                                        cur_nei_degree = sub_graph.degree[nei]
                                        line = [cre_nei, cur_nei_count, cur_nei_degree, node, node_count, node_degree]
                                        if edit_dis == 1:
                                            newline = self.err_type_classification(line)
                                        else:
                                            newline = line
                                        newline.insert(0, idx)
                                        ambiguous_df.loc[len(ambiguous_df)] = newline
                                        del newline   
                                idx += 1 
                        sub_graph.nodes[node]['flag'] = True        
                    else:
                        continue
        if self.config.verbose:
            genuine_df.to_csv(genuine_csv, index=False)  
            ambiguous_df.to_csv(ambiguous_csv, index=False)    
        return genuine_df, ambiguous_df
'''