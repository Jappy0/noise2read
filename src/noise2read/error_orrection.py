# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-09-05 21:32:02

import os
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from noise2read.isolates_correction import IsolatesErrorCorrection
from noise2read.classifier import MLClassifier 
import copy
from noise2read.utils import *
from noise2read.reads2vectors import Reads2Vectors
import sys
# import shutil
# from noise2read.utils import MemoryMonitor

class ErrorCorrection():
    """
    A class contains error correction related functions
    """
    def __init__(self, logger, config):
        """
        initialize the ErrorCorrection class

        Args:
            logger (class): customized logging
            config (class): parameters below setting using configparser
        """
        self.logger = logger
        self.config = config
        file_type = parse_file_type(config.input_file)
        if ".gz" in file_type:
            self.out_file_tye = file_type.split(".gz")[0]
        else:
            self.out_file_tye = file_type
        bases = config.input_file.split('/')[-1]
        self.base = bases.split('.')
        # Create an instance of the MemoryMonitor
        ##self.MM = MemoryMonitor(self.logger)

    def simplify_ambiguous_err_prediction(self, genuine_df, ambiguous_df):
        """
        ambiguous or high ambiguous errors prediction

        Args:
            total_reads (set): set contains all the reads
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction
            edit_dis (int): edit distance 1 or 2

        Returns:
            MultiVariables: genuine_df, new_negative_df, high_ambiguous_df
        """
        
        # merge genuine and ambiguous errors
        grouped_ambiguous_df = ambiguous_df.groupby("idx")
        # choose the error-free read for the ambiguous errors without using machine learning
        for name, group_df in grouped_ambiguous_df:
            if len(group_df) > 0:
                pred_index = group_df['StartReadCount'].idxmax()
                new_df = group_df.loc[:, ~group_df.columns.isin(['idx', 'predictand'])]
                entry = copy.deepcopy(new_df.loc[[pred_index]])
                genuine_df = pd.concat([genuine_df, entry], ignore_index=True)

        # if edit_dis == 1 and self.config.verbose:
        #     ambiguous_df.to_csv(self.config.result_dir + 'ambiguous_1nt_prediction.csv', index=False)
        # elif edit_dis == 2 and self.config.verbose:
        #     ambiguous_df.to_csv(self.config.result_dir + 'ambiguous_2nt_prediction.csv', index=False)
        del ambiguous_df, grouped_ambiguous_df
        return genuine_df   

    def simplify_2nt_correction(self, data_set, genuine_df, ambiguous_df):       
        """
        correcting ambiguous errors in 2nt-edit-distance-based read graph

        Args:
            data_set (str): raw data filename including path
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction

        Returns:
            str: corrected data filename including path
        """
        if not genuine_df.empty and not ambiguous_df.empty:
            genuine_ambi_errs_df = self.simplify_ambiguous_err_prediction(genuine_df, ambiguous_df)
            correct_file = self.all_in_one_2nt_correct_errors(data_set, genuine_ambi_errs_df)
            os.system("rm %s" % data_set)
            del genuine_ambi_errs_df
            self.logger.info("Error Correction finished.")
            return correct_file
        else:
            return data_set

    def simplify_correction(self, isolates_file, non_isolates_file, genuine_df, ambiguous_df):
        """
        Args:
            isolates_file (str):  The isolates' filename with path.
            non_isolates_file (str):  The non-isolates' filename with path
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction

        Returns:
            corrected_file (str): The corrected data filename including path.
        """
        ##self.MM.start()
        corrected_file = self.config.result_dir + self.base[0] + '_corrected.' + self.out_file_tye  

        if not genuine_df.empty and not ambiguous_df.empty:
            genuine_ambi_errs_df = self.simplify_ambiguous_err_prediction(genuine_df, ambiguous_df)
            ##self.MM.measure()
            # correct errors
            self.logger.info("Correcting 1nt-edit-distance based Errors")
            non_isolates_correct = self.correct_errors(non_isolates_file, genuine_ambi_errs_df)
            self.logger.info('1nt-edit-distance based Errors Correction Finished')
            del genuine_ambi_errs_df, ambiguous_df
        elif not genuine_df.empty:
            non_isolates_correct = self.correct_errors(non_isolates_file, genuine_df)
        else:
            self.logger.error("No genuine and ambiguous errors identified, failed to do error correction!")
            non_isolates_correct = non_isolates_file
            sys.exit(1)

        # # bcool correction
        IEC = IsolatesErrorCorrection(self.logger, self.config.num_workers, isolates_file, non_isolates_correct, self.config.result_dir, self.config.iso_change_detail)
        corrected_isolates = IEC.bcool_correct_isolates() 
        if corrected_isolates and non_isolates_correct:
            os.system("cat %s %s > %s" % (corrected_isolates, non_isolates_correct, corrected_file))
        else:
            self.logger.error("No corrected_isolates and/or non_isolates_correct, failed to do error correction")
        del IEC

        if os.path.exists(isolates_file):
            os.system("rm %s" % isolates_file)
        if os.path.exists(non_isolates_file):
            os.system("rm %s" % non_isolates_file)
        if os.path.exists(corrected_isolates):   
            os.system("rm %s" % corrected_isolates)
        if os.path.exists(non_isolates_correct):     
            os.system("rm %s" % non_isolates_correct)
        ##self.MM.measure()
        ##self.MM.stop()
        return corrected_file

    def all_in_one_correction(self, isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df):
        """
        the main function intergrates all_in_one_ambiguous_err_prediction, all_in_one_high_ambiguous_err_correction and IsolatesErrorCorrection

        Args:
            isolates_file (str):  The isolates' filename with path.
            non_isolates_file (str):  The non-isolates' filename with path
            unique_seqs (set): a set contains unique reads
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            negative_df (DataFrame): pandas dataframe containing negative samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction
            high_ambiguous_df (DataFrame): pandas dataframe containing high ambiguous samples for prediction

        Returns:
            MultiVariables: 
                corrected_file (str): The corrected data filename including path.
                tmp_correct: The filename including path of corrected data without prediction high ambiguous errors
                new_negative_df: pandas dataframe containing negative samples from the prediction result of ambiguous errors
        """
        ##self.MM.start()
        corrected_file = self.config.result_dir + self.base[0] + '_corrected.' + self.out_file_tye  
        if isinstance(high_ambiguous_df, pd.DataFrame) and high_ambiguous_df.empty:
            self.logger.info("No ambiguous between high-frequency reads identified!")
        if not genuine_df.empty and not ambiguous_df.empty and not negative_df.empty:
            if isinstance(high_ambiguous_df, pd.DataFrame) and not high_ambiguous_df.empty:
                genuine_ambi_errs_df, new_negative_df, high_ambi_df = self.all_in_one_ambiguous_err_prediction(unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df, edit_dis=1)
                ##self.MM.measure()
                # correct errors
                genuine_corrected_file = self.correct_errors(non_isolates_file, genuine_ambi_errs_df)
                self.logger.info("Genuine and ambiguous errors corrected.")
                ##self.MM.measure()
                del genuine_ambi_errs_df, genuine_df, negative_df, ambiguous_df, unique_seqs
                ###############################################################################################
                # IEC = IsolatesErrorCorrection(self.logger, self.config.num_workers, isolates_file, genuine_corrected_file, self.config.result_dir, self.config.iso_change_detail, self.config.min_iters)
                # tmp_corrected_isolates = IEC.bcool_correct_isolates() 
                # tmp_correct = self.config.result_dir + self.base[0] + '_correct_no_high.' + self.out_file_tye  
                # if tmp_corrected_isolates and genuine_corrected_file:
                #     os.system("cat %s %s > %s" % (tmp_corrected_isolates, genuine_corrected_file, tmp_correct))
                # # if os.path.exists(tmp_corrected_isolates):     
                # #     os.system("rm %s" % tmp_corrected_isolates)
                # del IEC
                #####################################################################
                non_isolates_correct = self.all_in_one_high_ambiguous_err_correction(genuine_corrected_file, high_ambi_df)
                del high_ambiguous_df
                self.logger.info("High ambiguous errors corrected.")
                self.logger.info('1nt-edit-distance based Errors Correction finished.')
                ##self.MM.measure()
                # self.logger.info("#############################################")
            else:
                genuine_ambi_errs_df, new_negative_df = self.all_in_one_ambiguous_err_prediction(unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df=None, edit_dis=1)
                ##self.MM.measure()
                # correct errors
                self.logger.info("Correcting 1nt-edit-distance based Errors")

                non_isolates_correct = self.correct_errors(non_isolates_file, genuine_ambi_errs_df)
                self.logger.info('1nt-edit-distance based Errors Correction Finished')
                del genuine_ambi_errs_df, negative_df, ambiguous_df, unique_seqs
                ##self.MM.measure()
        elif not genuine_df.empty:
            non_isolates_correct = self.correct_errors(non_isolates_file, genuine_df)
            new_negative_df = pd.DataFrame()
            ##self.MM.measure()
        else:
            self.logger.error("No genuine and ambiguous errors identified, failed to do error correction!")
            non_isolates_correct = non_isolates_file
            sys.exit(1)
        # # bcool correction
        IEC = IsolatesErrorCorrection(self.logger, self.config.num_workers, isolates_file, non_isolates_correct, self.config.result_dir, self.config.iso_change_detail)
        corrected_isolates = IEC.bcool_correct_isolates() 
        if corrected_isolates and non_isolates_correct:
            os.system("cat %s %s > %s" % (corrected_isolates, non_isolates_correct, corrected_file))
        else:
            self.logger.error("No corrected_isolates and/or non_isolates_correct, failed to do error correction")
        del IEC
        if os.path.exists(isolates_file):
            os.system("rm %s" % isolates_file)
        if os.path.exists(non_isolates_file):
            os.system("rm %s" % non_isolates_file)
        if os.path.exists(corrected_isolates):   
            os.system("rm %s" % corrected_isolates)
        if os.path.exists(non_isolates_correct):     
            os.system("rm %s" % non_isolates_correct)
        # if os.path.exists(self.config.result_dir + "bcool/"):
        #     shutil.rmtree(self.config.result_dir + "bcool/") 
        # if isinstance(high_ambiguous_df, pd.DataFrame):
        #     del high_ambiguous_df
        #     return corrected_file, tmp_correct, new_negative_df
        # else:
        ##self.MM.measure()
        ##self.MM.stop()
        return corrected_file, new_negative_df

    def all_in_one_ambiguous_err_prediction(self, total_reads, genuine_df, negative_df, ambiguous_df, high_ambiguous_df, edit_dis):
        """
        ambiguous or high ambiguous errors prediction

        Args:
            total_reads (set): set contains all the reads
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            negative_df (DataFrame): pandas dataframe containing negative samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction
            high_ambiguous_df (DataFrame): pandas dataframe containing high ambiguous samples for prediction
            edit_dis (int): edit distance 1 or 2

        Returns:
            MultiVariables: genuine_df, new_negative_df, high_ambiguous_df
        """
        
        RV = Reads2Vectors(self.logger, self.config, edit_dis)
        # RV = Reads2Vectors(self.logger, self.config.num_workers, self.config.result_dir, self.config.read_max_len, self.config.entropy_kmer, self.config.entropy_q, self.config.kmer_freq, self.config.read_type, edit_dis)
        if edit_dis == 1:
            study_name = 'ambiguous_1nt'
        elif edit_dis == 2:
            study_name = 'ambiguous_2nt'

        train_data, train_labels, ambiguous_data = RV.all_in_one_embedding(total_reads, genuine_df, negative_df, ambiguous_df, high_flag=False)  
        del total_reads
        TM = MLClassifier(self.logger, self.config, study_name, train_data, train_labels, ambiguous_data)
        predictions = TM.tunning(self.config.n_trials)
        ambiguous_df.insert(ambiguous_df.shape[1], 'predictand', predictions)

        # merge genuine and ambiguous errors
        grouped_ambiguous_df = ambiguous_df.groupby("idx")
        new_negative_df = pd.DataFrame()
        for name, group_df in grouped_ambiguous_df:
            self.logger.debug(group_df['predictand'])
            # if "N-X" or "X-N" in group_df['ErrorTye'].tolist():
            if len(group_df) > 0:
                pred_index = group_df['predictand'].idxmax()
                new_df = group_df.loc[:, ~group_df.columns.isin(['idx', 'predictand'])]
                entry = copy.deepcopy(new_df.loc[[pred_index]])
                bad_df = group_df.index.isin([pred_index])
                cur_negative = group_df[~bad_df]
                genuine_df = pd.concat([genuine_df, entry], ignore_index=True)
                new_negative_df = pd.concat([new_negative_df, cur_negative], ignore_index=True)
                
        if edit_dis == 1 and self.config.verbose:
            ambiguous_df.to_csv(self.config.result_dir + 'ambiguous_1nt_prediction.csv', index=False)
        elif edit_dis == 2 and self.config.verbose:
            ambiguous_df.to_csv(self.config.result_dir + 'ambiguous_2nt_prediction.csv', index=False)
        del TM, ambiguous_data, ambiguous_df, grouped_ambiguous_df
        #######################################################################################################
        if isinstance(high_ambiguous_df, pd.DataFrame):
            high_train_data, high_train_labels, high_ambiguous_data = RV.high_all_in_one_embedding(genuine_df, negative_df, new_negative_df, high_ambiguous_df)
            study_name = 'high_ambiguous_1nt'
            TM = MLClassifier(self.logger, self.config, study_name, high_train_data, high_train_labels, high_ambiguous_data)
            high_predictions = TM.tunning(self.config.n_trials)
            high_ambiguous_df.insert(high_ambiguous_df.shape[1], 'predictand', high_predictions)
            del RV, TM, high_train_data, high_train_labels, high_ambiguous_data, negative_df
            return genuine_df, new_negative_df, high_ambiguous_df           
        else:
            del negative_df
            return genuine_df, new_negative_df   

    def correct_errors(self, orginal_file, df_data):
        """
        correct errors in the raw to yield corrected dataset

        Args:
            orginal_file (str): raw data filename including path
            df_data (DataFrame): pandas dataframe save predicted result

        Returns:
            str: corrected data filename including path
        """
        record_dict, ori_file_type = parse_data_dict(orginal_file)
        seq2id = {}
        total_name_lst = []
        for record_id in record_dict:
            seq = str(record_dict[record_id].seq)
            seq2id.setdefault(seq, []).append(record_id)
            total_name_lst.append(record_id)
        err2cor_records = []
        err_name_lst = []
        for row_index, row in df_data.iterrows():
            err_count = row['EndReadCount']
            err_read = row['EndRead']
            cor_read = row['StartRead']
            err_pos = row['ErrorPosition']
            err_tye = row['ErrorTye']
            names = seq2id[err_read]
            endErrorTye = self.error_class(err_tye)
            err_name_lst.extend(names)
            if err_count == len(names):
                for name in names:
                    if ori_file_type == 'fastq' or ori_file_type == 'fq' or ori_file_type == 'fastq.gz' or ori_file_type == 'fq.gz':
                        qual = {}
                        val = record_dict[name].letter_annotations['phred_quality']
                        val.sort()
                        position_val = val[int(len(val)/2)]
                        if endErrorTye == '0':
                            val[err_pos] = position_val
                        elif endErrorTye == '2':
                            del val[err_pos]
                        elif endErrorTye == '1':
                            val.insert(err_pos, position_val)
                        qual['phred_quality'] = val
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description, letter_annotations=qual)  
                    else:
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description)  
                    err2cor_records.append(tmp_rec)

        bases = orginal_file.split('/')[-1]
        prefix = bases.split('.')
        err_free_name_lst = list(set(total_name_lst) - set(err_name_lst))
        corrected_err_file = self.config.result_dir + prefix[0] + '_err2cor.' + ori_file_type
        with open(corrected_err_file, "w") as handle:
            SeqIO.write(err2cor_records, handle, ori_file_type)

        err_free_reads_file = self.config.result_dir + prefix[0] + '_errfree.' + ori_file_type
        extract_records(self.config.result_dir, err_free_name_lst, orginal_file, err_free_reads_file)
        
        corrected_file = self.config.result_dir + prefix[0] + '_corrected.' + ori_file_type
        os.system("cat %s %s > %s" % (corrected_err_file, err_free_reads_file, corrected_file))
        if os.path.exists(corrected_err_file):
            os.system("rm %s" % corrected_err_file)
        if os.path.exists(err_free_reads_file):
            os.system("rm %s" % err_free_reads_file)
        return corrected_file

    def error_class(self, errorType):
        """
        indicate substution or indel error given an error type 

        Args:
            errorType (str): incidate an error type

        Returns:
            str: a flag to indicate substution or indel error
        """
        substitutions = ['A-G', 'A-C', 'A-T','G-A', 'G-C', 'G-T','C-A', 'C-G', 'C-T','T-A', 'T-C', 'T-G', 'A-N', 'T-N','G-N','C-N','N-A','N-T', 'N-C', 'N-G']
        deletions = ['A-X', 'C-X', 'T-X', 'G-X', 'N-X']
        insertions = ['X-A', 'X-G', 'X-C', 'X-T', 'X-N']
        if errorType in substitutions:
            return '0'
        elif errorType in deletions:
            return '1'
        elif errorType in insertions:
            return '2'
        else:
            return '3'

    def all_in_one_high_ambiguous_err_correction(self, orginal_file, ambi_prediction_result):
        """
        high ambiguous error correction

        Args:
            orginal_file (str): raw data filename including path
            ambi_prediction_result (DataFrame): pandas dataframe save predicted result

        Returns:
            str: corrected data filename including path
        """
        grouped_ambiguous_df = ambi_prediction_result.groupby("idx")
        pred_df = pd.DataFrame()
        for name, group_df in grouped_ambiguous_df:
            # self.logger.debug(group_df['predictand'])
            if len(group_df) == 2:
                preds = group_df['predictand'].tolist()
                pred1 = preds[0]
                pred2 = preds[1]
                if abs(pred1-pred2) >= self.config.proba_deviation:
                    pred_index = group_df['predictand'].idxmax()
                    # pred_index = group_df['predictand'].idxmin()
                    # new_df = group_df.loc[:, ~group_df.columns.isin(['idx', 'predictand'])]
                    entry = group_df.loc[[pred_index]]
                    pred_df = pd.concat([pred_df, entry], ignore_index=True)
        if self.config.verbose:
            pred_df.to_csv(self.config.result_dir + 'high_ambiguous_1nt_prediction.csv', index=False)

        if len(pred_df) > 0:
            corrected_file = self.correct_high_ambiguous_errors(orginal_file, pred_df)
            if os.path.exists(orginal_file):     
                os.system("rm %s" % orginal_file)
            return corrected_file
        else:
            return orginal_file

    ######################################################################################################
    def all_in_one_ed2_correction(self, data_set, total_reads, genuine_df, negative_df, ambiguous_df):       
        """
        predicting and correcting ambiguous errors in 2nt-edit-distance-based read graph

        Args:
            data_set (str): raw data filename including path
            total_reads (set): set contains all the reads
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            negative_df (DataFrame): pandas dataframe containing negative samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction

        Returns:
            str: corrected data filename including path
        """
        ##self.MM.start()
        if not genuine_df.empty and not ambiguous_df.empty and not negative_df.empty:
            genuine_ambi_prediction_result, new_negative_df = self.all_in_one_ambiguous_err_prediction(total_reads, genuine_df, negative_df, ambiguous_df, high_ambiguous_df=None, edit_dis=2)
            ##self.MM.measure()
            correct_file = self.all_in_one_2nt_correct_errors(data_set, genuine_ambi_prediction_result)
            os.system("rm %s" % data_set)
            del genuine_ambi_prediction_result, new_negative_df
            self.logger.info("Error Correction finished.")
            ##self.MM.measure()
            ##self.MM.stop()
            return correct_file
        else:
            return data_set


    def all_in_one_2nt_correct_errors(self, orginal_file, df_data):
        """
        correcting ambiguous errors in 2nt-edit-distance-based read graph

        Args:
            orginal_file (str): raw data filename including path
            df_data (DataFrame): pandas dataframe save predicted result

        Returns:
            str: corrected data filename including path
        """
        record_dict, ori_file_type = parse_data_dict(orginal_file)
        seq2id = {}
        total_name_lst = []
        for record_id in record_dict:
            seq = str(record_dict[record_id].seq)
            seq2id.setdefault(seq, []).append(record_id)
            total_name_lst.append(record_id)
        err2cor_records = []
        err_name_lst = []
        for row_index, row in df_data.iterrows():
            err_count = row['EndReadCount']
            err_read = row['EndRead']
            cor_read = row['StartRead']
            names = seq2id[err_read]
            err_name_lst.extend(names)
            if err_count == len(seq2id[err_read]):
                for name in names:
                    if ori_file_type == 'fastq' or ori_file_type == 'fq' or ori_file_type == 'fastq.gz' or ori_file_type == 'fq.gz':
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description, letter_annotations=record_dict[name].letter_annotations)  
                    else:
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description)  
                    err2cor_records.append(tmp_rec)

        bases = orginal_file.split('/')[-1]
        prefix = bases.split('.')
        err_free_name_lst = list(set(total_name_lst) - set(err_name_lst))
        corrected_err_file = self.config.result_dir + prefix[0] + '_err2cor.' + ori_file_type
        with open(corrected_err_file, "w") as handle:
            SeqIO.write(err2cor_records, handle, ori_file_type)

        err_free_reads_file = self.config.result_dir + prefix[0] + '_errfree.' + ori_file_type
        extract_records(self.config.result_dir, err_free_name_lst, orginal_file, err_free_reads_file)
        corrected_file = self.config.result_dir + prefix[0] + '_corrected.' + ori_file_type

        os.system("cat %s %s > %s" % (corrected_err_file, err_free_reads_file, corrected_file))
        if os.path.exists(corrected_err_file):
            os.system("rm %s" % corrected_err_file)
        if os.path.exists(err_free_reads_file):
            os.system("rm %s" % err_free_reads_file)
        return corrected_file

    ######################################################################################################
    def umi_correction(self, data, genuine_df):
        """
        umi correction

        Args:
            data (str): raw data filename including path
            genuine_df (DataFrame): pandas dataframe containing genuine errors

        Returns:
            str: corrected data filename including path
        """
        corrected_file = self.config.result_dir + self.base[0] + '_corrected.' + self.out_file_tye
        # errors_df = self.ambiguous_err_prediction(genuine_csv, ambiguous_csv, negative_csv, edit_dis=1)
        # correct errors
        self.logger.info("Correcting 1nt-edit-distance based Errors")
        # genuine_df = pd.read_csv(genuine_csv)
        corrected_file = self.correct_errors(data, genuine_df)
        self.logger.info('1nt-edit-distance based Errors Correction Finished')
        return corrected_file
    
    def correct_high_ambiguous_errors(self, orginal_file, df_data):
        """
        correcting high ambiguous errors from raw dataset to generate the corrected dataset

        Args:
            orginal_file (str): raw data filename including path
            df_data (DataFrame): pandas dataframe save predicted result

        Returns:
            str: corrected data filename including path
        """
        record_dict, ori_file_type = parse_data_dict(orginal_file)
        seq2id = {}
        total_name_lst = []
        for record_id in record_dict:
            seq = str(record_dict[record_id].seq)
            seq2id.setdefault(seq, []).append(record_id)
            total_name_lst.append(record_id)
            
        err2cor_records = []
        err_name_lst = []
        
        name2flag = {}
        for row_index, row in df_data.iterrows():
            # err_count = row['EndReadCount']
            err_read = row['EndRead']
            # cor_read = row['StartRead']
            names = seq2id[err_read]
            # to set a flag to prevent one id' read to be changed to diffterent true reads multi times
            for name in names:
                name2flag[name] = True

        for row_index, row in df_data.iterrows():
            err_count = row['EndReadCount']
            err_read = row['EndRead']
            cor_read = row['StartRead']
            err_pos = row['ErrorPosition']
            endErrorTye = self.error_class(row['ErrorTye'])
            names = seq2id[err_read]

            if len(names) >= 1:
                # name = names[0]
                name = 'name'
                for i in names:
                    if name2flag[i]:
                        name = copy.deepcopy(i)
                        name2flag[i] = False
                        break
                    else:
                        continue
                self.logger.debug(name)
                if name != 'name':
                    err_name_lst.append(name)
                    # names.remove(name)
                    # seq2id[err_read] = names
                    if ori_file_type == 'fastq' or ori_file_type == 'fq' or ori_file_type == 'fastq.gz' or ori_file_type == 'fq.gz':
                        qual = {}
                        val = record_dict[name].letter_annotations['phred_quality']
                        val.sort()
                        position_val = val[int(len(val)/2)]
                        if endErrorTye == '0':
                            val[err_pos] = position_val
                        elif endErrorTye == '2':
                            del val[err_pos]
                        elif endErrorTye == '1':
                            val.insert(err_pos, position_val)
                        qual['phred_quality'] = val
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description, letter_annotations=qual) 
                    else:
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description)  
                    self.logger.debug(f'{name, err_read}')
                    err2cor_records.append(tmp_rec)
                else:
                    self.logger.warning("This high frequency reads hasn't enough reads to be corrected to predicted countparts, and will ingnore this prediction.")
                
        bases = orginal_file.split('/')[-1]
        prefix = bases.split('.')
        err_free_name_lst = list(set(total_name_lst) - set(err_name_lst))
        corrected_err_file = self.config.result_dir + prefix[0] + '_err2cor.' + ori_file_type
        self.logger.debug(f'Number of high ambiguous errors has been corrected: {len(err2cor_records)}')
        
        with open(corrected_err_file, "w") as handle:
            SeqIO.write(err2cor_records, handle, ori_file_type)

        err_free_reads_file = self.config.result_dir + prefix[0] + '_errfree.' + ori_file_type
        extract_records(self.config.result_dir, err_free_name_lst, orginal_file, err_free_reads_file)
        corrected_file = self.config.result_dir + prefix[0] + '_correct.' + ori_file_type
        os.system("cat %s %s > %s" % (corrected_err_file, err_free_reads_file, corrected_file))
        if os.path.exists(corrected_err_file):
            os.system("rm %s" % corrected_err_file)
        if os.path.exists(err_free_reads_file):
            os.system("rm %s" % err_free_reads_file)
        return corrected_file

    ######################################################################################################
    def correct_amplicon_err(self, original_file, genuine_df, negative_df, new_negative_df, amplicon_df, predict_proba):
        """
        correcting amplicon errors from dataset to generate the corrected dataset

        Args:
            original_file (str): data filename including path
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            negative_df (DataFrame): pandas dataframe containing negative samples for training
            new_negative_df (DataFrame): pandas dataframe containing negative samples from ambiguous prediction result
            amplicon_df (DataFrame): pandas dataframe containing amplicon errors samples for prediction
            predict_proba (float): the threshold to indicate whetheru a sample error-prone or error-free

        Returns:
            _type_: _description_
        """
        ##self.MM.start()
        study_name = 'amplicon_1nt'
        RV = Reads2Vectors(self.logger, self.config, edit_dis=1)
        # RV = Reads2Vectors(self.logger, self.config.num_workers, self.config.result_dir, self.config.read_max_len, self.config.entropy_kmer, self.config.entropy_q, self.config.kmer_freq, self.config.read_type, edit_dis=1)
        # train_data, train_labels, amplicon_data = RV.all_in_one_embedding(total_reads, genuine_df, negative_df, amplicon_df, high_flag=False)  
        train_data, train_labels, amplicon_data = RV.high_all_in_one_embedding(genuine_df, negative_df, new_negative_df, amplicon_df)
        del RV
        ##self.MM.measure()
        TM = MLClassifier(self.logger, self.config, study_name, train_data, train_labels, amplicon_data)
        predictions = TM.tunning(self.config.n_trials)
        del TM
        ##self.MM.measure()
        amplicon_df.insert(amplicon_df.shape[1], 'predictand', predictions)
        if self.config.verbose:
            amplicon_df.to_csv(os.path.join(self.config.result_dir, 'amplicon_predition.csv'), index = False)
        grouped_amplicon_df = amplicon_df.groupby("idx")
        pred_df = pd.DataFrame()
        for name, group_df in grouped_amplicon_df:
            # self.logger.debug(group_df)
            if len(group_df) == 1:
                pred = group_df['predictand'].values[0]
                if pred >= predict_proba:
                    pred_df = pd.concat([pred_df, group_df], ignore_index=True)
            elif len(group_df) >= 2:
                # high_index = pred_df = pd.concat([pred_df, group_df], ignore_index=True)['StartReadCount'].idxmax()
                pred_index = group_df['predictand'].idxmax()
                entry = group_df.loc[[pred_index]]
                # self.logger.debug(len(entry))
                pred = entry['predictand'].values[0]
                if pred >= predict_proba:
                    pred_df = pd.concat([pred_df, entry], ignore_index=True)
        # corrected_file = self.config.result_dir + self.base[0] + '_corrected.' + self.out_file_tye 
        del train_data, train_labels, amplicon_data
        ##self.MM.measure()
        ##self.MM.stop()
        if len(pred_df) > 1:
            correct_file = self.correct_errors(original_file, pred_df)
            if os.path.exists(original_file):     
                os.system("rm %s" % original_file)
            return correct_file
        else:
            return original_file

    def get_deduplication(self, corrected_data):
        record_iterator, file_type = parse_data(corrected_data)
        seqs2id_dict = {}
        for item in record_iterator:
            seq = str(item.seq)
            seqs2id_dict.setdefault(seq, []).append(str(item.id))
        
        seq_records = []
        record_num = 1
        for read, ids in seqs2id_dict.items():
            frequency = len(ids)
            record_id = f"read {record_num}, abundance: {frequency}"

            seq_record = SeqRecord(Seq(read), id=record_id, description=" ".join(ids))
            seq_records.append(seq_record) 
            record_num += 1    
        deduplicated_file = self.config.result_dir + self.base[0] + ".deduplicated.reads.fasta"     
        SeqIO.write(seq_records, deduplicated_file, "fasta")

    def umi_read_correction(self, isolates_file, non_isolates_file, genuine_df):
        """
        Args:
            isolates_file (str):  The isolates' filename with path.
            non_isolates_file (str):  The non-isolates' filename with path
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
        Returns:
            corrected_file (str): The corrected data filename including path.
        """
        file_type = parse_file_type(isolates_file)
        if self.config.correcting_umi:
            corrected_file = self.config.result_dir + self.base[0] + '_umi_corrected.' + file_type
        else:
            corrected_file = self.config.result_dir + self.base[0] + '_read_corrected.' + file_type  
        if not genuine_df.empty:
            self.logger.info("Correcting 1nt-edit-distance based Errors")
            non_isolates_correct = self.correct_errors(non_isolates_file, genuine_df)
            self.logger.info('1nt-edit-distance based Errors Correction Finished')
        else:
            self.logger.error("No genuine and ambiguous errors identified, failed to do error correction!")
            non_isolates_correct = non_isolates_file
            sys.exit(1)

        # # bcool correction
        IEC = IsolatesErrorCorrection(self.logger, self.config.num_workers, isolates_file, non_isolates_correct, self.config.result_dir, self.config.iso_change_detail)
        corrected_isolates = IEC.bcool_correct_isolates() 
        if corrected_isolates and non_isolates_correct:
            os.system("cat %s %s > %s" % (corrected_isolates, non_isolates_correct, corrected_file))
        else:
            self.logger.error("No corrected_isolates and/or non_isolates_correct, failed to do error correction")
        del IEC

        if os.path.exists(isolates_file):
            os.system("rm %s" % isolates_file)
        if os.path.exists(non_isolates_file):
            os.system("rm %s" % non_isolates_file)
        if os.path.exists(corrected_isolates):   
            os.system("rm %s" % corrected_isolates)
        if os.path.exists(non_isolates_correct):     
            os.system("rm %s" % non_isolates_correct)

        return corrected_file

    def umi_read_2nt_correction(self, data_set, genuine_df):       
        """
        correcting ambiguous errors in 2nt-edit-distance-based read graph

        Args:
            data_set (str): raw data filename including path
            genuine_df (DataFrame): pandas dataframe containing positive samples for training
            ambiguous_df (DataFrame): pandas dataframe containing ambiguous samples for prediction

        Returns:
            str: corrected data filename including path
        """
        ##self.MM.start()
        if not genuine_df.empty:
            correct_file = self.all_in_one_2nt_correct_errors(data_set, genuine_df)
            os.system("rm %s" % data_set)
            del genuine_df
            self.logger.info("2nt-edit-distance based Errors Correction Finished")
            return correct_file
        else:
            return data_set

'''
    def correct_amplicon_errors(self, orginal_file, df_data):
        """
        correcting 2nt ambiguous errors from raw dataset to generate the corrected dataset

        Args:
            orginal_file (str): raw data filename including path
            df_data (DataFrame): pandas dataframe save predicted result

        Returns:
            str: corrected data filename including path
        """
        record_dict, ori_file_type = parse_data_dict(orginal_file)
        seq2id = {}
        total_name_lst = []
        for id in record_dict:
            seq = str(record_dict[id].seq)
            seq2id.setdefault(seq, []).append(id)
            total_name_lst.append(id)
        err2cor_records = []
        err_name_lst = []
        for row_index, row in df_data.iterrows():
            err_count = row['EndReadCount']
            err_read = row['EndRead']
            cor_read = row['StartRead']
 
            names = seq2id[err_read]
            err_name_lst.extend(names)
            if err_count == len(seq2id[err_read]):
                for name in names:
                    if ori_file_type == 'fastq' or ori_file_type == 'fq' or ori_file_type == 'fastq.gz' or ori_file_type == 'fq.gz':
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description, letter_annotations=record_dict[name].letter_annotations)  
                    else:
                        tmp_rec = SeqRecord(Seq(cor_read), id=name, description=record_dict[name].description)  
                    err2cor_records.append(tmp_rec)

        bases = orginal_file.split('/')[-1]
        prefix = bases.split('.')
        err_free_name_lst = list(set(total_name_lst) - set(err_name_lst))
        corrected_err_file = self.config.result_dir + prefix[0] + '_err2cor.' + self.out_file_tye
        with open(corrected_err_file, "w") as handle:
            SeqIO.write(err2cor_records, handle, self.out_file_tye)

        err_free_reads_file = self.config.result_dir + prefix[0] + '_errfree.' + self.out_file_tye
        extract_records(self.config.result_dir, err_free_name_lst, orginal_file, err_free_reads_file)
        corrected_file = self.config.result_dir + prefix[0] + '_correct.' + self.out_file_tye
        os.system("cat %s %s > %s" % (corrected_err_file, err_free_reads_file, corrected_file))
        return corrected_file
'''