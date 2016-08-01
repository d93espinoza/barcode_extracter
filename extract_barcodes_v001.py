__author__ = 'espinozada'
import pandas as pd
import numpy as np
import glob
import time
import copy
import cPickle
import multiprocessing as mp
from functools import partial


'''
this function will create a dictionary
of the barcodes in a fastq file with
reads over 100.
'''
def extractBarcodes(fastq, libID, threshold):
    print "Extracting " + str(fastq)
    file_stream = open(fastq,'r')
    nonthresh_dict = {}
    id_length = len(libID)
    reads = 0
    mapped = 0
    for i, line in enumerate(file_stream):
        if (i % 4 == 1):
            reads += 1
            if line[:id_length] == libID:
                mapped += 1
                if line[:50] in nonthresh_dict:
                    nonthresh_dict[line[:50]] += 1
                else:
                    nonthresh_dict[line[:50]] = 1
    file_stream.close()
    readmeinfo = (fastq, mapped, reads, (100*float(mapped)/reads),  threshold, libID)
    thresh_dict = {}
    for key in nonthresh_dict:
        if nonthresh_dict[key] > threshold:
            thresh_dict[key] = nonthresh_dict[key]
    return (thresh_dict, fastq, readmeinfo)


'''
this function takes two sequences and the start indices
and current number of errors, along with the number of corrects
and length of barcode and tries to match the next index...
uses recursion to return True if two barcodes match with X bp
mismatch/indels, false otherwise
'''
def zR3(seq1, seq2, seq1ind, seq2ind, n_err, totalmatched, barcodelength, maxmismatch):

    if n_err > n_err:
            exit("Critical Error: zR3 saw 3 mismatches as OK")


    if int(totalmatched) == int(barcodelength):
        return True

    if seq1[seq1ind] == seq2[seq2ind]:
        return zR3(seq1,seq2,seq1ind+1, seq2ind+1,n_err,totalmatched+1, barcodelength, maxmismatch)

    else:

        n_err += 1
        if n_err > maxmismatch:
            return False

        return zR3(seq1, seq2, seq1ind+1, seq2ind, n_err, totalmatched, barcodelength, maxmismatch) or \
               zR3(seq1, seq2, seq1ind, seq2ind+1, n_err, totalmatched, barcodelength, maxmismatch) or \
               zR3(seq1, seq2, seq1ind+1, seq2ind+1, n_err,totalmatched+1, barcodelength, maxmismatch)


'''
this function checks to see if the given sets are pairwise disjoint
'''
def all_disjoint(sets):
    all = set()
    for s in sets:
        for x in s:
            if x in all:
                return False
            all.add(x)
    return True


def barcode_edge_finder(index, barcodes):
    if index % 1000 == 0:
        print str(index) + " done"
    num_barcodes = len(barcodes)
    returnset = {index}
    for i in xrange(index+1, num_barcodes):
        if zR3(barcodes[index], barcodes[i], 6, 6, 0, 0, barcode_length, max_mismatches):
            returnset.add(i)
    return returnset


"""

|7~~~~~~~~~~~~~~~~~~~~O\_AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
|
 VVVVVVVVVVVVVVVVVVV\      |      |
|AAAAAAAAAAAAAAAAAAA/   ___\      /____
|~~~~~~~~~~~~~~~~~~~~~~/ /_/_/_/__>


"""
def combine_barcodes(initial_dataframe, libID_length, barc_length, m_mismatches, cores):
    all_barcodes = initial_dataframe.index.tolist()
    partial_func = partial(barcode_edge_finder, barcodes = all_barcodes)
    number_barcodes = len(all_barcodes)
    pool = mp.Pool(processes = cores)

    print("Building barcode graph structure")
    bccc_start = time.time()
    index_sets = pool.map(partial_func, xrange(0,len(all_barcodes)))
    bccc_end = time.time()
    print "Barcode graph structure built. Took: " + str(bccc_end-bccc_start) + " seconds"

    print "Finding connected components"
    cc_start = time.time()
    disjoint_index_sets = make_sets_disjoint(index_sets, xrange(0,number_barcodes))
    cc_end = time.time()
    print "Connected components found. Took: " + str(cc_end - cc_start) + " seconds"


    print "Now combining connected components"
    dfcr_start = time.time()
    pool2 = mp.Pool(processes = cores)
    partial_make_consensus = partial(make_consensus, your_data = initial_dataframe)
    list_of_series = pool2.map(partial_make_consensus, disjoint_index_sets)
    new_dataframe = pd.concat(list_of_series, axis = 1).transpose()
    dfcr_end = time.time()
    print "CC combined and dataframe built. Took: " + str(dfcr_end-dfcr_start) + " seconds"
    return new_dataframe


def make_consensus(your_setofbar, your_data):
    consensus_sequence = your_data.ix[your_setofbar].sum(axis = 1).order(ascending=False).index[0]
    temp_series = your_data.ix[your_setofbar].sum(axis = 0)
    temp_series.name = consensus_sequence
    return temp_series


'''
takes in a list of sets and the list of barcodes and combines
the list of sets into a set of disjoint sets... derived from
the theory of connected components in graphs
 /\   /\
/  \_/  \
| O   O |          __/\
\ = < = /         / __/
 <  ^  >__________\ \
  \       ) ) ) )    \
  /      ( ( ( (      |
  |  /____________\  /
  |  ||  |     |  |  |
 <<< |><_|>   <<<<<<_|>
this is O(n^2) so might want to fix in the future version
JK this is rather quick! not the bottleneck... bottleneck at the zR3 function
'''
def make_sets_disjoint(listofsets, listofindices):
    listy = copy.deepcopy(listofsets)
    for index in listofindices:
        disjointlistofsets = []
        newset = set()
        for setty in listy:
            if index in setty:
                newset = newset.union(setty)
            else:
                disjointlistofsets.append(setty)
        disjointlistofsets.append(newset)
        listy = disjointlistofsets
    assert all_disjoint(listy)
    return listy




if __name__ == "__main__":

    previous_boolean = raw_input("Does a previous file exist? (Y/N): ")
    if(previous_boolean == "Y"):
        old_file_name = raw_input("Enter previous pickled file name: ")
        old_file = open(old_file_name,'r')
        old_dictlist = cPickle.load(old_file)
        previous_files = [x[1] for x in old_dictlist]
        fileList = []
        for file in glob.glob("*.fastq"):
            if file not in previous_files:
                fileList.append(file)
    elif(previous_boolean == "N"):
        fileList = []
        for file in glob.glob("*.fastq"):
                fileList.append(file)
    else:
        exit("Must choose Y or N.")

    userlibID = raw_input("Please enter library ID sequence: ")
    barcode_length = int(raw_input("Please enter barcode length: "))
    input_thresh = int(raw_input("Please enter discrete threshold: "))
    newfile = raw_input("Please enter name of extracted barcode pickled output file: ")

    dictpool = mp.Pool(processes=4)
    partial_extract = partial(extractBarcodes, libID = userlibID, threshold = input_thresh)
    new_dictlist = dictpool.map(partial_extract, fileList)
    if(previous_boolean == "Y"):
        new_dictlist = new_dictlist + old_dictlist

    file_dump = open(newfile,'w')
    cPickle.dump(new_dictlist, file_dump)
    file_dump.close()



'''
    readme.write("Library ID used: " + str(userlibID) + "\n")
    readme.write("Threshold used: " + str(input_thresh) + "\n")
    readme.write("Barcode Length used: " + str(barcode_length) + "\n")
    readme.write("Max_mismatches used: " + str(max_mismatches) + "\n")
    readme.write("Total barcodes before combination: " + str(dataframe.shape[0]) + "\n")
    readme.write("FILENAME" + '\t' + "MAPPED" + '\t' + "READS" + '\t' + "MAP %" + '\n')
    for t in list_of_uncombined_dicts_runinfo:
        readme.write(str(t[0]) + '\t' + str(t[1]) + '\t' + str(t[2]) + '\t' + str(t[3]) + '\n')
    print "Total barcodes before combination: " + str(dataframe.shape[0]) + "\n"
    dataframe = combine_barcodes(dataframe, len(userlibID), barcode_length, max_mismatches, 4)
    readme.write("Total barcodes after combination: " + str(dataframe.shape[0]) + "\n")
    readme.close()
    print "Total barcodes after combination: " + str(dataframe.shape[0]) + "\n"

    #prints new dataframe to tab delimited file
    dataframe.to_csv(newfile, sep = '\t')

'''





