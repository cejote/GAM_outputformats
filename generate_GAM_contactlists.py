#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import os

#chrom_list_mouse = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
chrom_list_human = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', "chrY"]
chrom_list = chrom_list_human


def calculate_NPMI(seg_matrix_1, seg_matrix_2):
    #seg_matrix_1 will be on Y axis
    #seg_matrix_2 will be on X axis

    np.seterr(invalid='ignore',divide='ignore')

    # calculate M the number of samples
    M = len(seg_matrix_1[0])

    #calculate p(x,y) matrix i.e coseg matrix divided by M
    pxy = (seg_matrix_1.dot(seg_matrix_2.T))/M

    #p(x)p(y) matrix  - i.e detection frequency with the dot product of it's transposed partner
    pxpy = seg_matrix_1.sum(1).reshape(-1,1)/M * seg_matrix_2.sum(1)/M

    #define PMI matrix
    PMI = np.log2(pxy/pxpy)

    #bound values between -1 and 1
    NPMI = PMI/-np.log2(pxy)

    return NPMI


def calculate_cosegregation(seg_matrix_1, seg_matrix_2):
    # calculate M the number of samples
    M = len(seg_matrix_1[0])
    #calculate p(x,y) matrix i.e coseg matrix divided my M
    pxy = (seg_matrix_1.dot(seg_matrix_2.T))/M
    return pxy


def calculate_NPMI_cis(seg_matrix):
    # calculate M the number of samples
    M = len(seg_matrix[0])
    #calculate p(x,y) matrix i.e coseg matrix divided my M
    pxy = (seg_matrix.dot(seg_matrix.T))/M
    #p(x)p(y) matrix  - i.e detection frequency with the dot product of it's transposed self equals an N*N matrix
    pxpy = seg_matrix.sum(1).reshape(-1,1)/M * seg_matrix.sum(1)/M
    #define PMI matrix
    PMI = np.log2(pxy/pxpy)
    #bound values between -1 and 1
    NPMI = PMI/-np.log2(pxy)
    return NPMI

def calculate_cosegregation_cis(seg_matrix):
    # calculate M the number of samples
    M = len(seg_matrix[0])
    #calculate p(x,y) matrix i.e coseg matrix divided my M
    pxy = (seg_matrix.dot(seg_matrix.T))/M
    return pxy



def correlate(args):

    if len(args)==4:
        (segregation_table, matrix_type, chrom1, outpath) = args
        chrom2=chrom1
    if len(args)==5:
        (segregation_table, matrix_type, chrom1, chrom2, outpath) = args
    else:
        print ('You need to provide the path to segregation table, type of the matrices (npmi or coseg), chrom1  (cis) or chrom1 and chrom2 (trans), and the destination path')
        sys.exit(-1)

    print(args)

    #@todo: should set to lower case
    if not chrom1 in chrom_list:
        print ("Unrecognized chromosome 1: ", chrom1)
        sys.exit(-1)
    if not chrom2 in chrom_list:
        print ("Unrecognized chromosome 2: ", chrom2)
        sys.exit(-1)

    if matrix_type=='NPMI':
        matrix_function = calculate_NPMI
    elif matrix_type=='coseg':
        matrix_function = calculate_cosegregation
    else:
        print("This script does only support NMPI and coseg for correlation.")
        sys.exit(-1)


    #if organism == 'mouse' or organism == 'm':
    #    chromosomes = chrom_list_mouse
    #if organism == 'human' or organism == 'h':
    #    chromosomes = chrom_list_human
    
    seg_table_1 = pd.read_csv(segregation_table, sep='\t')
    seg_table_1.set_index(['chrom','start','stop'], inplace = True)

    seg_table_2 = pd.read_csv(segregation_table, sep='\t')
    seg_table_2.set_index(['chrom','start','stop'], inplace = True)
    #seg_table_2 = seg_table_1
    
    #for chrom1 in chromosomes:
    #get data
    subset_segregation_1 = seg_table_1.loc[chrom1]
    seg_matrix_1 = subset_segregation_1.values
    subset_segregation_1.reset_index(inplace = True)
    # Get coordinates for index and header
    coord_1 = chrom1 + ':' + subset_segregation_1['start'].astype(str) + '-' + subset_segregation_1['stop'].astype(str)

    #for chrom2 in chromosomes:
    #get data
    subset_segregation_2 = seg_table_2.loc[chrom2]
    seg_matrix_2 = subset_segregation_2.values
    subset_segregation_2.reset_index(inplace = True)
    # Get coordinates for index and header
    coord_2 = chrom2 + ':' + subset_segregation_2['start'].astype(str) + '-' + subset_segregation_2['stop'].astype(str)

    if matrix_type=='NPMI':
        matrix = matrix_function(seg_matrix_1, seg_matrix_2)
    elif matrix_type=='coseg':
        matrix = matrix_function(seg_matrix_1, seg_matrix_2)

    #print(coord_1)
    #print(coord_2)
    matrix_df = pd.DataFrame(matrix,index=coord_1, columns=coord_2) 
    #npmi_df.to_csv (outpath + matrix_type + '.' + str(segregation_table) + '.' + chrom + '.txt.gz', sep='\t', compression='gzip', na_rep='NA')
    ##### matrix_df.to_csv (outpath + 'NPMI.' + os.path.basename(segregation_table) + '.' + chrom1 +"_" + chrom2 + '.txt.gz', sep='\t', compression='gzip', na_rep='NA')

    #print(type(coord_1))
    #print(coord_1.values)
    #print(matrix_df)
    matrix_df["ids1"]=coord_1.values
    matrix_df_long = matrix_df.melt(id_vars=['ids1'])

    #exclude NA values    
    matrix_df_long= matrix_df_long[matrix_df_long['value'].notnull()]
    
    #split coordinate strings
    #print(matrix_df_long)
    #@TODO: should check for non-empty df here
    matrix_df_long[['chrom_1', 'start_1', 'end_1']] = matrix_df_long['ids1'].str.split('[:-]',expand=True)
    matrix_df_long[['chrom_2', 'start_2', 'end_2']] = matrix_df_long['variable'].str.split('[:-]',expand=True)
    
    matrix_df_long[['start_1', 'end_1', 'start_2', 'end_2']] = matrix_df_long[['start_1', 'end_1', 'start_2', 'end_2']].astype(int)
    #print(matrix_df_long.dtypes)
    #print(matrix_df_long.head())

    #remove diag and lower matrix
    matrix_df_long = matrix_df_long[(matrix_df_long.start_2) > (matrix_df_long.start_1)]

    #@todo: I should use sprintf here...
    fname=outpath +'/'+ matrix_type + '.' + os.path.basename(segregation_table) + '.' + chrom1 +"_" + chrom2 + '.long.gz'
    print(fname)

    matrix_df_long[['chrom_1', 'start_1', 'end_1', 'chrom_2', 'start_2', 'end_2', 'value']].to_csv (fname, sep='\t', compression='gzip', na_rep='NA', header=False, index=False)
    print("done "+ chrom1 + " " + chrom2)


if __name__ == "__main__":
    correlate(sys.argv[1:])



