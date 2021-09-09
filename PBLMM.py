'''
# PBLMM - Peptide based linear mixed models for differential expression analysis
# Copyright (C) 2021 Kevin Klann
#This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import sys
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import warnings
from statsmodels.stats.multitest import multipletests, local_fdr
import json
import time
import os

warnings.filterwarnings("ignore")

def normalization(input_file, channels):
    '''
    #Performs total intensity normalisation. Besides input file the function needs an array of all
    column names that contain the quantifications to be normalized (channels).
    '''
    # remove missing value rows
    #input_file = input_file.dropna(subset=channels)
    print("Normalization")
    # calculate summed intensity for each column and search minimum index
    minimum = np.argmin(input_file[channels].sum().values)
    summed = np.array(input_file[channels].sum().values)
    minimum = summed[minimum]
    # calculuate norm factors
    norm_factors = summed/minimum
    # normalize input with norm factors
    input_file[channels] = input_file[channels].divide(
        norm_factors, axis=1)
    print("Normalization done")
    return input_file

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

    
def sum_peptides_for_proteins(input_file, channels, mpa1):
    '''
    This function takes Peptide level (or PSM) dataframes and performs a sum based rollup to protein level.
    the channels variable takes an array of column names that contain the quantifictions. You can create such an
    array via this command:
    channels = [col for col in PSM.columns if 'Abundance:' in col]
    mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
    and use it for rollup.
    Returns Protein level DF.
    '''
   
    
    print('Calculate Protein quantifications from PSM')
    mpa = [col for col in input_file.columns if mpa1 in col]
    mpa = mpa[0]
    PSM_grouped = input_file.groupby(by=[mpa])
    result = {}
    for group in PSM_grouped.groups:
        temp = PSM_grouped.get_group(group)
        sums = temp[channels].sum()
        result[group] = sums
    protein_df = pd.DataFrame.from_dict(
        result, orient='index', columns=channels)
    print("Combination done")
    return protein_df


def tessa(source):
    result = []
    for p1 in range(len(source)):
            for p2 in range(p1+1,len(source)):
                    result.append([source[p1],source[p2]])
    return result

def peptide_based_lmm(self, input_file,labels, conditions,drop_missing=False, techreps=None, plexes=None, norm=normalization):

    columns = labels
    self.pair_names = []
    channels = [col for col in input_file.columns if columns[2] in col]
    if norm is not None:
        input_file = norm(input_file, channels)
    else:
        if drop_missing == True:
            input_file = input_file.dropna(subset=channels)
        else:
            pass
    
    # Protein level quantifications
    
    protein_data = sum_peptides_for_proteins(
        input_file=input_file, channels=channels)
    # Prepare Peptide data for LMM
    Peptides_for_LM = input_file[channels]
    sequence = [col for col in input_file.columns if columns[0] in col]
    sequence = sequence[0]
    Peptides_for_LM['Sequence'] = input_file[sequence]
    Acc = [col for col in input_file.columns if columns[1] in col]
    Acc = Acc[0]
    Peptides_for_LM['Accession'] = input_file[Acc]
    melted_Peptides = Peptides_for_LM.melt(
        id_vars=['Accession', 'Sequence'], value_vars=channels)
    # Replace column names with conditions
    
    
    if  techreps == None:
        pass
    else:
        melted_Peptides['Techreps']=melted_Peptides['variable']
        melted_Peptides['Techreps'].replace(to_replace=channels,
                            value=techreps, inplace=True)
    
    if  plexes == None:
        pass
    else:
        melted_Peptides['Multiplex']=melted_Peptides['variable']
        melted_Peptides['Multiplex'].replace(to_replace=channels,
                            value=plexes, inplace=True)
    melted_Peptides['variable'].replace(to_replace=channels,
                            value=conditions, inplace=True)    
    unique_conditions = list(set(conditions))
    
    pairs = tessa(unique_conditions)
    
    for pair in pairs:
        pair.sort()
        temp = melted_Peptides[(melted_Peptides['variable'].str.fullmatch(pair[0])) | (
            melted_Peptides['variable'].str.fullmatch(pair[1]))]
        temp['value'] = np.log2(temp['value'])
        temp = temp.dropna()
        grouped = temp.groupby(by=['Accession'])
        result_dict = {}
        fold_changes = []
        counter = 0
        for i in grouped.groups:
            temp2 = grouped.get_group(i)
            vc = {'Sequence': '0+Sequence'}
            #Base model
            model_form = "value ~ variable"
            #Extent model based on data
            if techreps is not None:
                
                vc['Techreps'] = '0+C(Techreps)'
            else:
                pass
            
        
            if plexes is not None:
                
                vc['Multiplex'] = '0+C(Multiplex)'
            else:
                pass
            
            model = smf.mixedlm(
                model_form, temp2, groups='Sequence', vc_formula=vc)
            try:
                result = model.fit()
 
                fc = result.params[1]
                pval = result.pvalues[1]
                fold_changes.append(fc)
                result_dict[i] = pval
            except:
                pass
        result_df_peptides_LMM = pd.DataFrame.from_dict(
            result_dict, orient='index', columns=['p_value'])
        result_df_peptides_LMM['fold_change'] = np.array(fold_changes)
        # Multiple testing correction:
        result_df_peptides_LMM['p_value'] = result_df_peptides_LMM['p_value'].fillna(
            value=1)
        pvals = result_df_peptides_LMM['p_value'].to_numpy()
        reject, pvals_corrected, a, b = multipletests(
            pvals, method='fdr_bh')
        result_df_peptides_LMM['q_value'] = pvals_corrected
        comparison = str(pair[0]) + '_' + str(pair[1])
        result_df_peptides_LMM = result_df_peptides_LMM.add_suffix(
            comparison)
        protein_data = protein_data.join(result_df_peptides_LMM)
   
    return protein_data

def main():
    settings_path = sys.argv[1]
    path = settings_path
    if not os.path.exists(os.path.join(path,"./Results")):
        os.mkdir(os.path.join(path,"./Results"))
    defaults =[
            'Annotated Sequence',
            'Master Protein Accessions',
            'Abundance:',
            ]
    
    with open(os.path.join(settings_path,'settings.json'), 'r') as infile:
        data = json.load(infile)
    
    
    input_file = pd.read_csv(data['file1'],sep='\t',header=0)
    conditions=[]
    design_matrix=pd.read_csv(data['file2'],sep='\t',header=0)
    conditions = list(design_matrix['conditions'])
    try:
        techreps = list(design_matrix['techreps'])
    except KeyError:
        techreps = None
    try:
        multiplex = list(design_matrix['multiplex'])
    except KeyError:
        multiplex = None
    
     


    labels = [data['seq_col'],data['acc_col'],data['abun_col']]
    
    for i in range(0,len(labels)):
        if labels[i] == '':
            labels[i] = defaults[i]
        else:
            pass
    if data['norm']=='Yes':
        normal = normalization
    else:
        normal = None
    result = peptide_based_lmm(input_file,conditions,labels=labels,norm=normal, plexes=multiplex,techreps=techreps)
    timestr=time.strftime("%Y%m%d-%H%M%S")
    result.to_csv(os.path.join(path,"./Results/")+timestr+'_LMM_Result.txt',sep='\t')
    return_data = {}
    unique_conditions = list(set(conditions))
    pairs = tessa(unique_conditions)

    #Parse data to return to app for visulaization.
    result=result.dropna()
    for pair in pairs:
        pair.sort()
        inner_dict = {}
        comparison = str(pair[0])+ '_' + str(pair[1])
        p_string = 'P value'+ comparison
        column_p = [col for col in result.columns if p_string in col]
        
        inner_dict['p_value']=result[column_p].to_numpy().flatten().tolist()
        q_string = 'corrected P value (q value)'+ comparison
        column_q = [col for col in result.columns if q_string in col]
        inner_dict['q_value']=result[column_q].to_numpy().flatten().tolist()
        fc_string = 'fold_change'+ comparison
        column_fc = [col for col in result.columns if fc_string in col]
        inner_dict['fold_change']=result[column_fc].to_numpy().flatten().tolist()
        return_data[comparison] = inner_dict
        
    print(json.dumps(return_data))

if __name__ == "__main__":
    main()
