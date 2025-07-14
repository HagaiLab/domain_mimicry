import pandas as pd
import os
import numpy as np

'''
This script matches CSU/Scannet results (surface, core or interface) with positive selection analysis at the residues level (PSR)- based on PAML results.
We take structural information of mimicked and target proteins and overlay them with PSR information, yielding percentages of positively selected residues (%PSR) in each each category- core, surface or interface. 
'''
output_path = r"output/"
convert_path = r"input/ensembl_ensg_uniprot_biomart.txt"
convert_df = pd.read_csv(convert_path, sep= '\t')
paml_stats_path = r"input/paml_stats_and residues_with_uniprot.csv"
mimicked_csu_path = r"input/csu_lists/mimicked/"
target_csu_path = r"input/csu_lists/target/"
overmasked_list_path = r"input/overmasked_ensembl_Ns_gaps.txt"
overmasked_genes = []
with open(overmasked_list_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        ensembl_id = line.strip('\n')
        overmasked_genes.append(ensembl_id)

def match_csu(csu_path):
    csu_results = os.listdir(csu_path)
    total_pos_sel_i_scanet_counter = 0
    total_scanet_i = 0
    total_s_counter = 0
    total_i_counter = 0
    total_c_counter = 0
    total_pos_sel_s_counter = 0
    total_pos_sel_i_counter = 0
    total_pos_sel_c_counter = 0
    total_pos_sel_counter = 0
    total_length = 0
    protein_count = 0

    for result in csu_results:
        ensembl = result.split('_')[0]
        if ensembl in overmasked_genes:
            print('overmasked')
            continue
        csu_file_path = os.path.join(csu_path, result)
        with open(csu_file_path) as file:
            csu_list = file.readlines()

        for csu in csu_list: #checking if residue is interface, surface, or core in CSU file
            if csu.strip('\n') == 'I':
                total_i_counter += 1
            elif csu.strip('\n') == 'S':
                total_s_counter += 1
            elif csu.strip('\n') == 'C':
                total_c_counter += 1

        protein_length = len(csu_list)
        total_length += protein_length

        # ScanNet analysis- adding general interface information
        if 'mimicked' in csu_path:
            scanet_data_path = r"input/scannet_interface_results/mimicked_0.7.csv"
        elif 'target' in csu_path:
            scanet_data_path = r"input/scannet_interface_results/target_0.7.csv"

        scanet_df = pd.read_csv(scanet_data_path)
        uniprot = convert_df.loc[
            convert_df['Gene stable ID'] == ensembl, 'UniProtKB/Swiss-Prot ID'].values[0]
        filtered_scanet_df = scanet_df[scanet_df['uniprot'] == uniprot]
        scanet_i_counter = len(filtered_scanet_df)
        total_scanet_i += scanet_i_counter

        # overlaying with positive selection information
        paml_stats_df = pd.read_csv(paml_stats_path)
        pos_sel_residues_uniprot_indices_list = paml_stats_df.loc[paml_stats_df['genes'] == ensembl, 'pos_sel_residues_uniprot'].values
        if len(pos_sel_residues_uniprot_indices_list) == 0: # no PSR
            continue

        pos_sel_residues_uniprot_indices = pos_sel_residues_uniprot_indices_list[0]
        if pos_sel_residues_uniprot_indices is None or (isinstance(pos_sel_residues_uniprot_indices, float) and np.isnan(pos_sel_residues_uniprot_indices)):
            continue

        protein_count += 1

        # at least one PSR identified
        pos_sel_residues_list = pos_sel_residues_uniprot_indices.split(',')
        total_pos_sel_counter += len(pos_sel_residues_list)

        for res in pos_sel_residues_list:
            pos_sel_index = int(res.split('_')[0])
            res_symbol = res.split('_')[1]
            probability = res.split('_')[2]


            matching_csu_result = csu_list[pos_sel_index-1].strip('\n')
            if matching_csu_result == 'I':
                total_pos_sel_i_counter += 1
                print(f'interface pos sel residue for {ensembl} in {pos_sel_index}, residue is {res_symbol} with probability of {probability}')

            elif matching_csu_result == 'S':
                total_pos_sel_s_counter += 1

            # Find positively selected residues in ScaNet
            if pos_sel_index in list(filtered_scanet_df['Residue Index']):
                print(f'SCannet interface pos sel residue for {ensembl} in {pos_sel_index}, residue is {res_symbol} with probability of {probability}')
                print(filtered_scanet_df)
                total_pos_sel_i_scanet_counter += 1




    # calculating PSR percentages for each group
    pos_sel_percent_total = (total_pos_sel_counter/total_length) * 100
    pos_sel_percent_c = (total_pos_sel_c_counter/total_c_counter) * 100
    pos_sel_percent_s = (total_pos_sel_s_counter/total_s_counter) * 100
    pos_sel_percent_i = (total_pos_sel_i_counter/total_i_counter) * 100
    pos_sel_percent_scanet_i = (total_pos_sel_i_scanet_counter/total_scanet_i) * 100


    column_names = ['total', 'csu_c' ,'csu_S', 'csu_I', 'scanet_I']
    total_count_list = [total_length, total_c_counter, total_s_counter, total_i_counter, total_scanet_i]
    pos_sel_count_list = [total_pos_sel_counter, total_pos_sel_c_counter, total_pos_sel_s_counter, total_pos_sel_i_counter, total_pos_sel_i_scanet_counter]
    pos_sel_percent_list = [pos_sel_percent_total, pos_sel_percent_c, pos_sel_percent_s, pos_sel_percent_i, pos_sel_percent_scanet_i]
    final_df = pd.DataFrame({
        'Category': column_names,
        'Total Count': total_count_list,
        'Pos Sel Count': pos_sel_count_list,
        '%PSR' : pos_sel_percent_list
    })


    print(f'there are {protein_count} proteins in this group')
    print(final_df)
    if 'mimicked' in csu_path:
        outpath = os.path.join(output_path, 'mimicked_PSR.csv')
    elif 'target' in csu_path:
        outpath = os.path.join(output_path, 'target_PSR.csv')


    final_df.to_csv(outpath, index= False)



print('mimicked')
match_csu(mimicked_csu_path)

print('target')
match_csu(target_csu_path)

