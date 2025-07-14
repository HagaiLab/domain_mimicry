import pandas as pd


'''
Reads the PAML statistics file and generates a file showing the positive selection percentages (%PSG) and numbers for each protein group,
using four different q-value thresholds for positive selection (0.0001, 0.001. 0.01, 0.05).
'''

paml_stats_path = r'input/paml_stats.csv'
paml_stats_df = pd.read_csv(paml_stats_path)

paml_stats_df_no_overmasking = paml_stats_df[paml_stats_df['overmasked_by_Ns_and_gaps'] == '-']
print(paml_stats_df_no_overmasking)

total_genes = len(paml_stats_df_no_overmasking)
total_mimicked = len(paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['mimicked'] == '+'])
total_target = len(paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['target'] == '+'])
total_viral_binding = len(paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['viral_binding'] == '+'])
total_func_mimicked =  len(paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['func_mimicked'] == '+'])
total_dmi_mimicked =  len(paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['dmi_mimicked'] == '+'])
total_dmi_target =  len(paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['dmi_target'] == '+'])


q_values = [0.0001, 0.001, 0.01, 0.05]
pos_sel_percentages = []
pos_sel_percentages_mimicked = []
pos_sel_percentages_target = []
pos_sel_genes_list = []
pos_sel_mimicked_list = []
pos_sel_target_list = []
pos_sel_viral_binding_list = []
pos_sel_percentages_viral_binding = []
pos_sel_func_mimicked_list = []
pos_sel_percentages_func_mimicked = []
pos_sel_dmi_mimicked_list = []
pos_sel_percentages_dmi_mimicked = []
pos_sel_dmi_target_list = []
pos_sel_percentages_dmi_target = []


for qvalue in q_values:
    filtered_df = paml_stats_df_no_overmasking[paml_stats_df_no_overmasking['q_values'] < qvalue]
    pos_sel_genes = len(filtered_df)
    pos_sel_genes_list.append(pos_sel_genes)
    pos_sel_mimicked = len(filtered_df[filtered_df['mimicked'] == '+'])
    pos_sel_mimicked_list.append(pos_sel_mimicked)
    pos_sel_target = len(filtered_df[filtered_df['target'] == '+'])
    pos_sel_target_list.append(pos_sel_target)
    pos_sel_viral_binding = len(filtered_df[filtered_df['viral_binding'] == '+'])
    pos_sel_viral_binding_list.append(pos_sel_viral_binding)
    pos_sel_func_mimicked = len(filtered_df[filtered_df['func_mimicked'] == '+'])
    pos_sel_func_mimicked_list.append(pos_sel_func_mimicked)
    pos_sel_dmi_mimicked = len(filtered_df[filtered_df['dmi_mimicked'] == '+'])
    pos_sel_dmi_mimicked_list.append(pos_sel_dmi_mimicked)
    pos_sel_dmi_target = len(filtered_df[filtered_df['dmi_target'] == '+'])
    pos_sel_dmi_target_list.append(pos_sel_dmi_target)

    pos_sel_percent = (pos_sel_genes/total_genes)*100
    pos_sel_percentages.append(pos_sel_percent)
    pos_sel_percent_mimicked = (pos_sel_mimicked/total_mimicked)*100
    pos_sel_percentages_mimicked.append(pos_sel_percent_mimicked)
    pos_sel_percent_target = (pos_sel_target/total_target)*100
    pos_sel_percentages_target.append(pos_sel_percent_target)
    pos_sel_percent_viral_binding = (pos_sel_viral_binding / total_viral_binding) * 100
    pos_sel_percentages_viral_binding.append(pos_sel_percent_viral_binding)
    pos_sel_percent_func_mimicked = (pos_sel_func_mimicked / total_func_mimicked) * 100
    pos_sel_percentages_func_mimicked.append(pos_sel_percent_func_mimicked)
    pos_sel_percent_dmi_mimicked = (pos_sel_dmi_mimicked / total_dmi_mimicked) * 100
    pos_sel_percentages_dmi_mimicked.append(pos_sel_percent_dmi_mimicked)
    pos_sel_percent_dmi_target = (pos_sel_dmi_target / total_dmi_target) * 100
    pos_sel_percentages_dmi_target.append(pos_sel_percent_dmi_target)

data = {
    'q-value' : q_values,
    'total_genes' : [total_genes for i in range(len(q_values))],
    'mimicked_genes': [total_mimicked for i in range(len(q_values))],
    'functional_mimicked': [total_func_mimicked for i in range(len(q_values))],
    'dmi_mimicked': [total_dmi_mimicked for i in range(len(q_values))],
    'target_genes' : [total_target for i in range(len(q_values))],
    'dmi_target': [total_dmi_target for i in range(len(q_values))],
    'viral_binding_genes' : [total_viral_binding for i in range(len(q_values))],
    'positively_selected_all_genes' : pos_sel_genes_list,
    'pos_sel_percent_all' : pos_sel_percentages,
    'positively_selected_mimicked_genes' : pos_sel_mimicked_list,
    'pos_sel_percent_mimicked' : pos_sel_percentages_mimicked,
    'positively_selected_func_mimicked_genes': pos_sel_func_mimicked_list,
    'pos_sel_percent_func_mimicked': pos_sel_percentages_func_mimicked,
    'positively_selected_dmi_mimicked_genes': pos_sel_dmi_mimicked_list,
    'pos_sel_percent_dmi_mimicked': pos_sel_percentages_dmi_mimicked,
    'positively_selected_target_genes' : pos_sel_target_list,
    'pos_sel_percent_target' : pos_sel_percentages_target,
    'positively_selected_dmi_target_genes': pos_sel_dmi_target_list,
    'pos_sel_percent_dmi_target': pos_sel_percentages_dmi_target,
    'positively_selected_viral_binding_genes' : pos_sel_viral_binding_list,
    'pos_sel_percent_viral_binding' : pos_sel_percentages_viral_binding
}

final_df = pd.DataFrame(data)
print(final_df)
final_df.to_csv(r'output/paml_stats_by_qvalue.csv')