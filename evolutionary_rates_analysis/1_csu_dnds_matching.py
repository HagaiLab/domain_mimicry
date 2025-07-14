import os
import pandas as pd
from uniprot_to_ensembel import get_ensg_id
from Bio import SeqIO

'''
 This script creates a file matching each residue to its spatial classification (interface, surface or core
  as identified by CSU), dN\dS value (evolutionary rates from Selecton), and matching residue by ensembl sequence.
'''

def read_fasta_biopython(file_path):
    sequences = []
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append({"id": record.id, "sequence": str(record.seq)})
    return sequences


def find_fasta(id, fasta_dictionary):
    for dictionary in fasta_dictionary:
        if id in dictionary['id']:
            return dictionary['sequence']


# for each file in the results folder, read txt file and save as a string
results_path = r"csu_dnds_input/BLAST_results"
results_files = os.listdir(results_path)
proteome_file = r'csu_dnds_input/human_proteome.fasta'
homology_dataframe_path = r"csu_dnds_input/homology_tables"
homology_results_files = os.listdir(homology_dataframe_path)
csu_results_path = r'csu_dnds_input/csu_lists'
csu_results_files = os.listdir(csu_results_path)
dnds_results_path = r"csu_dnds_input/dnds_lists"
dnds_results_files = os.listdir(dnds_results_path)


for homology_file_name in homology_results_files:
    homology_result_path = os.path.join(homology_dataframe_path, homology_file_name)
    mimicked_uniprot_id = homology_file_name.split("_")[1].split(".")[0]
    mimicking_uniprot_id = homology_file_name.split("_")[2].split(".")[0]
    target_uniprot_id = homology_file_name.split("_")[3].split(".")[0]
    homology_dataframe = pd.read_csv(homology_result_path, skiprows=1)


    #locate blast result of mimicked protein
    for filename in results_files:
        if mimicked_uniprot_id in filename:
            result_path = os.path.join(results_path, filename)
            uniprot_id = filename.split("_")[-1].split(".")[0]
            with open(result_path, 'r') as file:
                blast_alignment = file.read().replace('\n', '')
                blast_alignment_list = list(blast_alignment.split(","))
            qstart = blast_alignment_list[5]
            qend = blast_alignment_list[6]
            sstart = int(blast_alignment_list[7])
            send = int(blast_alignment_list[8])
            qseq = blast_alignment_list[9]
            qseq_list = [res for res in qseq]
            sseq = blast_alignment_list[10]
            sseq_list = [res for res in sseq]

            fasta_dictionary = read_fasta_biopython(proteome_file)
            uniprot_fasta = find_fasta(uniprot_id, fasta_dictionary)

            # creating empty lists for the uniprot and ensembel sequences after alignment
            # inset gaps before alignment using start values
            gaps_before_alignment_uniprot = int(qstart) - 1
            gaps_before_alignment_wnsembel = int(sstart) - 1
            uniprot_seq = []  # leave the entire uniprot seq, no gaps
            for i in range(gaps_before_alignment_uniprot):
                uniprot_seq.append(uniprot_fasta[i])
            ensembel_seq = ['None' for i in range(gaps_before_alignment_uniprot)]  # insert "None" before alignment part
            ensembel_index = ['None' for i in range(gaps_before_alignment_uniprot)]
            # add the aligned part
            for residue in qseq_list:
                uniprot_seq.append(residue)
            for i in range(int(qend), len(uniprot_fasta)):
                uniprot_seq.append(uniprot_fasta[i])

            for residue in sseq_list:
                ensembel_seq.append(residue)
            while len(ensembel_seq) < len(uniprot_seq):
                ensembel_seq.append('None')
            uniprot_index = []
            no_gaps_uniprot_index = 1
            for i in range(len(uniprot_seq)):
                if uniprot_seq[i] != '-':
                    uniprot_index.append(no_gaps_uniprot_index)
                    no_gaps_uniprot_index += 1
                else:
                    uniprot_index.append('')

            no_gaps_uniprot_ensembel = int(sstart)
            for i in range(len(sseq)):
                if sseq[i] != '-':
                    ensembel_index.append(no_gaps_uniprot_ensembel)
                    no_gaps_uniprot_ensembel += 1
                else:
                    ensembel_index.append('')

            while len(ensembel_index) < len(uniprot_seq):
                ensembel_index.append("None")

            uniprot_id_list = [uniprot_id for i in range(len(uniprot_seq))]
            target_uniprot_id_list = [target_uniprot_id for i in range(len(uniprot_seq))]

            # locate the correct csu list file
            for csu_filename in csu_results_files:
                if mimicked_uniprot_id in csu_filename and target_uniprot_id in csu_filename:
                    csu_result_path = os.path.join(csu_results_path, csu_filename)
                    with open(csu_result_path, 'r') as file:
                        # Read the content of the file line by line
                        lines = file.readlines()
                    # Strip newline characters from each line and store them in a list
                    csu_list = [line.strip() for line in lines]
                    #insert csu only when there os no gap
                    final_csu_list = []
                    i = 0
                    for char in uniprot_seq:
                        if char == '-':
                            final_csu_list.append('-')
                        else:
                            final_csu_list.append(csu_list[i])
                            i += 1
                    # locate dnds file
                    mimicked_ensembl_id = get_ensg_id(mimicked_uniprot_id)
                    for dnds_filename in dnds_results_files:
                        if mimicked_uniprot_id in dnds_filename:
                            dnds_result_path = os.path.join(dnds_results_path, dnds_filename)
                            with open(dnds_result_path, 'r') as file:
                                # Read the content of the file line by line
                                lines = file.readlines()
                            # Strip newline characters from each line and store them in a list
                            dnds_list = [line.strip() for line in lines]

                            dnds_list = dnds_list[sstart - 1:send]  # keeping the aligned part only
                            dnds_list_with_gaps = []
                            k = 0


                            for char in ensembel_seq:
                                if char == 'None':
                                    dnds_list_with_gaps.append('')
                                elif char == '-':
                                    dnds_list_with_gaps.append('')
                                else:
                                    dnds_list_with_gaps.append(dnds_list[k])
                                    k += 1

                            data = {
                                'mimicked_id': uniprot_id_list,
                                'uniprot_seq': uniprot_seq,
                                'uniprot_index': uniprot_index,
                                'ensembel_seq': ensembel_seq,
                                'ensembel_index': ensembel_index,
                                'dN/dS': dnds_list_with_gaps,
                                'csu_result': final_csu_list,
                                'target_id': target_uniprot_id_list
                            }

                            mimicked_target_df = pd.DataFrame(data)

                            final_table_path = f"csu_dnds_output/matching_results_{mimicked_uniprot_id}_{target_uniprot_id}.csv"
                            mimicked_target_df.to_csv(final_table_path, index=False)
