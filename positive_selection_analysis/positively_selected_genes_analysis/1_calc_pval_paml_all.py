import os
import re
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd


'''
Extracts LRT and q-values from PAML results, for positive selection analysis
'''


def read_codeml_file(filename):
    """Reads in the contents of a codeML output file and extracts the log likelihood value."""
    with open(filename, 'r') as f:
        contents = f.read()
        # Extract the log likelihood value from the output
        #m = re.search(r"log-likelihood:\s+([-0-9\.]+)", contents)
        m = re.search(r"lnL\(.+?\):\s+([-0-9\.]+)", contents)
        if m:
            return float(m.group(1))
        else:
            print(f"No 'log-likelihood' pattern match in file: {filename}")
    return None


def calculate_lrt(no_pos_file, pos_file):
    """Calculates the LRT between the two models."""
    lnL_no_pos = read_codeml_file(no_pos_file)
    lnL_pos = read_codeml_file(pos_file)
    print(lnL_no_pos)
    print(lnL_pos)
    if lnL_no_pos is None or lnL_pos is None:
        return None
    lrt = 2 * (lnL_pos - lnL_no_pos)
    return lrt


directory = r'input/paml_results'
p_vals = []
lrts = []
succeeded_prefix = []
# Files that pass the LRT test
passed_file_path = r'output/paml_results.csv'
with open(passed_file_path, 'w') as passed_file:
    # Keep track of the prefixes that have already been processed
    processed_prefixes = []
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Extract the first 3 letters from the filename
        prefix = filename.split("_")[0]
        #print(filename)
        # Check if the filename ends with '_no_positive.codeml'
        if filename.endswith('_np.codeml'):
            # Check if the prefix has already been processed
            if prefix not in processed_prefixes:
                # Mark the prefix as processed
                processed_prefixes.append(prefix)
                # Check if there is a corresponding positive file with the same prefix
                pos_filename = prefix + '_p.codeml'
                if os.path.exists(os.path.join(directory, pos_filename)):
                    # Calculate the LRT between the two files
                    lrt = calculate_lrt(os.path.join(directory, filename), os.path.join(directory, pos_filename))
                    if lrt is not None:
                        # Calculate the p-value
                        p_value = chi2.sf(lrt, df=1)
                        # Lookup the chi-square table for the critical value of 1 degree of freedom and alpha=0.05
                        lrts.append(lrt)
                        p_vals.append(p_value)
                        succeeded_prefix.append(prefix)
			#crit_value = chi2.ppf(0.95, df=1)
                        #if lrt > crit_value:
    failed_files = list(set(processed_prefixes) - set(succeeded_prefix))
    # with open(r"/tzachi_storage/ofirschor/positive_selection/10_primates_paml_results_failed.txt/", "w") as fd:
    #     fd.write("\n".join(failed_files))
    _, q_vals = fdrcorrection(p_vals)
    df = pd.DataFrame({"genes": succeeded_prefix,"lrt":lrts,"q_values":q_vals})
                            #ens_ids.append(ens_id)
                            # Save prefix and chi-square value to the output file
    df.to_csv(passed_file_path)
