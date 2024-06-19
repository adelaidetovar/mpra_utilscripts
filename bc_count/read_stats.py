import os
import re
import pandas as pd
import argparse

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument('--in_clip', required = True)
    opts.add_argument('--in_umi', default = None)
    opts.add_argument('--in_umibc', default = None)
    opts.add_argument('--in_clust', required = True)
    opts.add_argument('--out_tab', required = True)

    o = opts.parse_args()

    in_clip = o.in_clip
    in_umi = o.in_umi
    in_umibc = o.in_umibc
    in_clust = o.in_clust
    out_tab = o.out_tab

    # Initialize summary table
    summary_data = []

    # Process cutadapt logs
    for cutadapt_log in os.listdir(in_clip):
        if cutadapt_log.endswith(".clip.log"):
            libname = os.path.basename(cutadapt_log).replace(".clip.log", "")
            umitools_log = os.path.join(in_umi, f"{libname}.umi.log")
            deduplicated_barcodes = os.path.join(in_umibc, f"{libname}.starumi")
            clustered_barcodes = os.path.join(in_clust, f"{libname}.bc_cluster.txt")

            # Get the number of input read-pairs from cutadapt log
            with open(os.path.join(in_clip, cutadapt_log), 'r') as f:
                cutadapt_in = f.read() 
                input_read_pairs = int(re.search(r"Total read pairs processed:\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))

            # Get the number and percentage of pass-cutadapt read-pairs
                pass_cutadapt = int(re.search(r"Pairs written \(passing filters\):\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))
                pass_cutadapt_percentage = "{:.2f}%".format((pass_cutadapt / input_read_pairs) * 100)

            # Get the number and percentage of pass-umitools read-pairs
            if in_umi and os.path.exists(umitools_log):
                with open(umitools_log, 'r') as f:
                    pass_umitools = int(re.search(r"Reads output:\s+(\d+)", f.read()).group(1))
                pass_umitools_percentage = "{:.2f}%".format((pass_umitools / input_read_pairs) * 100)
            else:
                pass_umitools = 0
                pass_umitools_percentage = "N/A"

            if in_umibc and os.path.exists(deduplicated_barcodes):
                # Get the number of deduplicated barcodes
                with open(deduplicated_barcodes, 'r') as f:
                    deduplicated_barcode_count = sum(1 for _ in f)
            else:
                deduplicated_barcode_count = 0

            # Get the number of clustered barcodes
            with open(clustered_barcodes, 'r') as f:
                clustered_barcode_count = sum(1 for _ in f) - 1

            # Append to the summary table
            summary_data.append([libname, input_read_pairs, pass_cutadapt, pass_cutadapt_percentage,
                                 pass_umitools, pass_umitools_percentage, deduplicated_barcode_count, clustered_barcode_count])

    # Convert the summary data to a DataFrame and save to a file
    columns = ["Sample", "Input_Read_Pairs", "Pass_Cutadapt_Read_Pairs", "Pass_Cutadapt_Percentage",
               "Pass_Umitools_Read_Pairs", "Pass_Umitools_Percentage", "Deduplicated_Barcodes", "Clustered_Barcodes"]
    summary_df = pd.DataFrame(summary_data, columns=columns)
    summary_df.to_csv(out_tab, sep='\t', index=False)

    print(f"Summary table generated: {out_tab}")



if __name__ == '__main__':
    main()
