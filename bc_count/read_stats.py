import os
import re
import pandas as pd
import argparse

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument('--in_clip', required=True)
    opts.add_argument('--use_umi', default=True)
    opts.add_argument('--in_umi', default=None)
    opts.add_argument('--in_umibc', default=None)
    opts.add_argument('--in_clust', required=True)
    opts.add_argument('--out_tab', required=True)

    o = opts.parse_args()

    in_clip = o.in_clip
    use_umi = o.use_umi
    in_umi = o.in_umi
    in_umibc = o.in_umibc
    in_clust = o.in_clust
    out_tab = o.out_tab

    summary_data = []

    if use_umi is True:
        # get logs
        for umitools_log in os.listdir(in_umi):
            if umitools_log.endswith(".umi.log"):
                libname = os.path.basename(umitools_log).replace(".umi.log", "")
                cutadapt_log = os.path.join(in_clip, f"{libname}.clip.log")
                deduplicated_barcodes = os.path.join(in_umibc, f"{libname}.starumi") if in_umibc else None
                clustered_barcodes = os.path.join(in_clust, f"{libname}.bc_cluster.txt")

                # input read pairs, adapter read pairs from umitools
                pass_umitools = 0
                pass_umitools_percentage = "N/A"
                with open(os.path.join(in_umi, umitools_log), 'r') as f:
                    umitools_in = f.read()
                    try:
                        input_reads = int(re.search(r"Input Reads:\s+([\d]+)", umitools_in).group(1))
                        pass_umitools = int(re.search(r"Reads output:\s+(\d+)", umitools_in).group(1))
                        pass_umitools_percentage = "{:.2f}%".format((pass_umitools / input_reads) * 100)
                    except AttributeError:
                        print(f"Error: unable to parse {umitools_log}!")
                        continue

                # reads with barcode from cutadapt
                with open(os.path.join(in_clip, cutadapt_log), 'r') as f:
                    cutadapt_in = f.read()
                    try:
                        pass_cutadapt = int(re.search(r"Reads written \(passing filters\):\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))
                        pass_cutadapt_percentage = "{:.2f}%".format((pass_cutadapt / input_reads) * 100)
                    except AttributeError:
                        print(f"Error: unable to parse {cutadapt_log}!")
                        continue

                    # deduplicated barcodes
                deduplicated_barcode_count = 0
                if deduplicated_barcodes and os.path.exists(deduplicated_barcodes):
                    with open(deduplicated_barcodes, 'r') as f:
                        deduplicated_barcode_count = sum(1 for _ in f) - 1

                    # clustered barcodes
                clustered_barcode_count = 0
                if os.path.exists(clustered_barcodes):
                    with open(clustered_barcodes, 'r') as f:
                        clustered_barcode_count = sum(1 for _ in f) - 1

                    # add to summary table
                summary_data.append([libname, input_reads, pass_umitools, pass_umitools_percentage, pass_cutadapt, pass_cutadapt_percentage,
                                    deduplicated_barcode_count, clustered_barcode_count])

        # save to a file
        columns = ["Sample", "Input_Reads", "Pass_Umitools_Reads", "Pass_Umitools_Percentage", 
                    "Pass_Cutadapt_Reads", "Pass_Cutadapt_Percentage", "Deduplicated_Amplicons", "Clustered_Barcodes"]
        summary_df = pd.DataFrame(summary_data, columns=columns)
        summary_df.sort_values(by = "Sample", inplace=True)
        comma_cols = summary_df.columns.difference(["Sample", "Pass_Umitools_Percentage", "Pass_Cutadapt_Percentage"])
        summary_df[comma_cols] = summary_df[comma_cols].apply(lambda col: col.map(lambda x: f"{x:,}" if pd.notnull(x) else x))
        summary_df.to_csv(out_tab, sep='\t', index=False)

        print(f"Read stats table generated!")

    else:
        for cutadapt_log in os.listdir(in_clip):
            if cutadapt_log.endswith(".clip.log"):
                libname = os.path.basename(cutadapt_log).replace(".clip.log", "")
                clustered_barcodes = os.path.join(in_clust, f"{libname}.bc_cluster.txt")

                # input read pairs from cutadapt
                with open(os.path.join(in_clip, cutadapt_log), 'r') as f:
                    cutadapt_in = f.read()
                    try:
                        input_reads = int(re.search(r"Total reads processed:\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))
                    except AttributeError:
                        print(f"Error: Cannot find 'Total reads processed' in {cutadapt_log}!")
                        continue

                    # pass-cutadapt read-pairs
                    try:
                        pass_cutadapt = int(re.search(r"Reads written \(passing filters\):\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))
                        pass_cutadapt_percentage = "{:.2f}%".format((pass_cutadapt / input_reads) * 100)
                    except AttributeError:
                        print(f"Error: Cannot find 'Reads written (passing filters)' in {cutadapt_log}!")
                        continue

                clustered_barcode_count = 0
                if os.path.exists(clustered_barcodes):
                    with open(clustered_barcodes, 'r') as f:
                        clustered_barcode_count = sum(1 for _ in f) - 1

                summary_data.append([libname, input_reads, pass_cutadapt, pass_cutadapt_percentage, clustered_barcode_count])

        columns = ["Sample", "Input_Reads", "Pass_Cutadapt_Reads", "Pass_Cutadapt_Percentage", "Clustered_Barcodes"]
        summary_df = pd.DataFrame(summary_data, columns=columns)
        summary_df.sort_values(by="Sample", inplace=True)
        comma_cols = summary_df.columns.difference(["Sample", "Pass_Cutadapt_Percentage"])
        summary_df[comma_cols] = summary_df[comma_cols].apply(lambda col: col.map(lambda x: f"{x:,}" if pd.notnull(x) else x))
        summary_df.to_csv(out_tab, sep='\t', index=False)

        print(f"Read stats table generated!")

if __name__ == '__main__':
    main()

