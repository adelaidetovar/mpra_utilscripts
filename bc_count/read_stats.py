import os
import re
import pandas as pd
import argparse

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument('--in_clip', required=True)
    opts.add_argument('--in_umi', default=None)
    opts.add_argument('--in_umibc', default=None)
    opts.add_argument('--in_clust', required=True)
    opts.add_argument('--out_tab', required=True)

    o = opts.parse_args()

    in_clip = o.in_clip
    in_umi = o.in_umi
    in_umibc = o.in_umibc
    in_clust = o.in_clust
    out_tab = o.out_tab

    summary_data = []

    # get logs
    for cutadapt_log in os.listdir(in_clip):
        if cutadapt_log.endswith(".clip.log"):
            libname = os.path.basename(cutadapt_log).replace(".clip.log", "")
            umitools_log = os.path.join(in_umi, f"{libname}.umi.log") if in_umi else None
            deduplicated_barcodes = os.path.join(in_umibc, f"{libname}.starumi") if in_umibc else None
            clustered_barcodes = os.path.join(in_clust, f"{libname}.bc_cluster.txt")

            # input read pairs from cutadapt
            with open(os.path.join(in_clip, cutadapt_log), 'r') as f:
                cutadapt_in = f.read()
                try:
                    input_read_pairs = int(re.search(r"Total read pairs processed:\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))
                except AttributeError:
                    print(f"Error: Cannot find 'Total read pairs processed' in {cutadapt_log}")
                    continue

                # pass-cutadapt read-pairs
                try:
                    pass_cutadapt = int(re.search(r"Pairs written \(passing filters\):\s+([\d,]+)", cutadapt_in).group(1).replace(',', ''))
                    pass_cutadapt_percentage = "{:.2f}%".format((pass_cutadapt / input_read_pairs) * 100)
                except AttributeError:
                    print(f"Error: Cannot find 'Pairs written (passing filters)' in {cutadapt_log}")
                    continue

            # pass-umitools read-pairs
            pass_umitools = 0
            pass_umitools_percentage = "N/A"
            if umitools_log and os.path.exists(umitools_log):
                with open(umitools_log, 'r') as f:
                    try:
                        pass_umitools = int(re.search(r"Reads output:\s+(\d+)", f.read()).group(1))
                        pass_umitools_percentage = "{:.2f}%".format((pass_umitools / input_read_pairs) * 100)
                    except (AttributeError, ValueError):
                        print(f"Error: Unable to parse umitools log: {umitools_log}")
                        continue

            # deduplicated barcodes
            deduplicated_barcode_count = 0
            if deduplicated_barcodes and os.path.exists(deduplicated_barcodes):
                with open(deduplicated_barcodes, 'r') as f:
                    deduplicated_barcode_count = sum(1 for _ in f)

            # clustered barcodes
            clustered_barcode_count = 0
            if os.path.exists(clustered_barcodes):
                with open(clustered_barcodes, 'r') as f:
                    clustered_barcode_count = sum(1 for _ in f) - 1

            # add to summary table
            summary_data.append([libname, input_read_pairs, pass_cutadapt, pass_cutadapt_percentage,
                                 pass_umitools, pass_umitools_percentage, deduplicated_barcode_count, clustered_barcode_count])

    # save to a file
    columns = ["Sample", "Input_Read_Pairs", "Pass_Cutadapt_Read_Pairs", "Pass_Cutadapt_Percentage",
               "Pass_Umitools_Read_Pairs", "Pass_Umitools_Percentage", "Deduplicated_Barcodes", "Clustered_Barcodes"]
    summary_df = pd.DataFrame(summary_data, columns=columns)
    summary_df.to_csv(out_tab, sep='\t', index=False)

    print(f"Summary table generated: {out_tab}")


if __name__ == '__main__':
    main()

