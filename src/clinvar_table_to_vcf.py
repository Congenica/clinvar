import argparse
import collections
import gzip
import os
import re
import pandas as pd
import sys

from join_variant_summary_with_clinvar_alleles import FINAL_HEADER


def gzopen(path, mode='r', verbose=True):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)


def table_to_vcf(input_table_path, input_reference_genome):
    # validate args
    input_reference_genome_fai = input_reference_genome + ".fai"
    if not os.path.isfile(input_table_path):
        sys.exit("ERROR: %s not found" % input_table_path)
    if not os.path.isfile(input_reference_genome_fai):
        sys.exit("ERROR: %s (reference FASTA .fai) not found" %
                 input_reference_genome_fai)

    # read input table. low_memory allows dtypes to be inferred
    t = pd.read_table(gzopen(input_table_path), low_memory=False)

    missing_columns = {"chrom", "pos", "ref", "alt"} - set(t.columns)
    if missing_columns:
        sys.exit("ERROR: %s is missing columns: %s" % (input_table_path, str(missing_columns)))

    print("""##fileformat=VCFv4.1\n##source=clinvar""")

    descriptions = {
        'gold_stars': "Number of gold stars as shown on clinvar web pages to summarize review status. Lookup table described at http://www.ncbi.nlm.nih.gov/clinvar/docs/details/ was used to map the REVIEW_STATUS value to this number.",
    }
    not_required_in_header_or_info = ['chrom', 'pos', 'ref', 'alt', 'start', 'stop', 'strand',
                                      'clinical_significance_ordered', 'review_status_ordered',
                                      'dates_ordered', 'last_evaluated', 'submitters_ordered', 'scv', 'type',
                                      'age_of_onset', 'prevalence', 'disease_mechanism']
    single_value_fields = ['symbol', 'pathogenic', 'likely_pathogenic', 'uncertain_significance', 'likely_benign',
                           'benign', 'clnsig', 'gold_stars', 'conflicted']
    for key in FINAL_HEADER:
        if key not in not_required_in_header_or_info:
            number = '1' if key in single_value_fields else '.'
            print("""##INFO=<ID={},Number={},Type=String,Description="{}">""".format(
                key.upper(), number, descriptions.get(key, key.upper())))

    print('##INFO=<ID=TOTAL_VARIANTS,Number=1,Type=String,Description="{}"'.format(t.shape[0]))
    for pathogenicity in ('Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign'):
        pathogenicity_count = (t['clnsig'] == pathogenicity).sum()
        print('##INFO=<ID={}_VARIANTS,Number=1,Type=String,Description="{}"'.format(
            pathogenicity.upper().replace(' ', '_'), pathogenicity_count))

    with open(input_reference_genome_fai) as in_fai:
        for line in in_fai:
            chrom, length, _ = line.split("\t", 2)
            print("""##contig=<ID={},length={}>""".format(
                chrom.replace("chr", ""), length))
    print("""##reference={}""".format(input_reference_genome))

    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))
    for i, table_row in t.iterrows():
        vcf_row = []
        vcf_row.append(table_row["chrom"])
        vcf_row.append(table_row["pos"])
        vcf_row.append(table_row["rs"])  # ID = rsID
        vcf_row.append(table_row["ref"])
        vcf_row.append(table_row["alt"])
        vcf_row.append('.')  # QUAL
        vcf_row.append('.')  # FILTER

        info_field = collections.OrderedDict()

        # from VCF spec:
        #    INFO - additional information: (String, no white-space, semi-colons, or equals-signs permitted; commas are
        #    permitted only as delimiters for lists of values) INFO fields are encoded as a semicolon-separated series of short
        #    keys with optional values in the format: <key>=<data>[,data].
        for key in FINAL_HEADER:
            if key not in not_required_in_header_or_info:
                if pd.isnull(table_row[key]):
                    continue
                value = str(table_row[key])
                value = re.sub('\s*[,]\s*', '..', value)  # replace , with ..
                value = re.sub('\s*[;]\s*', ',', value)  # replace ; with |
                value = value.replace("=", " eq ").replace(" ", "_")
                
                info_field[key.upper()] = value
        vcf_row.append(";".join([key+"="+value for key, value in info_field.items()]))

        print("\t".join(map(str, vcf_row)))

    sys.stderr.write("Done\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_table_path', help="Tab-delimited input table")
    parser.add_argument('input_reference_genome', help="Reference FASTA used. The associated .fai file, e.g. b38.fa.fai, is necessary for the VCF header generation")
    args = parser.parse_args()

    table_to_vcf(args.input_table_path, args.input_reference_genome)
