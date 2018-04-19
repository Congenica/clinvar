import argparse
import collections
import gzip
import os
import re
import pandas as pd
import sys

def gzopen(path, mode='r', verbose=True):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)


def table_to_vcf(input_table_path,genome):
    # validate args
    if not os.path.isfile(input_table_path):
        sys.exit("ERROR: %s not found" % input_table_path)

    # read input table. low_memory allows dtypes to be inferred
    t = pd.read_table(gzopen(input_table_path), low_memory=False)

    missing_columns = {"chrom", "pos", "ref", "alt"} - set(t.columns)
    if missing_columns:
        sys.exit("ERROR: %s is missing columns: %s" % (input_table_path, str(missing_columns)))
    if genome == "b37":
        print("""
##fileformat=VCFv4.1
##source=clinvar
##INFO=<ID=RS,Number=1,Type=String,Description="RSID">
##INFO=<ID=MEASURESET_TYPE,Number=1,Type=String,Description="MEASURESET_TYPE">
##INFO=<ID=MEASURESET_ID,Number=1,Type=String,Description="MEASURESET_ID">
##INFO=<ID=CLNACC,Number=.,Type=String,Description="CLNACC">
##INFO=<ID=ALLELE_ID,Number=1,Type=String,Description="ALLELE_ID">
##INFO=<ID=SYMBOL,Number=1,Type=String,Description="SYMBOL">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="CLNHGVS">
##INFO=<ID=HGVS_P,Number=1,Type=String,Description="HGVS_P">
##INFO=<ID=MOLECULAR_CONSEQUENCE,Number=1,Type=String,Description="MOLECULAR_CONSEQUENCE">
##INFO=<ID=ORIGINAL_CLNSIG,Number=.,Type=String,Description="CLINICAL_SIGNIFICANCE">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="CLINICAL_SIGNIFICANCE_FROM_CLINVAR">
##INFO=<ID=PATHOGENIC,Number=1,Type=String,Description="PATHOGENIC">
##INFO=<ID=BENIGN,Number=1,Type=String,Description="BENIGN">
##INFO=<ID=CONFLICTED,Number=1,Type=String,Description="CONFLICTED">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="CLNREVSTAT">
##INFO=<ID=GOLD_STARS,Number=1,Type=String,Description="Number of gold stars as shown on clinvar web pages to summarize review status. Lookup table described at http://www.ncbi.nlm.nih.gov/clinvar/docs/details/ was used to map the CLNREVSTAT value to this number.">
##INFO=<ID=ALL_SUBMITTERS,Number=.,Type=String,Description="ALL_SUBMITTERS">
##INFO=<ID=CLNDBN,Number=.,Type=String,Description="CLNDBN">
##INFO=<ID=ALL_PMIDS,Number=.,Type=String,Description="ALL_PMIDS">
##INFO=<ID=INHERITANCE_MODES,Number=.,Type=String,Description="INHERITANCE_MODES">
##INFO=<ID=AGE_OF_ONSET,Number=1,Type=String,Description="AGE_OF_ONSET">
##INFO=<ID=PREVALENCE,Number=1,Type=String,Description="PREVALENCE">
##INFO=<ID=DISEASE_MECHANISM,Number=1,Type=String,Description="DISEASE_MECHANISM">
##INFO=<ID=CLNORIGIN,Number=.,Type=String,Description="CLNORIGIN">
##INFO=<ID=XREFS,Number=.,Type=String,Description="CROSS_REFERENCES">
##INFO=<ID=VC,Number=1,Type=String,Description="VARIANT CLASS">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=MT,length=16569>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##reference=Homo_sapiens_assembly19.fasta
        """.strip())

    if genome == "b38":
        print("""
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##source=clinvar
##INFO=<ID=RS,Number=1,Type=String,Description="RSID">
##INFO=<ID=MEASURESET_TYPE,Number=1,Type=String,Description="MEASURESET_TYPE">
##INFO=<ID=MEASURESET_ID,Number=1,Type=String,Description="MEASURESET_ID">
##INFO=<ID=CLNACC,Number=.,Type=String,Description="CLNACC">
##INFO=<ID=ALLELE_ID,Number=1,Type=String,Description="ALLELE_ID">
##INFO=<ID=SYMBOL,Number=1,Type=String,Description="SYMBOL">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="CLNHGVS">
##INFO=<ID=HGVS_P,Number=1,Type=String,Description="HGVS_P">
##INFO=<ID=MOLECULAR_CONSEQUENCE,Number=1,Type=String,Description="MOLECULAR_CONSEQUENCE">
##INFO=<ID=ORIGINAL_CLNSIG,Number=.,Type=String,Description="CLINICAL_SIGNIFICANCE">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="CLINICAL_SIGNIFICANCE_FROM_CLINVAR">
##INFO=<ID=PATHOGENIC,Number=1,Type=String,Description="PATHOGENIC">
##INFO=<ID=BENIGN,Number=1,Type=String,Description="BENIGN">
##INFO=<ID=CONFLICTED,Number=1,Type=String,Description="CONFLICTED">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="CLNREVSTAT">
##INFO=<ID=GOLD_STARS,Number=1,Type=String,Description="Number of gold stars as shown on clinvar web pages to summarize review status. Lookup table described at http://www.ncbi.nlm.nih.gov/clinvar/docs/details/ was used to map the CLNREVSTAT value to this number.">
##INFO=<ID=ALL_SUBMITTERS,Number=.,Type=String,Description="ALL_SUBMITTERS">
##INFO=<ID=CLNDBN,Number=.,Type=String,Description="CLNDBN">
##INFO=<ID=ALL_PMIDS,Number=.,Type=String,Description="ALL_PMIDS">
##INFO=<ID=INHERITANCE_MODES,Number=.,Type=String,Description="INHERITANCE_MODES">
##INFO=<ID=AGE_OF_ONSET,Number=1,Type=String,Description="AGE_OF_ONSET">
##INFO=<ID=PREVALENCE,Number=1,Type=String,Description="PREVALENCE">
##INFO=<ID=DISEASE_MECHANISM,Number=1,Type=String,Description="DISEASE_MECHANISM">
##INFO=<ID=CLNORIGIN,Number=.,Type=String,Description="CLNORIGIN">
##INFO=<ID=XREFS,Number=.,Type=String,Description="CROSS_REFERENCES">
##INFO=<ID=VC,Number=1,Type=String,Description="VARIANT CLASS">
##contig=<ID=1,assembly=b38,length=248956422>
##contig=<ID=2,assembly=b38,length=242193529>
##contig=<ID=3,assembly=b38,length=198295559>
##contig=<ID=4,assembly=b38,length=190214555>
##contig=<ID=5,assembly=b38,length=181538259>
##contig=<ID=6,assembly=b38,length=170805979>
##contig=<ID=7,assembly=b38,length=159345973>
##contig=<ID=8,assembly=b38,length=145138636>
##contig=<ID=9,assembly=b38,length=138394717>
##contig=<ID=10,assembly=b38,length=133797422>
##contig=<ID=11,assembly=b38,length=135086622>
##contig=<ID=12,assembly=b38,length=133275309>
##contig=<ID=13,assembly=b38,length=114364328>
##contig=<ID=14,assembly=b38,length=107043718>
##contig=<ID=15,assembly=b38,length=101991189>
##contig=<ID=16,assembly=b38,length=90338345>
##contig=<ID=17,assembly=b38,length=83257441>
##contig=<ID=18,assembly=b38,length=80373285>
##contig=<ID=19,assembly=b38,length=58617616>
##contig=<ID=20,assembly=b38,length=64444167>
##contig=<ID=21,assembly=b38,length=46709983>
##contig=<ID=22,assembly=b38,length=50818468>
##contig=<ID=X,assembly=b38,length=156040895>
##contig=<ID=Y,assembly=b38,length=57227415>
##contig=<ID=MT,assembly=b38,length=16571>
        """.strip())

    else:
        sys.exit("ERROR: Genome name is %s, not 'b37' or 'b38'" % genome)

    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]))
    for i, table_row in t.iterrows():
        vcf_row = []
        vcf_row.append(table_row["chrom"])
        vcf_row.append(table_row["pos"])
        vcf_row.append(table_row["clnhgvs"])  # ID
        vcf_row.append(table_row["ref"])
        vcf_row.append(table_row["alt"])
        vcf_row.append('.')  # QUAL
        vcf_row.append('.')  # FILTER

        info_field = collections.OrderedDict()

        # from VCF spec:
        #    INFO - additional information: (String, no white-space, semi-colons, or equals-signs permitted; commas are
        #    permitted only as delimiters for lists of values) INFO fields are encoded as a semicolon-separated series of short
        #    keys with optional values in the format: <key>=<data>[,data].

        for key in ['dbsnp', 'measureset_type', 'measureset_id', 'clnacc', 'allele_id', 
                    'symbol', 'clnhgvs', 'hgvs_p', 'molecular_consequence', 'original_clnsig',
                    'clnsig', 'pathogenic', 'benign', 'conflicted', 'clnrevstat', 'gold_stars',
                    'all_submitters', 'clndbn', 'all_pmids', 'inheritance_modes', 'age_of_onset', 
                    'prevalence', 'disease_mechanism', 'clnorigin', 'xrefs', 'type']:

            if pd.isnull(table_row[key]):
                continue
            value = str(table_row[key])
            value = re.sub('\s*[,]\s*', '|', value)  # replace , with ..
            value = re.sub('\s*[;]\s*', '|', value)  # replace ; with |
            value = value.replace("=", " eq ").replace(" ", "_")
            
            info_field[key.upper()] = value
        vcf_row.append(";".join([key+"="+value for key, value in info_field.items()]))

        print("\t".join(map(str, vcf_row)))

    sys.stderr.write("Done\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_table_path', help="Tab-delimited input table")
    args = parser.parse_args()

    table_to_vcf(args.input_table_path,args.genome)
