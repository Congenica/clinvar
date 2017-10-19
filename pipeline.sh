
#@XiaoleiZ has tracked down the cause of the error - which is that the “ClinicalSignificance" field
#is set to “-" for many variants (even though they have a clinical significance listed on the clinvar website). The "-" causes the join_data.R error.
#We're planning to use a previous version of the variant_summary file until the data is fixed.
#https://github.com/macarthur-lab/clinvar/pull/33
#JAMES NOTE: The one from March appears fine

source clinvar_pipeline.cfg

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/$variant_summary.txt.gz \
-P $clinvar_refdir

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/$xml \
-P /scratch/data/reference/ClinVar

#Do I want to have this auto-remove or rename old logs?
bsub -q docker -o $output_dir/Clinvar_XML_Parser_full.out \
-e $output_dir/Clinvar_XML_Parser_full.err 'cd clinvar/src/ && python2.7 master.py \
--b37-genome $b37_genome \
--b38-genome $b38_genome \
-X $clinvar_refdir/$xml \
-S $clinvar_refdir/$variant_summary \
--output-path $output_dir'

#Combine the multi and single-allele files into one big honking VCF
bcftools concat -a -d all $output_dir/b37/multi/clinvar_alleles.multi.b37.vcf.gz \
$output_dir/b37/single/clinvar_alleles.single.b37.vcf.gz \
-o $output_dir/clinvar_alleles.combined.b37.vcf.gz

#Remove the duplicate lines, since some alleles will apear in both lists:
cd $output_dir/clinvar_alleles.combined.b37.vcf.gz
grep "^#" clinvar_parsed_combined.b37.vcf > clinvar_parsed_combined.b37.unique.vcf && cat clinvar_parsed_combined.b37.vcf \
| egrep -v '^#' \
| sort \
| uniq \
>> clinvar_parsed_combined.b37.unique.vcf

bsub -P congenica -o $output_dir/ClinVar_Load_Log.out \
-e $output_dir/ClinVar_Load_Log.err \
python ~/sapientia-web/pipeline/ruffus/load_curated_snv_list.py --verbose 3 --jobs 4 \
--curated-variant-list $list_name --curated-variant-list-user $sapientia_user \
--no-qc --no-format-vcf --sample-vcf $output_dir/clinvar_alleles.combined.b37.vcf.gz
