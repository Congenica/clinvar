#!/usr/bin/env Rscript
#Script to join together the data fromt heClinVar tsv and from the parsed ClinVar XML

options(stringsAsFactors=F)
#options(warn=2) 
options(error = quote({
  dump.frames(to.file=T, dumpto='last.dump')
  load('last.dump.rda')
  print(last.dump)
  q()
}))

args = commandArgs(trailingOnly=TRUE)

convert_pathogenicity<-function(pathogenicity){
  #convert the list of pathogenicity terms into a single pathogenicity term that will fit in Sapientia.
  #From the list of terms, take the one with the greatest indication of pathogenicity.
  pathogenicity_delim=","
  pathogenicity_list<-unlist(strsplit(pathogenicity,","))
  
  pathogenic_terms<-c("Pathogenic/Likely pathogenic", "Pathogenic", 'pathogenic')
  likely_pathogenic_terms<-c('Likely pathogenic','likely pathogenic')
  vus_terms<-c('Uncertain significance','uncertain significance')
  likely_benign_terms<-c('Benign/Likely benign','Likely benign','likely benign')
  benign_terms<-c('Benign',"benign")
  unused_terms<-c('association not found','drug response','confers sensitivity','risk factor','other','association',
                  'protective','not provided','conflicting data from submitters', 'Conflicting interpretations of pathogenicity')
  
  matches<-pathogenicity_list %in% pathogenic_terms 
  if(TRUE %in% matches){
    return("Clearly pathogenic")
  }
  
  matches<-pathogenicity_list %in% likely_pathogenic_terms 
  if(TRUE %in% matches){
    return("Likely to be pathogenic")
  }
  
  matches<-pathogenicity_list %in% vus_terms 
  if(TRUE %in% matches){
    return("Unknown significance (VUS)")
  }
  
  matches<-pathogenicity_list %in% likely_benign_terms 
  if(TRUE %in% matches){
    return("Unlikely to be pathogenic")
  }
  
  matches<-pathogenicity_list %in% benign_terms 
  if(TRUE %in% matches){
    return("Clearly not pathogenic")
  }
  
  matches<-pathogenicity_list %in% unused_terms 
  if(TRUE %in% matches){
    return("")
  }
  
  else {
    warning("Unable to map \'",pathogenicity_list, "\'. Not a valid set of clinical significance terms")
    return("")
  }
}

variant_summary_table = 'variant_summary.txt.gz'
if (length(args) != 3) {
  print(paste("Requires 3 command line args: [variant_summary table] [clinvar_allele_trait_pairs table] [output_table_path]",multi))
  exit(-1)
}

#Read in the data files
variant_summary_table = gzfile(args[1])
clinvar_allele_trait_pairs_table = gzfile(args[2])
output_table = gzfile(args[3], 'w')


# load what we've extracted from the XML so far
xml_raw = read.table(clinvar_allele_trait_pairs_table, sep='\t', comment.char='', quote='', header=T, skipNul=T, check.names=F, stringsAsFactors=F)
print(dim(xml_raw))

# load the tab-delimited summary
txt_download = read.table(variant_summary_table, sep='\t', comment.char='', quote='', header=T, skipNul=T, check.names=F)
print(dim(txt_download))

# subset the tab-delimited summary to desired rows and cols
colnames(txt_download) = gsub('\\.','_',tolower(colnames(txt_download)))
colnames(txt_download) = replace(colnames(txt_download), 1, "allele_id")

desired_columns<-c('allele_id','clinicalsignificance','reviewstatus','lastevaluated')
txt_extract = subset(txt_download, assembly == 'GRCh37', select=desired_columns)
colnames(txt_extract)<-c('allele_id','clinical_significance','review_status','last_evaluated')
#drop the clinical_significance and review_status in clinvar_record.tsv 
#use the summary ones in variant_summary.txt
xml_extract = subset(xml_raw,select=-c(clinical_significance,review_status))

# join on allele id
combined = merge(xml_extract, txt_extract,by='allele_id',all.x=FALSE)

#Turn the list of pathogenicity terms into a string for the single most severe term, in a suitable form for Sapientia's
sapientia_clinsig<-sapply(as.character(combined$clinical_significance),convert_pathogenicity)
sapientia_clinsig<-sapply(sapientia_clinsig,as.character)

if(length(sapientia_clinsig)>0){
  combined<-cbind(combined,sapientia_clinsig)
}else{
  warning("Unable to bind clinsig terms to ClinVar data, no valid clinical significance terms in the TSV")
}
  
# lookup table based on http://www.ncbi.nlm.nih.gov/clinvar/docs/details/
gold_stars_table = list(
  'no assertion provided' = 0,
  'no assertion for the individual variant' = 0,
  'no assertion criteria provided' = 0,
  'criteria provided, single submitter' = 1,
  'criteria provided, conflicting interpretations' = 1, 
  'criteria provided, multiple submitters, no conflicts' = 2, 
  'reviewed by expert panel' = 3,
  'practice guideline' = 4
)

# Add some layers of interpretation on top of this
# Note: we are trying to get the "overall" interpretation that is displayed in the upper right of the clinvar web pages but
# it is not in any of the available FTP downloads, so this is a stopgap
combined$review_status<-sapply(combined$review_status,as.character)
combined$gold_stars = sapply(combined$review_status, function(k) { gold_stars_table[[k]] })

# pathogenic = 1 if at least one submission says path or likely path, 0 otherwise
combined$pathogenic = as.integer(grepl('athogenic',combined$clinical_significance))

# conflicted = 1 if at least one submission each of [likely] benign and [likely] pathogenic
combined$conflicted = as.integer(grepl('athogenic',combined$clinical_significance) & grepl('enign',combined$clinical_significance))

# benign = 1 if at least one submission says benign or likely benign, 0 otherwise
combined$benign = as.integer(grepl('enign',combined$clinical_significance))

# re-order the columns
combined = combined[,c('chrom','pos','ref','alt','dbsnp','measureset_type','measureset_id','rcv','allele_id','symbol', 'hgvs_c','hgvs_p','molecular_consequence','clinical_significance', 'sapientia_clinsig', 'pathogenic', 'benign', 'conflicted', 'review_status','last_evaluated', 'gold_stars','all_submitters','all_traits','all_pmids', 'inheritance_modes', 'age_of_onset','prevalence', 'disease_mechanism', 'origin', 'xrefs')]

write.table(combined, output_table, sep='\t', row.names=F, col.names=T, quote=F)

close(output_table)


