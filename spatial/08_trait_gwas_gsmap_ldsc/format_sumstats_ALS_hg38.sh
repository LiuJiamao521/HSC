# select columns
zcat 34873335-GCST90027163-MONDO_0004976.h.tsv.gz | \
awk -F'\t' -v OFS='\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      col[$i] = i
    }
    print "variant_id", "effect_allele", "other_allele", 
          "beta", "standard_error", "p_value", "N"
    next
  }
  
  {
    print $col["variant_id"], $col["effect_allele"], $col["other_allele"], 
          $col["beta"], $col["standard_error"], $col["p_value"], 80082
  }
' | gzip > ALS_Build38_pre.tsv.gz

# format sumstats
gsmap format_sumstats \
--sumstats ALS_Build38_pre.tsv.gz \
--out ALS \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--beta beta \
--p p_value \
--se standard_error \
--n N

