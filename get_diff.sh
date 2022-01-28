cat li|xargs -I {} perl diff_my_picard.pl picard_out/{}.alignment_summary_metrics my_out/{}.mt
