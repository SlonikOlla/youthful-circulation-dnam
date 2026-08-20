Reviewer-driven robustness analyses

1. replicate_aware_DMC_summary.csv
   Biological-replicate Welch sensitivity analysis of the >=10 percentage-point reversal sets.
   P values for aging and HPB were tested separately; BH correction used the full 1,008,133-CpG universe.
   Result: nominal support for 5,269 blood and 3,645 liver sites in both contrasts, but no sites passed q<0.05 in both contrasts.

2. blood_rep... / liver_rep...
   Per-CpG statistics for the effect-size-defined reversal sets.

3. genic_promoter_genebody_circular_null_summary.csv
   Conservative conditional null for promoter/gene-body convergence.
   Canonical chromosomes: observed 2,655 common genes; null mean 2,939; empirical P=1.0.
   Thus the 2,800 all-contig promoter/gene-body overlap should be treated descriptively.

4. promoter_window_reconstruction.csv
   Reverse-check of annotation labels against mm10 RefGene GTF.
   The strand-aware TSS -2000/+500-bp definition reproduced the stored promoter labels with >99.98% overall accuracy.

5. Scripts are included for transparency.
