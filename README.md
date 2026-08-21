# Youthful circulation and age-associated DNA methylation

Analysis repository supporting the manuscript:

**Youthful circulation reshapes age-associated DNA methylation across blood and liver**

This project analyzes processed reduced-representation bisulfite sequencing (RRBS) data from the prolonged heterochronic parabiosis study deposited as GEO accession **GSE224442**. The analysis asks how age-opposing methylation changes are organized across blood and liver, whether exact CpG responses persist after detachment, and whether tissue-specific responses converge at higher levels of genomic organization.

## Repository contents

- `data/DNAm_reanalysis_supporting_data.xlsx` — principal supporting tables.
- `results/` — selected overlap, null-model, annotation, persistence, and GO outputs.
- `reviewer_robustness/` — replicate-aware sensitivity analysis, promoter-window reconstruction, promoter/gene-body null analysis, and scripts used for those checks.

The raw/source RRBS data are not redistributed here. They should be obtained from GEO accession **GSE224442** and cited to the original study.

## Principal analysis choices

- Common analytical universe: 1,008,133 CpGs with >=10x coverage in every analyzed sample.
- Primary directional-reversal threshold: >=10 percentage-point aging effect and >=10 percentage-point HPB effect in the opposite direction.
- Cross-tissue and persistence tests: chromosome-aware circular-shift null models.
- Gene convergence: nearest-gene analysis with genomic null; promoter/gene-body analysis reported as a stricter sensitivity analysis.
- GO analysis: direct mouse GO annotations excluding IEA, with genomic-structure-aware null testing.

## Reproducibility

The derived tables in this repository correspond to the manuscript analyses. The original processed coverage files are available through the source dataset and are required to reconstruct the full CpG matrix from scratch.

## Citation

Archived on Zenodo: **DOI 10.5281/zenodo.22036635**

https://doi.org/10.5281/zenodo.22036635

## Authors

Igor Kovalchuk and Olga Kovalchuk.
