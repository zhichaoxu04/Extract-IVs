# Extract-IVs: Extract IVs for the summary-data-based Mendelian Randomization (MR)

## Summary-data-based Mendelian Randomization ([SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview))
[SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview) (Summary-data-based Mendelian Randomization) is a powerful method that allows the utilization of publicly available data for conducting Mendelian Randomization (MR) analysis. Instead of relying on individual-level data, SMR makes use of summary-level data from Genome-Wide Association Studies (GWAS) to identify potential causal relationships between a genetic trait (e.g., a genetic variant) and an outcome or phenotype.

- (eQTLs)[https://en.wikipedia.org/wiki/Expression_quantitative_trait_loci] are genomic loci that are associated with variation in gene expression levels across individuals. They provide insight into the genetic regulation of gene expression and are instrumental in interpreting GWAS findings, especially when the associated variants lie in non-coding regions.
- (mQTLs) are genomic loci that influence DNA methylation levels. DNA methylation is an epigenetic modification that can regulate gene expression.

### Query Expression Quantitative Trait Loci (eQTL) Summary Results
```bash
S:
cd S:\Your\path\to\SMR
smr --beqtl-summary Whole_Blood.lite --query 5.0e-8 --genes genelist.txt --out myquery_lite
```

- `--query` saves in text format a subset of the eQTL summary dataset based on the specified eQTL p-value threshold. The default value is 5.0e-8.
- `--genes` extracts a subset of probes which tag the genes in the list.
- `--beqtl-summary` reads summary-level data from a eQTL study in binary format. We store eQTL summary data in three separate files .esi (SNP information, in the same format as the PLINK .bim file), .epi (probe information) and .besd (eQTL summary statistics in binary format). See [Data Management](https://yanglab.westlake.edu.cn/software/smr/#DataManagement) for more information. 
