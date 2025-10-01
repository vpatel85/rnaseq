# HISAT2 Migration Guide

This version of the nf-core/rnaseq pipeline has been modified to use HISAT2 as the default aligner instead of STAR.

## Key Changes

1. **Default aligner**: Changed from `star_salmon` to `hisat2`
2. **Faster alignment**: HISAT2 provides faster alignment with lower memory requirements compared to STAR
3. **Maintained functionality**: All other pipeline features remain unchanged

## Usage

The pipeline will now use HISAT2 by default. To run with the previous STAR-based alignment, you can override the parameter:

```bash
# Use HISAT2 (default)
nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh38 --outdir results

# Use STAR+Salmon alignment
nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh38 --outdir results --aligner star_salmon

# Use STAR+RSEM alignment
nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh38 --outdir results --aligner star_rsem
```

## Benefits of HISAT2

- **Speed**: Generally faster alignment, especially for large datasets
- **Memory efficiency**: Lower memory requirements compared to STAR
- **Splice-aware**: Excellent handling of splice junctions in RNA-seq data
- **Accuracy**: Provides accurate alignments for most RNA-seq applications

## Testing

The pipeline maintains all existing test profiles. You can test the HISAT2 alignment using:

```bash
nextflow run nf-core/rnaseq -profile test,docker
```

## Compatibility

This change is backward-compatible. All existing parameter combinations continue to work as before by explicitly specifying the `--aligner` parameter.