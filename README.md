# g2bpy

A flexible GFF3-to-BED converter with dynamic attribute parsing, multi-criteria filtering, and intelligent column ordering.

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Why g2bpy?

While excellent GFF3-to-BED converters exist, `g2bpy` fills a specific niche for researchers who need:

- **Proper BED format output** with correctly ordered columns and 0-based start coordinates
- **Dynamic attribute parsing** without pre-specifying attribute names
- **Intelligent attribute ordering** with biologically relevant attributes prioritized
- **Complex multi-criteria filtering** with AND/OR logic and exclusion filters
- **Complete information preservation** from all GFF3 columns
- **Robust error handling** with informative diagnostics

## Quick Start

```bash
# Basic usage - extract protein-coding genes
python g2b.py annotation.gff3

# Extract all exons
python g2b.py annotation.gff3 -f column2=exon

# Complex filtering - exclude pseudogenes on specific chromosomes
python g2b.py annotation.gff3 -f column2=gene -f gene_biotype!=pseudogene -f column0=chr1,chr2,chrX

# Custom output file
python g2b.py annotation.gff3 -o filtered_features.bed
```

## Features

### 1. Proper BED Format Output
Converts GFF3 to standard BED format with correctly ordered columns:

| GFF3 Column           | --> | BED Column       | Description                           |
|-----------------------|-----|------------------|---------------------------------------|
| Column 0 (chrom)      | --> | Column 0 (chrom) | Chromosome/scaffold                   |
| Column 3 (start)      | --> | Column 1 (start) | Start position (converted to 0-based) |
| Column 4 (end)        | --> | Column 2 (end)   | End position                          |
| Column 1 (source)     | --> | Column 3 (source)| Feature source                        |
| Column 2 (type)       | --> | Column 4 (type)  | Feature type                          |
| Column 5 (score)      | --> | Column 5 (score) | Score value                           |
| Column 6 (strand)     | --> | Column 6 (strand)| Strand (+/-/.)                        |
| Column 7 (phase)      | --> | Column 7 (phase) | Reading frame phase                   |
| Column 8 (attributes) | --> | Columns 8+       | Parsed attributes                     |

### 2. Intelligent Attribute Ordering
Attributes are automatically ordered with biologically relevant fields first:

1. `gene` - Gene symbol/name
2. `description` - Gene description
3. `gene_name` - Official gene name
4. `gene_synonym` - Alternative names
5. `gene_biotype` - Gene classification
6. `ID` - Feature identifier
7. `db_xref` - Database cross-references
8. `extra_copy_number` - Copy number information
9. `copy_num_ID` - Copy identifier
10. *Additional attributes in alphabetical order*

### 3. Advanced Filtering System

#### Inclusion Filters (`=`)
```bash
# Include only genes
python g2b.py input.gff3 -f column2=gene

# Include multiple feature types (OR logic)
python g2b.py input.gff3 -f column2=gene,transcript,exon

# Filter by attribute values
python g2b.py input.gff3 -f gene_biotype=protein_coding
```

#### Exclusion Filters (`!=`)
```bash
# Exclude pseudogenes
python g2b.py input.gff3 -f column2=gene -f gene_biotype!=pseudogene

# Exclude multiple chromosomes
python g2b.py input.gff3 -f column0!=chrM,chrY

# Exclude specific feature types
python g2b.py input.gff3 -f column2!=exon,CDS
```

#### Complex Filter Combinations
```bash
# Protein-coding genes on autosomes only
python g2b.py input.gff3 \
    -f column2=gene \
    -f gene_biotype=protein_coding \
    -f column0!=chrX,chrY,chrM

# Positive strand features with high confidence
python g2b.py input.gff3 \
    -f column6=+ \
    -f column5!=. \
    -f confidence=high
```

### 4. Dynamic Attribute Discovery
Unlike tools requiring pre-specification of attributes, `g2bpy` automatically:
- Discovers all unique attributes in your GFF3 file
- Handles custom attributes from any annotation source
- Creates consistent column ordering across all output rows

### 5. Robust File Handling
- **Compressed file support**: Native handling of `.gz` files
- **URL decoding**: Automatically decodes encoded characters (e.g., `%2C` → `,`)
- **Error tolerance**: Continues processing despite malformed lines
- **Progress reporting**: Real-time feedback during processing

## 🔄 When to Use Alternative Tools

`g2bpy` prioritizes flexibility and correctness over raw speed. Consider these alternatives for specific use cases:

### For Maximum Speed
- **[gxf2bed](https://github.com/alejandrogzi/gxf2bed)** - Rust implementation, 3-5x faster for large files
- **[BEDOPS gff2bed](https://bedops.readthedocs.io/)** - C++ implementation with streaming efficiency

### For Malformed GFF3 Files
- **[AGAT](https://github.com/NBISweden/AGAT)** - Best-in-class error correction and GFF3 standardization

### For UCSC Genome Browser
- **[UCSC Kent Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)** - Ensures browser compatibility

### For R/Bioconductor Workflows
- **[rtracklayer](https://bioconductor.org/packages/rtracklayer/)** - Seamless R integration

## 📊 Feature Comparison

| Feature                   | g2bpy | gxf2bed | BEDOPS  | AGAT    | UCSC    |
|---------------------------|--------------|---------|---------|---------|---------|
| Proper BED column order   | Yes          | Yes     | Yes     | Partial | Yes     |
| Dynamic attribute parsing | Yes          | No      | Partial | Yes     | No      |
| Attribute prioritization  | Yes          | No      | No      | No      | No      |
| Multi-criteria filtering  | Yes          | No      | No      | Partial | No      |
| Exclusion filters         | Yes          | No      | No      | No      | No      |
| Retains all GFF3 columns  | Yes          | No      | Yes     | Yes     | Partial |
| Speed (relative)          | **           | ***     | ***     | *       | **      |
| Error tolerance           | ***          | **      | *       | ***     | *       |
| Compressed files          | Yes          | Yes     | Yes     | Yes     | Partial |
| Python native             | Yes          | No      | No      | No      | No      |

## Installation

### Requirements
- Python 3.6 or higher
- No external dependencies (uses Python standard library only)

### Download
```bash
# Clone the repository
git clone https://github.com/Thomas-X-Garcia/g2bpy.git
cd g2bpy

# Or download directly
wget https://raw.githubusercontent.com/Thomas-X-Garcia/g2bpy/main/g2b.py
chmod +x g2b.py
```

## Detailed Usage

### Command Line Options

```
usage: g2b.py [-h] [-o OUTPUT] [-f FILTERS] input

Convert GFF3 files to BED-like format with dynamic filtering and attribute parsing.

positional arguments:
  input                 Input GFF3 file (can be gzipped)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output BED file (default: input_filtered.bed)
  -f FILTERS, --filter FILTERS
                        Filter criteria (can specify multiple times)
```

### Input Formats
- Standard GFF3 files (`.gff3`, `.gff`)
- Compressed GFF3 files (`.gff3.gz`, `.gff.gz`)

### Output Format
Tab-delimited BED format with header:
```
#chrom  start   end     source  type    score   strand  phase   [prioritized_attributes...]
chr1    6046    13941   NCBI    gene    .       -       .       LOC124900618    ...
chr1    6049    13939   NCBI    mRNA    .       -       .       rna-XM_054983325.1  ...
```

### Filter Syntax

#### Column Filters
Use `column0` through `column8` to filter by GFF3 columns:
- `column0`: chromosome/scaffold
- `column1`: source
- `column2`: feature type
- `column3`: start position
- `column4`: end position
- `column5`: score
- `column6`: strand
- `column7`: phase
- `column8`: attributes (raw string)

#### Attribute Filters
Use any attribute name found in the GFF3 file:
- `gene_biotype=protein_coding`
- `ID=gene-ENSG00000223972`
- `coverage=1.0`

#### Filter Operators
- `=` : Include rows where field equals value(s)
- `!=` : Exclude rows where field equals value(s)

#### Multiple Values (OR Logic)
Use comma-separated values for OR logic within a single filter:
- `column0=chr1,chr2,chrX` : Include chr1 OR chr2 OR chrX
- `gene_biotype!=pseudogene,lncRNA` : Exclude pseudogenes AND lncRNAs

### Examples

#### Basic Gene Extraction
```bash
# Extract all genes (default includes only protein-coding)
python g2b.py annotation.gff3 -f column2=gene

# Extract all protein-coding genes explicitly
python g2b.py annotation.gff3 -f column2=gene -f gene_biotype=protein_coding
```

#### Feature Type Selection
```bash
# Extract multiple feature types
python g2b.py annotation.gff3 -f column2=gene,transcript,exon

# Extract CDS features on positive strand
python g2b.py annotation.gff3 -f column2=CDS -f column6=+
```

#### Chromosome Filtering
```bash
# Include specific chromosomes
python g2b.py annotation.gff3 -f column0=chr1,chr2,chr3

# Exclude mitochondrial and sex chromosomes
python g2b.py annotation.gff3 -f column0!=chrM,chrX,chrY
```

#### Complex Queries
```bash
# High-confidence protein-coding genes on autosomes
python g2b.py annotation.gff3 \
    -f column2=gene \
    -f gene_biotype=protein_coding \
    -f column0!=chrX,chrY,chrM \
    -f column5!=. \
    -o high_confidence_genes.bed

# Primary gene copies from RefSeq
python g2b.py refseq_annotation.gff3.gz \
    -f column2=gene \
    -f extra_copy_number=0 \
    -f column1=RefSeq \
    -o primary_refseq_genes.bed

# Multi-exon transcripts
python g2b.py annotation.gff3 \
    -f column2=transcript \
    -f exon_count!=1 \
    -o multi_exon_transcripts.bed
```

## Performance Considerations

`g2bpy` uses a two-pass approach:
1. **First pass**: Discovers all unique attributes in filtered rows
2. **Second pass**: Converts and writes output with consistent column ordering

This ensures proper column alignment but requires reading the file twice.

**Performance tips:**
- Use specific filters to reduce rows processed in both passes
- For simple conversions of very large files (>5GB), consider faster alternatives
- The overhead is negligible for files under 1GB
- Filtering can actually improve performance by reducing output size

**Benchmarks** (on 1GB GFF3 file):
- Full conversion: ~45 seconds
- With gene filter: ~12 seconds
- With complex filters: ~15-20 seconds

## Use Cases

### Ideal for:
- Converting annotations to proper BED format for genome browsers
- Research requiring complex feature selection criteria
- Working with custom annotations containing non-standard attributes
- Preparing input for BEDTools or other BED-based tools
- Creating consistent, analysis-ready datasets from diverse GFF3 sources

### Example Research Applications:

#### 1. Comparative Genomics
```bash
# Extract orthologous genes with specific properties
python g2b.py orthologs.gff3 \
    -f column2=gene \
    -f orthology_confidence=high \
    -f species_count=5
```

#### 2. Regulatory Analysis
```bash
# Extract promoter-associated features
python g2b.py regulatory.gff3 \
    -f feature_type=promoter,enhancer \
    -f validation_status=confirmed
```

#### 3. Structural Variation
```bash
# Extract high-confidence structural variants
python g2b.py structural_variants.gff3 \
    -f column2=insertion,deletion,duplication \
    -f size!=small \
    -f support_reads=10
```

#### 4. Alternative Splicing
```bash
# Extract alternatively spliced transcripts
python g2b.py annotation.gff3 \
    -f column2=transcript \
    -f transcript_biotype=protein_coding \
    -f alternative_splicing=yes
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Guidelines
1. Maintain backward compatibility
2. Add tests for new features
3. Update documentation
4. Follow PEP 8 style guidelines

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Inspired by the bioinformatics community's need for flexible annotation processing
- Built upon GFF3 specifications from the Sequence Ontology Consortium
- Thanks to the developers of existing GFF-to-BED tools for setting high standards

## Citation

If you use `g2bpy` in your research, please cite:
```
g2bpy: A flexible GFF3-to-BED converter with intelligent attribute handling
https://github.com/Thomas-X-Garcia/g2bpy
```

## Bug Reports

Please report bugs through the GitHub issue tracker with:
1. Your command line
2. Sample of input file (first 100 lines)
3. Error message or unexpected output
4. Python version

---

**Note**: While this tool was designed for clinical and research use, always validate output!
