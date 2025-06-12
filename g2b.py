#!/usr/bin/env python3
"""
g2bpy: Flexible GFF3 to BED Converter

This script converts GFF3 files to BED-like format with customizable filtering options
and dynamic column parsing from the attributes field. The output preserves all GFF3
information while using BED-style 0-based start coordinates.

Author: Thomas X. Garcia, PhD, HCLD
License: MIT
"""

import argparse
import gzip
import sys
import os
from collections import OrderedDict
from urllib.parse import unquote
import re


def parse_attributes(attributes_str):
    """
    Parse the attributes column (column 9 in GFF3) into a dictionary.
    
    Args:
        attributes_str: String containing semi-colon separated key=value pairs
        
    Returns:
        Dictionary of parsed attributes
    """
    attributes = OrderedDict()
    
    if not attributes_str or attributes_str == '.':
        return attributes
    
    # Split by semicolon and parse each key=value pair
    for item in attributes_str.split(';'):
        item = item.strip()
        if not item:
            continue
            
        # Handle key=value format
        if '=' in item:
            key, value = item.split('=', 1)
            # URL decode the value (handles %2C -> comma, etc.)
            value = unquote(value)
            attributes[key] = value
    
    return attributes


def parse_filter_argument(filter_str):
    """
    Parse filter argument in format 'key=value', 'key!=value', or 'key=value1,value2'.
    
    Args:
        filter_str: Filter string like:
            - 'column2=gene' (include genes)
            - 'column2!=exon' (exclude exons)
            - 'column0=chr1,chr2,chrX' (include chr1 OR chr2 OR chrX)
        
    Returns:
        Tuple of (filter_type, key, operator, values)
        - filter_type: 'column' or 'attribute'
        - key: column number (for columns) or attribute name
        - operator: '=' or '!='
        - values: list of values to match
    """
    # Check for != operator
    if '!=' in filter_str:
        key, value = filter_str.split('!=', 1)
        operator = '!='
    elif '=' in filter_str:
        key, value = filter_str.split('=', 1)
        operator = '='
    else:
        raise ValueError(f"Invalid filter format: {filter_str}. Expected format: 'key=value', 'key!=value', or 'key=value1,value2'")
    
    # Split comma-separated values for OR logic
    values = [v.strip() for v in value.split(',')]
    
    # Check if it's a column filter (column0-column8)
    if re.match(r'^column[0-8]$', key):
        col_num = int(key.replace('column', ''))
        return ('column', col_num, operator, values)
    else:
        # It's an attribute filter
        return ('attribute', key, operator, values)


def passes_filters(row, filters):
    """
    Check if a row passes all specified filters.
    AND logic between different filters, OR logic within comma-separated values.
    
    Args:
        row: List of column values
        filters: List of parsed filter tuples
        
    Returns:
        Boolean indicating if row passes all filters
    """
    for filter_type, key, operator, values in filters:
        if filter_type == 'column':
            # Column-based filter
            if key >= len(row):
                return False
            
            if operator == '=':
                # Must match at least one value (OR logic)
                if not any(row[key] == v for v in values):
                    return False
            elif operator == '!=':
                # Must not match any value
                if any(row[key] == v for v in values):
                    return False
                    
        elif filter_type == 'attribute':
            # Attribute-based filter
            attributes = parse_attributes(row[8] if len(row) > 8 else '')
            
            if operator == '=':
                # Must have the attribute and match at least one value
                if key not in attributes or not any(attributes[key] == v for v in values):
                    return False
            elif operator == '!=':
                # If attribute exists, must not match any value
                if key in attributes and any(attributes[key] == v for v in values):
                    return False
    
    return True


def order_attribute_keys(all_keys):
    """
    Order attribute keys with priority attributes first, then alphabetical.
    
    Args:
        all_keys: Set or list of all attribute keys
        
    Returns:
        Ordered list of attribute keys
    """
    # Priority order for specific attributes
    priority_attrs = [
        'gene',
        'description',
        'gene_name',
        'gene_synonym',
        'gene_biotype',
        'ID',
        'db_xref',
        'extra_copy_number',
        'copy_num_ID'
    ]
    
    # Create ordered list
    ordered_keys = []
    
    # First add priority attributes that exist
    for attr in priority_attrs:
        if attr in all_keys:
            ordered_keys.append(attr)
    
    # Then add remaining attributes in alphabetical order
    remaining_keys = sorted([k for k in all_keys if k not in priority_attrs])
    ordered_keys.extend(remaining_keys)
    
    return ordered_keys


def collect_all_attribute_keys(input_file, filters):
    """
    First pass through the file to collect all unique attribute keys
    that appear in rows matching the filters.
    
    Args:
        input_file: Path to input GFF3 file
        filters: List of parsed filter tuples
        
    Returns:
        Ordered list of attribute keys
    """
    all_keys = set()
    
    # Determine if file is gzipped
    open_func = gzip.open if input_file.endswith('.gz') else open
    
    try:
        with open_func(input_file, 'rt') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Split the line
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                
                # Check if row passes filters
                if not passes_filters(parts, filters):
                    continue
                
                # Parse attributes and collect keys
                attributes = parse_attributes(parts[8])
                all_keys.update(attributes.keys())
    
    except Exception as e:
        print(f"Error reading file for attribute collection: {e}", file=sys.stderr)
        raise
    
    # Order the keys according to priority
    return order_attribute_keys(all_keys)


def convert_gff_to_bed(input_file, output_file, filters):
    """
    Convert GFF3 file to BED-like format with filtering and attribute parsing.
    
    Args:
        input_file: Path to input GFF3 file
        output_file: Path to output BED file
        filters: List of parsed filter tuples
        
    Returns:
        Number of lines written (excluding header)
    """
    # First pass: collect all attribute keys
    print("Analyzing file structure...", file=sys.stderr)
    attribute_keys = collect_all_attribute_keys(input_file, filters)
    
    if not attribute_keys:
        print("Warning: No attribute keys found in filtered rows.", file=sys.stderr)
    
    # Determine if file is gzipped
    open_func = gzip.open if input_file.endswith('.gz') else open
    
    lines_written = 0
    
    try:
        with open_func(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
            # Write header with reordered columns
            header_parts = ['chrom', 'start', 'end', 'source', 'type', 'score', 'strand', 'phase']
            header_parts.extend(attribute_keys)
            outfile.write('#' + '\t'.join(header_parts) + '\n')
            
            # Process each line
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Split the line
                parts = line.split('\t')
                if len(parts) < 9:
                    print(f"Warning: Line {line_num} has fewer than 9 columns, skipping.", 
                          file=sys.stderr)
                    continue
                
                # Check if row passes filters
                if not passes_filters(parts, filters):
                    continue
                
                # Extract columns in original GFF3 order
                chrom = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = str(int(parts[3]) - 1)  # Convert to 0-based
                end = parts[4]
                score = parts[5]
                strand = parts[6]
                phase = parts[7]
                
                # Parse attributes
                attributes = parse_attributes(parts[8])
                
                # Build output row with reordered columns
                # New order: chrom, start, end, source, type, score, strand, phase
                output_parts = [chrom, start, end, source, feature_type, score, strand, phase]
                
                # Add attribute values in the order of attribute_keys
                for key in attribute_keys:
                    value = attributes.get(key, '')
                    output_parts.append(value)
                
                # Write the row
                outfile.write('\t'.join(output_parts) + '\n')
                lines_written += 1
    
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        raise
    
    return lines_written


def main():
    """Main function to handle command-line arguments and run conversion."""
    parser = argparse.ArgumentParser(
        description='Convert GFF3 files to BED-like format with dynamic filtering and attribute parsing.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filter Syntax:
  key=value         Include rows where key equals value
  key!=value        Exclude rows where key equals value
  key=val1,val2     Include rows where key equals val1 OR val2

Examples:
  # Default: filter for genes with protein_coding biotype
  %(prog)s input.gff3
  
  # Filter for exons
  %(prog)s input.gff3 -f column2=exon
  
  # Exclude pseudogenes
  %(prog)s input.gff3 -f column2=gene -f gene_biotype!=pseudogene
  
  # Multiple chromosomes (OR logic)
  %(prog)s input.gff3 -f column0=chr1,chr2,chrX
  
  # Complex filtering (AND between filters, OR within values)
  %(prog)s input.gff3 -f column2=gene,transcript -f column6=+ -f column0!=chrM
  
  # Custom output file
  %(prog)s input.gff3 -o output.bed
        """
    )
    
    parser.add_argument('input', 
                        help='Input GFF3 file (can be gzipped)')
    
    parser.add_argument('-o', '--output', 
                        help='Output BED file (default: input_filtered.bed)')
    
    parser.add_argument('-f', '--filter', 
                        action='append',
                        dest='filters',
                        help='Filter criteria: "key=value", "key!=value", or "key=val1,val2". '
                             'Can specify multiple times (AND logic between filters). '
                             'Use "column0" through "column8" for column filters, or attribute names. '
                             'Default: column2=gene and gene_biotype=protein_coding')
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Set default output file if not specified
    if not args.output:
        base_name = args.input
        # Remove .gz extension if present
        if base_name.endswith('.gz'):
            base_name = base_name[:-3]
        # Remove .gff3 extension if present
        if base_name.endswith('.gff3'):
            base_name = base_name[:-5]
        elif base_name.endswith('.gff'):
            base_name = base_name[:-4]
        args.output = base_name + '_filtered.bed'
    
    # Parse filters
    if not args.filters:
        # Default filters
        filter_strs = ['column2=gene', 'gene_biotype=protein_coding']
    else:
        filter_strs = args.filters
    
    # Parse filter strings
    filters = []
    try:
        for filter_str in filter_strs:
            filters.append(parse_filter_argument(filter_str))
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Display filters being used
    print(f"Input file: {args.input}", file=sys.stderr)
    print(f"Output file: {args.output}", file=sys.stderr)
    print("Filters:", file=sys.stderr)
    for filter_str in filter_strs:
        print(f"  - {filter_str}", file=sys.stderr)
    
    # Perform conversion
    try:
        lines_written = convert_gff_to_bed(args.input, args.output, filters)
        print(f"\nConversion complete. {lines_written} lines written to {args.output}")
    except Exception as e:
        print(f"\nError during conversion: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
