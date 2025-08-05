# CDS Processing Logging in BOFASA

## Overview

BOFASA now includes silent logging functionality for CDS features that are skipped or encounter issues during GenBank file processing. This logging system helps track and analyze problems with CDS feature processing without interrupting the main workflow.

## Features

- **Silent Operation**: Logging occurs in the background without affecting normal processing
- **Parallel-Safe**: Each sample writes to its own log file, avoiding write lock issues
- **Comprehensive Tracking**: Logs various types of CDS processing issues
- **Summary Reports**: Provides summary statistics at the end of processing

## Log File Format

Individual log files `cds_processing_issues_{sample_name}.log` are created in the output directory for each sample and contain tab-separated entries with the following columns:

1. **Timestamp**: When the issue occurred
2. **Scaffold**: Scaffold/contig identifier
3. **Feature_Info**: Information about the problematic feature
4. **Issue_Type**: Type of issue encountered
5. **Details**: Additional details about the issue

## Issue Types

The following types of CDS processing issues are logged:

- **no_cds_features**: No CDS features found in a record
- **complex_location**: Complex location formats (join{}, order{}) that are skipped
- **location_parse_error**: Location strings that cannot be parsed
- **location_parse_exception**: Exceptions during location parsing
- **no_translation**: CDS features without translation qualifier
- **short_protein**: Proteins shorter than the minimum length threshold

## Usage

The logging is automatically enabled when using the following functions:

- `process_genbank_file()`: Individual GenBank file processing
- `process_genomes_as_genbanks()`: Batch processing of GenBank files
- `process_annotation_directories()`: Processing annotation directories

### Example

```python
from bofasa.processing import process_genbank_file

# Process a GenBank file - logging happens automatically
process_genbank_file(
    input_file="sample.gbff",
    outdir="output",
    sample_name="sample1",
    locus_tag="SAMPLE1",
    min_length=20
)

# Check the log file
import os
log_file = os.path.join("output", "cds_processing_issues_sample1.log")
if os.path.exists(log_file):
    with open(log_file, 'r') as f:
        for line in f:
            print(line.strip())
```

## Output

### Individual Processing

When processing individual files, the function will print a message if issues were logged:

```
Processed 1234 proteins from 1 records in GenBank file.
CDS processing issues logged to: output/cds_processing_issues_sample1.log (5 issues)
```

### Batch Processing

When processing multiple genomes, a summary is provided at the end:

```
CDS processing summary: 15 total issues across 3 samples
  - complex_location: 8 issues
  - short_protein: 4 issues
  - translation_error: 3 issues
  - sample1: 5 issues
  - sample2: 7 issues
  - sample3: 3 issues
Detailed CDS processing issues logged to individual files: cds_processing_issues_*.log
```

## Testing

You can test the logging functionality using the provided test script:

```bash
python test_cds_logging.py
```

This script processes a test GenBank file and demonstrates the logging output.

## Implementation Details

### Thread Safety

The logging system uses:
- Individual log files per sample to avoid write locks
- Direct file appending for each sample
- No shared resources between parallel processes

### Error Handling

- Logging failures do not interrupt the main processing
- Each sample's logging is independent
- Graceful handling of file system errors

### Performance

- Minimal overhead during normal processing
- Efficient per-sample file management
- Non-blocking log file operations

## File Locations

- **Sample log files**: `{output_dir}/cds_processing_issues_{sample_name}.log`
- **One file per sample**: Each sample gets its own log file

## Troubleshooting

If you encounter issues with the logging system:

1. **Check file permissions**: Ensure write access to the output directory
2. **Monitor disk space**: Log files are created for each sample
3. **Review log files**: Check individual sample log files for detailed error information
4. **File organization**: Each sample has its own log file for easy identification

The logging system is designed to be robust and should not interfere with normal BOFASA operations. 