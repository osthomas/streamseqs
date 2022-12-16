# streamseqs
![pytest-badge](https://github.com/osthomas/streamseqs/actions/workflows/pytest.yaml/badge.svg)

`streamseqs` is a tiny package to quickly read in sequences from stdin or from
(compressed) FAST(A/Q) files from disk.

It only performs barebones error checking and parsing and expects its input to
be reasonably sane. If you need more sophisticated parsing and error checking,
you are probably better of using [Biopython](https://github.com/biopython/biopython).

## Installation

```bash
pip install streamseqs
```

## Usage

```python
from streamseqs import stream_records

# Stream from file
for record in stream_records("seqs.fasta"):
    process(record)

# Records are dicts:
# {"id": "sequence header", "seq": "ATGCT", "qual": "HHHHH"}
# For FASTA files, "qual" is `None`.

# Gzip compressed files are transparently handled
for record in stream_records("seqs.fastq.gz"):
    process(record)

# Can stream from stdin. Expects a stream, not a TTY!
import sys
for record in stream_records(sys.stdin):
    process(record)

# Format record for writing - as FASTQ if "qual" is present, as FASTA if not.
from streamseqs import format_record
format_record(record)
```
