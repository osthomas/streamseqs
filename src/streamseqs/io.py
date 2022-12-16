# io.py

"""
Utilities for input/output operations of (gzip compressed) FASTA/FASTQ files.
"""

import sys
import gzip
import logging
logger = logging.getLogger(__name__)


def stream(fpath):
    """
    Stream a file or stdin line by line.

    Parameters
    ----------
        fpath
            Path to an input file or `sys.stdin` to read from stdin.

    Yields
    ------
        line : str
            Lines from the input stream without trailing newline.
    """

    # No file given: read from stdin
    if fpath is sys.stdin:
        if sys.stdin.isatty():
            raise IOError("No input streamed to stdin!")
        else:
            logger.info("Reading from stdin.")
            for line in sys.stdin:
                yield line.strip()
    elif is_gzip(fpath):
        logger.info("Reading gzip compressed file.")
        with gzip.open(fpath, "r") as fh:
            for line in fh:
                yield line.strip()
    else:
        logger.info("Reading uncompressed file.")
        with open(fpath, "r") as fh:
            for line in fh:
                yield line.strip()


def is_gzip(fpath):
    """
    Check if a file is gzip compressed.
    """
    # gzip.BadGzipFile was added in version 3.8
    try:
        exc = gzip.BadGzipFile
    except AttributeError:
        exc = OSError
    with gzip.open(fpath, "r") as fh:
        try:
            fh.read(1)
            return True
        except exc:
            return False


def identify_fastx(stream):
    """
    Check if input is FASTA or FASTQ.

    The stream might be stdin - in that case, it is not possible to backtrack,
    therefore the identified first sequence id must also be returned to
    construct the first record.

    Returns
    -------
        tuple
            (sequence_type, first_id)
    """
    for i, line in enumerate(stream):
        if line == "":
            # Skip leading empty lines
            continue
        elif line.startswith(">"):
            seq_type = "FASTA"
            break
        elif line.startswith("@"):
            seq_type = "FASTQ"
            break
        else:
            raise IOError("Input is not a valid FASTA or FASTQ file!")
    else:
        # Empty file
        seq_type = "empty"
    logger.info(f"Identified sequence type: {seq_type}")

    return (seq_type, line[1:])


def new_record(seq_id, seq, qual):
    """
    Construct a sequence record as a dictionary.
    """
    return {"id": seq_id, "seq": seq, "qual": qual}


def fasta_records(stream, first_id):
    """
    Read lines until the next sequence is identified and construct a FASTA
    record.

    Yield records until the stream is exhausted.

    Parameters
    ----------
        stream
            Lines read from a FASTA file
        first_id
            Sequence ID in the first record of the file.
    """
    seq_id = first_id
    seq = []
    for line in stream:
        if not line.startswith(">"):
            seq.append(line)
        else:
            yield new_record(seq_id, "".join(seq), None)
            seq_id = line[1:]
            seq = []
    # Yield final record in the stream
    yield new_record(seq_id, "".join(seq), None)


def fastq_records(stream, first_id):
    """
    Read lines until the next sequence is identified and construct a FASTQ
    record.

    Yield records until the stream is exhausted.

    Parameters
    ----------
        stream
            Lines read from a FASTQ file
        first_id
            Sequence ID in the first record of the file.
    """
    first = True
    seq_id = first_id
    while True:
        if not first:
            try:
                seq_id = next(stream)[1:]
            except StopIteration:
                return
        try:
            seq = next(stream)
            plus = next(stream)
            qual = next(stream)
            first = False
        except StopIteration:
            raise IOError("Incomplete final record in FASTQ file!")
        if plus != "+" or len(seq) != len(qual):
            raise IOError("Invalid FASTQ record (ID: %s)" % seq_id)
        yield new_record(seq_id, seq, qual)


def stream_records(fpath):
    """
    Yield sequence records (FASTA or FASTQ) from a stream.

    Yields
    ------
        record : dict
            A dict with the following fields:
                "id": str, sequence ID / header
                "seq": str, the sequence
                "qual": str or None, the sequence quality.
    """
    _stream = stream(fpath)
    seq_type, first_id = identify_fastx(_stream)
    if seq_type == "FASTA":
        func = fasta_records
    elif seq_type == "FASTQ":
        func = fastq_records
    elif seq_type == "empty":
        return
    for record in func(_stream, first_id):
        yield record


def format_record(record):
    """
    Format record for writing.

    Parameters
    ----------
        record : dict
            A dict with the following fields:
                "id": str, sequence ID / header
                "seq": str, the sequence
                "qual": str (optional), the sequence quality.
            If "qual" is None, FASTA will be written, else FASTQ.
    """
    if record["qual"] is None:
        # FASTA
        return f">{record['id']}\n{record['seq']}\n"
    else:
        # FASTQ
        return f"@{record['id']}\n{record['seq']}\n+\n{record['qual']}\n"
