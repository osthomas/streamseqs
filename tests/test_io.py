import pytest
import pkg_resources

from streamseqs import io


def _testfile(fname):
    """
    Get path to a resource test file in directory "tests/resources".
    """
    return pkg_resources.resource_filename("tests.resources", fname)


def test_input_type_detected():
    fasta_stream = io.stream(_testfile("seqs.fasta"))
    assert io.identify_fastx(fasta_stream) == ("FASTA", "Seq1")
    fastq_stream = io.stream(_testfile("seqs.fastq"))
    assert io.identify_fastx(fastq_stream) == ("FASTQ", "Seq1")


def test_empty_file():
    records = list(io.stream_records(_testfile("empty.fasta")))
    assert len(records) == 0


def test_leading_empty_lines():
    fpath = _testfile("leading_empty_lines.fasta")
    expected = [
        io.new_record("Seq1", "ATGC", None),
        io.new_record("Seq2", "CCCCC", None)
    ]
    for r, e in zip(io.stream_records(fpath), expected):
        assert r == e


def test_parse_fasta_records():
    fpath = _testfile("seqs.fasta")
    expected = [
        io.new_record("Seq1", "ATGC", None),
        io.new_record("Seq2", "GATC", None),
        io.new_record("Seq3", "ATGCTGAGG", None)
    ]
    for r, e in zip(io.stream_records(fpath), expected):
        assert r == e


def test_parse_fastq_records():
    fpath = _testfile("seqs.fastq")
    expected = [
        io.new_record("Seq1", "ATGC", "JJJJ"),
        io.new_record("Seq2", "GATC", "HHHH")
    ]
    for r, e in zip(io.stream_records(fpath), expected):
        assert r == e


def test_incomplete_final_fastq_record_raises_error():
    fpath = _testfile("invalid2.fastq")
    with pytest.raises(IOError):
        for _ in io.stream_records(fpath):
            pass


def test_erroneous_fastq_raises_error():
    fpath = _testfile("invalid.fastq")
    with pytest.raises(IOError):
        for _ in io.stream_records(fpath):
            pass
