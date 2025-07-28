from pathlib import Path
from xml.etree.ElementTree import ParseError

# Standard BLAST6 format columns
_DEFAULT_BLAST6_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

def _get_blast6_defaults():
    """Get BLAST6 defaults from sunbeam, with fallback to standard format."""
    try:
        from sunbeam.bfx.parse import BLAST6_DEFAULTS
        return BLAST6_DEFAULTS
    except ImportError:
        try:
            # For compatibility with older versions of Sunbeam
            from sunbeamlib.parse import BLAST6_DEFAULTS
            return BLAST6_DEFAULTS
        except ImportError:
            # Fallback to standard BLAST6 format
            return _DEFAULT_BLAST6_COLUMNS


def parse_blast6(f, outfmt=None):
    if outfmt is None:
        outfmt = _get_blast6_defaults()
    for line in f.readlines():
        vals = line.strip().split("\t")
        if len(outfmt) == len(vals):
            yield dict(zip(outfmt, vals))


def blast_summary(blast_files):
    """Summarize BLAST results from an set of BLAST output files."""
    for infile in blast_files:
        sample = Path(infile).stem
        try:
            with open(infile) as f:
                for result in parse_blast6(f):
                    if len(result) > 0:
                        yield {
                            "sample": sample,
                            "query": result["qseqid"],
                            "hit": result["sseqid"],
                        }
        except ParseError:
            print("Skipping empty/malformed %s" % infile)
            continue
