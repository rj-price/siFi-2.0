"""Regenerate the committed inputs in tests/data/ (PLAN.md Phase 1).

Runs under the normal Python 3 `sifi2` env, since it needs bowtie and RNAplfold:

    conda run -n sifi2 python tests/py2_capture/make_inputs.py

You should not normally need this — everything it writes is committed. It exists
so the inputs are reproducible rather than magic, and so the reference set can be
changed deliberately. Changing anything here invalidates every golden fixture:
re-run tests/py2_capture/capture.py afterwards and review the diff.

Needs network access on the first run, to fetch the sequences from NCBI.

Produces:
    tests/data/reference.fasta          the beta-tubulin paralogue family
    tests/data/query.fasta              500 bp of TUB3
    tests/data/query_multi.fasta        3 records, for the Phase 4 batch test
    tests/data/bowtie_output_mm2.txt    real bowtie 1.3.1 output
    tests/data/TUB3_fragment_lunp       real RNAplfold 2.7.2 output
"""

import os
import shutil
import subprocess
import tempfile
import textwrap
import urllib.request

from Bio import SeqIO

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "tests", "data")

# Constants from legacy/sifi_pipeline.py, which must not drift from it.
SIRNA_SIZE = 21
MISMATCHES = 2
WINSIZE = 80
SPAN = 40
TEMPERATURE = 22

# Arabidopsis thaliana beta-tubulin mRNAs. A real paralogue family, so a query
# drawn from one member genuinely cross-hits the others: real main targets and
# real off-targets, which is the distinction siFi exists to make.
ACCESSIONS = [
    "NM_125665.4",      # TUB3 — the query's own gene, so the main target
    "NM_125664.4",      # TUB2 — near-identical to TUB3
    "NM_123801.2",      # TUB4
    "NM_122291.4",      # TUB8
    "NM_001203444.1",   # TUB8, second variant
]
QUERY_ACCESSION = "NM_125665.4"
QUERY_SLICE = (300, 800)
# Deposited antisense, so the fixtures contain minus-strand bowtie hits — see
# PLAN.md Phase 6 defect 1.
RC_ACCESSION = "NM_123801.2"


def fetch(path):
    """Fetch the accessions from NCBI, cached at `path`."""
    if os.path.exists(path):
        return path
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=%s&rettype=fasta&retmode=text" % ",".join(ACCESSIONS)
    )
    with urllib.request.urlopen(url) as response, open(path, "wb") as fh:
        fh.write(response.read())
    return path


def wrap(seq, width=70):
    return "\n".join(textwrap.wrap(str(seq), width))


def write_fastas(recs):
    with open(os.path.join(DATA, "reference.fasta"), "w") as fh:
        for acc in ACCESSIONS:
            rec = recs[acc]
            fh.write(">%s %s\n%s\n" % (rec.id, rec.description.split(" ", 1)[1], wrap(rec.seq)))
        rec = recs[RC_ACCESSION]
        fh.write(
            ">%s_rc %s, reverse complement\n%s\n"
            % (RC_ACCESSION, rec.description.split(" ", 1)[1], wrap(rec.seq.reverse_complement()))
        )

    start, stop = QUERY_SLICE
    with open(os.path.join(DATA, "query.fasta"), "w") as fh:
        fh.write(
            ">TUB3_fragment Arabidopsis thaliana tubulin beta chain 3 (TUB3) mRNA, bases %d-%d\n%s\n"
            % (start + 1, stop, wrap(recs[QUERY_ACCESSION].seq[start:stop]))
        )

    with open(os.path.join(DATA, "query_multi.fasta"), "w") as fh:
        for acc, name, lo, hi in [
            ("NM_125665.4", "TUB3_fragment", 300, 800),
            ("NM_123801.2", "TUB4_fragment", 300, 700),
            ("NM_122291.4", "TUB8_fragment", 300, 700),
        ]:
            fh.write(
                ">%s bases %d-%d of %s\n%s\n" % (name, lo + 1, hi, acc, wrap(recs[acc].seq[lo:hi]))
            )


def sirna_fasta(query_sequence, path):
    """Window into siRNAs byte-for-byte as legacy create_sirnas does — note the
    '> sirnaN' header really does have a space after the '>'."""
    start, end = 0, SIRNA_SIZE
    with open(path, "w") as fh:
        for _ in range(len(query_sequence)):
            if len(query_sequence[start:end]) == SIRNA_SIZE:
                fh.write("> sirna%d\n%s\n" % (start + 1, query_sequence[start:end]))
                start += 1
                end += 1


def run_bowtie(query_sequence, work):
    """Real bowtie 1.3.1, with legacy run_bowtie's exact argv."""
    index = os.path.join(work, "testdb")
    subprocess.check_call(
        ["bowtie-build", "-q", os.path.join(DATA, "reference.fasta"), index]
    )
    reads = os.path.join(work, "sirnas.fasta")
    sirna_fasta(query_sequence, reads)
    out = os.path.join(DATA, "bowtie_output_mm2.txt")
    subprocess.check_call(
        ["bowtie", "-a", "-v", str(MISMATCHES), "-y", index, "-f", reads, out]
    )
    return out


def run_rnaplfold(query_name, query_sequence, work):
    """Real RNAplfold 2.7.2, with legacy run_rnaplfold's exact argv and stdin."""
    prc = subprocess.Popen(
        [
            "RNAplfold",
            "-W", "%d" % WINSIZE,
            "-L", "%d" % SPAN,
            "-u", "%d" % SIRNA_SIZE,
            "-T", "%.2f" % TEMPERATURE,
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        cwd=work,
        text=True,
    )
    # legacy create_single_fasta_file writes '>' + id then the unwrapped sequence.
    prc.communicate(input=">%s\n%s\n\n" % (query_name, query_sequence))
    dest = os.path.join(DATA, query_name + "_lunp")
    shutil.copyfile(os.path.join(work, query_name + "_lunp"), dest)
    return dest


def main():
    work = tempfile.mkdtemp()
    recs = {
        r.id: r for r in SeqIO.parse(fetch(os.path.join(work, "tubulins.fasta")), "fasta")
    }
    missing = set(ACCESSIONS) - set(recs)
    if missing:
        raise SystemExit("NCBI did not return: %s" % ", ".join(sorted(missing)))

    write_fastas(recs)
    rec = next(SeqIO.parse(os.path.join(DATA, "query.fasta"), "fasta"))
    query_sequence = str(rec.seq)

    written = [
        os.path.join(DATA, name)
        for name in ("reference.fasta", "query.fasta", "query_multi.fasta")
    ]
    written.append(run_bowtie(query_sequence, work))
    written.append(run_rnaplfold(rec.id, query_sequence, work))

    print()
    for path in written:
        print("wrote %-40s %9d bytes" % (os.path.relpath(path, REPO), os.path.getsize(path)))
    print("\nNow re-run: conda run -n sifi2-py2 python tests/py2_capture/capture.py")


if __name__ == "__main__":
    main()
