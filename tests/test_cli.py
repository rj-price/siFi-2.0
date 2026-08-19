"""Phase 4: ``sifi2.cli`` — argument wiring, output files and batch mode.

Most of this exercises the CLI without touching bowtie or RNAplfold: the
question here is whether flags reach :class:`SifiConfig` intact and whether the
files land where they should. The two tests that do run the real pipeline are
marked ``needs_binaries`` and skip when the tools are not on PATH.

The load-bearing one is :func:`test_batch_mode_writes_one_result_set_per_record`
— the whole point of having a CLI is that the original's ``return`` inside the
FASTA loop meant only the first record was ever analysed.
"""

import json
import shutil

import pytest
from conftest import DATA

from sifi2 import bowtie
from sifi2.cli import (
    CliError,
    build_parser,
    config_from_args,
    default_db_location,
    main,
    make_main_target_selector,
    read_main_target_file,
    read_queries,
    safe_stem,
)
from sifi2.config import Mode

needs_binaries = pytest.mark.skipif(
    shutil.which("bowtie") is None or shutil.which("RNAplfold") is None,
    reason="needs bowtie and RNAplfold on PATH",
)


def parse(argv):
    return build_parser().parse_args(argv)


BASE = ["--query", "q.fasta", "--db", "mydb"]


# ----------------------------------------------------------------------
# Flags reaching the config
# ----------------------------------------------------------------------
def test_defaults_match_the_gui_design_settings():
    """``main.py``'s ``set_design_settings``: all four rules on, 21mers, 0 mismatches."""
    config = config_from_args(parse(["design", *BASE, "--all-targets-main"]), Mode.DESIGN)
    assert config.mode is Mode.DESIGN
    assert (config.sirna_size, config.mismatches) == (21, 0)
    assert (config.strand_check, config.terminal_check, config.end_check, config.accessibility_check) == (
        True,
        True,
        True,
        True,
    )
    assert (config.end_stability_threshold, config.accessibility_threshold, config.accessibility_window) == (
        1.0,
        0.1,
        8,
    )
    assert config.no_efficience is False


def test_defaults_match_the_gui_offtarget_settings():
    """``set_offtarget_settings`` turns end stability and accessibility off."""
    config = config_from_args(parse(["offtarget", *BASE]), Mode.OFFTARGET)
    assert config.mode is Mode.OFFTARGET
    assert (config.end_check, config.accessibility_check) == (False, False)
    assert (config.strand_check, config.terminal_check) == (True, True)
    # The thresholds keep their values even though the rules are off.
    assert (config.end_stability_threshold, config.accessibility_threshold) == (1.0, 0.1)


@pytest.mark.parametrize("mode", [Mode.DESIGN, Mode.OFFTARGET])
def test_mode_dependent_rules_are_overridable_in_both_directions(mode):
    command = "design" if mode is Mode.DESIGN else "offtarget"
    extra = ["--all-targets-main"] if mode is Mode.DESIGN else []
    on = config_from_args(parse([command, *BASE, *extra, "--end-stability", "--accessibility"]), mode)
    off = config_from_args(parse([command, *BASE, *extra, "--no-end-stability", "--no-accessibility"]), mode)
    assert (on.end_check, on.accessibility_check) == (True, True)
    assert (off.end_check, off.accessibility_check) == (False, False)


def test_mismatches_disable_efficiency_scoring_as_the_gui_did():
    """The GUI disabled every efficiency widget once mismatches > 0, and derived
    ``no_efficience`` from exactly that."""
    assert config_from_args(parse(["offtarget", *BASE]), Mode.OFFTARGET).no_efficience is False
    assert config_from_args(parse(["offtarget", *BASE, "--mismatches", "2"]), Mode.OFFTARGET).no_efficience is True
    # ...but it is still overridable either way.
    with_mm = parse(["offtarget", *BASE, "--mismatches", "2", "--efficiency"])
    assert config_from_args(with_mm, Mode.OFFTARGET).no_efficience is False
    assert config_from_args(parse(["offtarget", *BASE, "--no-efficiency"]), Mode.OFFTARGET).no_efficience is True


def test_previously_hidden_constants_are_now_flags():
    """The six values the original buried in the pipeline body (PLAN.md Phase 3)."""
    argv = [
        "offtarget", *BASE,
        "--winsize", "100", "--span", "50", "--temperature", "37",
        "--sirna-start-position", "5", "--overhang", "3", "--end-nucleotides", "4",
    ]  # fmt: skip
    config = config_from_args(parse(argv), Mode.OFFTARGET)
    assert (config.winsize, config.span, config.temperature) == (100, 50, 37)
    assert (config.sirna_start_position, config.overhang, config.end_nucleotides) == (5, 3, 4)


def test_db_location_defaults_to_platformdirs_not_a_windows_path():
    """``main.py`` did ``storageLocation(...).split('Local')[0] + '/Local/'``,
    which yields nonsense on any path with no ``Local`` component."""
    config = config_from_args(parse(["offtarget", *BASE]), Mode.OFFTARGET)
    assert config.db_location == default_db_location()
    assert "/Local/" not in config.db_location
    override = config_from_args(parse(["offtarget", *BASE, "--db-location", "/data/dbs"]), Mode.OFFTARGET)
    assert override.db_location == "/data/dbs"


def test_binary_overrides_are_passed_through():
    config = config_from_args(
        parse(["offtarget", *BASE, "--bowtie-path", "/opt/bowtie", "--rnaplfold-path", "/opt/RNAplfold"]),
        Mode.OFFTARGET,
    )
    assert (config.bowtie_path, config.rnaplfold_path) == ("/opt/bowtie", "/opt/RNAplfold")


# ----------------------------------------------------------------------
# Main target selection
# ----------------------------------------------------------------------
def test_design_requires_a_main_target_flag():
    """The pipeline raises ``PipelineCancelled`` without a selector, so the CLI
    must never reach it empty-handed."""
    with pytest.raises(SystemExit):
        parse(["design", *BASE])


def test_all_targets_main_selects_every_hit():
    selector = make_main_target_selector(parse(["design", *BASE, "--all-targets-main"]))
    assert selector([("TUB3", 300), ("TUB4", 120)]) == ["TUB3", "TUB4"]


def test_named_main_targets_are_returned_verbatim():
    argv = ["design", *BASE, "--main-target", "TUB3", "--main-target", "TUB4"]
    selector = make_main_target_selector(parse(argv))
    assert selector([("TUB3", 300), ("TUB8", 120)]) == ["TUB3", "TUB4"]


def test_a_main_target_that_matched_nothing_warns_but_does_not_fail(capsys):
    """Almost always a typo, and the consequence — every hit reported as an
    off-target — is silent otherwise."""
    selector = make_main_target_selector(parse(["design", *BASE, "--main-target", "TYPO"]))
    assert selector([("TUB3", 300)]) == ["TYPO"]
    assert "matched no hit" in capsys.readouterr().err


def test_main_target_file_skips_blanks_and_comments(tmp_path):
    target_file = tmp_path / "targets.txt"
    target_file.write_text("# my targets\nTUB3\n\n  TUB4  \n")
    assert read_main_target_file(target_file) == ["TUB3", "TUB4"]


def test_an_empty_main_target_file_is_an_error(tmp_path):
    target_file = tmp_path / "empty.txt"
    target_file.write_text("# nothing here\n")
    with pytest.raises(CliError, match="lists no target names"):
        read_main_target_file(target_file)


# ----------------------------------------------------------------------
# Query reading and filenames
# ----------------------------------------------------------------------
def test_read_queries_returns_every_record():
    queries = read_queries(str(DATA / "query_multi.fasta"))
    assert [name for name, _ in queries] == ["TUB3_fragment", "TUB4_fragment", "TUB8_fragment"]
    assert all(sequence for _, sequence in queries)


def test_duplicate_record_ids_are_rejected(tmp_path):
    """Two records with one id would write to one filename, silently losing the first."""
    query = tmp_path / "dup.fasta"
    query.write_text(">a\nACGT\n>a\nTGCA\n")
    with pytest.raises(CliError, match="duplicate record ids"):
        read_queries(str(query))


def test_an_empty_query_file_is_an_error(tmp_path):
    query = tmp_path / "empty.fasta"
    query.write_text("")
    with pytest.raises(CliError, match="no FASTA records"):
        read_queries(str(query))


@pytest.mark.parametrize(
    ("query_name", "expected"),
    [
        ("TUB3_fragment", "TUB3_fragment"),
        ("gi|123|ref|NM_001.1|", "gi_123_ref_NM_001.1"),
        ("seq with spaces", "seq_with_spaces"),
        ("...", "query"),
    ],
)
def test_safe_stem_makes_a_fasta_id_usable_as_a_filename(query_name, expected):
    assert safe_stem(query_name) == expected


# ----------------------------------------------------------------------
# End to end
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def test_db(tmp_path_factory):
    """A bowtie index built from the committed reference FASTA."""
    if shutil.which("bowtie-build") is None:
        pytest.skip("needs bowtie-build on PATH")
    db_location = tmp_path_factory.mktemp("dbs")
    bowtie.build_database("testdb", str(DATA / "reference.fasta"), str(db_location))
    return str(db_location)


@needs_binaries
def test_batch_mode_writes_one_result_set_per_record(test_db, tmp_path, capsys):
    """The reason a CLI is worth having: the original ``return``ed inside the
    ``SeqIO.parse`` loop, so records 2 and 3 were silently dropped."""
    outdir = tmp_path / "out"
    argv = [
        "offtarget", "--query", str(DATA / "query_multi.fasta"),
        "--db", "testdb", "--db-location", test_db, "--outdir", str(outdir),
    ]  # fmt: skip
    assert main(argv) == 0

    stems = ["TUB3_fragment", "TUB4_fragment", "TUB8_fragment"]
    assert sorted(path.name for path in outdir.iterdir()) == sorted(
        [f"{stem}.{extension}" for stem in stems for extension in ("json", "tsv", "png")]
    )

    # Record counts at --mismatches 0, after Phase 6 defect 1 (minus-strand hits
    # were dropped): 900/413/811 before the fix. TUB4 gains most, since the
    # reference holds it reverse-complemented and so it hits its own deposit.
    counts = [len(json.loads((outdir / f"{stem}.json").read_text())) for stem in stems]
    assert counts == [912, 793, 820]

    for stem, count in zip(stems, counts, strict=True):
        lines = (outdir / f"{stem}.tsv").read_text().splitlines()
        assert len(lines) == count + 1  # + the header
        assert lines[0].split("\t")[0] == "query_name"


@needs_binaries
def test_threaded_and_serial_runs_agree(test_db, tmp_path):
    outputs = {}
    for threads in ("1", "3"):
        outdir = tmp_path / f"threads{threads}"
        argv = [
            "design", "--query", str(DATA / "query_multi.fasta"), "--db", "testdb",
            "--db-location", test_db, "--outdir", str(outdir), "--all-targets-main",
            "--threads", threads, "--quiet", "--no-plot",
        ]  # fmt: skip
        assert main(argv) == 0
        outputs[threads] = {path.name: path.read_text() for path in sorted(outdir.iterdir())}
    assert outputs["1"] == outputs["3"]


@needs_binaries
def test_design_marks_unchosen_hits_as_off_targets(test_db, tmp_path):
    outdir = tmp_path / "out"
    argv = [
        "design", "--query", str(DATA / "query.fasta"), "--db", "testdb",
        "--db-location", test_db, "--outdir", str(outdir), "--main-target", "NM_125665.4", "--quiet",
    ]  # fmt: skip
    assert main(argv) == 0
    records = json.loads((outdir / "TUB3_fragment.json").read_text())
    assert {record["is_off_target"] for record in records} == {True, False}
    assert {record["hit_name"] for record in records if not record["is_off_target"]} == {"NM_125665.4"}


#: A random sequence, so no 21mer of it is in the tubulin reference.
NO_HIT_FASTA = """>nohits
ATGAACTGGAGTCTACGATGAGTGTACGAACGTCAGCTGGAACAGGCTTCCCACCAGGGT
TGCTACTTATCATTTATTGTACGTTCAAAGGCGTGGTTTGTTTCTTGTGGCTGGTTCGAT
"""


@needs_binaries
def test_design_without_database_hits_still_writes_efficiency_results(test_db, tmp_path, capsys):
    """PLAN.md Phase 6 defect 2: designing against a sequence with no database
    hits is a supported case, and now produces the efficiency-only output the
    original provided for but discarded."""
    query = tmp_path / "nohits.fasta"
    query.write_text(NO_HIT_FASTA)
    outdir = tmp_path / "out"
    argv = [
        "design", "--query", str(query), "--db", "testdb",
        "--db-location", test_db, "--outdir", str(outdir), "--all-targets-main",
    ]  # fmt: skip
    assert main(argv) == 0

    records = json.loads((outdir / "nohits.json").read_text())
    assert len(records) == 120 - 21 + 1
    assert {record["hit_name"] for record in records} == {None}
    assert sorted(path.name for path in outdir.iterdir()) == ["nohits.json", "nohits.png", "nohits.tsv"]
    # The run is still reported as having found no targets.
    assert "No targets found" in capsys.readouterr().out


@needs_binaries
def test_offtarget_without_database_hits_produces_nothing(test_db, tmp_path, capsys):
    """Off-target mode has no efficiency-only fallback: with no hits there is
    nothing to predict, so the run reports that and exits non-zero."""
    query = tmp_path / "nohits.fasta"
    query.write_text(NO_HIT_FASTA)
    outdir = tmp_path / "out"
    argv = [
        "offtarget", "--query", str(query), "--db", "testdb",
        "--db-location", test_db, "--outdir", str(outdir),
    ]  # fmt: skip
    assert main(argv) == 1
    assert "No targets found" in capsys.readouterr().out


# ----------------------------------------------------------------------
# Failures
# ----------------------------------------------------------------------
def test_a_missing_query_file_is_reported_not_raised(tmp_path, capsys):
    argv = ["offtarget", "--query", "nope.fasta", "--db", "d", "--db-location", str(tmp_path)]
    assert main(argv) == 1
    assert "query file not found" in capsys.readouterr().err


def test_a_missing_database_names_the_build_command(tmp_path, capsys):
    argv = ["offtarget", "--query", str(DATA / "query.fasta"), "--db", "nodb", "--db-location", str(tmp_path)]
    assert main(argv) == 1
    assert "sifi db build" in capsys.readouterr().err


def test_zero_threads_is_rejected(tmp_path, capsys):
    argv = [
        "offtarget", "--query", str(DATA / "query.fasta"), "--db", "d",
        "--db-location", str(tmp_path), "--threads", "0",
    ]  # fmt: skip
    assert main(argv) == 1
    assert "--threads" in capsys.readouterr().err


# ----------------------------------------------------------------------
# Database management
# ----------------------------------------------------------------------
def test_db_build_list_and_remove_round_trip(tmp_path, capsys):
    if shutil.which("bowtie-build") is None:
        pytest.skip("needs bowtie-build on PATH")
    db_location = str(tmp_path / "dbs")

    assert main(["db", "build", "--fasta", str(DATA / "reference.fasta"), "--name", "roundtrip",
                 "--db-location", db_location]) == 0  # fmt: skip
    assert bowtie.database_exists("roundtrip", db_location)

    assert main(["db", "list", "--db-location", db_location]) == 0
    assert "roundtrip" in capsys.readouterr().out

    assert main(["db", "remove", "roundtrip", "--db-location", db_location]) == 0
    assert not bowtie.database_exists("roundtrip", db_location)


def test_db_list_on_an_empty_location_exits_non_zero(tmp_path, capsys):
    assert main(["db", "list", "--db-location", str(tmp_path)]) == 1
    assert "no databases" in capsys.readouterr().out


def test_db_remove_reports_a_database_that_was_not_there(tmp_path, capsys):
    assert main(["db", "remove", "ghost", "--db-location", str(tmp_path)]) == 1
    assert "no such database" in capsys.readouterr().out
