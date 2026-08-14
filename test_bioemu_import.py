#!/usr/bin/env python3
"""Tests for bioemu_import.py.

Run with:  python -m pytest test_bioemu_import.py -q

The archive fixtures are built from real members pulled out of the published
Zenodo zips (a real GROMACS-written topology.pdb and real .cmprsd.xtc files), so
the header walking, atom-count checks and force-field splitting are exercised
against the actual byte layouts rather than synthetic stand-ins. Build them with
`python make_fixtures.py` — see FIXTURES below.
"""
from __future__ import annotations

import os
import re
import struct
import subprocess
import sys
import tomllib
import zipfile
from pathlib import Path

import pytest

import bioemu_import as bi

FIXTURES = Path(os.environ.get(
    "BIOEMU_FIXTURES",
    Path(__file__).with_name("fixtures"),
))

pytestmark = pytest.mark.skipif(
    not FIXTURES.is_dir(),
    reason=f"archive fixtures not found at {FIXTURES}; set BIOEMU_FIXTURES",
)


# --------------------------------------------------------------------------- #
# Identifier parsing
# --------------------------------------------------------------------------- #

def test_cath_system_ids():
    ids = bi.parse_system_ids(bi.DATASETS["cath1"], "cath1_1b43A02")
    assert (ids.pdb_id, ids.chain) == ("1b43", "A")
    ids = bi.parse_system_ids(bi.DATASETS["cath2"], "cath2_5kztB02")
    assert (ids.pdb_id, ids.chain) == ("5kzt", "B")


def test_cath_system_ids_reject_junk():
    with pytest.raises(ValueError):
        bi.parse_system_ids(bi.DATASETS["cath1"], "not_a_domain")


@pytest.mark.parametrize("system,pdb,mutation,variant,design", [
    ("1AOY", "1aoy", None, None, False),                 # bare wild-type
    ("1A0N_L7S", "1a0n", None, "L7S", False),            # wild-type variant
    ("2HBB_pross6", "2hbb", None, "pross6", False),      # PROSS redesign of a PDB entry
    ("1A0N_L7S__A12D", "1a0n", "A12D", "L7S", False),    # point mutant of the above
    ("EEHEE_rd3_0019", None, None, None, True),          # de novo design
    ("HEEH_KT_rd6_0007", None, None, None, True),        # 4 letters, but not a PDB code
    ("EA_run2_0290_0001", None, None, None, True),
])
def test_megasim_system_ids(system, pdb, mutation, variant, design):
    ids = bi.parse_system_ids(bi.DATASETS["megamerge"], system)
    assert ids.pdb_id == pdb
    assert ids.mutation == mutation
    assert ids.variant == variant
    assert ids.design is design


def test_megasim_mutant_records_its_parent():
    ids = bi.parse_system_ids(bi.DATASETS["megamut"], "1A0N_L7S__A12D")
    assert ids.parent == "1A0N_L7S"
    # A wild-type is nobody's mutant, so it has no parent to point at.
    assert bi.parse_system_ids(bi.DATASETS["megamerge"], "1A0N_L7S").parent is None


def test_octapeptides_have_no_identifiers():
    ids = bi.parse_system_ids(bi.DATASETS["opep"], "opep_0000")
    assert ids.pdb_id is None and ids.chain is None and not ids.design


# --------------------------------------------------------------------------- #
# Archive indexing / force-field grouping
# --------------------------------------------------------------------------- #

def _open(name: str):
    return zipfile.ZipFile(FIXTURES / name)


def test_index_archive_groups_members_by_system():
    ds = bi.DATASETS["cath1"]
    with _open("ONE_cath1.zip") as zf:
        index = bi.index_archive(zf, ds)
    assert set(index) == {"cath1_1b43A02", "cath1_9zzzA00"}
    sf = index["cath1_1b43A02"]
    assert sf.topology.endswith("topology.pdb")
    assert sf.dataset_json.endswith("dataset.json")
    assert len(sf.trajs) == 2
    assert sf.reference is None                      # CATH ships no reference.pdb


def test_index_archive_finds_reference_pdb_for_megasim():
    ds = bi.DATASETS["megamerge"]
    with _open("MSR_megasim_merge.zip") as zf:
        sf = bi.index_archive(zf, ds)["1AOY"]
    assert sf.reference.endswith("reference.pdb")


@pytest.mark.parametrize("member,expected", [
    ("trajs/run001_protein.cmprsd.xtc", "trajs/run001_protein.json"),
    ("trajs/aggr_trj_run0_clone0_folded_resim.cmprsd.xtc",
     "trajs/aggr_trj_run0_clone0_folded_resim.json"),
    ("trajs/trj_mutant_folded.xtc", "trajs/trj_mutant_folded.json"),
])
def test_sidecar_name(member, expected):
    assert bi.sidecar_name(member) == expected


def test_single_forcefield_system_is_one_group():
    ds = bi.DATASETS["cath1"]
    with _open("ONE_cath1.zip") as zf:
        groups = bi.discover_groups(zf, bi.index_archive(zf, ds)["cath1_1b43A02"])
    assert len(groups) == 1
    g = groups[0]
    assert g.slug == "ff99sb-ildn-300k"
    assert g.temperature_k == 300 and not g.overridden
    assert len(g.trajs) == 2


def test_sidecar_json_splits_a_system_into_two_groups():
    """The 77 MegaSim wild-types that mix ff14sb and a99SB-disp must not become
    one simulation: MDRepo records a single forcefield per simulation."""
    ds = bi.DATASETS["megamerge"]
    with _open("MSR_megasim_merge.zip") as zf:
        groups = bi.discover_groups(zf, bi.index_archive(zf, ds)["HEEH_KT_rd6_0007"])
    by_slug = {g.slug: g for g in groups}
    assert set(by_slug) == {"ff99sb-disp-295k", "ff14sb-295k"}
    # The overridden trajectory went to ff14sb, the inherited one to a99SB-disp.
    assert by_slug["ff14sb-295k"].overridden is True
    assert by_slug["ff99sb-disp-295k"].overridden is False
    assert [t.rsplit("/", 1)[-1] for t in by_slug["ff14sb-295k"].trajs] == \
        ["aggr_trj_run0_clone0_folded_resim.cmprsd.xtc"]


def test_system_without_sidecars_stays_single_group():
    ds = bi.DATASETS["megamerge"]
    with _open("MSR_megasim_merge.zip") as zf:
        groups = bi.discover_groups(zf, bi.index_archive(zf, ds)["1AOY"])
    assert [g.slug for g in groups] == ["ff99sb-disp-295k"]


def test_group_slug_is_stable_and_filename_safe():
    assert bi.group_slug("amber ff99sb-ildn", 300) == "ff99sb-ildn-300k"
    assert bi.group_slug("amber ff14sb", 295) == "ff14sb-295k"
    # An unknown force field still yields something usable as a path component.
    slug = bi.group_slug("Some/Weird FF v2.0", 310)
    assert "/" not in slug and " " not in slug and slug.endswith("-310k")


# --------------------------------------------------------------------------- #
# XTC header walking
# --------------------------------------------------------------------------- #

def _extract(zf, member, tmp_path) -> Path:
    dest = tmp_path / member.rsplit("/", 1)[-1]
    dest.write_bytes(zf.read(member))
    return dest


def _synth_xtc(path: Path, times, n_atoms: int = 3) -> Path:
    """A minimal but real xtc with the given frame times.

    Kept at <=9 atoms so coordinates take the uncompressed raw-float path,
    which lets the time series be dictated exactly without a compressor.
    """
    buf = bytearray()
    for step, t in enumerate(times):
        buf += struct.pack(">iiif", bi.XTC_MAGIC, n_atoms, step, t)
        buf += struct.pack(">9f", *([0.0] * 9))          # 3x3 box
        buf += struct.pack(">i", n_atoms)
        buf += struct.pack(f">{n_atoms * 3}f", *([0.0] * (n_atoms * 3)))
    path.write_bytes(bytes(buf))
    return path


def test_xtc_scan_reports_modal_not_mean_spacing(tmp_path):
    """An aggregated trajectory restarts its clock at each segment boundary.

    Every frame pair inside a segment is 10000 ps apart, so that is the file's
    real sampling interval; the mean over first..last would report 60 ps here,
    which is what the cath1 warnings were made of.
    """
    times = ([i * 10000.0 for i in range(300)]          # segment 1
             + [i * 10000.0 for i in range(198)]        # clock restarts
             + [i * 10000.0 for i in range(3)])         # short tail
    info = bi.xtc_scan(_synth_xtc(tmp_path / "aggregated.xtc", times))
    assert info.n_frames == 501
    assert info.ps_per_frame == 10000.0
    # The estimator this replaced; pinned so the regression can't come back.
    assert (times[-1] - times[0]) / (len(times) - 1) == pytest.approx(40.0)


def test_xtc_scan_ignores_duplicate_timestamps(tmp_path):
    times = [0.0, 10000.0, 10000.0, 20000.0, 30000.0, 30000.0, 40000.0]
    info = bi.xtc_scan(_synth_xtc(tmp_path / "dupes.xtc", times))
    assert info.ps_per_frame == 10000.0


def test_modal_dt_tolerates_float32_jitter():
    # float32 cannot hold 2 us exactly to the ps, so nominally equal gaps differ
    # slightly; they must still land in one bucket rather than fragmenting.
    gaps = [10000.0, 10000.25, 9999.75, 10000.125, 10000.0, 9999.875]
    assert bi.modal_dt(gaps) == pytest.approx(10000.0, abs=0.5)


def test_modal_dt_is_none_without_a_positive_gap():
    assert bi.modal_dt([]) is None                 # single-frame file
    assert bi.modal_dt([0.0, -5.0, 0.0]) is None   # nothing but restarts


def test_modal_dt_survives_a_minority_of_real_spacings():
    # A genuinely mis-declared file must still be caught: if most gaps are 500,
    # 500 is the answer, not the 10000 that a handful of frames happen to show.
    assert bi.modal_dt([500.0] * 20 + [10000.0] * 3) == 500.0


def test_scaled_stamp_is_recognised_only_at_exactly_1000x(tmp_path):
    """The octapeptide stamp, and nothing that merely resembles it.

    Every `.filtered.` file in ONE_octapeptides holds two frames at 100000 and
    10100000 ps against a declared 10000 ps. Widening this predicate would start
    swallowing real mis-declarations, so it is pinned to the exact shape.
    """
    stamp = bi.xtc_scan(_synth_xtc(tmp_path / "s.xtc", [100000.0, 10100000.0]))
    assert stamp.ps_per_frame == 10000000.0
    assert bi.is_scaled_stamp(stamp, 10000.0)

    # Two frames, but off by some other factor — a real disagreement.
    other = bi.xtc_scan(_synth_xtc(tmp_path / "o.xtc", [0.0, 500.0]))
    assert not bi.is_scaled_stamp(other, 10000.0)

    # The right factor but not two frames: an aggregate, not the stamp.
    many = bi.xtc_scan(
        _synth_xtc(tmp_path / "m.xtc", [i * 1e7 for i in range(6)]))
    assert not bi.is_scaled_stamp(many, 10000.0)

    # Single-frame files have no spacing at all.
    one = bi.xtc_scan(_synth_xtc(tmp_path / "1.xtc", [0.0]))
    assert one.ps_per_frame is None and not bi.is_scaled_stamp(one, 10000.0)


def test_real_octapeptide_filtered_file_carries_the_scaled_stamp(tmp_path):
    """Pinned against the published bytes, not a synthetic stand-in."""
    with _open("ONE_octapeptides.zip") as zf:
        p = _extract(zf, "ONE_octapeptides/opep_0000/trajs/"
                         "e10s1_e8s2p0f150-ADRIA_LARGEPEP_opep_0000-0-1-"
                         "RND0375_9.filtered.cmprsd.xtc", tmp_path)
        run = _extract(zf, "ONE_octapeptides/opep_0000/trajs/"
                           "run001_protein.cmprsd.xtc", tmp_path)
    info = bi.xtc_scan(p)
    assert info.n_frames == 2 and info.ps_per_frame == 10000000.0
    assert bi.is_scaled_stamp(info, 10000.0)
    # The runNNN files in the same system are the honest ones.
    assert bi.xtc_scan(run).ps_per_frame == 10000.0


def test_xtc_scan_reads_real_frame_headers(tmp_path):
    with _open("ONE_cath1.zip") as zf:
        p = _extract(zf, "ONE_cath1/cath1_1b43A02/trajs/run000_protein.cmprsd.xtc",
                     tmp_path)
    info = bi.xtc_scan(p)
    assert info.n_atoms == 1002
    assert info.n_frames == 201
    # The release stamps 10 ns per frame, matching dataset.json's save_traj_ns.
    assert info.ps_per_frame == pytest.approx(10000.0)


def test_xtc_scan_rejects_a_non_xtc(tmp_path):
    p = tmp_path / "garbage.xtc"
    p.write_bytes(b"not an xtc at all, no magic here....")
    with pytest.raises(bi.BadTrajectoryError, match="magic"):
        bi.xtc_scan(p)


def test_xtc_scan_rejects_a_truncated_file(tmp_path):
    with _open("ONE_cath1.zip") as zf:
        full = zf.read("ONE_cath1/cath1_1b43A02/trajs/run000_protein.cmprsd.xtc")
    p = tmp_path / "truncated.xtc"
    p.write_bytes(full[:5000])
    with pytest.raises(bi.BadTrajectoryError):
        bi.xtc_scan(p)


def test_xtc_scan_rejects_an_empty_file(tmp_path):
    p = tmp_path / "empty.xtc"
    p.write_bytes(b"")
    with pytest.raises(bi.BadTrajectoryError, match="no readable frames"):
        bi.xtc_scan(p)


def test_pdb_atom_count_matches_the_trajectory(tmp_path):
    with _open("ONE_cath1.zip") as zf:
        pdb = _extract(zf, "ONE_cath1/cath1_1b43A02/topology.pdb", tmp_path)
    assert bi.pdb_atom_count(pdb) == 1002


# --------------------------------------------------------------------------- #
# Metadata rendering
# --------------------------------------------------------------------------- #

def _render(ds_key: str, system: str, uniprot=(), **kw) -> dict:
    ds = bi.DATASETS[ds_key]
    group = bi.Group(slug="ff99sb-ildn-300k", force_field="amber ff99sb-ildn",
                     temperature_k=300, save_traj_ns=10.0,
                     trajs=["trajs/a.xtc", "trajs/b.xtc"])
    group.__dict__.update(kw.pop("group", {}))
    text = bi.render_metadata(
        ds, system, group, bi.parse_system_ids(ds, system), list(uniprot),
        traj_names=["a.xtc", "b.xtc"], pdb_name=f"{system}.pdb",
        psf_name=f"{system}.psf", n_atoms=1002, n_frames=251, **kw)
    return tomllib.loads(text)


def test_rendered_metadata_is_valid_toml_with_required_keys():
    meta = _render("cath1", "cath1_1b43A02", uniprot=["P12345"])
    for key in ("lead_contributor_orcid", "trajectory_file_names",
                "structure_file_name", "topology_file_name", "temperature_kelvin",
                "integration_timestep_fs", "short_description", "software_name",
                "software_version"):
        assert key in meta, f"missing required key {key}"
    assert meta["trajectory_file_names"] == ["a.xtc", "b.xtc"]
    assert meta["topology_file_name"].endswith(".psf")
    assert meta["structure_file_name"].endswith(".pdb")


def _orcid_checksum_ok(orcid: str) -> bool:
    """ISO 7064 MOD 11-2, the check-digit scheme ORCID uses."""
    digits = orcid.replace("-", "")
    total = 0
    for ch in digits[:15]:
        total = (total + int(ch)) * 2
    remainder = (12 - total % 11) % 11
    return ("X" if remainder == 10 else str(remainder)) == digits[15]


def test_contributors_are_the_paper_authors_in_order():
    # Credit in the deposition must match credit in the paper; a name that drifts
    # out of sync with the citation string is a provenance bug, not a typo.
    authors = [a.strip() for a in bi.PAPER_BIOEMU["authors"].split(",")]
    assert [c["name"] for c in bi.CONTRIBUTORS] == authors


def test_contributor_orcids_are_well_formed():
    with_orcid = [c for c in bi.CONTRIBUTORS if c.get("orcid")]
    assert len(with_orcid) == 8, "the 8 validated ORCIDs should all be present"
    for c in with_orcid:
        orcid = c["orcid"]
        assert re.fullmatch(r"\d{4}-\d{4}-\d{4}-\d{3}[\dX]", orcid), c["name"]
        assert _orcid_checksum_ok(orcid), f"{c['name']}: bad check digit {orcid}"
        assert orcid != bi.LEAD_CONTRIBUTOR_ORCID


def test_accented_contributor_names_survive_rendering():
    # The names are non-ASCII on purpose; this pins that they reach the TOML
    # intact rather than being stripped, mangled, or escaped into \u sequences.
    text = bi.render_metadata(
        bi.DATASETS["cath1"], "cath1_1b43A02",
        bi.Group(slug="ff99sb-ildn-300k", force_field="amber ff99sb-ildn",
                 temperature_k=300, save_traj_ns=10.0, trajs=["trajs/a.xtc"]),
        bi.parse_system_ids(bi.DATASETS["cath1"], "cath1_1b43A02"), [],
        traj_names=["a.xtc"], pdb_name="x.pdb", psf_name="x.psf",
        n_atoms=1002, n_frames=251)
    for name in ("José Jiménez-Luna", "Victor García Satorras", "Frank Noé",
                 "Freie Universität Berlin, Department of Physics"):
        assert name in text
    names = [c["name"] for c in tomllib.loads(text)["contributors"]]
    assert "Frank Noé" in names


def test_rendered_contributors_carry_names_and_orcids():
    meta = _render("cath1", "cath1_1b43A02")
    rendered = meta["contributors"]
    assert [c["name"] for c in rendered] == [c["name"] for c in bi.CONTRIBUTORS]
    for got, want in zip(rendered, bi.CONTRIBUTORS):
        assert got.get("orcid") == want.get("orcid"), want["name"]
        # Optional fields are omitted, never emitted empty or placeholdered.
        assert "orcid" in got or "orcid" not in want


def test_temperature_and_timestep_are_in_spec_range():
    meta = _render("cath1", "cath1_1b43A02")
    assert 275 <= meta["temperature_kelvin"] <= 700
    assert 1 <= meta["integration_timestep_fs"] <= 20


def test_short_description_respects_the_300_char_limit():
    meta = _render("megamut", "1A0N_L7S__A12D",
                   uniprot=[f"P{n:05d}" for n in range(40)])
    assert len(meta["short_description"]) <= bi.SHORT_DESCRIPTION_MAX


def test_truncate_short_only_trims_when_needed():
    assert bi.truncate_short("short text") == "short text"
    long = "word " * 200
    out = bi.truncate_short(long)
    assert len(out) <= bi.SHORT_DESCRIPTION_MAX
    assert out.endswith(".")


def test_software_is_the_closed_vocabulary_entry_not_openmm():
    """MDRepo's software_name vocabulary has no OpenMM; CUSTOM requires NA."""
    meta = _render("cath1", "cath1_1b43A02")
    assert meta["software_name"] == "CUSTOM"
    assert meta["software_version"] == "NA"
    assert "OpenMM" in meta["description"]


def test_cath_metadata_carries_pdb_and_uniprot():
    meta = _render("cath1", "cath1_1b43A02", uniprot=["P12345", "Q67890"])
    assert meta["pdb_id"] == "1b43"
    assert meta["uniprot_ids"] == ["P12345", "Q67890"]


def test_octapeptide_metadata_omits_absent_identifiers():
    meta = _render("opep", "opep_0000")
    assert "pdb_id" not in meta
    assert "uniprot_ids" not in meta


def test_mutant_metadata_names_the_mutation_and_parent():
    meta = _render("megamut", "1A0N_L7S__A12D")
    assert meta["pdb_id"] == "1a0n"
    assert "A12D" in meta["description"]
    assert "1A0N_L7S" in meta["description"]


def test_water_recorded_only_where_the_manuscript_states_it():
    """S.1.4 states TIP3P; S.1.5.4 restates neither water nor salt for MegaSim,
    so those tables must be absent rather than guessed."""
    standard = _render("cath1", "cath1_1b43A02")
    assert standard["water"]["model"] == "TIP3P"
    assert 900 <= standard["water"]["density_kg_m3"] <= 1100
    assert [s["name"] for s in standard["solutes"]] == ["Na+", "Cl-"]

    ds = bi.DATASETS["megamerge"]
    group = bi.Group(slug="ff99sb-disp-295k", force_field="amber ff99sb-disp",
                     temperature_k=295, save_traj_ns=10.0, trajs=["trajs/a.xtc"])
    meta = tomllib.loads(bi.render_metadata(
        ds, "1AOY", group, bi.parse_system_ids(ds, "1AOY"), [],
        traj_names=["a.xtc"], pdb_name="1AOY.pdb", psf_name="1AOY.psf",
        n_atoms=728, n_frames=100))
    assert "water" not in meta
    assert "solutes" not in meta


def test_megasim_cites_both_papers():
    ds = bi.DATASETS["megamerge"]
    group = bi.Group(slug="ff99sb-disp-295k", force_field="amber ff99sb-disp",
                     temperature_k=295, save_traj_ns=10.0, trajs=["trajs/a.xtc"])
    meta = tomllib.loads(bi.render_metadata(
        ds, "1AOY", group, bi.parse_system_ids(ds, "1AOY"), [],
        traj_names=["a.xtc"], pdb_name="p.pdb", psf_name="p.psf",
        n_atoms=728, n_frames=100))
    dois = {p["doi"] for p in meta["papers"]}
    assert dois == {"10.1126/science.adv9817", "10.1038/s41586-023-06328-6"}
    # CATH/opep have no MEGAscale ancestry, so they cite only BioEmu.
    assert len(_render("cath1", "cath1_1b43A02")["papers"]) == 1


def test_papers_carry_every_required_field():
    for paper in _render("megamut", "1A0N_L7S__A12D")["papers"]:
        for key in ("title", "authors", "journal", "year", "volume"):
            assert paper.get(key), f"{key} missing from {paper.get('title')}"
        assert isinstance(paper["year"], int)
        assert isinstance(paper["volume"], int)


def test_external_links_cite_source_model_and_processing_scripts():
    """Every entry gets three links: where the data came from, the model that
    motivated it, and the scripts that produced the deposited topology. The
    processing repo is cited bare in the description but must carry a scheme
    here, because the spec's url field is a URL rather than prose."""
    meta = _render("cath1", "cath1_1b43A02")
    by_url = {ln["url"]: ln["label"] for ln in meta["external_links"]}
    assert len(by_url) == 3
    assert all(u.startswith("https://") for u in by_url), by_url
    assert bi.PROCESSING_REPO_URL in by_url
    assert by_url[bi.PROCESSING_REPO_URL] == "MDRepo BioEmu processing scripts (GitHub)"
    assert bi.BIOEMU_REPO_URL in by_url
    # The Zenodo link is per-dataset, so it must track the record id.
    assert bi.ZENODO_DOI.format(record=bi.DATASETS["cath1"].record) in by_url
    assert bi.PROCESSING_REPO in meta["description"]


def test_split_group_says_so_in_the_description():
    ds = bi.DATASETS["megamerge"]
    group = bi.Group(slug="ff14sb-295k", force_field="amber ff14sb",
                     temperature_k=295, save_traj_ns=10.0,
                     trajs=["trajs/a.xtc"], overridden=True)
    meta = tomllib.loads(bi.render_metadata(
        ds, "HEEH_KT_rd6_0007", group,
        bi.parse_system_ids(ds, "HEEH_KT_rd6_0007"), [],
        traj_names=["a.xtc"], pdb_name="p.pdb", psf_name="p.psf",
        n_atoms=728, n_frames=100))
    assert meta["forcefield"] == "Amber ff14SB"
    assert "different force field" in meta["description"]


def test_toml_escaping_survives_quotes_and_backslashes():
    assert tomllib.loads("k = " + bi.toml_str('a "quoted" \\ path'))["k"] == \
        'a "quoted" \\ path'
    assert tomllib.loads("k = " + bi.toml_str("line\nbreak"))["k"] == "line\nbreak"
    # Control characters are dropped rather than emitted raw, which TOML forbids.
    assert "\x07" not in tomllib.loads("k = " + bi.toml_str("bell\x07here"))["k"]


def test_dropped_trajectories_are_disclosed_in_the_description():
    ds = bi.DATASETS["cath1"]
    group = bi.Group(slug="ff99sb-ildn-300k", force_field="amber ff99sb-ildn",
                     temperature_k=300, save_traj_ns=10.0, trajs=["trajs/a.xtc"])
    meta = tomllib.loads(bi.render_metadata(
        ds, "cath1_1b43A02", group,
        bi.parse_system_ids(ds, "cath1_1b43A02"), [], traj_names=["a.xtc"],
        pdb_name="p.pdb", psf_name="p.psf", n_atoms=1002, n_frames=201,
        dropped=[("bad.xtc", "728 atoms, topology has 1002")]))
    assert "excluded" in meta["description"] and "bad.xtc" in meta["description"]


# --------------------------------------------------------------------------- #
# Building an IN_DIR (requires parmed for the .psf)
# --------------------------------------------------------------------------- #

parmed = pytest.importorskip("parmed", reason="parmed is needed to write the .psf")


def _build(zip_name: str, ds_key: str, system: str, out: Path, which=0):
    ds = bi.DATASETS[ds_key]
    with _open(zip_name) as zf:
        sf = bi.index_archive(zf, ds)[system]
        group = bi.discover_groups(zf, sf)[which]
        ids = bi.parse_system_ids(ds, system)
        return bi.build_sim_dir(zf, ds, sf, group, ids, [], out / group.slug), group


def test_build_sim_dir_produces_a_complete_in_dir(tmp_path):
    built, _ = _build("ONE_cath1.zip", "cath1", "cath1_1b43A02", tmp_path)
    d = built.sim_dir
    names = {p.name for p in d.iterdir()}
    assert "mdrepo-metadata.toml" in names
    assert "cath1_1b43A02.pdb" in names
    assert "cath1_1b43A02.psf" in names          # generated, not shipped
    assert built.n_trajs == 2 and built.n_frames == 251

    meta = tomllib.loads((d / "mdrepo-metadata.toml").read_text(encoding="utf-8"))
    # Every file the metadata names must actually be in the directory.
    for key in ("structure_file_name", "topology_file_name"):
        assert (d / meta[key]).is_file()
    for traj in meta["trajectory_file_names"]:
        assert (d / traj).is_file()


def test_trajectories_are_deposited_byte_for_byte(tmp_path):
    """Nothing is re-encoded, so the published bytes must survive untouched."""
    built, _ = _build("ONE_cath1.zip", "cath1", "cath1_1b43A02", tmp_path)
    with _open("ONE_cath1.zip") as zf:
        original = zf.read(
            "ONE_cath1/cath1_1b43A02/trajs/run000_protein.cmprsd.xtc")
    assert (built.sim_dir / "run000_protein.cmprsd.xtc").read_bytes() == original


def test_generated_psf_atom_count_matches_the_pdb(tmp_path):
    built, _ = _build("ONE_cath1.zip", "cath1", "cath1_1b43A02", tmp_path)
    pdb_atoms = bi.pdb_atom_count(built.sim_dir / "cath1_1b43A02.pdb")
    psf_text = (built.sim_dir / "cath1_1b43A02.psf").read_text()
    assert f"{pdb_atoms:8d} !NATOM" in psf_text or f"{pdb_atoms} !NATOM" in psf_text


def test_reference_pdb_is_deposited_as_an_additional_file(tmp_path):
    built, _ = _build("MSR_megasim_merge.zip", "megamerge", "1AOY", tmp_path)
    meta = tomllib.loads((built.sim_dir / "mdrepo-metadata.toml").read_text(encoding="utf-8"))
    extra = meta["additional_files"]
    assert [f["file_name"] for f in extra] == ["1AOY_reference.pdb"]
    assert (built.sim_dir / "1AOY_reference.pdb").is_file()


def test_cath_in_dir_has_no_additional_files(tmp_path):
    built, _ = _build("ONE_cath1.zip", "cath1", "cath1_1b43A02", tmp_path)
    meta = tomllib.loads((built.sim_dir / "mdrepo-metadata.toml").read_text(encoding="utf-8"))
    assert "additional_files" not in meta


def test_metadata_is_written_utf8_under_an_ascii_locale(tmp_path):
    """The real write path, exercised under LANG=C.

    Contributor names carry diacritics, so a `write_text` that inherits the
    locale encoding raises UnicodeEncodeError on a POSIX-locale VM — which is
    where this importer runs. Forked rather than monkeypatched because the
    locale encoding is resolved below Python.
    """
    out = tmp_path / "out"
    script = (
        "import sys, zipfile, pathlib\n"
        "sys.path.insert(0, sys.argv[1])\n"
        "import bioemu_import as bi\n"
        "ds = bi.DATASETS['cath1']\n"
        "with zipfile.ZipFile(sys.argv[2]) as zf:\n"
        "    sf = bi.index_archive(zf, ds)['cath1_1b43A02']\n"
        "    g = bi.discover_groups(zf, sf)[0]\n"
        "    ids = bi.parse_system_ids(ds, 'cath1_1b43A02')\n"
        "    bi.build_sim_dir(zf, ds, sf, g, ids, [], pathlib.Path(sys.argv[3]))\n"
    )
    env = {**os.environ, "LANG": "C", "LC_ALL": "C", "PYTHONUTF8": "0",
           "PYTHONCOERCECLOCALE": "0"}
    proc = subprocess.run(
        [sys.executable, "-c", script, str(Path(bi.__file__).parent),
         str(FIXTURES / "ONE_cath1.zip"), str(out)],
        env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    raw = (out / "mdrepo-metadata.toml").read_bytes()
    assert "Frank Noé".encode("utf-8") in raw
    names = [c["name"] for c in tomllib.loads(raw.decode("utf-8"))["contributors"]]
    assert names == [c["name"] for c in bi.CONTRIBUTORS]


def test_each_forcefield_group_becomes_its_own_in_dir(tmp_path):
    ds = bi.DATASETS["megamerge"]
    with _open("MSR_megasim_merge.zip") as zf:
        sf = bi.index_archive(zf, ds)["HEEH_KT_rd6_0007"]
        groups = bi.discover_groups(zf, sf)
        ids = bi.parse_system_ids(ds, "HEEH_KT_rd6_0007")
        metas = {}
        for g in groups:
            built = bi.build_sim_dir(zf, ds, sf, g, ids, [], tmp_path / g.slug)
            metas[g.slug] = tomllib.loads(
                (built.sim_dir / "mdrepo-metadata.toml").read_text(encoding="utf-8"))
    assert metas["ff14sb-295k"]["forcefield"] == "Amber ff14SB"
    assert metas["ff99sb-disp-295k"]["forcefield"] == "Amber a99SB-disp"
    # No trajectory may appear in both simulations.
    a = set(metas["ff14sb-295k"]["trajectory_file_names"])
    b = set(metas["ff99sb-disp-295k"]["trajectory_file_names"])
    assert a and b and not (a & b)


def test_unreadable_and_mismatched_trajectories_are_dropped_not_fatal(tmp_path):
    notes = []
    ds = bi.DATASETS["cath1"]
    with _open("broken.zip") as zf:
        sf = bi.index_archive(zf, ds)["cath1_1b43A02"]
        group = bi.discover_groups(zf, sf)[0]
        built = bi.build_sim_dir(
            zf, ds, sf, group, bi.parse_system_ids(ds, "cath1_1b43A02"), [],
            tmp_path / "sim", on_note=lambda k, d: notes.append((k, d)))

    assert built.n_trajs == 1                       # only good.xtc survives
    assert {name for name, _why in built.dropped} == \
        {"garbage.xtc", "truncated.xtc", "wrong_atoms.xtc"}
    # Dropped files must not be left behind to confuse mdr-process.
    for name, _ in built.dropped:
        assert not (built.sim_dir / name).exists()
    kinds = {k for k, _ in notes}
    assert "unreadable-trajectory" in kinds and "atom-count-mismatch" in kinds

    meta = tomllib.loads((built.sim_dir / "mdrepo-metadata.toml").read_text(encoding="utf-8"))
    assert meta["trajectory_file_names"] == ["good.xtc"]


def test_scaled_stamp_is_summarised_once_not_noted_per_file(tmp_path):
    """One note per system, not one per trajectory.

    ~102 of every octapeptide system's trajectories carry the stamp — 112,756
    across the archive. A note each would bury every other observation in the
    manifest, so they collapse into a single `scaled-frame-stamp` row and no
    `frame-spacing-mismatch` is raised for them.
    """
    notes = []
    ds = bi.DATASETS["opep"]
    with _open("ONE_octapeptides.zip") as zf:
        sf = bi.index_archive(zf, ds)["opep_0000"]
        group = bi.discover_groups(zf, sf)[0]
        built = bi.build_sim_dir(
            zf, ds, sf, group, bi.parse_system_ids(ds, "opep_0000"), [],
            tmp_path / "sim", on_note=lambda k, d: notes.append((k, d)))

    assert built.n_trajs == 2                       # nothing is dropped for it
    assert built.n_frames == 203                    # 2 stamped + 201 honest
    kinds = [k for k, _ in notes]
    assert kinds == ["scaled-frame-stamp"]
    assert "1 of 2 trajectories" in notes[0][1]
    # The declared spacing is what reaches the deposition, not the stamp.
    meta = tomllib.loads(
        (built.sim_dir / "mdrepo-metadata.toml").read_text(encoding="utf-8"))
    assert "10 ns apart" in meta["description"]
    assert "10000000" not in meta["description"]


def test_a_two_frame_file_off_by_another_factor_is_still_reported(tmp_path):
    """The summary must not become a blanket amnesty for two-frame files."""
    notes = []
    ds = bi.DATASETS["opep"]
    doctored = tmp_path / "odd.zip"
    base = "ONE_octapeptides/opep_0000"
    with zipfile.ZipFile(FIXTURES / "ONE_octapeptides.zip") as z:
        filtered = next(n for n in z.namelist() if n.endswith(".filtered.cmprsd.xtc"))
        # The real stamped file with frame 0's time moved up, so the pair sits
        # 500 ps apart instead of 1000x the declared spacing. Everything else is
        # the published bytes, compressed coordinate blocks included.
        odd = bytearray(z.read(filtered))
        struct.pack_into(">f", odd, 12, 10099500.0)
        with zipfile.ZipFile(doctored, "w") as dst:
            for n in z.namelist():
                dst.writestr(n, bytes(odd) if n == filtered else z.read(n))
    with zipfile.ZipFile(doctored) as zf:
        sf = bi.index_archive(zf, ds)["opep_0000"]
        group = bi.discover_groups(zf, sf)[0]
        bi.build_sim_dir(zf, ds, sf, group,
                         bi.parse_system_ids(ds, "opep_0000"), [],
                         tmp_path / "sim",
                         on_note=lambda k, d: notes.append((k, d)))
    assert [k for k, _ in notes] == ["frame-spacing-mismatch"]
    assert "500 ps apart" in notes[0][1] and ".filtered." in notes[0][1]


def test_system_with_no_usable_trajectories_raises_skip(tmp_path):
    ds = bi.DATASETS["cath1"]
    broken = tmp_path / "all_bad.zip"
    with zipfile.ZipFile(FIXTURES / "ONE_cath1.zip") as src, \
            zipfile.ZipFile(broken, "w") as dst:
        base = "ONE_cath1/cath1_1b43A02"
        dst.writestr(f"{base}/dataset.json", src.read(f"{base}/dataset.json"))
        dst.writestr(f"{base}/topology.pdb", src.read(f"{base}/topology.pdb"))
        dst.writestr(f"{base}/trajs/bad.xtc", b"junk")
    with zipfile.ZipFile(broken) as zf:
        sf = bi.index_archive(zf, ds)["cath1_1b43A02"]
        group = bi.discover_groups(zf, sf)[0]
        with pytest.raises(bi.SkipSystem, match="no usable trajectories"):
            bi.build_sim_dir(zf, ds, sf, group,
                             bi.parse_system_ids(ds, "cath1_1b43A02"), [],
                             tmp_path / "sim")


def test_fresh_dir_refuses_to_wipe_protected_directories(tmp_path):
    for victim in (Path.cwd(), Path.home(), Path("/")):
        with pytest.raises(ValueError, match="refusing"):
            bi._fresh_dir(victim)


def test_fresh_dir_clears_stale_content(tmp_path):
    d = tmp_path / "sim"
    d.mkdir()
    (d / "stale.xtc").write_text("from a previous attempt")
    bi._fresh_dir(d)
    assert list(d.iterdir()) == []


# --------------------------------------------------------------------------- #
# Manifest
# --------------------------------------------------------------------------- #

@pytest.fixture()
def man(tmp_path):
    m = bi.Manifest(tmp_path / "manifest.sqlite")
    yield m
    m.close()


def test_add_systems_is_idempotent(man):
    assert man.add_systems("cath1", [("a", 2), ("b", 0)]) == 2
    man.add_systems("cath1", [("a", 2), ("c", 1)])
    rows = {r["system"] for r in man.conn.execute(
        "SELECT system FROM systems WHERE dataset='cath1'")}
    assert rows == {"a", "b", "c"}


def test_candidates_exclude_finished_and_exhausted_systems(man):
    man.add_systems("cath1", [("a", 1), ("b", 1), ("c", 1), ("d", 1)])
    man.mark_done("cath1", "a")
    man.mark_skipped("cath1", "b", "no trajectories")
    for _ in range(bi.MAX_ATTEMPTS):
        man.mark_failed("cath1", "c", "boom")
    assert man.candidate_systems("cath1") == ["d"]
    assert man.pending_count("cath1") == 1


def test_transient_requeue_preserves_the_attempt_budget(man):
    man.add_systems("cath1", [("a", 1)])
    man.mark_failed("cath1", "a", "boom")
    man.requeue("cath1", "a")
    row = man.system_row("cath1", "a")
    assert row["status"] == "pending" and row["attempts"] == 1
    man.reset_system("cath1", "a")
    assert man.system_row("cath1", "a")["attempts"] == 0


def test_import_records_are_per_group(man):
    man.add_systems("megamerge", [("HEEH_KT_rd6_0007", 3)])
    man.mark_imported("megamerge", "HEEH_KT_rd6_0007", "ff14sb-295k", "l")
    assert man.imported_groups("megamerge", "HEEH_KT_rd6_0007") == {"ff14sb-295k"}
    assert man.importing_groups("megamerge", "HEEH_KT_rd6_0007") == set()


def test_interrupted_push_is_left_ambiguous_until_resolved(man):
    man.add_systems("cath1", [("a", 1)])
    man.mark_importing("cath1", "a", "ff99sb-ildn-300k", "log")
    assert [dict(r) for r in man.ambiguous()][0]["grp"] == "ff99sb-ildn-300k"
    assert man.retry_import("cath1", "a", "ff99sb-ildn-300k") is True
    assert man.ambiguous() == []
    # Resolving the same record twice is a no-op, not an error.
    assert man.retry_import("cath1", "a", "ff99sb-ildn-300k") is False


def test_marking_imported_clears_the_ambiguous_state(man):
    man.add_systems("cath1", [("a", 1)])
    man.mark_importing("cath1", "a", "g", "log")
    man.mark_imported("cath1", "a", "g", "log")
    assert man.ambiguous() == []
    assert man.imported_groups("cath1", "a") == {"g"}


def test_uniprot_round_trip(man):
    man.set_uniprot([("cath1", "x", "P2"), ("cath1", "x", "P1"),
                     ("cath1", "y", "Q9")])
    assert man.uniprot_for("cath1", "x") == ["P1", "P2"]
    assert man.uniprot_for("cath1", "missing") == []


def test_reset_failed_scopes_to_a_dataset(man):
    man.add_systems("cath1", [("a", 1)])
    man.add_systems("opep", [("b", 1)])
    man.mark_failed("cath1", "a", "boom")
    man.mark_failed("opep", "b", "boom")
    assert man.reset_failed("cath1") == 1
    assert man.system_row("cath1", "a")["status"] == "pending"
    assert man.system_row("opep", "b")["status"] == "failed"


def test_notes_are_upserted_not_duplicated(man):
    man.note("cath1", "a", "atom-count-mismatch", "first")
    man.note("cath1", "a", "atom-count-mismatch", "second")
    rows = man.note_rows()
    assert len(rows) == 1 and rows[0]["detail"] == "second"


# --------------------------------------------------------------------------- #
# Failure classification
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("text", [
    "No space left on device", "connection reset by peer",
    "irods.exception.NetworkException: Could not connect to specified host",
    "HTTP 503 Service Unavailable", "Temporary failure in name resolution",
    "http.client.IncompleteRead: IncompleteRead(1019344 bytes read)",
])
def test_infrastructure_failures_are_transient(text):
    assert bi._is_transient_text(text)


@pytest.mark.parametrize("text", [
    "validate failed: topology_file_name has an unsupported extension",
    "psf has 500 atoms, pdb has 1002",
    "short_description exceeds 300 characters",
])
def test_input_failures_are_not_transient(text):
    assert not bi._is_transient_text(text)


def test_transient_exception_types():
    assert bi._is_transient_exc(bi.TransientError("x"))
    assert bi._is_transient_exc(ConnectionResetError("x"))
    assert bi._is_transient_exc(OSError(28, "No space left on device"))
    assert not bi._is_transient_exc(ValueError("bad domain id"))


def test_stalled_command_counts_as_transient():
    """A watchdog kill means the environment wedged, not that the input is bad —
    it must not spend one of the system's three attempts."""
    assert issubclass(bi.StalledCommandError, bi.TransientError)
    assert bi._is_transient_exc(bi.StalledCommandError("killed"))


# --------------------------------------------------------------------------- #
# Dataset table sanity
# --------------------------------------------------------------------------- #

def test_every_dataset_is_internally_consistent():
    for key, ds in bi.DATASETS.items():
        assert ds.key == key
        assert ds.protocol in bi.PROTOCOL_SUMMARY
        assert ds.id_style in ("cath", "megasim", "none")
        assert ds.ensemble in ("NPT", "NVT")
        assert ds.zip_name.endswith(".zip")


def test_dataset_order_is_smallest_archive_first():
    """`run` holds one archive at a time, so the order is the disk profile."""
    assert bi.DATASET_ORDER == ["opep", "megamerge", "cath1", "megamut", "cath2"]


def test_forcefields_cover_every_string_the_release_uses():
    assert set(bi.FORCEFIELDS) == {
        "amber ff99sb-ildn", "amber ff14sb", "amber ff99sb-disp"}


# --------------------------------------------------------------------------- #
# End-to-end: init + run through the real worker loop, with a stub mdr-process
# --------------------------------------------------------------------------- #

STUB_MDR = """#!/usr/bin/env python3
import sys, pathlib
log = pathlib.Path(__import__("os").environ["STUB_MDR_LOG"])
argv = sys.argv[1:]
with log.open("a") as fh:
    fh.write(" ".join(argv) + "\\n")
stage = next((a for a in argv if a in ("validate", "process")), "")
sim_dir = pathlib.Path(argv[argv.index(stage) + 1]) if stage else None
if stage == "validate" and sim_dir is not None:
    # A real `validate` reads the TOML and checks every file it names exists.
    import tomllib
    meta = tomllib.loads((sim_dir / "mdrepo-metadata.toml").read_bytes().decode())
    for key in ("structure_file_name", "topology_file_name"):
        assert (sim_dir / meta[key]).is_file(), meta[key]
    for t in meta["trajectory_file_names"]:
        assert (sim_dir / t).is_file(), t
if __import__("os").environ.get("STUB_MDR_FAIL") == stage:
    sys.stderr.write("stub failure in %s\\n" % stage)
    sys.exit(3)
sys.exit(0)
"""


@pytest.fixture()
def stub_mdr(tmp_path, monkeypatch):
    path = tmp_path / "stub-mdr-process"
    path.write_text(STUB_MDR)
    path.chmod(0o755)
    monkeypatch.setenv("STUB_MDR_LOG", str(tmp_path / "mdr-calls.log"))
    monkeypatch.delenv("STUB_MDR_FAIL", raising=False)
    return path


@pytest.fixture()
def root(tmp_path):
    """A --root with the fixture archives pre-staged, so nothing touches Zenodo."""
    lay = bi.Layout(tmp_path / "root")
    lay.ensure()
    for ds in bi.DATASETS.values():
        src = FIXTURES / ds.zip_name
        if src.exists():
            (lay.archives / ds.zip_name).write_bytes(src.read_bytes())
    return lay.root


def _cli(*argv):
    args = bi.build_parser().parse_args(list(argv))
    args.func(args)


def _init(root, *datasets):
    _cli("--root", str(root), "init", "--skip-sifts", "--datasets", ",".join(datasets))


def _run(root, *extra):
    _cli("--root", str(root), "run", "--no-verify", "--keep-archives",
         "-w", "1", *extra)


def test_init_builds_the_work_list_from_local_archives(root):
    _init(root, "cath1", "megamerge")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert set(man.candidate_systems("cath1")) == {
            "cath1_1b43A02", "cath1_9zzzA00"}
        assert man.pending_count("megamerge") == 3
        # A system the release ships without trajectories is still enrolled; it
        # is classified when a worker reaches it, not silently dropped here.
        row = man.system_row("cath1", "cath1_9zzzA00")
        assert row["n_trajs"] == 0 and row["status"] == "pending"
    finally:
        man.close()


def test_init_is_idempotent(root):
    _init(root, "cath1")
    _init(root, "cath1")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert man.pending_count("cath1") == 2
    finally:
        man.close()


def test_end_to_end_run_imports_every_system(root, stub_mdr):
    _init(root, "cath1", "megamerge", "megamut", "opep")
    _run(root, "--mdr-bin", str(stub_mdr))

    man = bi.Manifest(bi.Layout(root).db)
    try:
        statuses = {(r["dataset"], r["system"]): r["status"] for r in
                    man.conn.execute("SELECT dataset, system, status FROM systems")}
        # The trajectory-less system is skipped, not failed: it is a fact about
        # the release, not an error to retry.
        assert statuses[("cath1", "cath1_9zzzA00")] == "skipped"
        assert statuses[("cath1", "cath1_1b43A02")] == "done"
        assert all(s in ("done", "skipped") for s in statuses.values()), statuses

        sims = {(r["dataset"], r["system"], r["grp"]) for r in
                man.conn.execute("SELECT dataset, system, grp FROM sims "
                                 "WHERE state='imported'")}
        # The mixed-force-field wild-type became two simulations, everything
        # else exactly one.
        assert ("megamerge", "HEEH_KT_rd6_0007", "ff14sb-295k") in sims
        assert ("megamerge", "HEEH_KT_rd6_0007", "ff99sb-disp-295k") in sims
        assert len([s for s in sims if s[1] == "1AOY"]) == 1
        assert not man.ambiguous()
    finally:
        man.close()


def test_run_calls_validate_before_process_for_every_simulation(root, stub_mdr):
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr))
    calls = Path(os.environ["STUB_MDR_LOG"]).read_text().splitlines()
    stages = [c.split()[c.split().index("-l") + 2] for c in calls]
    assert stages == ["validate", "process"]
    assert "-s staging" in calls[1] and "-f" in calls[1].split()


def test_octapeptides_are_pushed_with_no_id(root, stub_mdr):
    """Octapeptides have neither a PDB nor a UniProt id, so mdr-process must be
    told not to look for one; CATH domains must NOT get that flag."""
    _init(root, "opep", "cath1")
    _run(root, "--mdr-bin", str(stub_mdr))
    calls = Path(os.environ["STUB_MDR_LOG"]).read_text().splitlines()
    process = [c for c in calls if " process " in c]
    opep = [c for c in process if "opep_0000" in c]
    cath = [c for c in process if "cath1_1b43A02" in c]
    assert opep and all("--no-id" in c for c in opep)
    assert cath and not any("--no-id" in c for c in cath)


def test_staging_is_reclaimed_after_a_successful_import(root, stub_mdr):
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr))
    staging = bi.Layout(root).staging
    leftover = [p for p in staging.rglob("*") if p.is_file()]
    assert leftover == [], leftover


def test_keep_retains_the_staged_in_dir(root, stub_mdr):
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--keep")
    sim = bi.Layout(root).sim_dir("cath1", "cath1_1b43A02", "ff99sb-ildn-300k")
    assert (sim / "mdrepo-metadata.toml").is_file()


def test_rerunning_a_completed_dataset_pushes_nothing_again(root, stub_mdr):
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr))
    first = Path(os.environ["STUB_MDR_LOG"]).read_text()
    _run(root, "--mdr-bin", str(stub_mdr))
    assert Path(os.environ["STUB_MDR_LOG"]).read_text() == first


def test_dry_run_leaves_the_system_pending_and_pushes_nothing(root, stub_mdr):
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--dry-run")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert man.system_row("cath1", "cath1_1b43A02")["status"] == "pending"
        assert man.conn.execute("SELECT COUNT(*) c FROM sims").fetchone()["c"] == 0
    finally:
        man.close()
    assert "-d" in Path(os.environ["STUB_MDR_LOG"]).read_text().split()


def test_a_failing_push_stops_at_max_attempts(root, stub_mdr, monkeypatch):
    monkeypatch.setenv("STUB_MDR_FAIL", "process")
    _init(root, "cath1")
    for _ in range(bi.MAX_ATTEMPTS + 2):
        _run(root, "--mdr-bin", str(stub_mdr), "--breaker-threshold", "0")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        row = man.system_row("cath1", "cath1_1b43A02")
        assert row["status"] == "failed"
        assert row["attempts"] == bi.MAX_ATTEMPTS      # budget respected, not exceeded
        assert man.pending_count("cath1") == 0
    finally:
        man.close()


def test_a_failed_push_is_left_ambiguous_rather_than_resubmitted(
        root, stub_mdr, monkeypatch):
    """`process` exiting non-zero does not prove nothing landed — it can fail
    part-way through the push. The system is therefore parked as ambiguous and
    never resubmitted automatically, even by reset-failed."""
    monkeypatch.setenv("STUB_MDR_FAIL", "process")
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--breaker-threshold", "0")

    man = bi.Manifest(bi.Layout(root).db)
    try:
        amb = man.ambiguous()
        assert [(r["system"], r["grp"]) for r in amb] == \
            [("cath1_1b43A02", "ff99sb-ildn-300k")]
    finally:
        man.close()

    # A working mdr-process plus reset-failed is deliberately NOT enough: the
    # operator has to say what happened in MDRepo first.
    monkeypatch.delenv("STUB_MDR_FAIL")
    _cli("--root", str(root), "reset-failed", "--datasets", "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--breaker-threshold", "0")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        row = man.system_row("cath1", "cath1_1b43A02")
        assert row["status"] == "failed"
        assert "resolve-import" in row["error"]
        assert man.ambiguous()
    finally:
        man.close()


def test_resolve_import_retry_unblocks_a_failed_push(root, stub_mdr, monkeypatch):
    monkeypatch.setenv("STUB_MDR_FAIL", "process")
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--breaker-threshold", "0")
    monkeypatch.delenv("STUB_MDR_FAIL")

    _cli("--root", str(root), "resolve-import", "cath1", "cath1_1b43A02",
         "ff99sb-ildn-300k", "--retry")
    _run(root, "--mdr-bin", str(stub_mdr))
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert man.system_row("cath1", "cath1_1b43A02")["status"] == "done"
        assert man.imported_groups("cath1", "cath1_1b43A02") == {"ff99sb-ildn-300k"}
        assert not man.ambiguous()
    finally:
        man.close()


def test_a_validate_failure_never_reaches_process_or_becomes_ambiguous(
        root, stub_mdr, monkeypatch):
    """validate runs before the import intent is recorded, so a metadata problem
    is an ordinary retryable failure — no manual resolution required."""
    monkeypatch.setenv("STUB_MDR_FAIL", "validate")
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--breaker-threshold", "0")
    calls = Path(os.environ["STUB_MDR_LOG"]).read_text()
    assert "validate" in calls and " process " not in calls

    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert not man.ambiguous()
        assert man.system_row("cath1", "cath1_1b43A02")["status"] == "failed"
    finally:
        man.close()

    monkeypatch.delenv("STUB_MDR_FAIL")
    _cli("--root", str(root), "reset-failed", "--datasets", "cath1")
    _run(root, "--mdr-bin", str(stub_mdr))
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert man.system_row("cath1", "cath1_1b43A02")["status"] == "done"
    finally:
        man.close()


def test_archive_is_deleted_once_its_systems_are_drained(root, stub_mdr):
    _init(root, "cath1")
    _cli("--root", str(root), "run", "--no-verify", "-w", "1",
         "--datasets", "cath1", "--mdr-bin", str(stub_mdr))
    assert not bi.Layout(root).archive(bi.DATASETS["cath1"]).exists()


def test_keep_archives_retains_the_zip(root, stub_mdr):
    _init(root, "cath1")
    _run(root, "--mdr-bin", str(stub_mdr), "--datasets", "cath1")
    assert bi.Layout(root).archive(bi.DATASETS["cath1"]).exists()


def test_scoped_run_touches_only_the_named_system(root, stub_mdr):
    _init(root, "megamerge")
    _run(root, "--mdr-bin", str(stub_mdr), "--datasets", "megamerge",
         "--systems", "1AOY")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert man.system_row("megamerge", "1AOY")["status"] == "done"
        assert man.system_row("megamerge", "1A0N_L7S")["status"] == "pending"
    finally:
        man.close()


def test_limit_caps_how_many_systems_are_attempted(root, stub_mdr):
    _init(root, "megamerge")
    _run(root, "--mdr-bin", str(stub_mdr), "--datasets", "megamerge", "--limit", "1")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        settled = man.conn.execute(
            "SELECT COUNT(*) c FROM systems WHERE dataset='megamerge' "
            "AND status != 'pending'").fetchone()["c"]
        assert settled == 1
    finally:
        man.close()


def test_extract_stages_a_system_without_importing(root, tmp_path):
    _init(root, "megamerge")
    out = tmp_path / "inspect"
    _cli("--root", str(root), "extract", "megamerge", "HEEH_KT_rd6_0007",
         "--out-dir", str(out), "--no-verify")
    dirs = sorted(p.name for p in (out / "megamerge" / "HEEH_KT_rd6_0007").iterdir())
    assert dirs == ["ff14sb-295k", "ff99sb-disp-295k"]
    man = bi.Manifest(bi.Layout(root).db)
    try:                                   # nothing was recorded as imported
        assert man.conn.execute("SELECT COUNT(*) c FROM sims").fetchone()["c"] == 0
    finally:
        man.close()


def test_resolve_import_clears_an_ambiguous_record(root):
    _init(root, "cath1")
    man = bi.Manifest(bi.Layout(root).db)
    man.mark_importing("cath1", "cath1_1b43A02", "ff99sb-ildn-300k", "log")
    man.close()

    _cli("--root", str(root), "resolve-import", "cath1", "cath1_1b43A02",
         "ff99sb-ildn-300k", "--imported")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        assert not man.ambiguous()
        assert man.imported_groups("cath1", "cath1_1b43A02") == {"ff99sb-ildn-300k"}
    finally:
        man.close()


def test_an_ambiguous_record_blocks_the_system_until_resolved(root, stub_mdr):
    _init(root, "cath1")
    man = bi.Manifest(bi.Layout(root).db)
    man.mark_importing("cath1", "cath1_1b43A02", "some-other-group", "log")
    man.close()
    _run(root, "--mdr-bin", str(stub_mdr), "--breaker-threshold", "0")
    man = bi.Manifest(bi.Layout(root).db)
    try:
        row = man.system_row("cath1", "cath1_1b43A02")
        assert row["status"] == "failed"
        assert "resolve-import" in row["error"]
        # Crucially: mdr-process was never invoked for that system while the
        # outcome of the earlier push was unknown.
        calls = Path(os.environ["STUB_MDR_LOG"])
        assert not calls.exists() or "cath1_1b43A02" not in calls.read_text()
    finally:
        man.close()


def test_unknown_dataset_name_is_rejected():
    with pytest.raises(SystemExit):
        _cli("--root", "/tmp/whatever", "init", "--datasets", "nope")


def test_several_workers_import_each_system_exactly_once(root, stub_mdr):
    """The per-system flock is the mutex; a double push would be a duplicate
    MDRepo record, which nothing downstream can undo."""
    _init(root, "megamerge", "cath1", "opep")
    _cli("--root", str(root), "run", "--no-verify", "--keep-archives",
         "-w", "4", "--mdr-bin", str(stub_mdr))

    calls = Path(os.environ["STUB_MDR_LOG"]).read_text().splitlines()
    sim_dirs = [c.split()[c.split().index("process") + 1] for c in calls
                if " process " in c]
    assert len(sim_dirs) == len(set(sim_dirs)), "a simulation was pushed twice"

    man = bi.Manifest(bi.Layout(root).db)
    try:
        rows = man.conn.execute(
            "SELECT dataset, system, status FROM systems").fetchall()
        assert all(r["status"] in ("done", "skipped") for r in rows)
        # 3 megamerge systems (one split in two) + 1 cath1 + 1 opep = 6 sims.
        assert man.conn.execute(
            "SELECT COUNT(*) c FROM sims WHERE state='imported'").fetchone()["c"] == 6
    finally:
        man.close()


def test_flock_self_test_passes_on_this_filesystem(tmp_path):
    """The guard that gates -w > 1. If this fails here, --root is on a
    filesystem where flock does not exclude and multi-worker runs are unsafe."""
    lay = bi.Layout(tmp_path / "root")
    lay.ensure()
    bi.verify_domain_locking(lay)          # raises RuntimeError if flock no-ops


def test_archive_is_kept_while_systems_remain_failed(root, stub_mdr, monkeypatch):
    """reset-failed must not force a 28 GB re-download."""
    monkeypatch.setenv("STUB_MDR_FAIL", "validate")
    _init(root, "cath1")
    _cli("--root", str(root), "run", "--no-verify", "-w", "1",
         "--datasets", "cath1", "--breaker-threshold", "0",
         "--mdr-bin", str(stub_mdr))
    assert bi.Layout(root).archive(bi.DATASETS["cath1"]).exists()


def test_scoped_run_never_deletes_the_archive(root, stub_mdr):
    """--systems drains only what was named; the rest of the dataset still needs
    the zip."""
    _init(root, "megamerge")
    _cli("--root", str(root), "run", "--no-verify", "-w", "1",
         "--datasets", "megamerge", "--systems", "1AOY",
         "--mdr-bin", str(stub_mdr))
    assert bi.Layout(root).archive(bi.DATASETS["megamerge"]).exists()
