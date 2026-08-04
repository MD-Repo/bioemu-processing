#!/usr/bin/env python3
"""
bioemu_import.py — download the BioEmu MD dataset release from Zenodo and import
it into MDRepo via `mdr-process`.

The release (Lewis et al., Science 2025; SI S.1.4/S.1.5) is published as five
zip archives across three Zenodo records:

    ONE_octapeptides.zip                    1100 systems   0.5 GB   record 15641199
    MSR_megasim_merge.zip                    271 systems   1.3 GB   record 15641184
    ONE_cath1.zip                             50 systems   2.1 GB   record 15629740
    MSR_megasim_mutants_disp_allatom.zip   21458 systems   8.4 GB   record 15641184
    MSR_cath2.zip                           1043 systems  27.7 GB   record 15629740

Every archive holds self-contained system directories:

    $ZIP_ROOT/$SYSTEM/dataset.json          force_field, temperature_K, save_traj_ns
                     topology.pdb           protein-only all-atom structure
                     reference.pdb          seed structure (MegaSim only)
                     trajs/*.xtc            coordinates, 10 ns/frame
                     trajs/*.json           per-trajectory force-field override

Pipeline (per system):
    claim (flock) -> read the system out of the local zip -> split its
    trajectories into force-field groups -> per group: build an IN_DIR
    (.pdb + generated .psf + .xtc + mdrepo-metadata.toml) -> `mdr-process
    validate` + `mdr-process process` -> delete that IN_DIR -> mark done.

An MDRepo **simulation = one (system, force-field group)**. Most systems yield
exactly one; the 77 MegaSim wildtypes that mix amber ff14sb (folded state) with
amber ff99sb-disp (unfolded state) yield two, because MDRepo records a single
`forcefield` per simulation. Expected totals:

    ONE_octapeptides    1100     MSR_megasim_merge      348  (271 + 77 split)
    ONE_cath1             50     MSR_megasim_mutants  21458
    MSR_cath2           1040  (3 of 1043 ship no trajectories)   TOTAL  23996

Design goals (shared with the sibling mdcath_import.py / dynamicpdb_import.py):
  * Bounded disk: one archive is held at a time, smallest first, and is deleted
    once drained. Peak is that zip plus one system's IN_DIR per worker.
  * Resumable: all state lives in a SQLite manifest; kill/restart at any time.
  * Safe with several worker processes: per-system flock() is the mutex, proven
    to exclude by a two-process self-test before any -w > 1 run starts.
  * Unattended-safe: infrastructure failures are classified as TransientError,
    requeue without spending an attempt, and trip a shared circuit breaker. A
    watchdog kills any mdr-process whose subtree stops making progress.

Two things differ from mdcath_import.py and are worth knowing:
  * Trajectories are deposited byte-for-byte as published. The XTCs already
    carry correct 10 ns frame times, so nothing is re-encoded — which also means
    mdtraj's silent xdrfile-overflow failure mode cannot occur here, and there
    is no coordinate-rewriting scan to run.
  * BioEmu ships no topology file, only topology.pdb, and MDRepo does not accept
    a PDB as topology. ParmEd generates a matching .psf per system.

Subcommands:
    init            build the work-list from the archives + the SIFTS UniProt map
    fetch           download (and verify) archives without importing
    run             launch N workers to import pending systems
    status          print progress
    extract         stage ONE system to a directory (debug; no import)
    reset-failed    requeue failed systems
    reset-sim       forget import records so a system's simulations are redone
    resolve-import  record the outcome of an interrupted push

Requires (on the processing VM): python3, parmed, and `mdr-process` on PATH.
See requirements.txt / setup_env.sh.
"""
from __future__ import annotations

import argparse
import contextlib
import fcntl
import gzip
import hashlib
import io
import json
import errno
import logging
import multiprocessing as mp
import os
import re
import shutil
import signal
import sqlite3
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# --------------------------------------------------------------------------- #
# Source archives
# --------------------------------------------------------------------------- #

ZENODO_API = "https://zenodo.org/api/records/{record}"
ZENODO_DOI = "https://doi.org/10.5281/zenodo.{record}"
SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
BIOEMU_REPO_URL = "https://github.com/microsoft/bioemu"


@dataclass(frozen=True)
class Dataset:
    """One published zip archive == one BioEmu sub-dataset."""
    key: str            # our short handle, also the manifest/CLI name
    record: str         # Zenodo record id
    zip_name: str       # file within that record
    zip_root: str       # top-level directory inside the zip
    label: str          # name used in the manuscript (Table S1)
    si_section: str     # SI section documenting the MD protocol
    ensemble: str       # production ensemble
    id_style: str       # how to derive pdb_id/uniprot_ids: cath | megasim | none
    protocol: str       # which provenance block applies: standard | megasim
    blurb: str          # one sentence for `description`, dataset-specific


# Ordered smallest-first: `run` drains one archive before fetching the next, so
# this is also the order in which disk is consumed.
DATASETS: dict[str, Dataset] = {
    "opep": Dataset(
        key="opep", record="15641199",
        zip_name="ONE_octapeptides.zip", zip_root="ONE_octapeptides",
        label="ONE-octapeptides", si_section="S.1.5.1", ensemble="NPT",
        id_style="none", protocol="standard",
        blurb=(
            "One of 1100 octapeptides sampled to represent peptide equilibrium "
            "ensembles. Trajectories combine the originally published adaptive-"
            "sampling set with five additional 1 us trajectories generated for "
            "BioEmu"
        ),
    ),
    "megamerge": Dataset(
        key="megamerge", record="15641184",
        zip_name="MSR_megasim_merge.zip", zip_root="MSR_megasim_merge",
        label="MSR-megasim", si_section="S.1.5.4", ensemble="NVT",
        id_style="megasim", protocol="megasim",
        blurb=(
            "One of 271 wild-type domains from the MEGAscale protein-stability "
            "set (Tsuboyama et al., Nature 2023), simulated to capture folding-"
            "unfolding transitions from both folded and thermally denatured "
            "starting structures"
        ),
    ),
    "cath1": Dataset(
        key="cath1", record="15629740",
        zip_name="ONE_cath1.zip", zip_root="ONE_cath1",
        label="ONE-cath1", si_section="S.1.5.2", ensemble="NPT",
        id_style="cath", protocol="standard",
        blurb=(
            "One of 50 CATH domains sampled to ~100 us cumulative simulation "
            "time each, by adaptive sampling over three epochs seeded from a "
            "reference PDB structure and then by minRMSD clustering of earlier "
            "epochs"
        ),
    ),
    "megamut": Dataset(
        key="megamut", record="15641184",
        zip_name="MSR_megasim_mutants_disp_allatom.zip",
        zip_root="MSR_megasim_mutants_disp_allatom",
        label="MSR-megasim-mutants", si_section="S.1.5.4", ensemble="NVT",
        id_style="megasim", protocol="megasim",
        blurb=(
            "One of 21458 point mutants of the MEGAscale wild-type domains "
            "(Tsuboyama et al., Nature 2023). Each mutant was seeded from its "
            "wild-type folded conformation with the side chain exchanged and "
            "energy-minimised, then simulated for 1 us in the folded state"
        ),
    ),
    "cath2": Dataset(
        key="cath2", record="15629740",
        zip_name="MSR_cath2.zip", zip_root="MSR_cath2",
        label="MSR-cath2", si_section="S.1.5.3", ensemble="NPT",
        id_style="cath", protocol="standard",
        blurb=(
            "One of ~1040 CATH v4.3.0 domains (50-200 residues, contiguous, no "
            "disulfides, coil fraction below 50%) selected for sequence "
            "coverage and sampled to ~39 us each by adaptive sampling over two "
            "reseeding epochs"
        ),
    ),
}
DATASET_ORDER = list(DATASETS)

# --------------------------------------------------------------------------- #
# MD provenance
#
# Everything here is transcribed from the BioEmu manuscript SI (S.1.4 general
# protocol, S.1.5.x per dataset) rather than inferred from the files. Items
# marked "NOT PUBLISHED" are the ones to confirm with the authors before a prod
# run — see the "Provenance notes" section of README.md.
# --------------------------------------------------------------------------- #

# S.1.4: "We internally developed code specifically tailored towards running
# large MD production campaigns on Azure compute resources. Our code is based on
# OpenMM as its compute engine, albeit setups are generated using OpenMM or
# GROMACS as a backend."  MDRepo's software_name vocabulary is closed
# (ACEMD/AMBER/CHARMM/CUSTOM/GROMACS/NAMD/SPONGE) and contains neither OpenMM nor
# "in-house harness", so CUSTOM is the honest entry; the spec then requires
# software_version = "NA". The real engine is preserved in the description.
SOFTWARE_NAME = "CUSTOM"
SOFTWARE_VERSION = "NA"
ENGINE_NAME = "OpenMM"

# S.1.4: production uses hydrogen mass repartitioning (H = 4 amu) with h-bond
# constraints at a 4 fs timestep. (Equilibration ran at 2 fs; MDRepo records the
# production value.) S.1.5.4 confirms the same 4 fs for MegaSim.
INTEGRATION_TIMESTEP_FS = 4

# Force fields, keyed by the exact string in the source dataset.json.
@dataclass(frozen=True)
class ForceField:
    label: str                    # MDRepo `forcefield`
    slug: str                     # filename-safe group id
    water: Optional[str]          # MDRepo [water].model; None omits the table
    comments: str


FORCEFIELDS: dict[str, ForceField] = {
    "amber ff99sb-ildn": ForceField(
        label="Amber ff99SB-ILDN", slug="ff99sb-ildn", water="TIP3P",
        comments=(
            "Explicit solvent, TIP3P water, cubic box with 1 nm padding and a "
            "0.1 M NaCl buffer. Equilibration: 0.1 ns NVT then 0.9 ns NPT with "
            "harmonic restraints on solute heavy atoms, released over a further "
            "0.1 ns, at a 2 fs timestep. Production: NPT, hydrogen mass "
            "repartitioning (hydrogen mass 4 amu) with hydrogen-bond "
            "constraints, 4 fs timestep, 300 K and 1 bar."
        ),
    ),
    "amber ff14sb": ForceField(
        label="Amber ff14SB", slug="ff14sb", water=None,
        comments=(
            "Explicit solvent in a rhombic dodecahedral box with 1.5 nm "
            "padding. Equilibration: 0.2 ns NVT then 0.6 ns NPT targeting 295 K "
            "and 1 bar with a Langevin integrator at a 4 fs timestep. "
            "Production: NVT at 295 K, hydrogen mass repartitioning (hydrogen "
            "mass 4 amu) with hydrogen-bond constraints, 4 fs timestep. Used "
            "for folded-state sampling of the wild-types whose native fold "
            "a99SB-disp destabilised. The water model is not stated in the "
            "manuscript for this dataset and is therefore not recorded here."
        ),
    ),
    "amber ff99sb-disp": ForceField(
        label="Amber a99SB-disp", slug="ff99sb-disp", water=None,
        comments=(
            "Explicit solvent in a rhombic dodecahedral box with 1.5 nm "
            "padding. Equilibration: 0.2 ns NVT then 0.6 ns NPT targeting 295 K "
            "and 1 bar with a Langevin integrator at a 4 fs timestep. "
            "Production: NVT at 295 K, hydrogen mass repartitioning (hydrogen "
            "mass 4 amu) with hydrogen-bond constraints, 4 fs timestep. "
            "a99SB-disp is parameterised for disordered states and is supplied "
            "with its own four-point water model; the manuscript does not state "
            "the water model used for this dataset, so it is not recorded here."
        ),
    ),
}

# S.1.4 states TIP3P and a 0.1 M NaCl buffer for the standard protocol. S.1.5.4
# re-specifies MegaSim's box and equilibration but restates neither, so solutes
# are recorded only for the standard-protocol datasets. NOT PUBLISHED: the
# solvent density, so [water].density_kg_m3 (a required field once [water] is
# present) is the nominal value below.
WATER_DENSITY_KG_M3 = 1000.0
STANDARD_SOLUTES = [("Na+", 0.1), ("Cl-", 0.1)]

# Ensemble/temperature prose per protocol, for `description`.
PROTOCOL_SUMMARY = {
    "standard": (
        "{engine}-based in-house simulation code, {ff} force field, "
        "explicit TIP3P solvent (cubic box, 1 nm padding, 0.1 M NaCl), "
        "NPT production at {temp} K and 1 bar with hydrogen mass "
        "repartitioning and a {dt} fs timestep"
    ),
    "megasim": (
        "{engine}-based in-house simulation code, {ff} force field, "
        "explicit solvent in a rhombic dodecahedral box with 1.5 nm padding, "
        "NVT production at {temp} K with hydrogen mass repartitioning and a "
        "{dt} fs timestep"
    ),
}

# [[papers]] — title/authors/journal/year/volume required; number/pages/doi optional.
PAPER_BIOEMU = {
    "title": "Scalable emulation of protein equilibrium ensembles with generative deep learning",
    "authors": (
        "Sarah Lewis, Tim Hempel, José Jiménez-Luna, Michael Gastegger, Yu Xie, "
        "Andrew Y. K. Foong, Victor García Satorras, Osama Abdin, Bastiaan S. Veeling, "
        "Iryna Zaporozhets, Yaoyi Chen, Soojung Yang, Adam E. Foster, Arne Schneuing, "
        "Jigyasa Nigam, Federico Barbero, Vincent Stimper, Andrew Campbell, Jason Yim, "
        "Marten Lienen, Yu Shi, Shuxin Zheng, Hannes Schulz, Usman Munir, "
        "Roberto Sordillo, Ryota Tomioka, Cecilia Clementi, Frank Noé"
    ),
    "journal": "Science",
    "year": 2025,
    "volume": 389,
    "number": "6761",
    "pages": "eadv9817",
    "doi": "10.1126/science.adv9817",
}
# MegaSim seed structures and the experimental stability measurements they came from.
PAPER_MEGASCALE = {
    "title": "Mega-scale experimental analysis of protein folding stability in biology and design",
    "authors": (
        "Kotaro Tsuboyama, Justas Dauparas, Jonathan Chen, Elodie Laine, "
        "Yasser Mohseni Behbahani, Jonathan J. Weinstein, Niall M. Mangan, "
        "Sergey Ovchinnikov, Gabriel J. Rocklin"
    ),
    "journal": "Nature",
    "year": 2023,
    "volume": 620,
    "number": "7973",
    "pages": "434-444",
    "doi": "10.1038/s41586-023-06328-6",
}

# `lead_contributor_orcid` is required by the spec; this placeholder matches the
# sibling importers and means "submitted locally as administrator". It is
# deliberately NOT one of the author ORCIDs below: the lead contributor is
# whoever submits the deposition, not an author of the paper.
LEAD_CONTRIBUTOR_ORCID = "0000-0000-0000-0000"
MSR_AI4SCIENCE = "Microsoft Research AI for Science"

# The full BioEmu author list, in manuscript order, so credit in the deposition
# matches credit in the paper — the names here are kept byte-identical to
# PAPER_BIOEMU["authors"] above, and a test enforces that.
#
# Names carry their diacritics (José Jiménez-Luna, Victor García Satorras,
# Frank Noé, Freie Universität), so this module emits non-ASCII text: every
# write of rendered metadata must name encoding="utf-8" rather than inherit the
# locale's, or a C/POSIX-locale VM raises UnicodeEncodeError. TOML is UTF-8 by
# spec, so this is also what a reader expects.
#
# `orcid` is optional and is present only for the 8 authors whose ORCID has been
# checked against their public profile; the release itself publishes no ORCIDs,
# so the remaining 20 carry a name alone rather than a guessed or placeholder id.
# `institution` is likewise recorded only where it is known from the dataset
# release (the four original contacts) — the manuscript's affiliations were not
# transcribed, and an unverified affiliation is worse than none.
CONTRIBUTORS = [
    {"name": "Sarah Lewis", "orcid": "0009-0009-6484-0352",
     "institution": MSR_AI4SCIENCE},
    {"name": "Tim Hempel", "orcid": "0000-0002-0073-9353",
     "email": "timhempel@microsoft.com", "institution": MSR_AI4SCIENCE},
    {"name": "José Jiménez-Luna", "orcid": "0000-0002-5335-7834"},
    {"name": "Michael Gastegger", "orcid": "0000-0001-7954-3275"},
    {"name": "Yu Xie", "orcid": "0000-0002-0088-5123"},
    {"name": "Andrew Y. K. Foong"},
    {"name": "Victor García Satorras"},
    {"name": "Osama Abdin", "orcid": "0000-0002-5471-2906"},
    {"name": "Bastiaan S. Veeling"},
    {"name": "Iryna Zaporozhets"},
    {"name": "Yaoyi Chen"},
    {"name": "Soojung Yang"},
    {"name": "Adam E. Foster"},
    {"name": "Arne Schneuing"},
    {"name": "Jigyasa Nigam"},
    {"name": "Federico Barbero"},
    {"name": "Vincent Stimper"},
    {"name": "Andrew Campbell"},
    {"name": "Jason Yim"},
    {"name": "Marten Lienen"},
    {"name": "Yu Shi"},
    {"name": "Shuxin Zheng"},
    {"name": "Hannes Schulz"},
    {"name": "Usman Munir"},
    {"name": "Roberto Sordillo"},
    {"name": "Ryota Tomioka"},
    {"name": "Cecilia Clementi", "orcid": "0000-0001-9221-2358",
     "institution": "Freie Universität Berlin, Department of Physics"},
    {"name": "Frank Noé", "orcid": "0000-0003-4169-9324",
     "email": "franknoe@microsoft.com", "institution": MSR_AI4SCIENCE},
]

# Deposit the MegaSim seed structure alongside the trajectory topology. Flip to
# False if `mdr-process validate` objects to [[additional_files]].
INCLUDE_REFERENCE_PDB = True

SHORT_DESCRIPTION_MAX = 300      # spec limit, enforced by truncate_short()

# --------------------------------------------------------------------------- #
# Operational guards
# --------------------------------------------------------------------------- #

MAX_ATTEMPTS = 3          # per-system retry budget for run
ERROR_TAIL_CHARS = 4000   # how much of a failed command's output to persist

MIN_FREE_GB = 25          # refuse to start an archive download/extraction below this
TRANSIENT_BACKOFF_SEC = 60
BREAKER_THRESHOLD = 10    # consecutive failures across all workers before halting
CONTENTION_ROUNDS = 60    # rounds to wait while every candidate is peer-locked

STALL_MINUTES = 10        # kill mdr-process after this long with zero subtree progress
MDR_MAX_HOURS = 4.0       # absolute per-command ceiling
WATCHDOG_POLL_SEC = 30
_CLK_TCK = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100

DOWNLOAD_CHUNK = 8 << 20
DOWNLOAD_RETRIES = 5

log = logging.getLogger("bioemu")


class TransientError(RuntimeError):
    """Infrastructure failure (disk, network, remote 5xx) rather than bad input.

    These do not consume the per-system attempt budget: a single outage would
    otherwise burn MAX_ATTEMPTS on thousands of systems at once, N times faster
    with N workers.
    """


class SkipSystem(RuntimeError):
    """The system holds nothing importable (e.g. no trajectories were released).

    Recorded and marked done rather than failed: 3 of MSR_cath2's 1043 system
    directories ship only dataset.json + topology.pdb, which is a fact about the
    release, not an error to retry.
    """


# mdr-process/HTTP output that indicates the environment failed, not the input.
# Matched against the last ERROR_TAIL_CHARS of a command's log, which also holds
# ordinary output, so every pattern must be one that cannot plausibly appear in a
# successful run: a false positive requeues forever without spending an attempt.
_TRANSIENT_PATTERNS = (
    "no space left", "disk quota exceeded", "connection refused",
    "connection reset", "connection aborted", "timed out", "timeout",
    "temporarily unavailable", "too many requests", "429", "500 server error",
    "502", "503", "504", "bad gateway", "service unavailable",
    "name or service not known", "network is unreachable",
    "incompleteread", "content-length",
    # iRODS/gocmd, reached through mdr-process's push step.
    "networkexception", "could not connect", "unable to connect",
    "no route to host", "host is unreachable",
    "temporary failure in name resolution", "broken pipe",
)


def _is_transient_text(text: str) -> bool:
    low = (text or "").lower()
    return any(p in low for p in _TRANSIENT_PATTERNS)


def _is_transient_exc(exc: BaseException) -> bool:
    if isinstance(exc, TransientError):
        return True
    if isinstance(exc, (ConnectionError, TimeoutError, urllib.error.URLError)):
        return True
    if isinstance(exc, OSError) and exc.errno in (errno.ENOSPC, errno.EDQUOT, errno.EIO):
        return True
    return _is_transient_text(f"{type(exc).__name__}: {exc}")


def free_gb(path: Path) -> float:
    return shutil.disk_usage(str(path)).free / 1024**3


def require_free_space(path: Path, min_gb: float, what: str) -> None:
    """Raise TransientError rather than fail a system when the volume is full."""
    have = free_gb(path)
    if have < min_gb:
        raise TransientError(
            f"only {have:.1f} GB free on {path} (need {min_gb:.1f} GB for {what})"
        )


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# --------------------------------------------------------------------------- #
# Layout / paths
# --------------------------------------------------------------------------- #

@dataclass
class Layout:
    root: Path

    @property
    def db(self) -> Path: return self.root / "manifest.sqlite"
    @property
    def archives(self) -> Path: return self.root / "archives"
    @property
    def staging(self) -> Path: return self.root / "staging"
    @property
    def locks(self) -> Path: return self.root / "locks"
    @property
    def logs(self) -> Path: return self.root / "logs"

    def ensure(self) -> None:
        for p in (self.root, self.archives, self.staging, self.locks, self.logs):
            p.mkdir(parents=True, exist_ok=True)

    def archive(self, ds: Dataset) -> Path:
        return self.archives / ds.zip_name

    def system_stage(self, dataset: str, system: str) -> Path:
        return self.staging / dataset / system

    def sim_dir(self, dataset: str, system: str, grp: str) -> Path:
        return self.staging / dataset / system / grp

    def lock_file(self, dataset: str, system: str) -> Path:
        return self.locks / f"{dataset}__{system}.lock"

    def archive_lock(self, ds: Dataset) -> Path:
        return self.locks / f"archive__{ds.key}.lock"

    def log_file(self, dataset: str, system: str, grp: str) -> Path:
        return self.logs / f"{dataset}__{system}__{grp}.log"


def _fresh_dir(in_dir: Path) -> Path:
    """Make `in_dir` exist and be empty — but ONLY if it is safe to wipe.

    IN_DIRs are always a `<dataset>/<system>/<group>` leaf the tool creates and
    owns. This guard refuses to rmtree the current directory, $HOME, the
    filesystem root, or any ancestor of the cwd, so a mistaken `--out-dir .`
    (which would otherwise delete sibling files) is rejected instead of obeyed.
    """
    resolved = in_dir.resolve()
    forbidden = {Path("/").resolve(), Path.home().resolve(), Path.cwd().resolve()}
    if resolved in forbidden or resolved in Path.cwd().resolve().parents:
        raise ValueError(
            f"refusing to use {in_dir} as an IN_DIR (would clear a protected "
            f"directory). Pass an --out-dir under which per-system subdirs are created."
        )
    if in_dir.exists():
        shutil.rmtree(in_dir)
    in_dir.mkdir(parents=True)
    return in_dir


# --------------------------------------------------------------------------- #
# Identifiers
#
# cath1/cath2 system names embed a CATH domain id: cath1_1b43A02 -> PDB 1b43,
# chain A. MegaSim names are MEGAscale entry names: a bare PDB code (1AOY), a
# PDB code with a variant suffix (1A0N_L7S, 2HBB_pross6), a de novo design
# (EEHEE_rd3_0019, HEEH_KT_rd6_0007), or a mutant of any of those
# (1A0N_L7S__A12D). Octapeptides are synthetic and carry no identifier at all.
# --------------------------------------------------------------------------- #

CATH_SYSTEM_RE = re.compile(r"^cath[12]_(?P<domain>[0-9A-Za-z]{4}[0-9A-Za-z][0-9]{2})$")
PDB_CODE_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")


@dataclass(frozen=True)
class SystemIds:
    pdb_id: Optional[str] = None        # lowercase 4-char accession
    chain: Optional[str] = None         # None => union across all chains
    parent: Optional[str] = None        # MegaSim wild-type this derives from
    mutation: Optional[str] = None      # MegaSim point mutation, e.g. A12D
    variant: Optional[str] = None       # MegaSim wild-type variant tag, e.g. L7S
    design: bool = False                # de novo design, no PDB ancestry


def parse_system_ids(ds: Dataset, system: str) -> SystemIds:
    """Derive PDB provenance from a system directory name."""
    if ds.id_style == "cath":
        m = CATH_SYSTEM_RE.match(system)
        if not m:
            raise ValueError(f"unexpected CATH system name {system!r}")
        dom = m.group("domain")
        return SystemIds(pdb_id=dom[:4].lower(), chain=dom[4])

    if ds.id_style == "megasim":
        # `__` separates a mutant from its wild-type; a single `_` separates a
        # wild-type's variant tag from the PDB code it derives from.
        parent, _, mutation = system.partition("__")
        head, _, variant = parent.partition("_")
        if PDB_CODE_RE.match(head):
            # Chainless: MEGAscale entries name a domain, not a chain, so UniProt
            # accessions are taken as the union across the entry's chains.
            return SystemIds(pdb_id=head.lower(),
                             parent=parent if mutation else None,
                             mutation=mutation or None,
                             variant=variant or None)
        return SystemIds(parent=parent if mutation else None,
                         mutation=mutation or None, design=True)

    return SystemIds()      # octapeptides


# --------------------------------------------------------------------------- #
# Manifest (SQLite). The rollback journal works cleanly for the intended
# single-VM deployment; busy_timeout also permits concurrent status reads.
# --------------------------------------------------------------------------- #

class Manifest:
    def __init__(self, path: Path):
        self.path = path
        self.conn = sqlite3.connect(str(path), timeout=60)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=DELETE")
        self.conn.execute("PRAGMA busy_timeout=60000")
        # The import-intent commit must be durable before mdr-process can push.
        self.conn.execute("PRAGMA synchronous=FULL")
        self._init_schema()

    def _init_schema(self) -> None:
        self.conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS systems (
                dataset    TEXT NOT NULL,
                system     TEXT NOT NULL,
                status     TEXT NOT NULL DEFAULT 'pending',  -- pending|done|failed|skipped
                attempts   INTEGER NOT NULL DEFAULT 0,
                n_trajs    INTEGER,
                error      TEXT,
                updated_at TEXT,
                PRIMARY KEY (dataset, system)
            );
            CREATE TABLE IF NOT EXISTS sims (      -- import lifecycle per (system, ff group)
                dataset     TEXT NOT NULL,
                system      TEXT NOT NULL,
                grp         TEXT NOT NULL,
                state       TEXT NOT NULL DEFAULT 'imported',  -- importing|imported
                log_path    TEXT,
                started_at  TEXT,
                imported_at TEXT,
                PRIMARY KEY (dataset, system, grp)
            );
            CREATE TABLE IF NOT EXISTS uniprot (
                dataset TEXT NOT NULL,
                system  TEXT NOT NULL,
                acc     TEXT NOT NULL,
                PRIMARY KEY (dataset, system, acc)
            );
            CREATE TABLE IF NOT EXISTS archives (
                dataset     TEXT PRIMARY KEY,
                size        INTEGER,
                md5         TEXT,
                url         TEXT,
                state       TEXT NOT NULL DEFAULT 'absent',   -- absent|present|drained
                updated_at  TEXT
            );
            CREATE TABLE IF NOT EXISTS notes (     -- non-fatal observations worth keeping
                dataset  TEXT NOT NULL,
                system   TEXT NOT NULL,
                kind     TEXT NOT NULL,
                detail   TEXT,
                noted_at TEXT,
                PRIMARY KEY (dataset, system, kind)
            );
            CREATE INDEX IF NOT EXISTS idx_systems_status ON systems(dataset, status);
            """
        )
        self.conn.commit()

    # ---- population ----
    def add_systems(self, dataset: str, rows: list[tuple[str, int]]) -> int:
        cur = self.conn.executemany(
            "INSERT OR IGNORE INTO systems(dataset, system, n_trajs, updated_at) "
            "VALUES (?, ?, ?, ?)",
            [(dataset, s, n, _now()) for s, n in rows],
        )
        self.conn.commit()
        return cur.rowcount

    def set_uniprot(self, rows: list[tuple[str, str, str]]) -> None:
        self.conn.executemany(
            "INSERT OR IGNORE INTO uniprot(dataset, system, acc) VALUES (?, ?, ?)", rows
        )
        self.conn.commit()

    def clear_uniprot(self) -> None:
        self.conn.execute("DELETE FROM uniprot")
        self.conn.commit()

    def uniprot_for(self, dataset: str, system: str) -> list[str]:
        rows = self.conn.execute(
            "SELECT acc FROM uniprot WHERE dataset=? AND system=? ORDER BY acc",
            (dataset, system),
        ).fetchall()
        return [r["acc"] for r in rows]

    def systems_needing_ids(self, dataset: str) -> list[str]:
        rows = self.conn.execute(
            "SELECT system FROM systems WHERE dataset=?", (dataset,)
        ).fetchall()
        return [r["system"] for r in rows]

    # ---- archives ----
    def set_archive(self, dataset: str, size: int, md5: Optional[str],
                    url: str, state: str) -> None:
        self.conn.execute(
            """INSERT INTO archives(dataset, size, md5, url, state, updated_at)
               VALUES (?, ?, ?, ?, ?, ?)
               ON CONFLICT(dataset) DO UPDATE SET
                   size=excluded.size, md5=excluded.md5, url=excluded.url,
                   state=excluded.state, updated_at=excluded.updated_at""",
            (dataset, size, md5, url, state, _now()),
        )
        self.conn.commit()

    def archive_row(self, dataset: str) -> Optional[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM archives WHERE dataset=?", (dataset,)
        ).fetchone()

    def set_archive_state(self, dataset: str, state: str) -> None:
        self.conn.execute(
            "UPDATE archives SET state=?, updated_at=? WHERE dataset=?",
            (state, _now(), dataset),
        )
        self.conn.commit()

    # ---- claiming / status ----
    def candidate_systems(self, dataset: str, limit: int = 256,
                          only: Optional[tuple] = None) -> list[str]:
        sql = [
            "SELECT system FROM systems",
            "WHERE dataset = ?",
            "AND (status = 'pending' OR (status = 'failed' AND attempts < ?))",
        ]
        params: list = [dataset, MAX_ATTEMPTS]
        if only:
            sql.append("AND system IN (%s)" % ",".join("?" for _ in only))
            params.extend(only)
        sql.append("ORDER BY RANDOM() LIMIT ?")
        params.append(limit)
        return [r["system"] for r in self.conn.execute("\n".join(sql), params).fetchall()]

    def pending_count(self, dataset: str, only: Optional[tuple] = None) -> int:
        sql = ["SELECT COUNT(*) c FROM systems WHERE dataset = ?",
               "AND (status = 'pending' OR (status = 'failed' AND attempts < ?))"]
        params: list = [dataset, MAX_ATTEMPTS]
        if only:
            sql.append("AND system IN (%s)" % ",".join("?" for _ in only))
            params.extend(only)
        return self.conn.execute("\n".join(sql), params).fetchone()["c"]

    def system_row(self, dataset: str, system: str) -> Optional[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM systems WHERE dataset=? AND system=?", (dataset, system)
        ).fetchone()

    def imported_groups(self, dataset: str, system: str) -> set[str]:
        rows = self.conn.execute(
            "SELECT grp FROM sims WHERE dataset=? AND system=? AND state='imported'",
            (dataset, system),
        ).fetchall()
        return {r["grp"] for r in rows}

    def importing_groups(self, dataset: str, system: str) -> set[str]:
        rows = self.conn.execute(
            "SELECT grp FROM sims WHERE dataset=? AND system=? AND state='importing'",
            (dataset, system),
        ).fetchall()
        return {r["grp"] for r in rows}

    def mark_importing(self, dataset: str, system: str, grp: str, log_path: str) -> None:
        self.conn.execute(
            """INSERT INTO sims(dataset, system, grp, state, log_path, started_at, imported_at)
               VALUES (?, ?, ?, 'importing', ?, ?, NULL)
               ON CONFLICT(dataset, system, grp) DO UPDATE SET
                   state='importing', log_path=excluded.log_path,
                   started_at=excluded.started_at, imported_at=NULL""",
            (dataset, system, grp, log_path, _now()),
        )
        self.conn.commit()

    def mark_imported(self, dataset: str, system: str, grp: str, log_path: str) -> None:
        self.conn.execute(
            """INSERT INTO sims(dataset, system, grp, state, log_path, started_at, imported_at)
               VALUES (?, ?, ?, 'imported', ?, NULL, ?)
               ON CONFLICT(dataset, system, grp) DO UPDATE SET
                   state='imported', log_path=excluded.log_path,
                   imported_at=excluded.imported_at""",
            (dataset, system, grp, log_path, _now()),
        )
        self.conn.commit()

    def retry_import(self, dataset: str, system: str, grp: str) -> bool:
        cur = self.conn.execute(
            "DELETE FROM sims WHERE dataset=? AND system=? AND grp=? AND state='importing'",
            (dataset, system, grp),
        )
        self.conn.commit()
        return bool(cur.rowcount)

    def reset_sims(self, dataset: str, system: str,
                   groups: Optional[list[str]] = None) -> int:
        sql = "DELETE FROM sims WHERE dataset=? AND system=?"
        params: list = [dataset, system]
        if groups:
            sql += " AND grp IN (%s)" % ",".join("?" for _ in groups)
            params.extend(groups)
        cur = self.conn.execute(sql, params)
        self.conn.commit()
        return cur.rowcount

    def ambiguous(self) -> list[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM sims WHERE state='importing' ORDER BY dataset, system, grp"
        ).fetchall()

    # ---- state transitions ----
    def requeue(self, dataset: str, system: str) -> None:
        self.conn.execute(
            "UPDATE systems SET status='pending', error=NULL, updated_at=? "
            "WHERE dataset=? AND system=?", (_now(), dataset, system),
        )
        self.conn.commit()

    def reset_system(self, dataset: str, system: str) -> None:
        """Requeue one system *and* clear its attempt budget (explicit CLI ask)."""
        self.conn.execute(
            "UPDATE systems SET status='pending', attempts=0, error=NULL, updated_at=? "
            "WHERE dataset=? AND system=?", (_now(), dataset, system),
        )
        self.conn.commit()

    def mark_done(self, dataset: str, system: str) -> None:
        self.conn.execute(
            "UPDATE systems SET status='done', error=NULL, updated_at=? "
            "WHERE dataset=? AND system=?", (_now(), dataset, system),
        )
        self.conn.commit()

    def mark_skipped(self, dataset: str, system: str, reason: str) -> None:
        self.conn.execute(
            "UPDATE systems SET status='skipped', error=?, updated_at=? "
            "WHERE dataset=? AND system=?", (reason, _now(), dataset, system),
        )
        self.conn.commit()

    def mark_failed(self, dataset: str, system: str, error: str) -> None:
        self.conn.execute(
            "UPDATE systems SET status='failed', attempts=attempts+1, error=?, "
            "updated_at=? WHERE dataset=? AND system=?",
            (error[-ERROR_TAIL_CHARS:], _now(), dataset, system),
        )
        self.conn.commit()

    def reset_failed(self, dataset: Optional[str] = None) -> int:
        sql = ("UPDATE systems SET status='pending', attempts=0, error=NULL, "
               "updated_at=? WHERE status='failed'")
        params: list = [_now()]
        if dataset:
            sql += " AND dataset=?"
            params.append(dataset)
        cur = self.conn.execute(sql, params)
        self.conn.commit()
        return cur.rowcount

    def note(self, dataset: str, system: str, kind: str, detail: str) -> None:
        self.conn.execute(
            """INSERT INTO notes(dataset, system, kind, detail, noted_at)
               VALUES (?, ?, ?, ?, ?)
               ON CONFLICT(dataset, system, kind) DO UPDATE SET
                   detail=excluded.detail, noted_at=excluded.noted_at""",
            (dataset, system, kind, detail[-ERROR_TAIL_CHARS:], _now()),
        )
        self.conn.commit()

    def note_rows(self, limit: int = 0) -> list[sqlite3.Row]:
        sql = "SELECT * FROM notes ORDER BY kind, dataset, system"
        if limit:
            sql += " LIMIT ?"
            return self.conn.execute(sql, (limit,)).fetchall()
        return self.conn.execute(sql).fetchall()

    # ---- reporting ----
    def counts(self) -> dict[str, dict[str, int]]:
        out: dict[str, dict[str, int]] = {}
        for r in self.conn.execute(
            "SELECT dataset, status, COUNT(*) c FROM systems GROUP BY dataset, status"
        ):
            out.setdefault(r["dataset"], {})[r["status"]] = r["c"]
        for r in self.conn.execute(
            "SELECT dataset, COUNT(*) c FROM sims WHERE state='imported' GROUP BY dataset"
        ):
            out.setdefault(r["dataset"], {})["sims_imported"] = r["c"]
        return out

    def failures(self, limit: int) -> list[sqlite3.Row]:
        return self.conn.execute(
            "SELECT dataset, system, attempts, error FROM systems WHERE status='failed' "
            "ORDER BY updated_at DESC LIMIT ?", (limit,),
        ).fetchall()

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self.conn.close()


# --------------------------------------------------------------------------- #
# Zenodo: archive metadata + resumable, checksummed download
# --------------------------------------------------------------------------- #

def _urlopen(url: str, headers: Optional[dict] = None, timeout: int = 120):
    req = urllib.request.Request(url, headers=headers or {})
    return urllib.request.urlopen(req, timeout=timeout)


def _download_text(url: str) -> str:
    with _urlopen(url) as r:
        return r.read().decode("utf-8", "replace")


def zenodo_file_info(ds: Dataset) -> tuple[str, int, Optional[str]]:
    """(download url, size, md5) for a dataset's archive, from the Zenodo API."""
    try:
        meta = json.loads(_download_text(ZENODO_API.format(record=ds.record)))
    except Exception as e:  # noqa: BLE001 - network/API failures are transient
        raise TransientError(f"could not read Zenodo record {ds.record}: {e}") from e
    for f in meta.get("files", []):
        if f.get("key") == ds.zip_name:
            checksum = f.get("checksum") or ""
            md5 = checksum.split("md5:", 1)[1] if checksum.startswith("md5:") else None
            url = (f.get("links") or {}).get("self")
            if not url:
                url = (f"https://zenodo.org/records/{ds.record}/files/"
                       f"{ds.zip_name}?download=1")
            return url, int(f.get("size") or 0), md5
    raise RuntimeError(f"{ds.zip_name} not present in Zenodo record {ds.record}")


def file_md5(path: Path, chunk: int = DOWNLOAD_CHUNK) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def download_archive(lay: Layout, ds: Dataset, man: Manifest,
                     min_free_gb: float = MIN_FREE_GB,
                     verify: bool = True) -> Path:
    """Fetch ds's zip into <root>/archives, resuming a partial file if present.

    Zenodo serves ranged requests, so an interrupted multi-GB transfer resumes
    instead of restarting; the published md5 is then checked before the archive
    is used. A mismatch deletes the file and raises a plain RuntimeError so a
    persistently corrupt mirror cannot requeue forever.
    """
    dest = lay.archive(ds)
    part = dest.with_suffix(dest.suffix + ".part")
    if dest.exists() and not verify:
        # An archive already on disk is taken as-is when verification is waived.
        # This is the escape hatch for a VM that cannot reach Zenodo, or one
        # where the zips were staged by hand — without it, `--no-verify` would
        # still need the API round-trip just to learn the expected size.
        log.info("[%s] using the archive already in %s (unverified)",
                 ds.key, lay.archives)
        man.set_archive(ds.key, dest.stat().st_size, None, "", "present")
        return dest
    url, size, md5 = zenodo_file_info(ds)
    man.set_archive(ds.key, size, md5, url, "absent")

    if dest.exists() and dest.stat().st_size == size:
        if verify and md5 and file_md5(dest) != md5:
            log.warning("[%s] cached archive failed md5; re-downloading", ds.key)
            dest.unlink()
        else:
            man.set_archive(ds.key, size, md5, url, "present")
            return dest

    need = (size - (part.stat().st_size if part.exists() else 0)) / 1024**3
    require_free_space(lay.archives, max(min_free_gb, need + 2), f"{ds.zip_name}")

    for attempt in range(1, DOWNLOAD_RETRIES + 1):
        have = part.stat().st_size if part.exists() else 0
        if have > size:                      # stale/corrupt partial
            part.unlink()
            have = 0
        if have == size:
            break
        headers = {"Range": f"bytes={have}-"} if have else {}
        log.info("[%s] downloading %s (%.2f/%.2f GB, attempt %d/%d)",
                 ds.key, ds.zip_name, have / 1024**3, size / 1024**3,
                 attempt, DOWNLOAD_RETRIES)
        try:
            with _urlopen(url, headers) as r, open(part, "ab" if have else "wb") as fh:
                if have and r.status != 206:
                    # Server ignored the range: start over rather than concatenate.
                    fh.close()
                    part.unlink()
                    continue
                while True:
                    block = r.read(DOWNLOAD_CHUNK)
                    if not block:
                        break
                    fh.write(block)
        except Exception as e:  # noqa: BLE001 - retry any transport failure
            if attempt == DOWNLOAD_RETRIES:
                raise TransientError(f"{ds.zip_name}: download failed: {e}") from e
            log.warning("[%s] download interrupted (%s); retrying", ds.key, e)
            time.sleep(5 * attempt)

    actual = part.stat().st_size if part.exists() else 0
    if actual != size:
        raise TransientError(
            f"{ds.zip_name}: got {actual} bytes, expected {size}"
        )
    if verify and md5:
        got = file_md5(part)
        if got != md5:
            part.unlink()
            raise RuntimeError(
                f"{ds.zip_name}: md5 {got} != published {md5}; deleted, "
                f"re-download on the next attempt"
            )
        log.info("[%s] md5 verified", ds.key)
    part.replace(dest)
    man.set_archive(ds.key, size, md5, url, "present")
    return dest


@contextlib.contextmanager
def archive_ready(lay: Layout, ds: Dataset, man: Manifest, **kw):
    """Ensure ds's archive is on disk, serialising concurrent fetches via flock."""
    lock_path = lay.archive_lock(ds)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)      # blocking: one downloader, others wait
        path = download_archive(lay, ds, man, **kw)
    yield path


# --------------------------------------------------------------------------- #
# Reading systems out of an archive
# --------------------------------------------------------------------------- #

TRAJ_SUFFIX = ".xtc"
CMPRSD = ".cmprsd"


def sidecar_name(traj_member: str) -> str:
    """trajs/run001_protein.cmprsd.xtc -> trajs/run001_protein.json

    The release names its per-trajectory force-field overrides after the
    trajectory with `.cmprsd.xtc` replaced by `.json` (dataset card example).
    """
    base = traj_member[: -len(TRAJ_SUFFIX)]
    if base.endswith(CMPRSD):
        base = base[: -len(CMPRSD)]
    return base + ".json"


@dataclass
class SystemFiles:
    system: str
    dataset_json: str
    topology: Optional[str]
    reference: Optional[str]
    trajs: list[str] = field(default_factory=list)
    sidecars: dict[str, str] = field(default_factory=dict)   # traj member -> json member


def index_archive(zf: zipfile.ZipFile, ds: Dataset) -> dict[str, SystemFiles]:
    """Group an archive's namelist into per-system file sets."""
    prefix = ds.zip_root + "/"
    out: dict[str, SystemFiles] = {}
    members = [n for n in zf.namelist() if n.startswith(prefix) and not n.endswith("/")]
    by_system: dict[str, list[str]] = {}
    for n in members:
        rest = n[len(prefix):]
        system, sep, tail = rest.partition("/")
        if not sep or not system:
            continue
        by_system.setdefault(system, []).append(n)

    for system, paths in by_system.items():
        sf = SystemFiles(system=system, dataset_json="", topology=None, reference=None)
        jsons = set()
        for p in paths:
            base = p.rsplit("/", 1)[-1]
            if base == "dataset.json":
                sf.dataset_json = p
            elif base == "topology.pdb":
                sf.topology = p
            elif base == "reference.pdb":
                sf.reference = p
            elif p.endswith(TRAJ_SUFFIX):
                sf.trajs.append(p)
            elif p.endswith(".json"):
                jsons.add(p)
        sf.trajs.sort()
        for t in sf.trajs:
            s = sidecar_name(t)
            if s in jsons:
                sf.sidecars[t] = s
        out[system] = sf
    return out


@dataclass
class Group:
    """Trajectories of one system sharing a force field and temperature."""
    slug: str
    force_field: str          # raw string from the json
    temperature_k: int
    save_traj_ns: float
    trajs: list[str] = field(default_factory=list)
    overridden: bool = False  # at least one traj came from a sidecar json


def _read_json_member(zf: zipfile.ZipFile, member: str) -> dict:
    return json.loads(zf.read(member).decode("utf-8"))


def group_slug(force_field: str, temperature_k: int) -> str:
    ff = FORCEFIELDS.get(force_field)
    base = ff.slug if ff else re.sub(r"[^a-z0-9]+", "-", force_field.lower()).strip("-")
    return f"{base}-{temperature_k}k"


def discover_groups(zf: zipfile.ZipFile, sf: SystemFiles) -> list[Group]:
    """Split a system's trajectories by (force field, temperature).

    MDRepo records a single `forcefield` per simulation, so a system whose
    trajectories were produced with more than one force field has to become more
    than one simulation. Trajectories inherit dataset.json unless a sidecar json
    overrides it — which is how the 77 MegaSim wild-types that mix ff14sb (folded
    state) with a99SB-disp (unfolded state) are marked.
    """
    base = _read_json_member(zf, sf.dataset_json)
    groups: dict[tuple[str, int], Group] = {}
    for traj in sf.trajs:
        meta = base
        overridden = False
        side = sf.sidecars.get(traj)
        if side:
            meta = _read_json_member(zf, side)
            overridden = True
        ff = str(meta.get("force_field", base.get("force_field", "")))
        temp = int(round(float(meta.get("temperature_K",
                                        base.get("temperature_K", 0)))))
        save_ns = float(meta.get("save_traj_ns", base.get("save_traj_ns", 0.0)))
        key = (ff, temp)
        g = groups.get(key)
        if g is None:
            g = groups[key] = Group(slug=group_slug(ff, temp), force_field=ff,
                                    temperature_k=temp, save_traj_ns=save_ns)
        g.trajs.append(traj)
        g.overridden = g.overridden or overridden
    return [groups[k] for k in sorted(groups)]


# --------------------------------------------------------------------------- #
# XTC inspection
#
# Frame headers are walked directly rather than via mdtraj: we only need the atom
# count (to prove the trajectory matches the topology) and the frame count and
# spacing (for the description). Walking skips over each compressed block by its
# stored length, so no coordinates are decompressed and the cost is a seek per
# frame. Nothing is re-encoded, so mdtraj's silent xdrfile-overflow failure mode
# — the one mdcath_import.py has to guard against — cannot arise here.
# --------------------------------------------------------------------------- #

XTC_MAGIC = 1995
_XTC_HEADER = struct.Struct(">iiif")     # magic, natoms, step, time
_XTC_INT = struct.Struct(">i")


class BadTrajectoryError(RuntimeError):
    """An xtc is unreadable or does not match its topology."""


@dataclass
class XtcInfo:
    n_atoms: int
    n_frames: int
    ps_per_frame: Optional[float]     # None when the file holds a single frame


def xtc_scan(path: Path) -> XtcInfo:
    """Walk every frame header of an xtc. Raises BadTrajectoryError on garbage."""
    n_atoms = None
    n_frames = 0
    first_t = last_t = None
    size = path.stat().st_size
    with open(path, "rb") as fh:
        while True:
            head = fh.read(_XTC_HEADER.size)
            if len(head) < _XTC_HEADER.size:
                break
            magic, natoms, _step, t = _XTC_HEADER.unpack(head)
            if magic != XTC_MAGIC:
                raise BadTrajectoryError(
                    f"{path.name}: bad frame magic {magic} at byte "
                    f"{fh.tell() - _XTC_HEADER.size} (frame {n_frames})"
                )
            if n_atoms is None:
                n_atoms = natoms
            elif natoms != n_atoms:
                raise BadTrajectoryError(
                    f"{path.name}: atom count changes mid-file "
                    f"({n_atoms} -> {natoms} at frame {n_frames})"
                )
            fh.seek(36, os.SEEK_CUR)                       # 3x3 box
            raw = fh.read(_XTC_INT.size)
            if len(raw) < _XTC_INT.size:
                raise BadTrajectoryError(f"{path.name}: truncated at frame {n_frames}")
            (lsize,) = _XTC_INT.unpack(raw)
            if lsize != natoms:
                raise BadTrajectoryError(
                    f"{path.name}: coordinate block declares {lsize} atoms, "
                    f"header says {natoms} (frame {n_frames})"
                )
            if lsize <= 9:
                fh.seek(lsize * 3 * 4, os.SEEK_CUR)        # stored as raw floats
            else:
                fh.seek(4 + 4 * 6 + 4, os.SEEK_CUR)        # precision, min/max int, smallidx
                raw = fh.read(_XTC_INT.size)
                if len(raw) < _XTC_INT.size:
                    raise BadTrajectoryError(
                        f"{path.name}: truncated at frame {n_frames}")
                (nbytes,) = _XTC_INT.unpack(raw)
                if nbytes < 0:
                    raise BadTrajectoryError(
                        f"{path.name}: negative block length at frame {n_frames}")
                fh.seek(nbytes + ((4 - nbytes % 4) % 4), os.SEEK_CUR)
            if fh.tell() > size:
                raise BadTrajectoryError(f"{path.name}: truncated at frame {n_frames}")
            if first_t is None:
                first_t = t
            last_t = t
            n_frames += 1

    if not n_frames or n_atoms is None:
        raise BadTrajectoryError(f"{path.name}: no readable frames")
    dt = ((last_t - first_t) / (n_frames - 1)) if n_frames > 1 else None
    return XtcInfo(n_atoms=n_atoms, n_frames=n_frames, ps_per_frame=dt)


def pdb_atom_count(path: Path) -> int:
    n = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                n += 1
            elif line.startswith("ENDMDL"):
                break            # topology.pdb holds a single MODEL
    return n


# --------------------------------------------------------------------------- #
# Topology generation
# --------------------------------------------------------------------------- #

def write_psf_from_pdb(pdb_path: Path, dest_psf: Path) -> int:
    """Generate a PSF topology from the protein-only PDB, preserving exact atoms.

    BioEmu ships no topology file and MDRepo does not accept a .pdb as topology.
    ParmEd reads the (all-atom, hydrogen-containing) PDB and writes a PSF with
    the identical atom set and order, so it lines up with the xtc frame by frame.

    NOTE: a PSF built from a bare PDB carries connectivity and atom names but not
    force-field charges/types. If `mdr-process validate` rejects it, the fallback
    is an OpenMM-typed topology — see README.md ("verify on the VM").
    """
    try:
        import parmed as pmd
    except ImportError as e:      # pragma: no cover - environment problem
        raise RuntimeError(
            "parmed is required to generate the .psf topology (pip install parmed)"
        ) from e
    st = pmd.load_file(str(pdb_path))
    st.save(str(dest_psf), format="psf", overwrite=True)
    return len(st.atoms)


# --------------------------------------------------------------------------- #
# mdrepo-metadata.toml
# --------------------------------------------------------------------------- #

def toml_str(value: str) -> str:
    """A TOML basic string. Escapes what the spec requires; strips control chars."""
    out = []
    for ch in str(value):
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        elif ord(ch) < 0x20 or ord(ch) == 0x7F:
            continue
        else:
            out.append(ch)
    return '"' + "".join(out) + '"'


def toml_str_list(items) -> str:
    return "[" + ", ".join(toml_str(i) for i in items) + "]"


def truncate_short(text: str, limit: int = SHORT_DESCRIPTION_MAX) -> str:
    """Keep short_description inside the spec's character limit, at a word boundary."""
    if len(text) <= limit:
        return text
    cut = text[: limit - 1]
    space = cut.rfind(" ")
    if space > limit * 0.6:
        cut = cut[:space]
    return cut.rstrip(" ,;:.") + "."


def format_ns(ns: float) -> str:
    if ns >= 1000:
        return f"{ns / 1000:g} us"
    return f"{ns:g} ns"


def describe_system(ds: Dataset, system: str, ids: SystemIds,
                    uniprot_ids: list[str]) -> str:
    """The identity clause shared by short_description and description."""
    bits = []
    if ids.pdb_id:
        bits.append(f"PDB entry {ids.pdb_id}"
                    + (f" chain {ids.chain}" if ids.chain else ""))
    if uniprot_ids:
        bits.append("Uniprot entr" + ("y " if len(uniprot_ids) == 1 else "ies ")
                    + ", ".join(uniprot_ids))
    if ids.mutation:
        bits.append(f"point mutation {ids.mutation} of wild-type {ids.parent}")
    elif ids.variant:
        bits.append(f"variant {ids.variant}")
    elif ids.design:
        bits.append("a de novo designed sequence with no PDB entry")
    return "; ".join(bits)


def render_metadata(ds: Dataset, system: str, group: Group, ids: SystemIds,
                    uniprot_ids: list[str], traj_names: list[str],
                    pdb_name: str, psf_name: str, n_atoms: int, n_frames: int,
                    reference_name: Optional[str] = None,
                    dropped: Optional[list[tuple[str, str]]] = None) -> str:
    ff = FORCEFIELDS.get(group.force_field)
    ff_label = ff.label if ff else group.force_field
    ff_comments = ff.comments if ff else ""
    ident = describe_system(ds, system, ids, uniprot_ids)
    sampled_ns = n_frames * group.save_traj_ns

    short = (f"BioEmu {ds.label} {system}: {format_ns(sampled_ns)} of all-atom MD "
             f"in {len(traj_names)} trajector"
             f"{'y' if len(traj_names) == 1 else 'ies'} at {group.temperature_k} K, "
             f"{ff_label}")
    if ident:
        short += f" ({ident})"
    short += "."

    desc_parts = [
        f"All-atom molecular dynamics from the BioEmu training data release "
        f"(Lewis et al., Science 2025), dataset {ds.label}, system "
        f"{ds.zip_root}/{system}.",
        f"{ds.blurb}.",
    ]
    if ident:
        desc_parts.append(f"This system corresponds to {ident}.")
    summary = PROTOCOL_SUMMARY[ds.protocol].format(
        engine=ENGINE_NAME, ff=ff_label, temp=group.temperature_k,
        dt=INTEGRATION_TIMESTEP_FS,
    )
    # Only the first character: str.capitalize() would also lower-case the rest,
    # turning OpenMM into Openmm and a99SB-disp into A99sb-disp.
    desc_parts.append(summary[:1].upper() + summary[1:] + ".")
    desc_parts.append(
        f"This deposition holds the {len(traj_names)} trajector"
        f"{'y' if len(traj_names) == 1 else 'ies'} produced with {ff_label}: "
        f"{n_frames} frames total, {format_ns(group.save_traj_ns)} apart, "
        f"{format_ns(sampled_ns)} of sampled time, covering the "
        f"{n_atoms} solute (protein) atoms only (solvent stripped)."
    )
    if group.overridden:
        desc_parts.append(
            "The source system directory also contains trajectories produced "
            "with a different force field; those are deposited as a separate "
            "simulation, because a force field is a property of a simulation "
            "rather than of a file."
        )
    if ds.protocol == "megasim":
        desc_parts.append(
            "Production trajectories were run for 1.5 us per starting structure "
            "and the first 500 ns of each was discarded upstream, so the frames "
            "here are the retained portion."
            if ds.key == "megamerge" else
            "Each mutant was simulated for 1 us in the folded state; an upstream "
            "burn-in period was removed, so the frames here are the retained "
            "portion."
        )
    desc_parts.append(
        "Trajectory files carry their published names and content unchanged. "
        "The topology file is a PSF generated from the released topology.pdb "
        "(connectivity and atom naming only, no force-field parameters), since "
        "the release ships no topology in a format MDRepo accepts."
    )
    if reference_name:
        desc_parts.append(
            f"{reference_name} is the seed structure the simulations were "
            f"started from, deposited alongside the trajectory topology."
        )
    if dropped:
        desc_parts.append(
            "Note: " + str(len(dropped)) + " released trajector"
            + ("y was" if len(dropped) == 1 else "ies were")
            + " excluded because "
            + ("it does" if len(dropped) == 1 else "they do")
            + " not match the system topology or could not be read ("
            + "; ".join(f"{name}: {why}" for name, why in dropped) + ")."
        )
    desc = " ".join(desc_parts)

    # All bare (top-level) keys MUST precede any [table]/[[table]] header in TOML.
    lines = [
        f"lead_contributor_orcid = {toml_str(LEAD_CONTRIBUTOR_ORCID)}",
        f"trajectory_file_names = {toml_str_list(traj_names)}",
        f"structure_file_name = {toml_str(pdb_name)}",
        f"topology_file_name = {toml_str(psf_name)}",
        f"temperature_kelvin = {group.temperature_k}",
        f"integration_timestep_fs = {INTEGRATION_TIMESTEP_FS}",
        f"short_description = {toml_str(truncate_short(short))}",
        f"software_name = {toml_str(SOFTWARE_NAME)}",
        f"software_version = {toml_str(SOFTWARE_VERSION)}",
        f"description = {toml_str(desc)}",
        f"forcefield = {toml_str(ff_label)}",
    ]
    if ff_comments:
        lines.append(f"forcefield_comments = {toml_str(ff_comments)}")
    if ids.pdb_id:
        lines.append(f"pdb_id = {toml_str(ids.pdb_id)}")
    if uniprot_ids:
        lines.append(f"uniprot_ids = {toml_str_list(uniprot_ids)}")
    lines.append("")

    # [water] — model + density both required once the table is present. Recorded
    # only where the manuscript states the model (see FORCEFIELDS).
    if ff and ff.water:
        lines += [
            "[water]",
            f"model = {toml_str(ff.water)}",
            f"density_kg_m3 = {WATER_DENSITY_KG_M3}",
            "",
        ]

    if ds.protocol == "standard":
        for name, conc in STANDARD_SOLUTES:
            lines += ["[[solutes]]", f"name = {toml_str(name)}",
                      f"concentration_mol_liter = {conc}", ""]

    if reference_name:
        lines += [
            "[[additional_files]]",
            f"file_name = {toml_str(reference_name)}",
            'file_type = "Structure"',
            'file_description = "Seed structure the production runs were started from."',
            "",
        ]

    papers = [PAPER_BIOEMU]
    if ds.protocol == "megasim":
        papers.append(PAPER_MEGASCALE)
    for p in papers:
        lines.append("[[papers]]")
        lines.append(f"title = {toml_str(p['title'])}")
        lines.append(f"authors = {toml_str(p['authors'])}")
        lines.append(f"journal = {toml_str(p['journal'])}")
        lines.append(f"year = {p['year']}")
        lines.append(f"volume = {p['volume']}")
        if p.get("number"):
            lines.append(f"number = {toml_str(p['number'])}")
        if p.get("pages"):
            lines.append(f"pages = {toml_str(p['pages'])}")
        if p.get("doi"):
            lines.append(f"doi = {toml_str(p['doi'])}")
        lines.append("")

    lines += [
        "[[external_links]]",
        f"url = {toml_str(ZENODO_DOI.format(record=ds.record))}",
        f"label = {toml_str('BioEmu ' + ds.label + ' source archive (Zenodo)')}",
        "",
        "[[external_links]]",
        f"url = {toml_str(BIOEMU_REPO_URL)}",
        'label = "BioEmu model and inference code (GitHub)"',
        "",
    ]

    for c in CONTRIBUTORS:
        lines.append("[[contributors]]")
        lines.append(f"name = {toml_str(c['name'])}")
        if c.get("orcid"):
            lines.append(f"orcid = {toml_str(c['orcid'])}")
        if c.get("email"):
            lines.append(f"email = {toml_str(c['email'])}")
        if c.get("institution"):
            lines.append(f"institution = {toml_str(c['institution'])}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Building an IN_DIR
# --------------------------------------------------------------------------- #

@dataclass
class BuiltSim:
    sim_dir: Path
    group: Group
    n_trajs: int
    n_frames: int
    dropped: list[tuple[str, str]]


def _extract_member(zf: zipfile.ZipFile, member: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zf.open(member) as src, open(dest, "wb") as out:
        shutil.copyfileobj(src, out, DOWNLOAD_CHUNK)


def build_sim_dir(zf: zipfile.ZipFile, ds: Dataset, sf: SystemFiles, group: Group,
                  ids: SystemIds, uniprot_ids: list[str], sim_dir: Path,
                  on_note=None) -> BuiltSim:
    """Materialise one (system, force-field group) as an mdr-process IN_DIR."""
    system = sf.system
    _fresh_dir(sim_dir)

    if not sf.topology:
        raise SkipSystem(f"{system}: no topology.pdb in the archive")

    pdb_name = f"{system}.pdb"
    psf_name = f"{system}.psf"
    pdb_path = sim_dir / pdb_name
    _extract_member(zf, sf.topology, pdb_path)
    n_pdb_atoms = pdb_atom_count(pdb_path)
    if not n_pdb_atoms:
        raise RuntimeError(f"{system}: topology.pdb holds no ATOM records")

    n_psf_atoms = write_psf_from_pdb(pdb_path, sim_dir / psf_name)
    if n_psf_atoms != n_pdb_atoms:
        raise RuntimeError(
            f"{system}: generated psf has {n_psf_atoms} atoms, pdb has {n_pdb_atoms}"
        )

    reference_name = None
    if INCLUDE_REFERENCE_PDB and sf.reference:
        reference_name = f"{system}_reference.pdb"
        _extract_member(zf, sf.reference, sim_dir / reference_name)

    traj_names: list[str] = []
    dropped: list[tuple[str, str]] = []
    total_frames = 0
    for member in group.trajs:
        name = member.rsplit("/", 1)[-1]
        dest = sim_dir / name
        _extract_member(zf, member, dest)
        try:
            info = xtc_scan(dest)
        except BadTrajectoryError as e:
            log.error("[%s/%s] dropping %s — %s", ds.key, system, name, e)
            dropped.append((name, str(e).split(": ", 1)[-1]))
            dest.unlink(missing_ok=True)
            if on_note:
                on_note("unreadable-trajectory", f"{name}: {e}")
            continue
        if info.n_atoms != n_pdb_atoms:
            why = f"{info.n_atoms} atoms, topology has {n_pdb_atoms}"
            log.error("[%s/%s] dropping %s — %s", ds.key, system, name, why)
            dropped.append((name, why))
            dest.unlink(missing_ok=True)
            if on_note:
                on_note("atom-count-mismatch", f"{name}: {why}")
            continue
        # The published frame spacing should match dataset.json. MDRepo derives
        # sampling from the file itself, so a disagreement is recorded rather
        # than corrected — the file is the authority, dataset.json is the claim.
        if info.ps_per_frame is not None and group.save_traj_ns:
            expected = group.save_traj_ns * 1000.0
            if abs(info.ps_per_frame - expected) > max(1.0, 0.01 * expected):
                detail = (f"{name}: frames {info.ps_per_frame:g} ps apart, "
                          f"dataset.json declares {expected:g} ps")
                log.warning("[%s/%s] %s", ds.key, system, detail)
                if on_note:
                    on_note("frame-spacing-mismatch", detail)
        traj_names.append(name)
        total_frames += info.n_frames

    if not traj_names:
        raise SkipSystem(
            f"{system} [{group.slug}]: no usable trajectories"
            + (f" ({len(dropped)} dropped)" if dropped else "")
        )

    meta = render_metadata(
        ds, system, group, ids, uniprot_ids, traj_names,
        pdb_name=pdb_name, psf_name=psf_name, n_atoms=n_pdb_atoms,
        n_frames=total_frames, reference_name=reference_name, dropped=dropped,
    )
    # encoding is explicit: contributor names carry diacritics, and the locale
    # default is ASCII under LANG=C. TOML is UTF-8 by spec regardless.
    (sim_dir / "mdrepo-metadata.toml").write_text(meta, encoding="utf-8")
    return BuiltSim(sim_dir=sim_dir, group=group, n_trajs=len(traj_names),
                    n_frames=total_frames, dropped=dropped)


# --------------------------------------------------------------------------- #
# mdr-process invocation (watchdogged)
# --------------------------------------------------------------------------- #

class StalledCommandError(TransientError):
    """mdr-process made no progress and was killed by the watchdog."""


def _descendants(pid: int) -> list[int]:
    """Every pid in the subtree rooted at `pid`, via /proc (Linux only)."""
    children: dict[int, list[int]] = {}
    try:
        for entry in os.listdir("/proc"):
            if not entry.isdigit():
                continue
            try:
                with open(f"/proc/{entry}/stat", "rb") as fh:
                    data = fh.read().decode("utf-8", "replace")
                ppid = int(data[data.rindex(")") + 2:].split()[1])
            except (OSError, ValueError, IndexError):
                continue
            children.setdefault(ppid, []).append(int(entry))
    except OSError:
        return [pid]
    out, stack = [], [pid]
    while stack:
        cur = stack.pop()
        out.append(cur)
        stack.extend(children.get(cur, []))
    return out


def _subtree_progress(pid: int) -> int:
    """A monotonic counter of work done by the subtree: I/O bytes plus CPU ticks.

    A wedged gocmd push burns neither, so a constant value means "stuck" — where
    a value that keeps climbing means the command is merely slow.
    """
    total = 0
    for p in _descendants(pid):
        try:
            with open(f"/proc/{p}/io", "rb") as fh:
                for line in fh:
                    if line.startswith((b"read_bytes:", b"write_bytes:")):
                        total += int(line.split()[1])
        except (OSError, ValueError):
            pass
        try:
            with open(f"/proc/{p}/stat", "rb") as fh:
                fields = fh.read().decode("utf-8", "replace")
            rest = fields[fields.rindex(")") + 2:].split()
            total += (int(rest[11]) + int(rest[12])) * (1000 // max(_CLK_TCK, 1))
        except (OSError, ValueError, IndexError):
            pass
    return total


def _signal_pids(pids: list[int], sig: int) -> None:
    for p in pids:
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.kill(p, sig)


def _kill_subtree(proc: "subprocess.Popen") -> None:
    pids = _descendants(proc.pid)
    _signal_pids(pids, signal.SIGTERM)
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline and proc.poll() is None:
        time.sleep(0.5)
    if proc.poll() is None:
        _signal_pids(_descendants(proc.pid) or pids, signal.SIGKILL)
    with contextlib.suppress(Exception):
        proc.wait(timeout=15)


def _read_tail(log_path: Path) -> str:
    with contextlib.suppress(OSError):
        with open(log_path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            fh.seek(max(0, fh.tell() - ERROR_TAIL_CHARS))
            return fh.read().decode("utf-8", "replace")
    return ""


def run_mdr(argv: list[str], log_path: Path,
            stall_sec: float = STALL_MINUTES * 60,
            max_sec: float = MDR_MAX_HOURS * 3600,
            poll_sec: float = WATCHDOG_POLL_SEC) -> tuple[int, str]:
    """Run an mdr-process command, streaming output to log_path, watchdogged.

    Output is written live (not buffered until exit) so `tail -f` shows progress.
    A watchdog polls the process subtree: if it makes no progress for
    `stall_sec`, or runs past `max_sec`, the whole subtree is killed and
    StalledCommandError is raised. Set either to 0 to disable that guard.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    killed_reason = None
    with open(log_path, "ab", buffering=0) as fh:
        fh.write((">>> " + " ".join(argv) + "\n").encode())
        # No start_new_session: the child stays in our process group so a manual
        # Ctrl-C still reaches it and its gocmd descendants. The watchdog kills
        # the subtree explicitly rather than by group signal.
        proc = subprocess.Popen(argv, stdout=fh, stderr=subprocess.STDOUT)
        start = last_progress = time.monotonic()
        last_counter = _subtree_progress(proc.pid)
        while True:
            try:
                proc.wait(timeout=poll_sec)
                break
            except subprocess.TimeoutExpired:
                pass
            now = time.monotonic()
            counter = _subtree_progress(proc.pid)
            if counter != last_counter:
                last_counter, last_progress = counter, now
            if max_sec and (now - start) >= max_sec:
                killed_reason = f"exceeded {max_sec / 3600:.1f}h ceiling"
            elif stall_sec and (now - last_progress) >= stall_sec:
                killed_reason = (f"no subtree progress for "
                                 f"{(now - last_progress) / 60:.0f} min")
            if killed_reason:
                fh.write(f"\n>>> WATCHDOG: killing mdr-process subtree — "
                         f"{killed_reason}\n".encode())
                _kill_subtree(proc)
                break

    tail = _read_tail(log_path)
    if killed_reason is not None:
        raise StalledCommandError(
            f"mdr-process killed by watchdog ({killed_reason}):\n{tail}")
    return proc.returncode, tail


def _mdr_error(stage: str, rc: int, tail: str) -> RuntimeError:
    """A full disk or an unreachable MDRepo is not the system's fault."""
    cls = TransientError if _is_transient_text(tail) else RuntimeError
    return cls(f"{stage} failed (rc={rc}):\n{tail}")


def import_sim(cfg: "RunConfig", sim_dir: Path, log_path: Path, has_ids: bool,
               before_process=None) -> None:
    """validate + process one IN_DIR. Raises RuntimeError on failure."""
    base = [cfg.mdr_bin, "-l", cfg.log_level]
    if cfg.num_threads:
        base += ["-t", str(cfg.num_threads)]
    wd = dict(stall_sec=cfg.stall_sec, max_sec=cfg.max_sec)

    rc, tail = run_mdr(base + ["validate", str(sim_dir)], log_path, **wd)
    if rc != 0:
        raise _mdr_error("validate", rc, tail)

    proc_cmd = base + ["process", str(sim_dir), "-s", cfg.server, "-f"]
    if not has_ids:              # neither PDB nor UniProt id (octapeptides, designs)
        proc_cmd += ["--no-id"]
    if cfg.dry_run:
        proc_cmd += ["-d"]
    if cfg.work_dir:
        proc_cmd += ["-w", str(cfg.work_dir)]
    if cfg.out_dir:
        proc_cmd += ["-o", str(cfg.out_dir)]
    # blastp threads are per-process, so N workers each spawn their own search:
    # the real CPU load is workers x blast_num_threads, not this number alone.
    if cfg.blast_num_threads:
        proc_cmd += ["--blast-num-threads", str(cfg.blast_num_threads)]
    # Persist the intent immediately before the command that can mutate MDRepo.
    if before_process is not None:
        before_process()
    rc, tail = run_mdr(proc_cmd, log_path, **wd)
    if rc != 0:
        raise _mdr_error("process", rc, tail)


# --------------------------------------------------------------------------- #
# Per-system pipeline
# --------------------------------------------------------------------------- #

@dataclass
class RunConfig:
    root: Path
    server: str = "staging"
    mdr_bin: str = "mdr-process"
    num_threads: Optional[int] = None
    blast_num_threads: Optional[int] = None
    log_level: str = "info"
    dry_run: bool = False
    keep_on_success: bool = False
    keep_archives: bool = False
    work_dir: Optional[Path] = None
    out_dir: Optional[Path] = None
    datasets: tuple = ()                 # restrict to these dataset keys
    only_systems: Optional[tuple] = None
    workers: int = 1
    limit: int = 0
    min_free_gb: float = MIN_FREE_GB
    breaker_threshold: int = BREAKER_THRESHOLD
    verify_archives: bool = True
    stall_sec: float = STALL_MINUTES * 60
    max_sec: float = MDR_MAX_HOURS * 3600


def process_system(cfg: RunConfig, lay: Layout, man: Manifest, ds: Dataset,
                   zf: zipfile.ZipFile, index: dict[str, SystemFiles],
                   system: str) -> None:
    importing = man.importing_groups(ds.key, system)
    if importing:
        raise RuntimeError(
            f"ambiguous prior import for group(s) {', '.join(sorted(importing))}; "
            f"verify MDRepo, then run `resolve-import {ds.key} {system} GROUP "
            f"--imported` or `--retry`"
        )

    sf = index.get(system)
    if sf is None:
        raise RuntimeError(f"{system}: not present in {ds.zip_name}")
    if not sf.trajs:
        raise SkipSystem(f"{system}: the release ships no trajectories for this system")

    ids = parse_system_ids(ds, system)
    uniprot_ids = man.uniprot_for(ds.key, system)
    has_ids = bool(ids.pdb_id or uniprot_ids)

    groups = discover_groups(zf, sf)
    done = man.imported_groups(ds.key, system)
    pending = [g for g in groups if g.slug not in done]
    if not pending:
        _finalize(cfg, lay, man, ds, system)
        return

    note = lambda kind, detail: man.note(ds.key, system, kind, detail)  # noqa: E731
    failures, transient = [], []
    for group in pending:
        require_free_space(lay.staging, cfg.min_free_gb, f"{system} [{group.slug}]")
        sim_dir = lay.sim_dir(ds.key, system, group.slug)
        log_path = lay.log_file(ds.key, system, group.slug)
        try:
            built = build_sim_dir(zf, ds, sf, group, ids, uniprot_ids, sim_dir,
                                  on_note=note)
        except SkipSystem as e:
            log.warning("[%s/%s] %s", ds.key, system, e)
            man.note(ds.key, system, "empty-group", str(e))
            shutil.rmtree(sim_dir, ignore_errors=True)
            continue
        except Exception as e:  # noqa: BLE001 - continue with the other groups
            failures.append(f"{group.slug}: {e}")
            if _is_transient_exc(e):
                transient.append(group.slug)
            log.warning("[%s/%s] %s BUILD FAILED: %s", ds.key, system, group.slug, e)
            continue

        try:
            log.info("[%s/%s] importing %s (%d trajs, %d frames)",
                     ds.key, system, group.slug, built.n_trajs, built.n_frames)
            before = None
            if not cfg.dry_run:
                before = lambda g=group.slug, lp=log_path: man.mark_importing(
                    ds.key, system, g, str(lp))
            import_sim(cfg, sim_dir, log_path, has_ids, before_process=before)
            if not cfg.dry_run:
                man.mark_imported(ds.key, system, group.slug, str(log_path))
                if not cfg.keep_on_success:
                    # Reclaim disk immediately: with 21458 mutant systems the
                    # staging tree would otherwise grow for the whole run.
                    shutil.rmtree(sim_dir, ignore_errors=True)
        except Exception as e:  # noqa: BLE001
            failures.append(f"{group.slug}: {e}")
            if _is_transient_exc(e):
                transient.append(group.slug)
            first = str(e).splitlines()[0] if str(e) else e
            log.warning("[%s/%s] %s FAILED: %s", ds.key, system, group.slug, first)

    if failures:
        # One transient failure taints the system: it is retried as a whole, and
        # already-imported groups are skipped by imported_groups().
        cls = TransientError if transient else RuntimeError
        raise cls("; ".join(failures))
    if cfg.dry_run:
        log.info("[%s/%s] dry-run complete; manifest and staged inputs retained",
                 ds.key, system)
        return
    if not man.imported_groups(ds.key, system):
        raise SkipSystem(f"{system}: no group yielded a usable simulation")
    _finalize(cfg, lay, man, ds, system)


def _finalize(cfg: RunConfig, lay: Layout, man: Manifest, ds: Dataset,
              system: str) -> None:
    if not cfg.keep_on_success:
        shutil.rmtree(lay.system_stage(ds.key, system), ignore_errors=True)
    man.mark_done(ds.key, system)
    log.info("[%s/%s] DONE", ds.key, system)


# --------------------------------------------------------------------------- #
# Locking / workers
# --------------------------------------------------------------------------- #

@contextlib.contextmanager
def try_system_lock(lay: Layout, dataset: str, system: str):
    """Non-blocking per-system flock. Yields True if this process claimed it."""
    path = lay.lock_file(dataset, system)
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(path, "w")
    try:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            yield False
            return
        yield True
    finally:
        with contextlib.suppress(Exception):
            fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()


class Breaker:
    """Consecutive-failure counter shared by every worker.

    A single outage (MDRepo down, volume full, Zenodo throttling) fails whichever
    system each worker happens to hold, over and over. Tripping stop_evt halts
    the run instead of grinding the outage into the manifest.

    `counter` is an mp.Value rather than a Manager proxy: Manager objects are not
    picklable, and a worker is started with the spawn context, which pickles
    everything passed to Process(args=...).
    """

    def __init__(self, counter, stop_evt, threshold: int):
        self.counter = counter          # mp.Value('i'), has its own lock
        self.stop_evt = stop_evt
        self.threshold = threshold

    def success(self) -> None:
        if self.threshold:
            with self.counter.get_lock():
                self.counter.value = 0

    def failure(self, what: str = "") -> bool:
        """Returns True if this failure tripped the breaker."""
        if not self.threshold:
            return False
        with self.counter.get_lock():
            self.counter.value += 1
            n = self.counter.value
        if n >= self.threshold:
            log.error(
                "circuit breaker: %d consecutive failures (last: %s) — stopping all "
                "workers. Investigate, then re-run (transient failures kept their "
                "attempt budget).", n, what,
            )
            self.stop_evt.set()
            return True
        return False


class Limiter:
    """Run-wide cap on how many systems are attempted, shared by every worker.

    A worker reserves a slot *before* it commits to a system; once `limit` slots
    are taken the run stops claiming work. Every attempt counts — success, hard
    failure, or transient requeue alike — so a bounded first-batch run touches
    exactly `limit` systems rather than overshooting by up to N.
    """

    def __init__(self, counter, stop_evt, limit: int):
        self.counter = counter          # mp.Value('i'), has its own lock
        self.stop_evt = stop_evt
        self.limit = limit

    def reserve(self) -> bool:
        if self.limit <= 0:
            return True
        with self.counter.get_lock():
            if self.counter.value >= self.limit:
                self.stop_evt.set()
                return False
            self.counter.value += 1
            reached = self.counter.value >= self.limit
        if reached:
            log.info("--limit reached (%d systems attempted); workers will stop",
                     self.limit)
            self.stop_evt.set()
        return True


def filesystem_type(path: Path) -> str:
    """Best-effort filesystem name for `path` (used in the flock refusal message)."""
    try:
        target = str(Path(path).resolve())
        best, best_type = "", "unknown"
        with open("/proc/self/mounts", "r", errors="replace") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 3:
                    continue
                mount, fstype = parts[1], parts[2]
                if (target == mount or target.startswith(mount.rstrip("/") + "/")) \
                        and len(mount) > len(best):
                    best, best_type = mount, fstype
        return best_type
    except OSError:
        return "unknown"


def _hold_lock_probe(path: str, acquired, release_evt) -> None:
    """Child half of the flock self-test: take the lock, hold it, then release."""
    with open(path, "w") as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return
        acquired.set()
        release_evt.wait(30)


def verify_domain_locking(lay: Layout) -> None:
    """Prove flock() actually excludes a second process before starting workers.

    Some FUSE filesystems accept flock() and silently no-op; NFS mounted `nolock`
    degrades it the same way. Either would let two workers import the same system
    twice. Rather than guess from the mount type, contend for real.
    """
    probe = lay.locks / ".locktest"
    ctx = mp.get_context("spawn")
    acquired, release_evt = ctx.Event(), ctx.Event()
    child = ctx.Process(target=_hold_lock_probe,
                        args=(str(probe), acquired, release_evt))
    child.start()
    try:
        if not acquired.wait(30):
            raise RuntimeError("lock self-test: child never acquired the lock")
        with open(probe, "w") as fh:
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                return          # correct: the child's lock excluded us
            fcntl.flock(fh, fcntl.LOCK_UN)
        raise RuntimeError(
            f"flock() on {lay.locks} does NOT exclude a second process "
            f"(filesystem: {filesystem_type(lay.locks)}). Two workers could "
            f"import the same system twice. Move --root to a local filesystem "
            f"(ext4/xfs/btrfs), or run with -w 1."
        )
    finally:
        release_evt.set()
        child.join(30)
        if child.is_alive():
            child.terminate()
        with contextlib.suppress(OSError):
            probe.unlink()


def worker_loop(cfg: RunConfig, ds_key: str, worker_id: int, stop_evt,
                breaker=None, limiter=None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s w{worker_id} %(levelname)s %(message)s",
    )
    lay = Layout(cfg.root)
    ds = DATASETS[ds_key]
    man = Manifest(lay.db)
    zf = zipfile.ZipFile(lay.archive(ds))
    index = index_archive(zf, ds)
    if breaker is None:   # direct/testing invocation without a shared counter
        breaker = Breaker(mp.Value("i", 0), stop_evt, cfg.breaker_threshold)
    if limiter is None:
        limiter = Limiter(mp.Value("i", 0), stop_evt, cfg.limit)
    idle_rounds = 0
    contended_rounds = 0
    # A successful dry-run deliberately leaves the system pending, so selecting
    # from the manifest would hand back the same system forever. Snapshot the
    # work list instead, so each is checked exactly once.
    dry_run_systems = iter(
        man.candidate_systems(ds.key, limit=2_147_483_647, only=cfg.only_systems)
    ) if cfg.dry_run else None
    try:
        while not stop_evt.is_set():
            claimed = contended = backoff = False
            if dry_run_systems is not None:
                system = next(dry_run_systems, None)
                if system is None:
                    log.info("[%s] dry-run work list complete; worker exiting", ds.key)
                    break
                candidates = (system,)
            else:
                candidates = man.candidate_systems(ds.key, limit=64,
                                                   only=cfg.only_systems)
            for system in candidates:
                if stop_evt.is_set():
                    break
                with try_system_lock(lay, ds.key, system) as got:
                    if not got:
                        contended = True     # a peer owns it; try the next one
                        continue
                    # Re-read under the lock: a peer may have finished it since
                    # the candidate query ran.
                    row = man.system_row(ds.key, system)
                    if row is None or row["status"] in ("done", "skipped"):
                        continue
                    if row["status"] == "failed" and row["attempts"] >= MAX_ATTEMPTS:
                        continue
                    # Reserve a --limit slot only once we hold the lock and know
                    # the system is really ours, so the cap counts work attempted
                    # rather than candidates looked at.
                    if not limiter.reserve():
                        break
                    claimed = True
                    try:
                        process_system(cfg, lay, man, ds, zf, index, system)
                        breaker.success()
                    except SkipSystem as e:
                        # Nothing importable was released for this system. That
                        # is a property of the data, not a failure to retry, and
                        # it must not trip the breaker.
                        log.info("[%s/%s] SKIP: %s", ds.key, system, e)
                        man.mark_skipped(ds.key, system, str(e))
                        shutil.rmtree(lay.system_stage(ds.key, system),
                                      ignore_errors=True)
                        breaker.success()
                    except Exception as e:  # noqa: BLE001
                        shutil.rmtree(lay.system_stage(ds.key, system),
                                      ignore_errors=True)
                        first = str(e).splitlines()[0] if str(e) else repr(e)
                        if _is_transient_exc(e):
                            # Infrastructure, not this system: requeue without
                            # spending an attempt, then back off so a full disk
                            # or a down server isn't hammered in a tight loop.
                            man.requeue(ds.key, system)
                            log.warning("[%s/%s] TRANSIENT (not counted): %s",
                                        ds.key, system, first)
                            backoff = True
                        else:
                            man.mark_failed(ds.key, system,
                                            f"{type(e).__name__}: {e}")
                            log.warning("[%s/%s] FAILED: %s", ds.key, system, first)
                        breaker.failure(f"{ds.key}/{system}")
                    break     # re-query candidates after finishing one system
            if backoff:
                # Outside the flock: don't hold a requeued system hostage while
                # waiting for the disk or the server to come back.
                stop_evt.wait(TRANSIENT_BACKOFF_SEC)
            if claimed:
                idle_rounds = contended_rounds = 0
            elif contended:
                # Candidates exist but peers hold every lock. Waiting is correct
                # while those peers are alive — systems leave the candidate set
                # as they complete — but don't wait forever on a stuck lock.
                contended_rounds += 1
                if contended_rounds >= CONTENTION_ROUNDS:
                    log.info("[%s] all candidates locked by peers for too long; "
                             "worker exiting", ds.key)
                    break
                stop_evt.wait(5)
            else:
                # Nothing claimable. With peers running this may be temporary —
                # a peer's transient failure requeues a system — so wait a few
                # rounds before giving up. Alone, no new work can appear, so
                # waiting would only stall the next dataset.
                idle_rounds += 1
                if cfg.workers <= 1 or idle_rounds >= 3:
                    log.info("[%s] no claimable systems; worker exiting", ds.key)
                    break
                stop_evt.wait(5)
    finally:
        with contextlib.suppress(Exception):
            zf.close()
        man.close()


# --------------------------------------------------------------------------- #
# init: work-list + SIFTS UniProt map
# --------------------------------------------------------------------------- #

class _RemoteFile(io.RawIOBase):
    """Range-request-backed file object, so a zip's central directory can be read
    without downloading the whole archive. Used only by `init`."""

    def __init__(self, url: str, size: int):
        self.url, self.size, self.pos = url, size, 0

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, off, whence=0):
        self.pos = (off if whence == 0 else
                    self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n <= 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        last = None
        for attempt in range(DOWNLOAD_RETRIES):
            try:
                with _urlopen(self.url, {"Range": f"bytes={self.pos}-{end}"}) as r:
                    data = r.read()
                break
            except Exception as e:  # noqa: BLE001 - Zenodo throttles ranged reads
                last = e
                time.sleep(2 * (attempt + 1))
        else:
            raise TransientError(f"ranged read failed: {last}")
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def open_archive_index(lay: Layout, ds: Dataset) -> zipfile.ZipFile:
    """Open ds's zip — locally if present, else over HTTP range requests.

    `init` only needs the central directory (a few MB), so a work-list can be
    built before any of the 41 GB is downloaded.
    """
    local = lay.archive(ds)
    if local.exists():
        return zipfile.ZipFile(local)
    url, size, _md5 = zenodo_file_info(ds)
    log.info("[%s] reading the central directory of %s over HTTP", ds.key, ds.zip_name)
    return zipfile.ZipFile(io.BufferedReader(_RemoteFile(url, size),
                                             buffer_size=1 << 16))


def build_sifts_map(wanted: dict[str, set[Optional[str]]]) -> dict[tuple[str, Optional[str]], set[str]]:
    """(pdb, chain) -> UniProt accessions, from the SIFTS flat file.

    `wanted` maps a lowercase PDB accession to the chains asked for; a None chain
    means "union across every chain of this entry", which is what the MEGAscale
    entry names (a domain, not a chain) call for.
    """
    log.info("downloading SIFTS pdb_chain_uniprot (%d PDB entries wanted)", len(wanted))
    out: dict[tuple[str, Optional[str]], set[str]] = {}
    with _urlopen(SIFTS_URL, timeout=300) as resp, gzip.GzipFile(fileobj=resp) as gz:
        for raw in io.TextIOWrapper(gz, encoding="utf-8", errors="replace"):
            if raw.startswith("#") or raw.startswith("PDB"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            pdb, chain, acc = parts[0].lower(), parts[1], parts[2]
            chains = wanted.get(pdb)
            if chains is None:
                continue
            if None in chains:
                out.setdefault((pdb, None), set()).add(acc)
            if chain in chains:
                out.setdefault((pdb, chain), set()).add(acc)
    return out


def cmd_init(args) -> None:
    lay = Layout(Path(args.root))
    lay.ensure()
    man = Manifest(lay.db)
    keys = _selected_datasets(args)

    wanted: dict[str, set[Optional[str]]] = {}
    system_ids: list[tuple[str, str, SystemIds]] = []
    for key in keys:
        ds = DATASETS[key]
        zf = open_archive_index(lay, ds)
        try:
            index = index_archive(zf, ds)
        finally:
            with contextlib.suppress(Exception):
                zf.close()
        rows = [(s, len(sf.trajs)) for s, sf in sorted(index.items())]
        added = man.add_systems(key, rows)
        empty = sum(1 for _s, n in rows if n == 0)
        log.info("[%s] %d systems in %s (%d new, %d ship no trajectories)",
                 key, len(rows), ds.zip_name, added, empty)
        for system, _n in rows:
            try:
                ids = parse_system_ids(ds, system)
            except ValueError as e:
                log.warning("[%s] %s", key, e)
                continue
            system_ids.append((key, system, ids))
            if ids.pdb_id:
                wanted.setdefault(ids.pdb_id, set()).add(ids.chain)

    if args.skip_sifts or not wanted:
        log.info("SIFTS map skipped")
    else:
        sifts = build_sifts_map(wanted)
        rows = []
        for key, system, ids in system_ids:
            if not ids.pdb_id:
                continue
            for acc in sorted(sifts.get((ids.pdb_id, ids.chain), ())):
                rows.append((key, system, acc))
        man.clear_uniprot()
        man.set_uniprot(rows)
        log.info("SIFTS: %d (system, uniprot) pairs for %d PDB entries",
                 len(rows), len(wanted))

    man.close()
    print("init complete. Next: run a single-system dry-run, e.g.")
    first = keys[0]
    print(f"  python {Path(sys.argv[0]).name} --root {args.root} run "
          f"--datasets {first} --systems <SYSTEM> -w 1 --dry-run")


# --------------------------------------------------------------------------- #
# Subcommands
# --------------------------------------------------------------------------- #

def _selected_datasets(args) -> list[str]:
    raw = getattr(args, "datasets", None)
    if not raw:
        return list(DATASET_ORDER)
    keys: list[str] = []
    bad: list[str] = []
    for value in raw:
        for tok in value.split(","):
            tok = tok.strip()
            if not tok:
                continue
            if tok not in DATASETS:
                bad.append(tok)
            elif tok not in keys:
                keys.append(tok)
    if bad:
        sys.exit("unknown dataset(s): %s (choose from %s)"
                 % (", ".join(bad), ", ".join(DATASET_ORDER)))
    if not keys:
        sys.exit("--datasets given with no usable names")
    return [k for k in DATASET_ORDER if k in keys]


def _normalize_systems(values: Optional[list[str]]) -> Optional[tuple]:
    if not values:
        return None
    out: list[str] = []
    for value in values:
        for tok in value.split(","):
            tok = tok.strip()
            if tok and tok not in out:
                out.append(tok)
    if not out:
        sys.exit("--systems given with no usable names")
    return tuple(out)


def _config_from_args(args) -> RunConfig:
    return RunConfig(
        root=Path(args.root),
        server=args.server,
        mdr_bin=args.mdr_bin,
        num_threads=args.num_threads,
        blast_num_threads=args.blast_num_threads,
        log_level=args.mdr_log_level,
        dry_run=args.dry_run,
        keep_on_success=args.keep,
        keep_archives=args.keep_archives,
        work_dir=Path(args.work_dir) if args.work_dir else None,
        out_dir=Path(args.out_dir) if args.out_dir else None,
        datasets=tuple(_selected_datasets(args)),
        only_systems=_normalize_systems(getattr(args, "systems", None)),
        workers=getattr(args, "workers", 1),
        limit=args.limit,
        min_free_gb=args.min_free_gb,
        breaker_threshold=args.breaker_threshold,
        verify_archives=not args.no_verify,
        stall_sec=args.stall_minutes * 60,
        max_sec=args.mdr_max_hours * 3600,
    )


def cmd_run(args) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    lay = Layout(Path(args.root))
    lay.ensure()
    cfg = _config_from_args(args)

    if cfg.dry_run and args.workers != 1:
        log.warning("dry-run forces a single worker (a dry-run leaves systems "
                    "pending, so peers cannot tell what has been checked)")
        args.workers = cfg.workers = 1
    if args.workers > 1 and not args.skip_lock_check:
        try:
            verify_domain_locking(lay)
        except RuntimeError as e:
            sys.exit(f"refusing to start {args.workers} workers: {e}")

    man = Manifest(lay.db)
    amb = man.ambiguous()
    if amb:
        print(f"⚠ {len(amb)} ambiguous import(s) awaiting resolve-import "
              f"(verify in MDRepo, then --imported or --retry):")
        for r in amb[:10]:
            print(f"  {r['dataset']}/{r['system']} [{r['grp']}]  "
                  f"(since {r['started_at']})")

    if cfg.only_systems:
        for key in cfg.datasets:
            for s in cfg.only_systems:
                if man.system_row(key, s):
                    man.reset_system(key, s)

    ctx = mp.get_context("spawn")
    requested = list(cfg.datasets)
    for key in requested:
        ds = DATASETS[key]
        if not man.pending_count(key, cfg.only_systems):
            log.info("[%s] nothing pending", key)
            continue
        try:
            with archive_ready(lay, ds, man, min_free_gb=cfg.min_free_gb,
                               verify=cfg.verify_archives):
                pass
        except Exception as e:  # noqa: BLE001
            log.error("[%s] archive unavailable: %s", key, e)
            if not _is_transient_exc(e):
                raise
            continue

        stop_evt = ctx.Event()
        breaker = Breaker(ctx.Value("i", 0), stop_evt, cfg.breaker_threshold)
        limiter = Limiter(ctx.Value("i", 0), stop_evt, cfg.limit)
        procs = [ctx.Process(target=worker_loop,
                             args=(cfg, key, i, stop_evt, breaker, limiter))
                 for i in range(args.workers)]
        log.info("[%s] starting %d worker(s) on %s", key, args.workers, ds.zip_name)
        try:
            for p in procs:
                p.start()
            for p in procs:
                p.join()
        except KeyboardInterrupt:
            log.warning("interrupted — stopping workers")
            stop_evt.set()
            for p in procs:
                p.join(60)
            raise

        remaining = man.pending_count(key, cfg.only_systems)
        stuck = man.counts().get(key, {}).get("failed", 0)
        if remaining:
            log.warning("[%s] %d system(s) still pending", key, remaining)
        elif stuck:
            # Systems that burned their attempt budget are exactly the ones a
            # later `reset-failed` will retry. Deleting the archive now would
            # make that retry re-download up to 28 GB.
            log.warning("[%s] %d system(s) failed; keeping %s so `reset-failed` "
                        "can retry without re-downloading", key, stuck, ds.zip_name)
        elif not cfg.keep_archives and not cfg.dry_run and not cfg.only_systems:
            man.set_archive_state(key, "drained")
            with contextlib.suppress(OSError):
                lay.archive(ds).unlink()
                log.info("[%s] archive drained and deleted", key)

    _print_status(man, show_failures=0)
    man.close()


def cmd_fetch(args) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    lay = Layout(Path(args.root))
    lay.ensure()
    man = Manifest(lay.db)
    for key in _selected_datasets(args):
        ds = DATASETS[key]
        with archive_ready(lay, ds, man, min_free_gb=args.min_free_gb,
                           verify=not args.no_verify) as path:
            print(f"{key}: {path} ({path.stat().st_size / 1024**3:.2f} GB)")
    man.close()


def _print_status(man: Manifest, show_failures: int) -> None:
    counts = man.counts()
    header = f"{'dataset':<12}{'total':>8}{'done':>8}{'skipped':>9}{'failed':>8}{'pending':>9}{'sims':>8}"
    print(header)
    print("-" * len(header))
    tot = dict(total=0, done=0, skipped=0, failed=0, pending=0, sims_imported=0)
    for key in DATASET_ORDER:
        c = counts.get(key)
        if not c:
            continue
        total = sum(v for k, v in c.items() if k != "sims_imported")
        row = dict(total=total, done=c.get("done", 0), skipped=c.get("skipped", 0),
                   failed=c.get("failed", 0), pending=c.get("pending", 0),
                   sims_imported=c.get("sims_imported", 0))
        for k in tot:
            tot[k] += row[k]
        print(f"{key:<12}{row['total']:>8}{row['done']:>8}{row['skipped']:>9}"
              f"{row['failed']:>8}{row['pending']:>9}{row['sims_imported']:>8}")
    print("-" * len(header))
    print(f"{'TOTAL':<12}{tot['total']:>8}{tot['done']:>8}{tot['skipped']:>9}"
          f"{tot['failed']:>8}{tot['pending']:>9}{tot['sims_imported']:>8}")

    amb = man.ambiguous()
    if amb:
        print(f"\n⚠ {len(amb)} ambiguous import(s) awaiting resolve-import "
              f"(verify in MDRepo, then --imported or --retry):")
        for r in amb[:20]:
            print(f"  {r['dataset']}/{r['system']} [{r['grp']}]  "
                  f"(since {r['started_at']})")

    notes = man.note_rows(limit=0)
    if notes:
        kinds: dict[str, int] = {}
        for r in notes:
            kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
        print("\ndata notes: " + ", ".join(f"{k}={v}" for k, v in sorted(kinds.items()))
              + "   (see `status --show-notes N`)")

    if show_failures:
        rows = man.failures(show_failures)
        if rows:
            print(f"\nlast {len(rows)} failure(s):")
            for r in rows:
                first = (r["error"] or "").strip().splitlines()
                print(f"  {r['dataset']}/{r['system']} (attempts {r['attempts']}): "
                      f"{first[0] if first else ''}")


def cmd_status(args) -> None:
    lay = Layout(Path(args.root))
    man = Manifest(lay.db)
    _print_status(man, args.show_failures)
    if args.show_notes:
        rows = man.note_rows(limit=args.show_notes)
        if rows:
            print(f"\nfirst {len(rows)} data note(s):")
            for r in rows:
                print(f"  [{r['kind']}] {r['dataset']}/{r['system']}: {r['detail']}")
    man.close()


def cmd_extract(args) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    lay = Layout(Path(args.root))
    lay.ensure()
    man = Manifest(lay.db)
    ds = DATASETS[args.dataset] if args.dataset in DATASETS else None
    if ds is None:
        sys.exit(f"unknown dataset {args.dataset!r} (choose from "
                 f"{', '.join(DATASET_ORDER)})")
    with archive_ready(lay, ds, man, min_free_gb=args.min_free_gb,
                       verify=not args.no_verify) as archive:
        zf = zipfile.ZipFile(archive)
        try:
            index = index_archive(zf, ds)
            sf = index.get(args.system)
            if sf is None:
                sys.exit(f"{args.system!r} is not in {ds.zip_name}")
            if not sf.trajs:
                sys.exit(f"{args.system}: the release ships no trajectories")
            ids = parse_system_ids(ds, args.system)
            uniprot_ids = man.uniprot_for(ds.key, args.system)
            out_root = Path(args.out_dir)
            for group in discover_groups(zf, sf):
                sim_dir = out_root / ds.key / args.system / group.slug
                built = build_sim_dir(zf, ds, sf, group, ids, uniprot_ids, sim_dir)
                print(f"{sim_dir}  ({built.n_trajs} trajectories, "
                      f"{built.n_frames} frames"
                      + (f", {len(built.dropped)} dropped" if built.dropped else "")
                      + ")")
        finally:
            zf.close()
    man.close()


def cmd_reset_failed(args) -> None:
    lay = Layout(Path(args.root))
    man = Manifest(lay.db)
    keys = _selected_datasets(args)
    n = sum(man.reset_failed(k) for k in keys) if args.datasets else man.reset_failed()
    print(f"requeued {n} failed system(s)")
    man.close()


def cmd_reset_sim(args) -> None:
    lay = Layout(Path(args.root))
    man = Manifest(lay.db)
    n = man.reset_sims(args.dataset, args.system, args.groups or None)
    man.reset_system(args.dataset, args.system)
    print(f"forgot {n} import record(s) for {args.dataset}/{args.system}; "
          f"system requeued")
    man.close()


def cmd_resolve_import(args) -> None:
    lay = Layout(Path(args.root))
    man = Manifest(lay.db)
    if args.imported:
        man.mark_imported(args.dataset, args.system, args.group, "")
        print(f"recorded {args.dataset}/{args.system} [{args.group}] as imported")
    else:
        if man.retry_import(args.dataset, args.system, args.group):
            print(f"cleared the ambiguous record for {args.dataset}/{args.system} "
                  f"[{args.group}]; it will be redone")
        else:
            print("no ambiguous record found (already resolved?)")
    man.reset_system(args.dataset, args.system)
    man.close()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _add_mdr_args(p) -> None:
    p.add_argument("-s", "--server", default="staging",
                   choices=["staging", "prod"], help="MDRepo target")
    p.add_argument("--mdr-bin", default="mdr-process")
    p.add_argument("--mdr-log-level", default="info")
    # No short flags here: -t/-w/-o belong to mdr-process, and -w is this tool's
    # worker count (as in the sibling importers). These are forwarded, not ours.
    p.add_argument("--num-threads", type=int, default=None,
                   help="forwarded to mdr-process -t")
    p.add_argument("--blast-num-threads", type=int, default=None,
                   help="forwarded to mdr-process process --blast-num-threads")
    p.add_argument("--work-dir", default=None,
                   help="forwarded to mdr-process -w")
    p.add_argument("--out-dir", default=None,
                   help="forwarded to mdr-process -o")
    p.add_argument("-d", "--dry-run", action="store_true",
                   help="pass -d to mdr-process (build import JSON, no push)")
    p.add_argument("--stall-minutes", type=float, default=STALL_MINUTES,
                   help="kill mdr-process after this long with no subtree progress (0 disables)")
    p.add_argument("--mdr-max-hours", type=float, default=MDR_MAX_HOURS,
                   help="absolute per-command ceiling in hours (0 disables)")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bioemu_import.py",
        description="Import the BioEmu MD dataset release (Zenodo) into MDRepo.",
    )
    p.add_argument("--root", required=True,
                   help="working directory for the manifest, archives and staging")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("init", help="build the work-list + SIFTS UniProt map")
    q.add_argument("--datasets", action="append",
                   help="restrict to these datasets (comma-separated; repeatable)")
    q.add_argument("--skip-sifts", action="store_true",
                   help="do not download the SIFTS map (no uniprot_ids)")
    q.set_defaults(func=cmd_init)

    q = sub.add_parser("fetch", help="download archives without importing")
    q.add_argument("--datasets", action="append")
    q.add_argument("--min-free-gb", type=float, default=MIN_FREE_GB)
    q.add_argument("--no-verify", action="store_true",
                   help="skip the md5 check against Zenodo")
    q.set_defaults(func=cmd_fetch)

    q = sub.add_parser("run", help="import pending systems")
    q.add_argument("--datasets", action="append",
                   help="restrict to these datasets (comma-separated; repeatable). "
                        "Archives are processed smallest-first, one at a time.")
    q.add_argument("--systems", action="append",
                   help="restrict to these system names (comma-separated; repeatable)")
    q.add_argument("-w", "--workers", type=int, default=1, dest="workers")
    q.add_argument("--limit", type=int, default=0,
                   help="stop after this many systems are attempted per dataset")
    q.add_argument("--keep", action="store_true",
                   help="do not delete staged files after success")
    q.add_argument("--keep-archives", action="store_true",
                   help="do not delete an archive once its systems are drained")
    q.add_argument("--min-free-gb", type=float, default=MIN_FREE_GB)
    q.add_argument("--breaker-threshold", type=int, default=BREAKER_THRESHOLD,
                   help="halt all workers after N consecutive failures (0 disables)")
    q.add_argument("--no-verify", action="store_true",
                   help="skip the md5 check against Zenodo")
    q.add_argument("--skip-lock-check", action="store_true",
                   help="skip the flock exclusion self-test that guards -w > 1. Don't.")
    _add_mdr_args(q)
    q.set_defaults(func=cmd_run)

    q = sub.add_parser("status", help="print progress")
    q.add_argument("--show-failures", type=int, default=0)
    q.add_argument("--show-notes", type=int, default=0)
    q.set_defaults(func=cmd_status)

    q = sub.add_parser("extract", help="stage ONE system for inspection (no import)")
    q.add_argument("dataset", help=", ".join(DATASET_ORDER))
    q.add_argument("system")
    q.add_argument("--out-dir", required=True)
    q.add_argument("--min-free-gb", type=float, default=MIN_FREE_GB)
    q.add_argument("--no-verify", action="store_true")
    q.set_defaults(func=cmd_extract)

    q = sub.add_parser("reset-failed", help="requeue failed systems")
    q.add_argument("--datasets", action="append")
    q.set_defaults(func=cmd_reset_failed)

    q = sub.add_parser("reset-sim",
                       help="forget import records so a system is redone")
    q.add_argument("dataset")
    q.add_argument("system")
    q.add_argument("groups", nargs="*",
                   help="force-field group slugs; omit for all")
    q.set_defaults(func=cmd_reset_sim)

    q = sub.add_parser("resolve-import",
                       help="record the outcome of an interrupted push")
    q.add_argument("dataset")
    q.add_argument("system")
    q.add_argument("group")
    g = q.add_mutually_exclusive_group(required=True)
    g.add_argument("--imported", action="store_true",
                   help="it landed in MDRepo; record it as imported")
    g.add_argument("--retry", action="store_true",
                   help="it is absent from MDRepo; clear the record so it is redone")
    q.set_defaults(func=cmd_resolve_import)
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    args.func(args)


if __name__ == "__main__":
    main()
