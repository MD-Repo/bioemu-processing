# BioEmu → MDRepo import

Download the [BioEmu MD dataset release](https://github.com/microsoft/bioemu)
(Lewis et al., *Science* **389**, eadv9817, 2025) from Zenodo and import it into
MDRepo via `mdr-process`. Four of the release's five archives are in scope;
**`opep` is deliberately not imported** — see [Scope](#scope).

This is the sibling of `mdcath_import.py` (whose unit was a 5-replica
`(domain, temperature)` group) and `dynamicpdb_import.py` (one trajectory per
protein). BioEmu sits in between: **one MDRepo simulation = one
(system, force-field group)**, where a group is all of a system's trajectories
sharing a force field and temperature.

## What gets imported

Four zip archives across three Zenodo records, **39.5 GB compressed**:

| Dataset | Archive | Size | Systems | Trajectories | Force field / T | Simulations |
|---|---|---:|---:|---:|---|---:|
| `megamerge` | `MSR_megasim_merge.zip` | 1.3 GB | 271 | 3,789 | **mixed** / 295 K | **348** |
| `cath1` | `ONE_cath1.zip` | 2.1 GB | 50 | 1,248 | ff99SB-ILDN / 300 K | 50 |
| `megamut` | `MSR_megasim_mutants_disp_allatom.zip` | 8.4 GB | 21,458 | 21,458 | a99SB-disp / 295 K | 21,458 |
| `cath2` | `MSR_cath2.zip` | 27.7 GB | 1,043 | 40,984 | ff99SB-ILDN / 300 K | **1,040** |
| | | **39.5 GB** | **22,822** | **67,479** | | **22,896** |
| ~~`opep`~~ | ~~`ONE_octapeptides.zip`~~ | 0.5 GB | 1100 | 118,252 | ff99SB-ILDN / 300 K | *out of scope* |

Two counts differ from the system count, and both are properties of the release:

- **`cath2`: 1,043 → 1,040.** Three system directories
  (`cath2_3bpqD00`, `cath2_3dh3A01`, `cath2_3bdlA01`) ship `dataset.json` and
  `topology.pdb` but no trajectories. They are recorded as `skipped`, not
  failed. 1,040 is also the number the manuscript reports (Table S1).
- **`megamerge`: 271 → 348.** 77 wild-types mix two force fields (see below).

## Scope

**`opep` (ONE-octapeptides) is not imported.** The decision is about what the
data is worth to us, not about any defect in it:

- The systems are 8-residue peptides — 113 to 146 atoms, solvent stripped.
- Of its 118,252 trajectories, **112,756 (95.4%) hold exactly two frames**,
  a single 10 ns gap each. They contribute 2.26 ms of the archive's 8.06 ms but
  carry almost no continuous dynamics.
- The continuous sampling is only the 5,496 `runNNN_protein.cmprsd.xtc` files,
  five per system at 1-2 us each.

Short peptides sampled in us-scale fragments are not what this deposition is
for, so the archive is left out rather than deposited for completeness.

**The tool still supports it in full** — `opep` remains in `DATASETS`, and
`run --datasets opep` imports it correctly, including the frame-stamp handling
described below. Nothing is disabled; the dataset is simply not part of the
planned import. If it is ever wanted, everything needed is already here and
tested.

Practical consequences:

- `init` still builds a work-list covering all 23,922 systems. The 1,100 `opep`
  rows stay `pending` for good and show in `status` totals. Pass
  `init --datasets megamerge,cath1,megamut,cath2` to leave them out entirely.
- A plain `run` with no `--datasets` **would** import `opep` first, since it is
  the smallest archive. Always name the four, or run `init` without `opep` so
  there is nothing pending for it.

## Per-system source layout

```
$ZIP_ROOT/$SYSTEM/
  dataset.json                 {force_field, save_traj_ns, system_id, temperature_K}
  topology.pdb                 protein-only, all-atom, solvent stripped
  reference.pdb                seed structure                     (MegaSim only)
  trajs/*.xtc                  coordinates, 10 ns/frame
  trajs/*.json                 per-trajectory force-field override (MegaSim only)
```

## How it works

`run` processes **one archive at a time, smallest first**, and deletes it once
its systems are drained — so peak disk is one zip (28 GB worst case, for
`cath2`) plus one IN_DIR per worker, never the whole 39.5 GB.

For each system a worker:

1. **claims** it with a non-blocking `flock` (per-system lockfile);
2. **splits** its trajectories into force-field groups (see below);
3. per group, **builds an IN_DIR** — `<system>.pdb`, a generated `<system>.psf`,
   the group's `.xtc` files, and `mdrepo-metadata.toml`;
4. **imports** it: `mdr-process validate` then `mdr-process process`;
5. **reclaims disk** — the IN_DIR goes right after its import, then the next
   group is built; the whole staging dir goes when the system completes.

State lives in a SQLite manifest (`<root>/manifest.sqlite`); the job is fully
**resumable** — kill and restart any time. Already-imported
`(system, group)` pairs are skipped on restart.

### Why a system can become more than one simulation

MDRepo records a single `forcefield` per simulation. In `MSR_megasim_merge`, a
trajectory may carry a sidecar json (`aggr_trj_run0_clone0_folded_resim.json`
next to `…_folded_resim.cmprsd.xtc`) that overrides `dataset.json`'s force
field. Per SI S.1.5.4, that is how the release marks the **77 of 271**
wild-types whose *folded* state was sampled with amber ff14sb because a99SB-disp
destabilised the native fold, while their *unfolded* state stayed on a99SB-disp.

Trajectories are therefore grouped by `(force_field, temperature_K)` and each
group becomes its own simulation, tagged `ff14sb-295k` / `ff99sb-disp-295k`.
Depositing them as one simulation would attach a force field to files it did not
produce; dropping the minority would discard released data.

### Coordinates are deposited byte-for-byte

Unlike the mdCATH and Dynamic PDB importers, **no coordinate block is ever
re-encoded**. In `cath1`, `cath2` and `megamerge` the published XTCs are
deposited byte-for-byte, carrying correct frame times — 10,000 ps between
consecutive frames, matching `save_traj_ns: 10.0`.

**⚠ `megamut` is the exception: its frame times are rewritten before deposit.**
Every trajectory in `MSR_megasim_mutants_disp_allatom` stamps its frame times in
*nanoseconds*, in a header field the XTC format defines as picoseconds, so the
file understates its own sampling interval by 1000x. The importer multiplies
every frame stamp by 1000 on the staged copy; see
[the mutant frame times](#megamut-frame-times-are-stamped-in-ns-not-ps) below.
Only the 4-byte time field of each frame header changes — no frame changes
length, and the coordinate blocks stay byte-identical to the release.

**⚠ Some released trajectories are aggregates whose clock restarts.** A file's
frame times are not necessarily monotonic: several cath1 trajectories are
concatenations of shorter segments, and the time field returns to ~0 at each
boundary. The spacing *within* a segment is still 10,000 ps, so the sampling
interval is what `dataset.json` claims — but any estimator of the form
`(last_t - first_t) / (n_frames - 1)` reads the restarts as a much shorter
interval (`cath1_1b43A02/run006` yields 60 ps that way). This importer therefore
takes the **modal positive gap between consecutive frames**, which ignores the
boundaries; see `modal_dt()`. **If `mdr-process` derives sampling the naive way,
it will mis-measure these files** — worth checking on the VM against a system
that logs no `frame-spacing-mismatch` note but is known to be an aggregate.

**Every octapeptide `*.filtered.cmprsd.xtc` carries a timestamp scaled by
1000.** This is handled, but `opep` is [out of scope](#scope), so it only
matters if that archive is ever imported. All 112,756 of them — 95% of the
`opep` trajectories, ~102 per system —
hold exactly two frames stamped `100000 ps` and `10100000 ps`, the *same* pair
in every file of every system. That gap is 1000x the 10 ns `dataset.json`
declares, which is what an ns→ps conversion using `1e6` instead of `1e3`
produces (0.1 ns and 10 ns becoming 100000 and 10100000 ps rather than 100 and
10000). A field identical across 112,756 files is not a per-trajectory clock, so
**`dataset.json`'s 10 ns is what this importer records.**

The manuscript decides it rather than leaving it to inference. The release holds
805,608 frames — 5,246 `runNNN` files of 101 frames, 250 of 201, and 112,756
filtered pairs — and at the declared 10 ns those sum to **8.06 ms, the
octapeptide total the paper reports**. Read literally, the filtered files alone
would be 11.3 seconds, ~1400x the published figure. That arithmetic also pins
the convention `sampled_ns` uses: frames × dt gives 8.06 ms and matches, while
elapsed spans give 6.87 ms and do not — so the per-simulation figures deposited
across the 1100 systems sum to the published total.

The 1000x case is recognised narrowly — exactly two frames, and a gap exactly
that factor off (`is_scaled_stamp()`) — and summarised in a single
`scaled-frame-stamp` note per system rather than one `frame-spacing-mismatch`
per file, which would otherwise write 112,756 rows saying the same thing. Any
*other* disagreement, including a 2-frame file off by some different factor, is
still reported per file.

#### `megamut` frame times are stamped in ns, not ps

This is the mirror image of the octapeptide case, and the one place this importer
alters a released file. Every trajectory in `MSR_megasim_mutants_disp_allatom`
stamps consecutive frames **10 ps** apart against the **10 ns** its
`dataset.json` declares — 1000x too *short*, where the octapeptide stamp is 1000x
too long. Three systems spanning the archive were range-read from Zenodo to
confirm it:

| System | index | frames | first *t* | last *t* | observed spacing | declared |
|---|---:|---:|---:|---:|---:|---:|
| `1A0N_L7S__A12D` | 0 | 89 | 119.0 | 999.0 | 10 ps | 10,000 ps |
| `2LYQ__L56C` | 10,729 | 97 | 30.0 | 990.0 | 10 ps | 10,000 ps |
| `v2_6IVS__Y46T` | 21,457 | 97 | 30.0 | 990.0 | 10 ps | 10,000 ps |

The ratio is exactly 0.001 in every case. What settles the direction is the
*start* and *end*: each file ends at t ≈ 1000 from a per-system start below it.
Read as ns those are the paper's 1 μs folded-state runs with the upstream burn-in
removed — which is exactly what the deposited description already says. Read as
ps they are 1 ns runs that all inexplicably stop at 1000 ps. The generating
script left corroborating fingerprints: `topology.pdb` carries `REMARK 1 CREATED
WITH MDTraj`, the `step` field counts frames rather than integration steps
(GROMACS would write 2,500,000 per 10 ns at 4 fs), and these are the only
trajectories in the release not named `*.cmprsd.xtc`.

**Why this one is corrected rather than noted.** MDRepo derives sampling from the
trajectory itself, so depositing these bytes unaltered would publish 21,458
mutants as ~1 ns of sampling and hand anyone who opens a file a time axis 1000x
short. The octapeptide stamp could be left alone because that archive is out of
scope; this one cannot. The correction multiplies the absolute stamp rather than
rebasing to zero — a nonzero start records how much burn-in was dropped for that
mutant, and 119.0 → 119,000 ps keeps it.

Only the staged copy under `staging/` is touched; the archive under `archives/`
is never modified, so a restage always re-derives from the untouched release and
the operation is idempotent. The trigger is `is_ns_stamped_as_ps()`, which fires
only when the observed spacing sits within 1% of exactly `save_traj_ns / 1000` —
it keys on the measured spacing, not the dataset name, so a partly-affected
archive is still handled and any other disagreement still falls through to a
per-file `frame-spacing-mismatch`. Each affected system records one
`rescaled-frame-times` note (`status --show-notes`), and the deposited
`description` states that the times were multiplied and that coordinates are
unchanged. `cath1`, `cath2` and `megamerge` were sampled the same way and stamp
correctly, so nothing in them is rewritten.

mdtraj's silent xdrfile-overflow failure mode (which `mdcath_import.py` has to
guard against, and re-audit with a `scan` subcommand) **cannot arise here**: no
`save_xtc` call exists to corrupt anything, and the rescale seeks to a fixed
offset in each frame header without touching a coordinate block. There is
consequently no coordinate-rewriting scan in this tool.

Each trajectory is still checked before deposit — its frame headers are walked
directly (no mdtraj needed) to confirm it parses, is not truncated, and holds
the same atom count as `topology.pdb`. A trajectory that fails is dropped from
the simulation, the drop is disclosed in the deposited `description`, and the
reason is recorded in the manifest (`status --show-notes`). If a group loses
*every* trajectory, that group is skipped rather than deposited empty.

### Topology is generated, not shipped

BioEmu ships only `topology.pdb`, and MDRepo does not accept a PDB as topology.
ParmEd reads the (hydrogen-containing) PDB and writes a matching `<system>.psf`
with the identical atom set and order. The deposited `description` says so and
points at this repo (`PROCESSING_REPO`) for the scripts that did it. See [bioemu_import.py](bioemu_import.py).

## Setup (on the VM, as exouser)

```bash
./setup_env.sh ~/bioemu-venv          # venv with parmed + numpy
source ~/bioemu-venv/bin/activate     # mdr-process must already be on PATH
```

## Usage

> `opep` is [out of scope](#scope). Every command below names the four datasets
> that are in scope; a bare `init`/`run` would pull it in.

```bash
# 0. the four in-scope datasets, in smallest-first order
IN_SCOPE=megamerge,cath1,megamut,cath2

# 1. one-time: build the 22,822-system work-list + the SIFTS UniProt map.
#    Only each zip's central directory is read (a few MB over HTTP range
#    requests) — no archive is downloaded yet.
python bioemu_import.py --root /opt/bioemu init --datasets $IN_SCOPE

# 2. dry-run a single system end-to-end first (no push to MDRepo)
python bioemu_import.py --root /opt/bioemu run \
    --datasets cath1 --systems cath1_1b43A02 -w 1 --dry-run --keep

# 3. real run against staging, smallest archive first
python bioemu_import.py --root /opt/bioemu run --datasets $IN_SCOPE -s staging -w 4

# just one dataset (archives are still processed smallest-first within the set)
python bioemu_import.py --root /opt/bioemu run --datasets megamerge,cath1 -s prod

# progress / failures / data notes
python bioemu_import.py --root /opt/bioemu status --show-failures 20
python bioemu_import.py --root /opt/bioemu status --show-notes 50

# requeue anything that failed, then re-run
python bioemu_import.py --root /opt/bioemu reset-failed

# stage ONE system for manual inspection (no import) — one dir per group
python bioemu_import.py --root /opt/bioemu extract megamerge HEEH_KT_rd6_0007 \
    --out-dir /tmp/insp

# download archives without importing (e.g. overnight, before a run)
python bioemu_import.py --root /opt/bioemu fetch --datasets cath2

# after an interrupted import, verify whether that simulation exists in MDRepo,
# then explicitly record the outcome (automatic resubmission is blocked)
python bioemu_import.py --root /opt/bioemu resolve-import cath1 cath1_1b43A02 \
    ff99sb-ildn-300k --imported
# or, only when it is absent from MDRepo:
python bioemu_import.py --root /opt/bioemu resolve-import cath1 cath1_1b43A02 \
    ff99sb-ildn-300k --retry

# redo a system's simulations from scratch
python bioemu_import.py --root /opt/bioemu reset-sim cath1 cath1_1b43A02
```

### Key `run` options

- `-w/--workers N` — worker processes (default 1). Workers share one archive.
- `-s/--server prod|staging` — MDRepo target (default `staging`).
- `-d/--dry-run` — pass `-d` to `mdr-process` (build import JSON, no push).
  Forced to one worker: a dry-run deliberately leaves the system `pending`, so
  independent workers cannot tell which they have already checked. The work
  list is snapshotted at start so each system is checked exactly once.
- `--datasets`, `--systems` — restrict the run (comma-separated, repeatable).
- `--keep` — don't delete the staged IN_DIR after success (debugging).
- `--keep-archives` — don't delete an archive once drained (re-runs need no
  re-download; costs 39.5 GB steady-state).
- `--no-verify` — skip the md5 check against Zenodo. Also makes an archive
  *already present* under `<root>/archives/` be used as-is without contacting
  Zenodo at all, which is the escape hatch for a VM that cannot reach it.
- `--num-threads`, `--work-dir`, `--out-dir` — forwarded to `mdr-process`
  (`-t/-w/-o`). Note `-w N` (workers) × `--num-threads T` is your total CPU load.
- `--blast-num-threads N` — forwarded to `mdr-process process`. Each worker runs
  its own `blastp`, so tune it *with* the worker count: `-w 8
  --blast-num-threads 2` and `-w 2 --blast-num-threads 8` both ask for 16 cores.
- `--min-free-gb` — pause instead of failing below this much free space
  (default 25). Must exceed the next archive's size for `run` to make progress.
- `--breaker-threshold` — halt all workers after N consecutive failures
  (default 10; `0` disables).
- `--skip-lock-check` — skip the `flock` exclusion self-test that guards
  multi-worker starts. Don't.

## Running multiple workers

Workers are **separate processes** (`mp.Process`), each with its own SQLite
connection and its own read-only handle on the archive — never threads sharing
one. That distinction is what makes the manifest safe:

- **The manifest uses a rollback journal (`journal_mode=DELETE`), not WAL.** WAL
  needs a shared-memory `-shm` mapping that network filesystems can't provide,
  so `DELETE` keeps the root relocatable. With `busy_timeout=60000` and
  `synchronous=FULL`, concurrent writers are correct.
- **`flock` is the real mutex, not the database.** A worker takes a non-blocking
  per-system lock, then *re-reads the system row under that lock* before doing
  anything. Two workers that pick the same candidate cannot both process it.
- **This requires a filesystem where `flock` actually excludes.** Some FUSE
  filesystems accept `flock()` and silently no-op; NFS mounted `nolock` degrades
  it. `run` performs a real two-process contention self-test before starting any
  multi-worker run and refuses to start if exclusion can't be demonstrated.

Choosing `N`:

- **Disk** is the binding constraint for the archive, not the workers: one zip
  (up to 28 GB) is shared, and each worker adds only one IN_DIR (a few hundred
  MB at most — the largest single system is well under 1 GB).
- **`megamut` is 21,458 tiny simulations**, so its cost is dominated by
  `mdr-process` round-trips rather than I/O. It benefits most from more workers.
- **Remote limits:** how many parallel pushes MDRepo tolerates may cap you below
  the CPU/disk limit.

Failure handling is built for unattended runs:

- **Transient failures don't consume the retry budget.** Full disk, connection
  errors, Zenodo throttling and remote 5xx requeue the system as `pending`
  without incrementing `attempts`, then back off 60 s.
- **A watchdog kills a wedged `mdr-process`.** If its whole process subtree makes
  no progress (zero I/O **and** zero CPU) for `--stall-minutes` (default 10), or
  it runs past `--mdr-max-hours` (default 4), the entire subtree is killed
  (`gocmd` grandchildren included) and it's treated as transient. Output streams
  to `logs/<dataset>__<system>__<group>.log` live, so `tail -f` works.
- **A circuit breaker stops the run** after `--breaker-threshold` consecutive
  failures across all workers.
- **A system with nothing importable is `skipped`, not `failed`** — it doesn't
  consume attempts and doesn't trip the breaker.

**Ambiguous imports.** `mdr-process process` can fail *part-way through* a push,
leaving the simulation partially landed under an allocated id. A non-zero exit
therefore does **not** prove nothing landed. Rather than risk a duplicate, the
`(system, group)` is left in state `importing` and `status` flags it:

```
⚠ 2 ambiguous import(s) awaiting resolve-import (verify in MDRepo, then --imported or --retry):
  cath1/cath1_1b43A02 [ff99sb-ildn-300k]  (since 2026-08-03T18:00:08Z)
```

Such a system is **never resubmitted automatically** — not even by
`reset-failed`. Check MDRepo, then record the outcome explicitly with
`resolve-import`. This is the one manual step, and it never blocks the rest of
the run. A `validate` failure creates no ambiguity (validation runs before the
import intent is recorded), so those are ordinary retryable failures.

## Provenance notes — verify before a prod run

MD settings are transcribed from the manuscript SI, not inferred from the files.
Items marked **⚠** should be confirmed with the authors.

**S.1.4, the standard protocol** (`cath1`, `cath2`, and `opep` if ever
imported): TIP3P water, cubic
box with 1 nm padding, 0.1 M NaCl. Equilibration 0.1 ns NVT + 0.9 ns NPT with
restraints on solute heavy atoms, released over 0.1 ns, at 2 fs. **Production
NPT at 300 K / 1 bar, hydrogen mass repartitioning (4 amu) with h-bond
constraints, 4 fs timestep.**

**S.1.5.4, MegaSim** (`megamerge`, `megamut`): rhombic dodecahedral box, 1.5 nm
padding. Equilibration 0.2 ns NVT + 0.6 ns NPT at 295 K / 1 bar, Langevin, 4 fs.
**Production NVT at 295 K, 4 fs**, same constraints and HMR as S.1.4.

- **`software_name = "CUSTOM"`, `software_version = "NA"`** — S.1.4 says "We
  internally developed code … based on OpenMM as its compute engine". MDRepo's
  `software_name` vocabulary is closed (ACEMD/AMBER/CHARMM/CUSTOM/GROMACS/
  NAMD/SPONGE) and holds neither OpenMM nor "in-house harness", so `CUSTOM` is
  the honest entry and the spec then requires `software_version = "NA"`. OpenMM
  is preserved in the `description`. (Do **not** use AMBER/GROMACS — those name
  programs that did not run these simulations; the force field is Amber, the
  engine is not.)
- **`[water]` is recorded only for the standard-protocol datasets.** S.1.4 states
  TIP3P; S.1.5.4 re-specifies MegaSim's box and equilibration but restates
  neither the water model nor the salt, so MegaSim depositions carry no `[water]`
  and no `[[solutes]]` rather than a guessed value. a99SB-disp is supplied with
  its own four-point water model, which is noted in `forcefield_comments`.
  **⚠ If the authors confirm the MegaSim water model and ion concentration, add
  them to `FORCEFIELDS` / `STANDARD_SOLUTES`.**
- **⚠ `[water].density_kg_m3 = 1000.0`** is nominal — the release publishes no
  per-system solvent density, and the field is required once `[water]` is
  present.
- **⚠ `lead_contributor_orcid = "0000-0000-0000-0000"`** — the placeholder used
  by the sibling importers, meaning "submitted locally as administrator". It is
  deliberately not one of the author ORCIDs: the lead contributor is whoever
  submits the deposition, not an author of the paper.
- **`[[contributors]]` is the full 28-author list, in manuscript order**, so
  credit in the deposition matches credit in the paper. `orcid` is optional and
  is recorded for the 8 authors whose ORCID has been checked against their
  public profile (Lewis, Hempel, Jimenez-Luna, Gastegger, Xie, Abdin, Clementi,
  Noe); the release itself publishes no ORCIDs, so the other 20 carry a name
  alone rather than a guessed id. `institution` is likewise kept only for the
  four original dataset contacts — the manuscript affiliations were not
  transcribed, and an unverified affiliation is worse than none.
- **Names carry their diacritics** (`José Jiménez-Luna`, `Victor García
  Satorras`, `Frank Noé`, `Freie Universität`), in `[[contributors]]` and in the
  citation string alike; a test pins the two to each other. This means the
  emitted TOML is non-ASCII, so `mdrepo-metadata.toml` is written with an
  explicit `encoding="utf-8"` — inheriting the locale's encoding would raise
  `UnicodeEncodeError` on a VM running under `LANG=C`. TOML is UTF-8 by spec, so
  `mdr-process` should be reading it that way regardless. **⚠ If `mdr-process
  validate` turns out to choke on non-ASCII contributor names, transliterate
  `CONTRIBUTORS` and `PAPER_BIOEMU["authors"]` together** — they are asserted
  equal, so changing one alone fails the suite.
- **IDs.** `cath1`/`cath2` system names embed a CATH domain
  (`cath1_1b43A02` → PDB `1b43`, chain `A`), so UniProt comes from SIFTS for
  that exact chain. (`opep` peptides are synthetic and get neither field —
  `mdr-process` is invoked with `--no-id` for those — but that path is unused
  while the dataset is [out of scope](#scope).) MegaSim names are MEGAscale
  entry names: `pdb_id` is taken from the 4-character prefix wherever one exists
  (`1AOY`, `1A0N_L7S`, `2HBB_pross6`, and mutants like `1A0N_L7S__A12D` → `1a0n`),
  with the variant or point mutation stated in the `description`; UniProt is the
  union across the entry's chains, since a MEGAscale entry names a domain rather
  than a chain. De novo designs (`EEHEE_rd3_0019`, `HEEH_KT_rd6_0007`,
  `EA_run2_*`, `r10_572_TrROS_Hall`) get neither. **Note this deliberately
  attaches a parent `pdb_id` to ~21.5k mutants whose sequence differs from the
  deposited entry** — the fold is the parent's, the sequence is not.

  Eight entries invert that order: a version-tagged redesign of a natural
  structure leads with its tag and puts the code **last** — `v2_6IVS`,
  `v2K43S_2KVV`, and `v2R14S_R16S_2L3X`, where the code is the third token.
  Parsed as a prefix these look like designs, which cost all 455 systems built
  on them (447 mutants plus the 8 wild-types in `megamerge`) both their `pdb_id`
  and, through it, their UniProt accessions. They resolve to `2hdz`, `2kvv`,
  `2l3x`, `2lc2`, `2ldm`, `2lxe`, `4uzx` and `6ivs`. Recognition requires *both*
  a leading `v<digit>` tag and a trailing PDB code, which keeps it off the
  `r6`/`r7`/`r10`/`r11` TrRosetta designs — they lead with a similar token but
  end in `_Hall`.

## ⚠ Two things that can only be verified on the VM

1. **PSF acceptance.** A PSF built by ParmEd from a bare PDB carries
   connectivity and atom names but not force-field charges/types. If
   `mdr-process validate` rejects it, switch `write_psf_from_pdb()` to an
   OpenMM-typed topology
   (`openmm.app.ForceField('amber99sbildn.xml').createSystem(...)` →
   `parmed.openmm.load_topology(...).save(psf)`). Test with `extract` plus a
   manual `mdr-process validate` on one system before the full run.
2. **Large `trajectory_file_names` lists.** With `opep` out of scope this is
   much reduced — the ~107-file octapeptide simulations were the stress case.
   `cath2` still averages 39 files per simulation and tops out at 53
   (`cath2_2vhvA02`), so confirm `mdr-process` handles that on one such system
   before running the dataset — it has to stat and read every file to sum the
   duration.

Two smaller ones worth a glance on the first dry-run:

3. **`[[additional_files]]`** carries MegaSim's `reference.pdb`. If `validate`
   objects, set `INCLUDE_REFERENCE_PDB = False` and it is described in prose
   instead.
4. **`.cmprsd.xtc` double extensions** (`run001_protein.cmprsd.xtc`) — the name
   ends in `.xtc`, so it should satisfy the spec, but it is unusual.
5. **Non-ASCII `μ`.** Durations are written `1 μs`, not `1 us`, in both
   `description` and `short_description`. `SHORT_DESCRIPTION_MAX` counts
   characters, so the 300-char limit is unaffected; confirm `mdr-process`
   round-trips the UTF-8 on the first dry-run.

## Layout under `--root`

```
manifest.sqlite               SQLite state (systems, sims, uniprot, archives, notes)
archives/<name>.zip           one at a time; deleted once drained
archives/<name>.zip.part      resumable partial download
staging/<dataset>/<system>/<group>/    transient IN_DIR (deleted on success)
locks/<dataset>__<system>.lock         per-system flock files
locks/archive__<dataset>.lock          serialises concurrent archive fetches
logs/<dataset>__<system>__<group>.log  mdr-process output per simulation
```

## Tests

```bash
python make_fixtures.py                      # ~3 MB of real members from Zenodo
python -m pytest test_bioemu_import.py -q     # 125 tests, ~1 min
```

The fixtures are small zips assembled from **real** archive members — a real
GROMACS-written `topology.pdb` and real `.cmprsd.xtc` files — so header walking,
atom-count checks and force-field splitting run against the actual byte layouts.
Each dataset's members come from *its own* archive, which matters more than it
looks, and has now cost twice. The octapeptide fixture originally stood in cath1
content, and that alone kept the 1000x-scaled two-frame stamp — 95% of that
archive — out of the suite until the first real run hit it. The mutant fixture
then stood in `megamerge` content, which stamps its frame times *correctly*, so
the suite actively asserted the opposite of what `megamut` does and hid the
ns-in-a-ps-field bug the same way. Both now use members from their own archive.
The end-to-end tests drive the real worker loop, manifest and locking against a
stub `mdr-process` that asserts every file the TOML names is present.

## Citation

Lewis, S. *et al.* Scalable emulation of protein equilibrium ensembles with
generative deep learning. *Science* **389**, eadv9817 (2025).
[doi:10.1126/science.adv9817](https://doi.org/10.1126/science.adv9817)

MegaSim seed structures and stability measurements: Tsuboyama, K. *et al.*
Mega-scale experimental analysis of protein folding stability in biology and
design. *Nature* **620**, 434–444 (2023).

Data: [10.5281/zenodo.15629740](https://doi.org/10.5281/zenodo.15629740) (CATH),
[10.5281/zenodo.15641184](https://doi.org/10.5281/zenodo.15641184) (MegaSim),
[10.5281/zenodo.15641199](https://doi.org/10.5281/zenodo.15641199) (Octapeptides).
Released under CDLA-Permissive-2.0.
