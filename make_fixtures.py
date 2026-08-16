#!/usr/bin/env python3
"""Build the small archive fixtures the test suite runs against.

The fixtures mirror the real BioEmu archive layouts, but are assembled from a
handful of *real* members pulled out of the published Zenodo zips — a real
GROMACS-written topology.pdb and real .cmprsd.xtc files — so the tests exercise
the actual byte layouts rather than synthetic stand-ins.

Only the few members needed are fetched, using HTTP range requests against
Zenodo (the full archives are 41 GB), so this costs a few MB and a minute.

    python make_fixtures.py [--out fixtures]

Members are cached under <out>/.members so re-running is free.
"""
from __future__ import annotations

import argparse
import io
import json
import time
import urllib.request
import zipfile
from pathlib import Path

ZIP_URL = "https://zenodo.org/records/{record}/files/{name}?download=1"

# (cache name, record, archive, member) — the smallest set that covers every
# shape the importer has to handle.
MEMBERS = [
    ("cath1_topology.pdb", "15629740", "ONE_cath1.zip",
     "ONE_cath1/cath1_1b43A02/topology.pdb"),
    ("cath1_run000.xtc", "15629740", "ONE_cath1.zip",
     "ONE_cath1/cath1_1b43A02/trajs/run000_protein.cmprsd.xtc"),
    ("cath1_rep0.xtc", "15629740", "ONE_cath1.zip",
     "ONE_cath1/cath1_1b43A02/trajs/1b43A02_0.cmprsd.xtc"),
    ("megasim_topology.pdb", "15641184", "MSR_megasim_merge.zip",
     "MSR_megasim_merge/HEEH_KT_rd6_0007/topology.pdb"),
    ("megasim_reference.pdb", "15641184", "MSR_megasim_merge.zip",
     "MSR_megasim_merge/HEEH_KT_rd6_0007/reference.pdb"),
    ("megasim_disp.xtc", "15641184", "MSR_megasim_merge.zip",
     "MSR_megasim_merge/HEEH_KT_rd6_0007/trajs/aggr_trj_reseed_run5.cmprsd.xtc"),
    ("megasim_ff14sb.xtc", "15641184", "MSR_megasim_merge.zip",
     "MSR_megasim_merge/HEEH_KT_rd6_0007/trajs/"
     "aggr_trj_run0_clone0_folded_resim.cmprsd.xtc"),
    # Real mutant members. These matter for the same reason the octapeptide ones
    # below do: every trajectory in this archive stamps its frame times in ns
    # inside the ps field, and standing in megamerge content here hid that from
    # the suite entirely — megamerge stamps correctly, so the substitution
    # silently asserted the opposite of what the archive does.
    ("megamut_topology.pdb", "15641184", "MSR_megasim_mutants_disp_allatom.zip",
     "MSR_megasim_mutants_disp_allatom/1A0N_L7S__A12D/topology.pdb"),
    ("megamut_mutant.xtc", "15641184", "MSR_megasim_mutants_disp_allatom.zip",
     "MSR_megasim_mutants_disp_allatom/1A0N_L7S__A12D/trajs/trj_mutant_folded.xtc"),
    # Real octapeptide members. The filtered one matters: 95% of that archive is
    # two-frame files carrying a timestamp 1000x dataset.json (see the README),
    # and standing in cath1 content here would hide that shape from the suite
    # entirely — which is exactly what it did until the first real opep run.
    ("opep_topology.pdb", "15641199", "ONE_octapeptides.zip",
     "ONE_octapeptides/opep_0000/topology.pdb"),
    ("opep_filtered.xtc", "15641199", "ONE_octapeptides.zip",
     "ONE_octapeptides/opep_0000/trajs/e10s1_e8s2p0f150-ADRIA_LARGEPEP_"
     "opep_0000-0-1-RND0375_9.filtered.cmprsd.xtc"),
    ("opep_run001.xtc", "15641199", "ONE_octapeptides.zip",
     "ONE_octapeptides/opep_0000/trajs/run001_protein.cmprsd.xtc"),
]


class RemoteFile(io.RawIOBase):
    """Range-request-backed file object, so zipfile can read a remote archive."""

    def __init__(self, url: str):
        self.url = url
        self.pos = 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=120) as r:
            self.size = int(r.headers["Content-Length"])

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
        for attempt in range(5):          # Zenodo throttles ranged reads
            try:
                req = urllib.request.Request(
                    self.url, headers={"Range": f"bytes={self.pos}-{end}"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 * (attempt + 1))
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def fetch_members(cache: Path) -> dict[str, bytes]:
    cache.mkdir(parents=True, exist_ok=True)
    out: dict[str, bytes] = {}
    by_archive: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for name, record, archive, member in MEMBERS:
        if (cache / name).exists():
            out[name] = (cache / name).read_bytes()
        else:
            by_archive.setdefault((record, archive), []).append((name, member))
    for (record, archive), wanted in by_archive.items():
        url = ZIP_URL.format(record=record, name=archive)
        print(f"reading {archive} over HTTP for {len(wanted)} member(s)…")
        zf = zipfile.ZipFile(io.BufferedReader(RemoteFile(url), buffer_size=1 << 16))
        try:
            for name, member in wanted:
                data = zf.read(member)
                (cache / name).write_bytes(data)
                out[name] = data
                print(f"  {name:26s} {len(data):>9,} bytes")
        finally:
            zf.close()
    return out


def dataset_json(system_id: str, ff: str, temp: float) -> str:
    return json.dumps({"force_field": ff, "save_traj_ns": 10.0,
                       "system_id": system_id, "temperature_K": temp}, indent=4)


def build(out: Path, m: dict[str, bytes]) -> None:
    out.mkdir(parents=True, exist_ok=True)

    # ONE_cath1: one ordinary system, plus one that ships no trajectories (the
    # real MSR_cath2 has 3 such directories).
    with zipfile.ZipFile(out / "ONE_cath1.zip", "w", zipfile.ZIP_DEFLATED) as z:
        base = "ONE_cath1/cath1_1b43A02"
        z.writestr(f"{base}/dataset.json",
                   dataset_json(base, "amber ff99sb-ildn", 300.0))
        z.writestr(f"{base}/topology.pdb", m["cath1_topology.pdb"])
        z.writestr(f"{base}/trajs/run000_protein.cmprsd.xtc", m["cath1_run000.xtc"])
        z.writestr(f"{base}/trajs/1b43A02_0.cmprsd.xtc", m["cath1_rep0.xtc"])
        base = "ONE_cath1/cath1_9zzzA00"
        z.writestr(f"{base}/dataset.json",
                   dataset_json(base, "amber ff99sb-ildn", 300.0))
        z.writestr(f"{base}/topology.pdb", m["cath1_topology.pdb"])

    # MSR_megasim_merge: one system with a per-trajectory force-field override
    # (the 77 real ones), two without.
    with zipfile.ZipFile(out / "MSR_megasim_merge.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for system in ("HEEH_KT_rd6_0007", "1A0N_L7S", "1AOY"):
            base = f"MSR_megasim_merge/{system}"
            z.writestr(f"{base}/dataset.json",
                       dataset_json(base, "amber ff99sb-disp", 295.0))
            z.writestr(f"{base}/topology.pdb", m["megasim_topology.pdb"])
            z.writestr(f"{base}/reference.pdb", m["megasim_reference.pdb"])
            z.writestr(f"{base}/trajs/aggr_trj_reseed_run5.cmprsd.xtc",
                       m["megasim_disp.xtc"])
            if system == "HEEH_KT_rd6_0007":
                stem = "aggr_trj_run0_clone0_folded_resim"
                z.writestr(f"{base}/trajs/{stem}.cmprsd.xtc", m["megasim_ff14sb.xtc"])
                z.writestr(f"{base}/trajs/{stem}.json",
                           dataset_json(base, "amber ff14sb", 295.0))

    # MSR_megasim_mutants: exactly one trajectory per system, 21458 of them, and
    # every one of them stamped in ns rather than ps — real bytes, so the suite
    # sees the 10-ps-against-a-declared-10-ns spacing the importer has to correct.
    with zipfile.ZipFile(out / "MSR_megasim_mutants_disp_allatom.zip", "w",
                         zipfile.ZIP_DEFLATED) as z:
        base = "MSR_megasim_mutants_disp_allatom/1A0N_L7S__A12D"
        z.writestr(f"{base}/dataset.json",
                   dataset_json(base, "amber ff99sb-disp", 295))
        z.writestr(f"{base}/topology.pdb", m["megamut_topology.pdb"])
        z.writestr(f"{base}/trajs/trj_mutant_folded.xtc", m["megamut_mutant.xtc"])

    # ONE_octapeptides: no PDB or UniProt identifier anywhere, and both shapes
    # the archive actually holds — one `runNNN` file whose spacing matches
    # dataset.json, and one `.filtered.` file with the 1000x-scaled two-frame
    # stamp that 112,756 of its trajectories carry.
    with zipfile.ZipFile(out / "ONE_octapeptides.zip", "w", zipfile.ZIP_DEFLATED) as z:
        base = "ONE_octapeptides/opep_0000"
        z.writestr(f"{base}/dataset.json",
                   dataset_json(base, "amber ff99sb-ildn", 300.0))
        z.writestr(f"{base}/topology.pdb", m["opep_topology.pdb"])
        z.writestr(f"{base}/trajs/e10s1_e8s2p0f150-ADRIA_LARGEPEP_opep_0000-0-1-"
                   f"RND0375_9.filtered.cmprsd.xtc", m["opep_filtered.xtc"])
        z.writestr(f"{base}/trajs/run001_protein.cmprsd.xtc", m["opep_run001.xtc"])

    # Trajectories that must be dropped rather than deposited: not an xtc at
    # all, truncated mid-frame, and one whose atom count does not match the
    # topology (728 vs 1002).
    with zipfile.ZipFile(out / "broken.zip", "w", zipfile.ZIP_DEFLATED) as z:
        base = "ONE_cath1/cath1_1b43A02"
        z.writestr(f"{base}/dataset.json",
                   dataset_json(base, "amber ff99sb-ildn", 300.0))
        z.writestr(f"{base}/topology.pdb", m["cath1_topology.pdb"])
        z.writestr(f"{base}/trajs/good.xtc", m["cath1_run000.xtc"])
        z.writestr(f"{base}/trajs/garbage.xtc", b"not an xtc at all, no magic here....")
        z.writestr(f"{base}/trajs/truncated.xtc", m["cath1_run000.xtc"][:5000])
        z.writestr(f"{base}/trajs/wrong_atoms.xtc", m["megasim_disp.xtc"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(Path(__file__).with_name("fixtures")))
    args = ap.parse_args()
    out = Path(args.out)
    build(out, fetch_members(out / ".members"))
    print()
    for p in sorted(out.glob("*.zip")):
        print(f"{p.name:45s} {p.stat().st_size:>9,} bytes")
    print(f"\nfixtures ready in {out}\nrun:  python -m pytest test_bioemu_import.py -q")


if __name__ == "__main__":
    main()
