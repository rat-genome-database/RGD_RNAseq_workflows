#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path
from datetime import datetime

"""
Build a JBrowse2 session JSON for one BIOProject.

Usage:
  python make_jbrowse_session_for_bioproject.py BIOProjectID PRJdir baseDir [passAccList]

What it does:
  - Finds per-sample RNAseq track JSON files under PRJdir (RNAseq_*.json),
    excluding *geneTPMbed.json and *TXTPMbed.json
  - Includes ONLY samples that passed STARQC. The set of PASS samples is read
    from the PASS accession list that the orchestrator builds and drives RSEM /
    ComputeSex / BWjson from (single source of truth):
        <baseDir>/log_files/STARQC/<BIOProjectID>_Unique_AccList_PASS.txt
    (override by passing the PASS AccList path as the optional 4th argument).
    geo_accession is read from column 2 (header row skipped). A track is kept
    only if its sample (metadata["Sample Accession ID"], or the geo_accession
    parsed from the trackId) is in that list.
    If the PASS AccList cannot be found the script aborts rather than silently
    including failed samples; set ALLOW_MISSING_QC=1 to fall back to no filtering.
  - Forces BigWig URLs to the public RGD location:
#     6 March 2026 change location to https://download.rgd.mcw.edu/expression/<BIOProjectID>/ExpressionProfileBigWig/<trackId>.bigwig
      https://download.rgd.mcw.edu/expression/<BIOProjectID>/Genome-wide_read_coverage_BigWig_files/<trackId>.bigwig

  - Assigns colors deterministically based on composite grouping:
      (Tissue, Strain, Sex, Sample Characteristic)
  - Writes a combined session JSON to:
      <baseDir>/<BIOProjectID>_jbrowse_session_GRCr8.json

Important fix:
  - JBrowse often uses the display config in views[0].tracks to render initially-open tracks.
    So we also set the view track display color (displays[].color),
    in addition to renderer.color1 in sessionTracks.
5 Feb 2026 need to update so that this works in GRCr8 make sure no mRatBN7.2 in this file.
"""


def pick_color_for_key(key, palette, assigned):
    """Deterministic color assignment by composite key."""
    if key in assigned:
        return assigned[key]
    color = palette[len(assigned) % len(palette)]
    assigned[key] = color
    return color


def normalize_str(x) -> str:
    """Normalize metadata strings (strip, collapse None)."""
    if x is None:
        return ""
    return str(x).strip()


def load_pass_samples(pass_acc_list: Path):
    """Return the set of PASS sample IDs from the PASS accession list.

    This is the same PASS-only, deduplicated list the orchestrator builds from
    the STARQC report and uses to drive RSEM / ComputeSex / BWjson:
        <baseDir>/log_files/STARQC/<BIOProjectID>_Unique_AccList_PASS.txt
    It is tab-separated with a header row; geo_accession is in column 2:
        Run  geo_accession  Tissue  Strain  Sex  PMID  GEOpath  Title  ...
    Returns None if the list cannot be read.
    """
    if not pass_acc_list.is_file():
        return None

    pass_ids = set()
    with pass_acc_list.open() as fh:
        header = fh.readline()  # skip header row
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            geo_accession = cols[1].strip()
            if geo_accession:
                pass_ids.add(geo_accession)
    return pass_ids


def sample_id_for_track(track, track_id) -> str:
    """Best-effort sample (geo_accession) for a track.

    Prefer metadata["Sample Accession ID"] (set by BWjson_v7.sh); fall back to
    the last underscore-delimited token of the trackId, which is built as
    RNAseq_<Tissue>_<Strain>_<Sex>_<geo_accession>.
    """
    metadata = track.get("metadata", {}) or {}
    sid = normalize_str(metadata.get("Sample Accession ID", ""))
    if sid:
        return sid
    return track_id.split("_")[-1] if track_id else ""


def main():
    if len(sys.argv) not in (4, 5):
        print("Usage: make_jbrowse_session_for_bioproject.py BIOProjectID PRJdir baseDir [passAccList]", file=sys.stderr)
        sys.exit(1)

    BIOProjectID = sys.argv[1]
    PRJdir = Path(sys.argv[2])
    baseDir = Path(sys.argv[3])

    if not PRJdir.is_dir():
        print(f"Error: PRJdir does not exist or is not a directory: {PRJdir}", file=sys.stderr)
        sys.exit(1)

    # ---- Resolve PASS accession list and the set of PASS samples ----
    # Default location matches the orchestrator: Logdir=<baseDir>/log_files
    if len(sys.argv) == 5:
        pass_acc_list = Path(sys.argv[4])
    else:
        pass_acc_list = baseDir / "log_files" / "STARQC" / f"{BIOProjectID}_Unique_AccList_PASS.txt"

    pass_samples = load_pass_samples(pass_acc_list)

    allow_missing_qc = os.environ.get("ALLOW_MISSING_QC", "") == "1"
    if pass_samples is None:
        if allow_missing_qc:
            print(
                f"WARNING: PASS AccList not found ({pass_acc_list}); "
                "ALLOW_MISSING_QC=1 set, proceeding WITHOUT QC filtering.",
                file=sys.stderr,
            )
        else:
            print(
                f"Error: PASS AccList not found: {pass_acc_list}\n"
                "Refusing to build a session that may include QC-failed samples.\n"
                "Pass the PASS AccList path as the 4th argument, or set "
                "ALLOW_MISSING_QC=1 to bypass QC filtering.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        print(f"Loaded {len(pass_samples)} PASS sample(s) from {pass_acc_list}", file=sys.stderr)

    # Find only RNAseq_*.json, excluding helper JSONs
    track_files = []
    for p in PRJdir.rglob("RNAseq_*.json"):
        name = p.name
        if "geneTPMbed" in name or "TXTPMbed" in name:
            continue
        track_files.append(p)

    track_files = sorted(track_files)

    if not track_files:
        print(f"No RNAseq_*.json track files found under {PRJdir}", file=sys.stderr)
        sys.exit(1)

    print("Found track JSON files:", file=sys.stderr)
    for p in track_files:
        print(f"  {p}", file=sys.stderr)

    # Palette (cycles if more groups than colors)
    palette = [
        "#1f77b4",  # blue
        "#ff7f0e",  # orange
        "#2ca02c",  # green
        "#d62728",  # red
        "#9467bd",  # purple
        "#8c564b",  # brown
        "#e377c2",  # pink
        "#7f7f7f",  # gray
        "#bcbd22",  # olive
        "#17becf",  # teal
    ]

    combo_to_color = {}   # composite key -> color
    trackid_to_color = {} # trackId -> color (used for view config)

    session_tracks = []

    # Debug: show which group got what color
    print("\n--- Color assignment debug (trackId -> composite key -> color) ---", file=sys.stderr)

    for p in track_files:
        with p.open() as fh:
            track = json.load(fh)

        track_id = track.get("trackId")
        if not track_id:
            print(f"Warning: track JSON {p} has no 'trackId'; skipping", file=sys.stderr)
            continue

        # ---- QC gate: skip any sample that did not PASS STARQC ----
        if pass_samples is not None:
            sid = sample_id_for_track(track, track_id)
            if sid not in pass_samples:
                print(
                    f"Skipping QC-failed/unknown sample: trackId={track_id} "
                    f"sample={sid or '?'} ({p})",
                    file=sys.stderr,
                )
                continue

        metadata = track.get("metadata", {})

        # Force correct track type
        track["type"] = "QuantitativeTrack"

        # Force correct public BigWig URL
        new_uri = (
            f"https://download.rgd.mcw.edu/expression/"
            f"{BIOProjectID}/Genome-wide_read_coverage_BigWig_files/"
            f"{track_id}.bigwig"
        )
        adapter = track.setdefault("adapter", {})
        adapter["type"] = "BigWigAdapter"
        bw = adapter.setdefault("bigWigLocation", {})
        bw["locationType"] = "UriLocation"
        bw["uri"] = new_uri

        # Composite grouping key: Tissue + Strain + Sex + Sample Characteristic
        tissue = normalize_str(metadata.get("Tissue", ""))

        # Use only strain name for grouping (strip URL/extra after comma)
        strain_raw = normalize_str(metadata.get("Strain", ""))
        strain_name = strain_raw.split(",")[0].strip() if strain_raw else ""

        sex = normalize_str(metadata.get("Sex", ""))

        sample_char = normalize_str(metadata.get("Sample Characteristic", ""))
        sample_char = " ".join(sample_char.split())  # collapse whitespace

        if any([tissue, strain_name, sex, sample_char]):
            combo_key = f"Tissue={tissue};Strain={strain_name};Sex={sex};SampleChar={sample_char}"
        else:
            combo_key = "Unknown"

        color = pick_color_for_key(combo_key, palette, combo_to_color)
        trackid_to_color[track_id] = color

        # Debug line
        print(f"{track_id}\t{combo_key}\t{color}", file=sys.stderr)

        # Ensure displays exist and apply color
        displays = track.get("displays")
        if not isinstance(displays, list) or not displays:
            displays = [{
                "type": "LinearWiggleDisplay",
                "displayId": f"{track_id}-LinearWiggleDisplay",
            }]
            track["displays"] = displays

        for disp in displays:
            disp["type"] = "LinearWiggleDisplay"
            disp.setdefault("displayId", f"{track_id}-LinearWiggleDisplay")

            # SessionTracks: keep renderer.color1 (fine)
            renderer = disp.get("renderer", {})
            renderer["type"] = "XYPlotRenderer"
            renderer["color1"] = color
            disp["renderer"] = renderer

            renderers = disp.get("renderers", {})
            renderers["XYPlotRenderer"] = {"type": "XYPlotRenderer", "color1": color}
            disp["renderers"] = renderers

            disp["defaultRendering"] = "xyplot"

        session_tracks.append(track)

    if not session_tracks:
        print("No usable RNAseq tracks found after filtering; exiting.", file=sys.stderr)
        sys.exit(1)

    # Build a LinearGenomeView with tracks opened by default
    # CRITICAL FOR YOUR COLOR ISSUE:
    # Put displays[].color directly in views[0].tracks for each track.
    view_tracks = []
    for t in session_tracks:
        tid = t.get("trackId")
        if not tid:
            continue

        color = trackid_to_color.get(tid, "#1f77b4")

        view_tracks.append({
            "type": "QuantitativeTrack",
            "configuration": tid,
            "displays": [{
                "type": "LinearWiggleDisplay",
                "displayId": f"{tid}-LinearWiggleDisplay",
                "color": color,                 # THIS is what exported config uses
                "defaultRendering": "xyplot"
            }]
        })

    # ---- Gene track first (top) ----
    GENE_TRACK_BLOCK = {
        "id": "F-8qwRhumS",
        "type": "FeatureTrack",
        "configuration": "Rat GRCr8 (rn8) Genes and Transcripts-GRCr8",
        "minimized": False,
        "displays": [
            {
                "id": "uZq89S4_XC",
                "type": "LinearBasicDisplay",
                "heightPreConfig": 152,
                "configuration": "Rat GRCr8 (rn8) Genes and Transcripts-GRCr8-LinearBasicDisplay",
            }
        ],
    }
    view_tracks = [GENE_TRACK_BLOCK] + view_tracks
    # -------------------------------

    # Default view: Chr4:12,315,273..14,497,135 (1-based)
    TARGET_REFNAME = "Chr4"
    TARGET_START_1 = 12_315_273
    TARGET_END_1 = 14_497_135

    # Broad span so the session is not limited
    WHOLE_CHR4_END = 1_000_000_000

    assumed_viewport_px = 2000
    window_bp = max(1, TARGET_END_1 - TARGET_START_1 + 1)
    bp_per_px = max(1.0, window_bp / float(assumed_viewport_px))

    target_start_0 = TARGET_START_1 - 1
    offset_px = float(target_start_0) / float(bp_per_px)

    linear_view = {
        "id": "lgv1",
        "type": "LinearGenomeView",
        "tracks": view_tracks,
        "displayedRegions": [{
            "assemblyName": "GRCr8",
            "refName": TARGET_REFNAME,
            "start": 0,
            "end": WHOLE_CHR4_END
        }],
        "bpPerPx": bp_per_px,
        "offsetPx": offset_px
    }

    session_name = f"{BIOProjectID}_RNAseq_expression"
    timestamp = datetime.now().isoformat(timespec="seconds")

    session_root = {
        "session": {
            "name": session_name,
            "description": f"Auto-generated session for {BIOProjectID} on {timestamp}",
            "views": [linear_view],
            "sessionTracks": session_tracks
        }
    }

    baseDir.mkdir(parents=True, exist_ok=True)
    out_path = baseDir / f"{BIOProjectID}_jbrowse_session_GRCr8.json"

    with out_path.open("w") as out_fh:
        json.dump(session_root, out_fh, indent=2)

    print(f"\nSession JSON written to: {out_path}", file=sys.stderr)
    print(f"Default view should open near: {TARGET_REFNAME}:{TARGET_START_1:,}..{TARGET_END_1:,}", file=sys.stderr)
    print("Genes/Transcripts track inserted as first (top) view track.", file=sys.stderr)
    print("View track colors are set via displays[].color (matches exported example).", file=sys.stderr)


if __name__ == "__main__":
    main()
