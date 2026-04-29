#!/usr/bin/env python3
"""
make_jbrowse_session_for_combined_bioproject_v2.py

Build a JBrowse2 session JSON for a combined BioProject directory.

Usage:
  python make_jbrowse_session_for_combined_bioproject_v2.py COMBINED_ID PRJdir baseDir

Arguments:
  COMBINED_ID   The merged BioProject/GEO accession ID (e.g. GSE70012).
                Used to update public-facing metadata links and to name output files.
  PRJdir        Path to the combined project directory that contains per-sample
                reads_fastq/ subdirectories, each of which holds RNAseq_*.json
                track definition files.
  baseDir       Directory where the output session JSON will be written.
                Typically the same as PRJdir.

Example:
  python make_jbrowse_session_for_combined_bioproject_v2.py GSE70012 \\
    /path/to/data/expression/GEO/GSE70012 \\
    /path/to/data/expression/GEO/GSE70012

Behavior:
  - Finds RNAseq_*.json files recursively under PRJdir.
  - Excludes *geneTPMbed.json and *TXTPMbed.json.
  - Updates only public-facing metadata links:
      metadata["RGD Metadata Report"]      -> geoAcc=COMBINED_ID
      metadata["Project Repository Link"]  -> acc=COMBINED_ID
  - Leaves metadata["Project Accession ID"] unchanged for traceability.
  - Forces BigWig URIs to point to the combined project download directory.
  - Assigns colors by biological group (Tissue / Strain / Sex / Sample Characteristic).
  - Writes:
      <baseDir>/<COMBINED_ID>_jbrowse_session_GRCr8.json

Configuration:
  Two values near the top of main() may need to be adjusted for your deployment:
    BIGWIG_BASE_URL  — base URL from which BigWig files are publicly served.
    GENE_TRACK_CFG   — JBrowse2 configuration key for the reference gene/transcript track.
"""

import json
import re
import sys
from pathlib import Path
from datetime import datetime


def normalize_str(x) -> str:
    if x is None:
        return ""
    return str(x).strip()


def update_project_links(track, combined_id):
    """
    Update ONLY:
      - metadata['RGD Metadata Report']
      - metadata['Project Repository Link']

    Do NOT modify:
      - metadata['Project Accession ID']
    """
    metadata = track.get("metadata", {})

    rgd_url = metadata.get("RGD Metadata Report")
    if isinstance(rgd_url, str):
        metadata["RGD Metadata Report"] = re.sub(
            r"geoAcc=[^&]+",
            f"geoAcc={combined_id}",
            rgd_url,
        )

    repo_url = metadata.get("Project Repository Link")
    if isinstance(repo_url, str):
        metadata["Project Repository Link"] = re.sub(
            r"acc=[^&]+",
            f"acc={combined_id}",
            repo_url,
        )

    track["metadata"] = metadata
    return track


def pick_color_for_key(key, palette, assigned):
    if key in assigned:
        return assigned[key]
    color = palette[len(assigned) % len(palette)]
    assigned[key] = color
    return color


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: make_jbrowse_session_for_combined_bioproject_v2.py "
            "COMBINED_ID PRJdir baseDir",
            file=sys.stderr,
        )
        sys.exit(1)

    combined_id = sys.argv[1]
    PRJdir = Path(sys.argv[2])
    baseDir = Path(sys.argv[3])

    # ------------------------------------------------------------------
    # CONFIGURATION — adjust these for your deployment environment
    # ------------------------------------------------------------------
    # Base URL from which BigWig files are publicly served.
    # The full URI per track will be:
    #   {BIGWIG_BASE_URL}/{combined_id}/Genome-wide_read_coverage_BigWig_files/{track_id}.bigwig
    BIGWIG_BASE_URL = "https://your.download.server/expression"   # <-- replace

    # JBrowse2 configuration key for the reference gene/transcript annotation track.
    # This is the value of the "configuration" field in the gene FeatureTrack block.
    GENE_TRACK_CFG = "Your_Species_Assembly_Genes_and_Transcripts-AssemblyName"  # <-- replace
    GENE_TRACK_DISPLAY_CFG = f"{GENE_TRACK_CFG}-LinearBasicDisplay"              # auto-derived
    # ------------------------------------------------------------------

    if not PRJdir.is_dir():
        print(f"Error: PRJdir does not exist or is not a directory: {PRJdir}", file=sys.stderr)
        sys.exit(1)

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

    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    combo_to_color = {}
    trackid_to_color = {}
    session_tracks = []

    print("\n--- Color assignment debug ---", file=sys.stderr)

    for p in track_files:
        with p.open() as fh:
            track = json.load(fh)

        track = update_project_links(track, combined_id)

        track_id = track.get("trackId")
        if not track_id:
            print(f"Warning: track JSON {p} has no trackId; skipping", file=sys.stderr)
            continue

        metadata = track.get("metadata", {})

        track["type"] = "QuantitativeTrack"

        new_uri = (
            f"{BIGWIG_BASE_URL}/"
            f"{combined_id}/Genome-wide_read_coverage_BigWig_files/"
            f"{track_id}.bigwig"
        )

        adapter = track.setdefault("adapter", {})
        adapter["type"] = "BigWigAdapter"

        bw = adapter.setdefault("bigWigLocation", {})
        bw["locationType"] = "UriLocation"
        bw["uri"] = new_uri

        tissue = normalize_str(metadata.get("Tissue", ""))

        strain_raw = normalize_str(metadata.get("Strain", ""))
        strain_name = strain_raw.split(",")[0].strip() if strain_raw else ""

        sex = normalize_str(metadata.get("Sex", ""))

        sample_char = normalize_str(metadata.get("Sample Characteristic", ""))
        sample_char = " ".join(sample_char.split())

        if any([tissue, strain_name, sex, sample_char]):
            combo_key = (
                f"Tissue={tissue};"
                f"Strain={strain_name};"
                f"Sex={sex};"
                f"SampleChar={sample_char}"
            )
        else:
            combo_key = "Unknown"

        color = pick_color_for_key(combo_key, palette, combo_to_color)
        trackid_to_color[track_id] = color

        print(f"{track_id}\t{combo_key}\t{color}", file=sys.stderr)

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

            renderer = disp.get("renderer", {})
            renderer["type"] = "XYPlotRenderer"
            renderer["color1"] = color
            disp["renderer"] = renderer

            renderers = disp.get("renderers", {})
            renderers["XYPlotRenderer"] = {
                "type": "XYPlotRenderer",
                "color1": color,
            }
            disp["renderers"] = renderers

            disp["defaultRendering"] = "xyplot"

        session_tracks.append(track)

    if not session_tracks:
        print("No usable RNAseq tracks found after filtering; exiting.", file=sys.stderr)
        sys.exit(1)

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
                "color": color,
                "defaultRendering": "xyplot",
            }],
        })

    gene_track_block = {
        "id": "gene-track-1",
        "type": "FeatureTrack",
        "configuration": GENE_TRACK_CFG,
        "minimized": False,
        "displays": [
            {
                "id": "gene-track-display-1",
                "type": "LinearBasicDisplay",
                "heightPreConfig": 152,
                "configuration": GENE_TRACK_DISPLAY_CFG,
            }
        ],
    }

    view_tracks = [gene_track_block] + view_tracks

    # ------------------------------------------------------------------
    # Default viewport — update refName, coordinates, and assembly name
    # to match your reference genome.
    # ------------------------------------------------------------------
    target_refname = "Chr1"          # <-- replace with a representative chromosome
    target_start_1 = 1               # 1-based start of the default view window
    target_end_1   = 10_000_000      # 1-based end of the default view window
    whole_chr_end  = 1_000_000_000   # large sentinel; JBrowse clips to actual chrom length
    assembly_name  = "YourAssembly"  # <-- replace (e.g. "GRCr8", "GRCh38")
    # ------------------------------------------------------------------

    assumed_viewport_px = 2000
    window_bp = max(1, target_end_1 - target_start_1 + 1)
    bp_per_px = max(1.0, window_bp / float(assumed_viewport_px))

    target_start_0 = target_start_1 - 1
    offset_px = float(target_start_0) / float(bp_per_px)

    linear_view = {
        "id": "lgv1",
        "type": "LinearGenomeView",
        "tracks": view_tracks,
        "displayedRegions": [{
            "assemblyName": assembly_name,
            "refName": target_refname,
            "start": 0,
            "end": whole_chr_end,
        }],
        "bpPerPx": bp_per_px,
        "offsetPx": offset_px,
    }

    timestamp = datetime.now().isoformat(timespec="seconds")

    session_root = {
        "session": {
            "name": f"{combined_id}_RNAseq_expression",
            "description": f"Auto-generated combined session for {combined_id} on {timestamp}",
            "views": [linear_view],
            "sessionTracks": session_tracks,
        }
    }

    baseDir.mkdir(parents=True, exist_ok=True)
    out_path = baseDir / f"{combined_id}_jbrowse_session_{assembly_name}.json"

    with out_path.open("w") as out_fh:
        json.dump(session_root, out_fh, indent=2)

    print(f"\nSession JSON written to: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
