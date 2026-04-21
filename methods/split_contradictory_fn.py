#!/usr/bin/env python3
"""
make_fn_pair_tables.py

Extract contradictory-only rows from a gene-level FN TSV and write:
  1) *_c_pairs.tsv  : Focus_Gene -> Cactus_Ortholog_Gene
  2) *_n_pairs.tsv  : Focus_Gene -> ProtSyn_Ortholog_Gene
  3) *_focus_candidate_map.tsv : mapping used later to merge X/Y scores for plotting

IMPORTANT:
- Output files are written with Unix newlines to avoid CRLF (^M) issues with awk.
"""

import argparse
import csv
import os
import sys


REQUIRED_COLS = [
    "FN_ID",
    "Focus_Gene",
    "Focus_Species",
    "Cactus_Ortholog_Gene",
    "Cactus_Ortholog_Species",
    "ProtSyn_Ortholog_Gene",
    "FN_Class",
]

OUT_HEADERS_FOR_BUILD_PAIR = ["Gene1", "Species1", "Gene2_prot_syn", "Species2_prot_syn"]


def die(msg: str, code: int = 2) -> None:
    """Print an error message to stderr and exit with a non-zero code."""
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def norm(s: str) -> str:
    """Normalize a TSV cell value: None -> '', and strip whitespace."""
    return (s or "").strip()


def is_missing(v: str) -> bool:
    """Return True if a TSV cell should be treated as missing."""
    v = norm(v)
    return v == "" or v.lower() in {"na", "nan", "none", "null", "-", "no_result"}


def open_tsv(path: str):
    """Open a TSV file and return (file_handle, csv.DictReader)."""
    f = open(path, "r", newline="")
    reader = csv.DictReader(f, delimiter="\t")
    return f, reader


def ensure_cols(reader: csv.DictReader, required) -> None:
    """Exit with error if the TSV is missing any required columns."""
    hdr = reader.fieldnames or []
    missing = [c for c in required if c not in hdr]
    if missing:
        die(
            "Missing required columns in input FN table:\n"
            + "\n".join(f"  - {m}" for m in missing)
            + "\n\nFound columns:\n  "
            + ", ".join(hdr)
        )


def write_tsv(path: str, header, rows_iter) -> None:
    """
    Write rows (dicts) to TSV at `path` with a fixed `header` order.

    IMPORTANT:
    - csv.DictWriter defaults to CRLF (\r\n). We force Unix newlines via
      lineterminator="\n" to avoid ^M issues in awk.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    # Use newline="" as recommended by Python csv docs, then control line endings
    # via lineterminator in the writer.
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",   # <-- THIS is the real fix
        )
        w.writeheader()
        for r in rows_iter:
            w.writerow(r)



def main() -> None:
    ap = argparse.ArgumentParser(
        description="Extract contradictory-only c-pairs and n-pairs from gene-level FN table."
    )
    ap.add_argument("--fn", required=True, help="Path to gene-level FN table (TSV).")
    ap.add_argument(
        "--out_prefix",
        required=True,
        help="Output prefix, e.g. results/all_species_contradictory (no extension).",
    )
    ap.add_argument(
        "--allow_missing_protsyn",
        action="store_true",
        help="Keep contradictory rows even if ProtSyn ortholog gene is missing (not recommended).",
    )
    args = ap.parse_args()

    out_c_pairs = args.out_prefix + "_c_pairs.tsv"
    out_n_pairs = args.out_prefix + "_n_pairs.tsv"
    out_focus_map = args.out_prefix + "_focus_candidate_map.tsv"

    fin, reader = open_tsv(args.fn)
    try:
        ensure_cols(reader, REQUIRED_COLS)

        has_ps_species_col = "ProtSyn_Ortholog_Species" in (reader.fieldnames or [])

        c_rows = []
        n_rows = []
        focus_map_rows = []

        kept = 0
        skipped_non_contra = 0
        skipped_missing = 0
        row_id = 0

        for row in reader:
            fn_class = norm(row.get("FN_Class", ""))
            if fn_class != "contradictory":
                skipped_non_contra += 1
                continue

            fn_id = norm(row.get("FN_ID", ""))
            fg = norm(row.get("Focus_Gene", ""))
            fs = norm(row.get("Focus_Species", ""))

            cg = norm(row.get("Cactus_Ortholog_Gene", ""))
            cs = norm(row.get("Cactus_Ortholog_Species", ""))

            pg = norm(row.get("ProtSyn_Ortholog_Gene", ""))

            ps = norm(row.get("ProtSyn_Ortholog_Species", "")) if has_ps_species_col else ""
            if is_missing(ps):
                ps = cs  # infer ProtSyn species from cactus ortholog species

            # Require cactus side
            if is_missing(fn_id) or is_missing(fg) or is_missing(fs) or is_missing(cg) or is_missing(cs):
                skipped_missing += 1
                continue

            # Require ProtSyn gene by default
            if not args.allow_missing_protsyn and is_missing(pg):
                skipped_missing += 1
                continue

            row_id += 1

            # C-pairs (X)
            c_rows.append(
                {"Gene1": fg, "Species1": fs, "Gene2_prot_syn": cg, "Species2_prot_syn": cs}
            )

            # N-pairs (Y)
            if not is_missing(pg):
                n_rows.append(
                    {"Gene1": fg, "Species1": fs, "Gene2_prot_syn": pg, "Species2_prot_syn": ps}
                )

            # Focus-candidate mapping
            focus_map_rows.append(
                {
                    "Row_ID": str(row_id),
                    "FN_ID": fn_id,
                    "Focus_Gene": fg,
                    "Focus_Species": fs,
                    "Cactus_Ortholog_Gene": cg,
                    "Cactus_Ortholog_Species": cs,
                    "ProtSyn_Ortholog_Gene": pg,
                    "ProtSyn_Ortholog_Species": ps,
                    "FN_Class": fn_class,
                }
            )

            kept += 1

        write_tsv(out_c_pairs, OUT_HEADERS_FOR_BUILD_PAIR, c_rows)
        write_tsv(out_n_pairs, OUT_HEADERS_FOR_BUILD_PAIR, n_rows)
        write_tsv(
            out_focus_map,
            [
                "Row_ID",
                "FN_ID",
                "Focus_Gene",
                "Focus_Species",
                "Cactus_Ortholog_Gene",
                "Cactus_Ortholog_Species",
                "ProtSyn_Ortholog_Gene",
                "ProtSyn_Ortholog_Species",
                "FN_Class",
            ],
            focus_map_rows,
        )

        print("Done.")
        print(f"Input: {args.fn}")
        print(f"Kept contradictory rows (points): {kept}")
        print(f"Skipped non-contradictory rows: {skipped_non_contra}")
        print(f"Skipped due to missing required fields: {skipped_missing}")
        print(f"Wrote: {out_c_pairs}")
        print(f"Wrote: {out_n_pairs}")
        print(f"Wrote: {out_focus_map}")

        if kept == 0:
            die("No contradictory rows kept. Check FN_Class values and column names.", code=3)

    finally:
        fin.close()


if __name__ == "__main__":
    main()
