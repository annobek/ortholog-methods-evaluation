"""
plot_cactus_score_contradictory_fn.py

Goal
-----
Create a scatter plot for contradictory false negatives (gene-level plotting unit):

- X-axis: CactusSupportScore for the *C-pair*  (Focus_Gene -> Cactus_Ortholog_Gene)
- Y-axis: CactusSupportScore for the *N-pair*  (Focus_Gene -> ProtSyn_Ortholog_Gene)

To do this, the script merges:
1) focus_candidate_map.tsv produced by make_fn_pair_tables.py
2) get_cactus_score.sh output for c-pairs
3) get_cactus_score.sh output for n-pairs

Input files
-----------
--focus_map  : TSV with one row per contradictory focus gene (one plotted point).
               Required columns:
                 Row_ID, FN_ID, Focus_Gene, Focus_Species,
                 Cactus_Ortholog_Gene, Cactus_Ortholog_Species,
                 ProtSyn_Ortholog_Gene, ProtSyn_Ortholog_Species,
                 FN_Class

--c_scores   : TSV from get_cactus_score.sh run on *_c_pairs_with_coords.tsv
--n_scores   : TSV from get_cactus_score.sh run on *_n_pairs_with_coords.tsv
               Required columns in both score TSVs:
                 gene1, species1, gene2, species2, Score

Outputs
-------
--out_tsv : merged TSV with two new columns:
              CactusSupportScore_cactus
              CactusSupportScore_protsyn
--out_png : scatter plot PNG

Data types / behavior notes
---------------------------
- All table values are read as strings and stripped.
- Score values are parsed as float; missing/invalid scores become empty in out_tsv.
- Points with missing X or Y are not plotted by default; use --plot_missing_as_zero
  if you prefer plotting missing as 0.0 (not recommended for contradictory-only).
- Negative scores are allowed and will be plotted (linear axes).
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")  # required for Slurm/headless nodes
import matplotlib.pyplot as plt


def die(msg: str, code: int = 2) -> None:
    """
    Print an error message to stderr and exit the program with a non-zero code.

    Parameters
    ----------
    msg : str
        Human-readable error message.
    code : int
        Exit code (non-zero means failure).

    Returns
    -------
    None (the program exits).
    """
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def norm(s: str) -> str:
    """
    Normalize a value from TSV:
    - Convert None to empty string
    - Strip surrounding whitespace

    Parameters
    ----------
    s : str

    Returns
    -------
    str
    """
    return (s or "").strip()


def is_missing(v: str) -> bool:
    """
    Decide whether a TSV cell should be treated as missing.

    Treats any of these (case-insensitive) as missing:
      - empty string
      - na / nan
      - none / null
      - -
      - no_result

    Parameters
    ----------
    v : str

    Returns
    -------
    bool
    """
    v = norm(v)
    return v == "" or v.lower() in {"na", "nan", "none", "null", "-", "no_result"}


def parse_float(v: str):
    """
    Parse a float from a TSV cell.

    Parameters
    ----------
    v : str

    Returns
    -------
    float or None
        None if missing or not parseable.
    """
    v = norm(v)
    if is_missing(v):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def ensure_outdir(path: str) -> None:
    """
    Ensure the parent directory of `path` exists.

    Parameters
    ----------
    path : str

    Returns
    -------
    None
    """
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def read_scores(path: str):
    """
    Read get_cactus_score.sh output (TSV) and build a lookup:
      (gene1, species1, gene2, species2) -> Score(float)

    Notes
    -----
    - If the same key appears multiple times and has different scores, we keep the first
      and print a warning.

    Parameters
    ----------
    path : str
        Path to score TSV.

    Returns
    -------
    dict[tuple[str,str,str,str], float]
    """
    scores = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            die(f"Score file has no header: {path}")

        needed = {"gene1", "gene2", "species1", "species2", "Score"}
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            die(
                f"Score file missing columns {missing}: {path}\n"
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            key = (
                norm(row.get("gene1", "")),
                norm(row.get("species1", "")),
                norm(row.get("gene2", "")),
                norm(row.get("species2", "")),
            )
            sc = parse_float(row.get("Score", ""))
            if sc is None:
                continue

            if key in scores and abs(scores[key] - sc) > 1e-12:
                print(
                    f"WARNING: Duplicate key with different Score in {path}: {key} "
                    f"{scores[key]} vs {sc} (keeping first)",
                    file=sys.stderr,
                )
                continue

            scores.setdefault(key, sc)

    return scores


def read_focus_candidate_map(path: str):
    """
    Read the focus-candidate mapping table produced by make_fn_pair_tables.py.

    Parameters
    ----------
    path : str
        Path to *_focus_candidate_map.tsv

    Returns
    -------
    list[dict[str,str]]
        Each dict contains the required columns as strings.
    """
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            die(f"Focus-candidate map has no header: {path}")

        needed = {
            "Row_ID",
            "FN_ID",
            "Focus_Gene",
            "Focus_Species",
            "Cactus_Ortholog_Gene",
            "Cactus_Ortholog_Species",
            "ProtSyn_Ortholog_Gene",
            "ProtSyn_Ortholog_Species",
            "FN_Class",
        }
        missing = [c for c in needed if c not in (reader.fieldnames or [])]
        if missing:
            die(
                f"Focus-candidate map missing columns {missing}: {path}\n"
                f"Found: {reader.fieldnames}"
            )

        for r in reader:
            rows.append({k: norm(r.get(k, "")) for k in needed})

    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Merge CactusSupportScore for c-pairs and n-pairs and plot scatter (X=cactus, Y=prot-syn)."
    )
    ap.add_argument(
        "--focus_map",
        required=True,
        help="*_focus_candidate_map.tsv from make_fn_pair_tables.py",
    )
    ap.add_argument(
        "--c_scores",
        required=True,
        help="get_cactus_score.sh output for c-pairs (TSV).",
    )
    ap.add_argument(
        "--n_scores",
        required=True,
        help="get_cactus_score.sh output for n-pairs (TSV).",
    )
    ap.add_argument("--out_tsv", required=True, help="Output merged TSV path.")
    ap.add_argument("--out_png", required=True, help="Output plot PNG path.")
    ap.add_argument("--title", default="Cactus support: C-pair vs N-pair", help="Plot title.")
    ap.add_argument("--point_size", type=float, default=12.0, help="Scatter marker size.")
    ap.add_argument("--alpha", type=float, default=0.6, help="Point alpha (0..1).")
    ap.add_argument(
        "--plot_missing_as_zero",
        action="store_true",
        help="If set, missing scores are plotted as 0.0 instead of being skipped.",
    )
    args = ap.parse_args()

    focus_rows = read_focus_candidate_map(args.focus_map)
    c_map = read_scores(args.c_scores)
    n_map = read_scores(args.n_scores)

    merged_rows = []
    xs = []
    ys = []
    missing_x = 0
    missing_y = 0

    for r in focus_rows:
        # We expect all rows to be contradictory, but keep this guard.
        if r.get("FN_Class") != "contradictory":
            continue

        fg = r["Focus_Gene"]
        fs = r["Focus_Species"]

        cg = r["Cactus_Ortholog_Gene"]
        cs = r["Cactus_Ortholog_Species"]

        pg = r["ProtSyn_Ortholog_Gene"]
        ps = r["ProtSyn_Ortholog_Species"]

        key_c = (fg, fs, cg, cs)
        key_n = (fg, fs, pg, ps)

        x = c_map.get(key_c)
        y = n_map.get(key_n)

        if x is None:
            missing_x += 1
        if y is None:
            missing_y += 1

        out = dict(r)
        out["CactusSupportScore_cactus"] = "" if x is None else str(x)
        out["CactusSupportScore_protsyn"] = "" if y is None else str(y)
        merged_rows.append(out)

        # Decide what gets plotted
        if args.plot_missing_as_zero:
            xs.append(0.0 if x is None else x)
            ys.append(0.0 if y is None else y)
        else:
            if x is None or y is None:
                continue
            xs.append(x)
            ys.append(y)

    if len(merged_rows) == 0:
        die("No rows found in focus_map (or none were contradictory).", code=3)

    # Write merged TSV (force Unix newlines to avoid ^M problems)
    ensure_outdir(args.out_tsv)
    with open(args.out_tsv, "w", newline="") as f:
        header = [
            "Row_ID",
            "FN_ID",
            "Focus_Gene",
            "Focus_Species",
            "Cactus_Ortholog_Gene",
            "Cactus_Ortholog_Species",
            "ProtSyn_Ortholog_Gene",
            "ProtSyn_Ortholog_Species",
            "FN_Class",
            "CactusSupportScore_cactus",
            "CactusSupportScore_protsyn",
        ]
        w = csv.DictWriter(
            f,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",  # avoid CRLF (^M)
        )
        w.writeheader()
        for row in merged_rows:
            w.writerow({h: row.get(h, "") for h in header})

    # Plot
    ensure_outdir(args.out_png)
    plt.figure()
    if xs:
        plt.scatter(xs, ys, s=args.point_size, alpha=args.alpha)

    # Reference lines
    plt.axhline(0.0)
    plt.axvline(0.0)

    # y=x line for reference (use symmetric range if possible)
    if xs:
        minv = min(min(xs), min(ys))
        maxv = max(max(xs), max(ys))
        pad = (maxv - minv) * 0.05 if maxv > minv else 0.1
        lo = minv - pad
        hi = maxv + pad
        plt.plot([lo, hi], [lo, hi])
        plt.xlim(lo, hi)
        plt.ylim(lo, hi)

    plt.xlabel("CactusSupportScore (C-pair: Focus → Cactus ortholog)")
    plt.ylabel("CactusSupportScore (N-pair: Focus → ProtSyn ortholog)")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=300)
    plt.close()

    print("Done.")
    print(f"Merged rows written: {args.out_tsv}")
    print(f"Plot written: {args.out_png}")
    print(f"Points plotted: {len(xs)}")
    print(f"Missing X (CactusSupportScore_cactus): {missing_x}")
    print(f"Missing Y (CactusSupportScore_protsyn): {missing_y}")


if __name__ == "__main__":
    main()
