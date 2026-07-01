import pandas as pd
import os

"""
Automate tracing of Prot-Syn anomaly cases (ProtSyn candidate has Nscore_ProtSyn = 0)
through Prot-Syn intermediate tables, using the key:

    (focus_gene, target_species)

Interpretation:
  - focus_gene is the gene from the anomaly row (Focus_Gene)
  - target_species is the "other species in the comparison" (e.g., rat in a dog-rat pair)

Because gene IDs are globally unique across species, this key avoids false matches across
different comparisons while remaining simple and intuitive.

Orientation handling:
  In Prot-Syn files, a focus_gene can appear in either Gene1 or Gene2 column.
  This script always treats the row as an edge and records the partner gene on the other side,
  but only if the partner species matches target_species.

Canonicalization:
  Species names vary ("Canis lupus familiaris" vs "Canis_Lupus" vs "Mus").
  By default we use genus only as species_key:
    "Canis lupus familiaris" -> "canis"
    "Canis_Lupus"            -> "canis"
    "Mus musculus"           -> "mus"
    "Mus"                    -> "mus"

Inputs (TSV)
------------
Required:
- anomaly.tsv
    Must include: Focus_Gene, Focus_Species, Cactus_Ortholog_Gene, Cactus_Ortholog_Species,
                  ProtSyn_Ortholog_Gene
- pairs_original.tsv
    Columns: Gene1 Species1 Protein1 Gene2 Species2 Protein2 SharedNeighbors GeometricMean
- pairs_transitive_local.tsv
    Same structure as pairs_original
- global_pairs.tsv
    Columns: Gene1 Protein1 Species1 Gene2 Protein2 Species2
- one_to_one.tsv
    Same structure as global_pairs

Optional:
- false_negatives.tsv
    If given and includes FN_ID, used for extra context merge (not required for tracing).

Output (TSV)
------------
anomaly_trace_report.tsv:
  Original anomaly columns plus:
    - TargetSpeciesKey
    - for Focus->ProtSyn and Focus->Cactus:
        * present_in_pairs_original, nscore_pairs_original
        * present_in_pairs_transitive_local, nscore_pairs_transitive_local
        * present_in_global_pairs
        * present_in_one_to_one
    - top partners for focus / protsyn genes within the target species context
    - FLAG_final_pair_unscored:
        True if Focus->ProtSyn is in one_to_one but absent from both scored tables
"""



def read_tsv(path: str) -> pd.DataFrame:
    """
    Read TSV into a DataFrame; force all columns to string.
    """
    return pd.read_csv(path, sep="\t", dtype=str)


def normalize_gene(g):
    """
    Normalize gene identifier to stripped string, or None.
    """
    if pd.isna(g):
        return None
    g = str(g).strip()
    return g if g else None


def species_key(s):
    """
    Canonicalize species label to a stable key (default: genus-only, lowercased).
    """
    if pd.isna(s):
        return None
    s = str(s).strip()
    if not s:
        return None
    genus = s.replace("_", " ").split()[0]
    return genus.lower()


def index_scored_by_target(df: pd.DataFrame,
                           gene1_col="Gene1", sp1_col="Species1",
                           gene2_col="Gene2", sp2_col="Species2",
                           nscore_col="SharedNeighbors"):
    """
    Build an index for scored tables (pairs_original / pairs_transitive_local):

        key = (focus_gene, target_species_key)
        value = dict {partner_gene: max_nscore}

    Also builds an index for retrieving top partners easily:

        partners_index[(focus_gene, target_species_key)] = list of (partner_gene, nscore) sorted desc.

    Orientation handling:
      If focus_gene is Gene1, partner is Gene2 and partner species must match target.
      If focus_gene is Gene2, partner is Gene1 and partner species must match target.

    Returns
    -------
    pair_nscore : dict
        (focus_gene, target_species_key, partner_gene) -> max_nscore
        (this is convenient for direct existence checks)
    partners_index : dict
        (focus_gene, target_species_key) -> [(partner_gene, nscore), ...] sorted desc
    """
    pair_nscore = {}
    partners_tmp = {}

    for _, r in df.iterrows():
        g1 = normalize_gene(r[gene1_col]); s1 = species_key(r[sp1_col])
        g2 = normalize_gene(r[gene2_col]); s2 = species_key(r[sp2_col])
        if g1 is None or g2 is None or s1 is None or s2 is None:
            continue

        try:
            ns = float(r[nscore_col])
        except Exception:
            continue

        # record edge both ways: (g1 -> species of g2) partner g2, and (g2 -> species of g1) partner g1
        k12 = (g1, s2, g2)
        k21 = (g2, s1, g1)

        if k12 not in pair_nscore or ns > pair_nscore[k12]:
            pair_nscore[k12] = ns
        if k21 not in pair_nscore or ns > pair_nscore[k21]:
            pair_nscore[k21] = ns

        partners_tmp.setdefault((g1, s2), []).append((g2, ns))
        partners_tmp.setdefault((g2, s1), []).append((g1, ns))

    # sort and deduplicate partners by max score
    partners_index = {}
    for k, lst in partners_tmp.items():
        best = {}
        for p, ns in lst:
            if p not in best or ns > best[p]:
                best[p] = ns
        partners_index[k] = sorted(best.items(), key=lambda x: x[1], reverse=True)

    return pair_nscore, partners_index


def index_unscored_by_target(df: pd.DataFrame,
                             gene1_col="Gene1", sp1_col="Species1",
                             gene2_col="Gene2", sp2_col="Species2"):
    """
    Build an index for unscored tables (global_pairs / one_to_one):

        key = (focus_gene, target_species_key)
        value = set(partner_genes)

    Returns
    -------
    partners_set : dict
        (focus_gene, target_species_key) -> set(partner_gene)
    """
    partners_set = {}

    for _, r in df.iterrows():
        g1 = normalize_gene(r[gene1_col]); s1 = species_key(r[sp1_col])
        g2 = normalize_gene(r[gene2_col]); s2 = species_key(r[sp2_col])
        if g1 is None or g2 is None or s1 is None or s2 is None:
            continue

        partners_set.setdefault((g1, s2), set()).add(g2)
        partners_set.setdefault((g2, s1), set()).add(g1)

    return partners_set


def top_partners(partners_index, focus_gene, target_species_key, k=5):
    """
    Format top-k partners for (focus_gene, target_species_key).
    """
    focus_gene = normalize_gene(focus_gene)
    if focus_gene is None or target_species_key is None:
        return ""
    lst = partners_index.get((focus_gene, target_species_key), [])[:k]
    return ";".join(f"{p}({ns:g})" for p, ns in lst)


def main(anomaly_tsv, fn_tsv,
         pairs_original_tsv, pairs_transitive_local_tsv,
         global_pairs_tsv, one_to_one_tsv,
         outdir,
         out_name="anomaly_trace_report.tsv"):


    anomaly = read_tsv(anomaly_tsv)

    # Optional FN merge for context only
    if fn_tsv:
        fn = read_tsv(fn_tsv)
        if "FN_ID" in anomaly.columns and "FN_ID" in fn.columns:
            anomaly = anomaly.merge(
                fn[["FN_ID", "Gene1", "Gene2", "Pair", "Swapped"]].drop_duplicates(),
                on="FN_ID",
                how="left",
                suffixes=("", "_FN")
            )

    po = read_tsv(pairs_original_tsv)
    ptl = read_tsv(pairs_transitive_local_tsv)
    gp = read_tsv(global_pairs_tsv)
    oto = read_tsv(one_to_one_tsv)

    po_pair_nscore, po_partners = index_scored_by_target(po)
    ptl_pair_nscore, ptl_partners = index_scored_by_target(ptl)
    gp_partners_set = index_unscored_by_target(gp)
    oto_partners_set = index_unscored_by_target(oto)

    out_rows = []

    for _, r in anomaly.iterrows():
        focus_gene = normalize_gene(r.get("Focus_Gene"))
        target_species = species_key(r.get("Cactus_Ortholog_Species"))  # "other species" in your comparison
        ps_gene = normalize_gene(r.get("ProtSyn_Ortholog_Gene"))
        cactus_gene = normalize_gene(r.get("Cactus_Ortholog_Gene"))

        row = dict(r)
        row["TargetSpeciesKey"] = "" if target_species is None else target_species

        # --- Focus -> ProtSyn candidate checks ---
        po_ns_ps = po_pair_nscore.get((focus_gene, target_species, ps_gene))
        ptl_ns_ps = ptl_pair_nscore.get((focus_gene, target_species, ps_gene))

        row["in_pairs_original_focus_ps"] = po_ns_ps is not None
        row["nscore_pairs_original_focus_ps"] = 0.0 if po_ns_ps is None else po_ns_ps

        row["in_pairs_transitive_local_focus_ps"] = ptl_ns_ps is not None
        row["nscore_pairs_transitive_local_focus_ps"] = 0.0 if ptl_ns_ps is None else ptl_ns_ps

        row["in_global_pairs_focus_ps"] = ps_gene in gp_partners_set.get((focus_gene, target_species), set())
        row["in_one_to_one_focus_ps"] = ps_gene in oto_partners_set.get((focus_gene, target_species), set())

        # --- Focus -> Cactus expected ortholog checks ---
        po_ns_c = po_pair_nscore.get((focus_gene, target_species, cactus_gene))
        ptl_ns_c = ptl_pair_nscore.get((focus_gene, target_species, cactus_gene))

        row["in_pairs_original_focus_cactus"] = po_ns_c is not None
        row["nscore_pairs_original_focus_cactus"] = 0.0 if po_ns_c is None else po_ns_c

        row["in_pairs_transitive_local_focus_cactus"] = ptl_ns_c is not None
        row["nscore_pairs_transitive_local_focus_cactus"] = 0.0 if ptl_ns_c is None else ptl_ns_c

        row["in_global_pairs_focus_cactus"] = cactus_gene in gp_partners_set.get((focus_gene, target_species), set())
        row["in_one_to_one_focus_cactus"] = cactus_gene in oto_partners_set.get((focus_gene, target_species), set())

        # --- Context: top partners within target species ---
        row["top_partners_pairs_original_focus"] = top_partners(po_partners, focus_gene, target_species, k=5)
        row["top_partners_pairs_original_protsyn_gene"] = top_partners(po_partners, ps_gene, species_key(r.get("Focus_Species")), k=5)
        row["top_partners_pairs_transitive_local_focus"] = top_partners(ptl_partners, focus_gene, target_species, k=5)
        row["top_partners_pairs_transitive_local_protsyn_gene"] = top_partners(ptl_partners, ps_gene, species_key(r.get("Focus_Species")), k=5)

        # Flag: final pair exists, but never appeared in scored tables (=> explains nscore=0)
        row["FLAG_final_pair_unscored"] = bool(
            row["in_one_to_one_focus_ps"]
            and row["nscore_pairs_original_focus_ps"] == 0.0
            and row["nscore_pairs_transitive_local_focus_ps"] == 0.0
        )

    out_rows.append(row)

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, out_name)

    pd.DataFrame(out_rows).to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path} ({len(out_rows)} rows)")



if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--anomaly", required=True)
    ap.add_argument("--fn", default=None)
    ap.add_argument("--pairs_original", required=True)
    ap.add_argument("--pairs_transitive_local", required=True)
    ap.add_argument("--global_pairs", required=True)
    ap.add_argument("--one_to_one", required=True)
    ap.add_argument(
        "--outdir",
        required=True,
        help="Output directory for anomaly trace report"
    )
    ap.add_argument(
        "--out-name",
        default="anomaly_trace_report.tsv",
        help="Output filename (default: anomaly_trace_report.tsv)"
    )

    args = ap.parse_args()

    main(
      anomaly_tsv=args.anomaly,
      fn_tsv=args.fn,
      pairs_original_tsv=args.pairs_original,
      pairs_transitive_local_tsv=args.pairs_transitive_local,
      global_pairs_tsv=args.global_pairs,
      one_to_one_tsv=args.one_to_one,
      outdir=args.outdir,
      out_name=args.out_name
    )

