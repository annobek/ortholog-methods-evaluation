import pandas as pd
import numpy as np
from pathlib import Path

# Extract pairs that are in forward and reverse files
# Resolve ambiguous pairs with max best GRIS score forward and reverse


# Full names of the species (for search for genes in syntenic gene tables)
full_names = {
    "can": "Canis lupus familiaris",
    "mac": "Macaca fascicularis",
    "rat": "Rattus norvegicus",
    "mus": "Mus musculus", 
    "ar": "Arabidopsis thaliana",
    "card": "Cardamine hirsuta",
    "dmel": "Drosophila melanogaster",
    "dsim": "Drosophila simulans",
    "dsec": "Drosophila sechellia"
}


def _clean_cols(df, cols):
    # Processes DataFrame, converting the columns to string to avoid errors from string-numeric mismatch
    # Parameters: DataFrame; column names
    # Return: Processed DataFrame

    for c in cols:
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip()
        else:
            raise KeyError(f"Missing column '{c}'")
    return df



def extract_from_blast(blast_file):
    # Extracts from BLAST/DIAMOND following information: Protein_ID for species1; Protein_ID for species2; geomean
    # Parameters: BLAST/DIAMOND file path (str)
    # Output: Number of valid DIAMOND rows after cleanup; first three lines
    # Return: Protein_ID for species1; Protein_ID for species2; geomean (DataFrame)
    colnames = ["qseqid","sseqid","qstart","qend","qlen","sstart","send","slen","pident","evalue","bitscore","overlap","geometric_mean"]
    blast_df = pd.read_csv(blast_file, sep=",", names=colnames, engine="python")
    blast_df = blast_df[["qseqid","sseqid","geometric_mean"]]
    blast_df.columns = blast_df.columns.str.strip()
    blast_df = blast_df.dropna(subset=["qseqid", "sseqid", "geometric_mean"])
    blast_df = blast_df[blast_df["qseqid"] != blast_df["sseqid"]].copy()
    blast_df["geometric_mean"] = pd.to_numeric(blast_df["geometric_mean"], errors="coerce")
    blast_df = blast_df.dropna(subset=["geometric_mean"])

    return blast_df



def extract_from_map(map_file):
    # Extracts from the map following information: Protein_ID; Gene_ID
    # Parameters: map file path (str)
    # Return: Protein_ID; Gene_ID (dict)
  
    map_df = pd.read_csv(map_file, sep="\t")

    # Address different column order
    if map_df.columns[0].lower().startswith("protein"):
        map_df = map_df[[map_df.columns[1], map_df.columns[0]]]  # swap columns


    # For mammals and plants
    #map_df = map_df[["GeneID", "ProteinID"]]
    #map_information = map_df[["GeneID", "ProteinID"]]

    # For drosophila
    map_df = map_df[["ProteinID", "GeneID"]]
    map_information = map_df[["ProteinID", "GeneID"]]

    # Create dictionary for faster lookup: ProteinID -> GeneID
    protein_to_gene_map = dict(zip(map_information['ProteinID'], map_information['GeneID']))

    return protein_to_gene_map



def map_proteins_to_genes(geomean_df, prot2gene):
     # Maps protein IDs to gene IDs from the obtained DataFrame with extract_from_blast(),
     # but keep original protein IDs (qseqid,sseqid) in output
    # Parameters: 
    # - Protein_ID for species1; Protein_ID for species2; geomean (DataFrame) (from extract_from_blast())
    # - Protein_ID; Gene_ID (dict) (from extract_from_map())
    # Output: Number of mapped gene pairs
    # Return: Query_Gene, Target_Gene, Query_Protein, Target_Protein, geometric_mean (DataFrame)

    gene_pairs = geomean_df.copy()
    gene_pairs["Query_Gene"]   = gene_pairs["qseqid"].map(prot2gene)
    gene_pairs["Target_Gene"]  = gene_pairs["sseqid"].map(prot2gene)
    gene_pairs = gene_pairs.dropna(subset=["Query_Gene", "Target_Gene"])

    gene_pairs = gene_pairs.rename(columns={"qseqid": "Query_Protein", "sseqid": "Target_Protein"})
    gene_pairs = gene_pairs[["Query_Gene", "Target_Gene", "Query_Protein", "Target_Protein", "geometric_mean"]]
    # normalize dtypes
    gene_pairs = _clean_cols(gene_pairs, ["Query_Gene", "Target_Gene", "Query_Protein", "Target_Protein"])
    return gene_pairs




def load_gene_table(genes_table, species1, species2):
    # Loads synteny gene table, filters to the two species and return two dicts:
    # gene->region for species1_full and gene->region for species2_full
    
    sp1 = full_names[species1]
    sp2 = full_names[species2]

    df = pd.read_csv(genes_table, sep="\t", dtype=str)
    df = _clean_cols(df, ["Region_ID", "Species", "Gene-ID"])
    df = df[df["Species"].isin([sp1, sp2])].copy()

    gene2region_sp1 = dict(zip(df.loc[df["Species"] == sp1, "Gene-ID"], df.loc[df["Species"] == sp1, "Region_ID"]))
    gene2region_sp2 = dict(zip(df.loc[df["Species"] == sp2, "Gene-ID"], df.loc[df["Species"] == sp2, "Region_ID"]))
    return gene2region_sp1, gene2region_sp2



def load_region_gris(gris_tsv):
    # Load region GRIS file with columns including: Region_ID, reg_GRIS
    # Returns dict {Region_ID -> reg_GRIS (float)} (missing -> empty dict)
    
    path = Path(gris_tsv)
    if not path.exists():
        return {}
    
    gris_df = pd.read_csv(path, sep="\t", dtype=str)

    if "Region_ID" not in gris_df.columns or "reg-GRIS" not in gris_df.columns:
        # also accept 'reg_GRIS' if the header is underscore
        if "reg_GRIS" in gris_df.columns:
            gris_df = gris_df.rename(columns={"reg_GRIS": "reg-GRIS"})
        else:
            return {}
        
    gris_df = _clean_cols(gris_df, ["Region_ID", "reg-GRIS"])
    gris_df["reg-GRIS"] = pd.to_numeric(gris_df["reg-GRIS"], errors="coerce")

    return dict(zip(gris_df["Region_ID"], gris_df["reg-GRIS"]))



def reciprocal(forward_df, reverse_df):
    # Keeps only pairs that are in forward and reverse.
    # Input columns for forward and reverse: Query_Gene, Target_Gene, Query_Protein, Target_Protein, geometric_mean
    # Output: reciprocal pairs df


    # Normalize ID types to string for safe comparison
    forward_df["Query_Gene"] = forward_df["Query_Gene"].astype(str)
    forward_df["Target_Gene"] = forward_df["Target_Gene"].astype(str)
    reverse_df["Query_Gene"] = reverse_df["Query_Gene"].astype(str)
    reverse_df["Target_Gene"] = reverse_df["Target_Gene"].astype(str)

    # Add a helper tuple key to identify pairs
    forward_df["pair"] = list(zip(forward_df["Query_Gene"], forward_df["Target_Gene"]))
    reverse_df["pair"] = list(zip(reverse_df["Target_Gene"], reverse_df["Query_Gene"]))  # reversed direction

    # Intersect based on these tuple keys
    reciprocal_pairs = set(forward_df["pair"]) & set(reverse_df["pair"])

    reciprocal_df = forward_df[forward_df["pair"].isin(reciprocal_pairs)].copy()
    reciprocal_df.drop(columns=["pair"], inplace=True)

    print(f"Reciprocal pairs retained: {len(reciprocal_df):,}")
    return reciprocal_df



def annotate_regions_and_gris(rbh_df,
                              genes_fwd_tsv, genes_rev_tsv,
                              gris_fwd_tsv=None, gris_rev_tsv=None,
                              sp1="can", sp2="mac"):
   
    # Add Region IDs from forward & reverse gene tables and region GRIS values (if provided).
    
    # region mappings (forward)
    gene2region_fwd_sp1, gene2region_fwd_sp2 = load_gene_table(genes_fwd_tsv, sp1, sp2)

    # region mappings (reverse)
    gene2region_rev_sp2, gene2region_rev_sp1 = load_gene_table(genes_rev_tsv, sp2, sp1)  # note: order flipped

    out = rbh_df.copy()
    # Regions from forward table context
    out["Region_Fwd_Query"]  = out["Query_Gene"].map(gene2region_fwd_sp1)
    out["Region_Fwd_Target"] = out["Target_Gene"].map(gene2region_fwd_sp2)
    # Regions from reverse table context
    out["Region_Rev_Query"]  = out["Query_Gene"].map(gene2region_rev_sp1)
    out["Region_Rev_Target"] = out["Target_Gene"].map(gene2region_rev_sp2)

    # If both genes co-locate to the same region in a given context, use that region id; else NaN
    out["Region_Fwd"] = np.where(out["Region_Fwd_Query"] == out["Region_Fwd_Target"],
                                 out["Region_Fwd_Query"], np.nan)
    out["Region_Rev"] = np.where(out["Region_Rev_Query"] == out["Region_Rev_Target"],
                                 out["Region_Rev_Query"], np.nan)

    # region GRIS lookup
    reg_gris_fwd = load_region_gris(gris_fwd_tsv) if gris_fwd_tsv else {}
    reg_gris_rev = load_region_gris(gris_rev_tsv) if gris_rev_tsv else {}

    out["region_gris_fwd"] = out["Region_Fwd"].map(reg_gris_fwd) if reg_gris_fwd else np.nan
    out["region_gris_rev"] = out["Region_Rev"].map(reg_gris_rev) if reg_gris_rev else np.nan

    # total region support (used only as tiebreaker)
    out["region_gris_total"] = out[["region_gris_fwd", "region_gris_rev"]].sum(axis=1, min_count=1)

    return out



def pick_one_to_one(sorted_candidates):
    
    # Greedy 1:1 selection:
    #  - Input must be sorted in desired priority (desc).
    #  - Walk top->down, keep a pair if both genes unused.
    
    used_q, used_t = set(), set()
    kept, discarded = [], []
    for _, r in sorted_candidates.iterrows():
        q, t = r["Query_Gene"], r["Target_Gene"]
        if (q not in used_q) and (t not in used_t):
            kept.append(r)
            used_q.add(q)
            used_t.add(t)
        else:
            discarded.append(r)
            
    return pd.DataFrame(kept), pd.DataFrame(discarded)




def build_one_to_one_for_pair(
    sp1, sp2,
    genes_fwd_tsv,
    genes_rev_tsv,
    gris_fwd_tsv,
    gris_rev_tsv,
    diamond_pairs_gene_df,
    out_tsv
):
    
    # Create 1:1 ortholog table for (sp1, sp2)
   
    # Build a base view of DIAMOND gene pairs for this pair
    base_df = diamond_pairs_gene_df.rename(columns={"Query": "Query_Gene", "Target": "Target_Gene"}).copy()

    # Restrict DIAMOND pairs to this species pair using the gene tables
    gene2region_fwd_sp1, gene2region_fwd_sp2 = load_gene_table(genes_fwd_tsv, sp1, sp2)
    gene2region_rev_sp2, gene2region_rev_sp1 = load_gene_table(genes_rev_tsv, sp2, sp1)

    sp1_genes_fwd = set(gene2region_fwd_sp1.keys())
    sp2_genes_fwd = set(gene2region_fwd_sp2.keys())
    sp2_genes_rev = set(gene2region_rev_sp2.keys())
    sp1_genes_rev = set(gene2region_rev_sp1.keys())

    # Forward = real DIAMOND rows with (sp1 -> sp2)
    fwd_df = base_df[
        base_df["Query_Gene"].isin(sp1_genes_fwd) &
        base_df["Target_Gene"].isin(sp2_genes_fwd)
    ].copy()

    # Reverse = real DIAMOND rows with (sp2 -> sp1) 
    rev_df = base_df[
        base_df["Query_Gene"].isin(sp2_genes_rev) &
        base_df["Target_Gene"].isin(sp1_genes_rev)
    ].copy()


    # For reporting: forward-only and reverse-only (discarded) pairs
    fwd_pairs = set(zip(fwd_df["Query_Gene"], fwd_df["Target_Gene"]))
    rev_pairs = set(zip(rev_df["Target_Gene"], rev_df["Query_Gene"]))  # normalized to (sp1, sp2)
    forward_only = fwd_pairs - rev_pairs
    reverse_only = rev_pairs - fwd_pairs

    print(f"{sp1}-{sp2}: forward-only (discarded) = {len(forward_only):,}, "
          f"reverse-only (discarded) = {len(reverse_only):,}")

    # Reciprocal filtering
    rbh = reciprocal(fwd_df, rev_df)   
    print(f"{sp1}-{sp2}: reciprocal pairs = {len(rbh):,}")

    # Annotate regions and GRIS and resolve 1:many 
    rbh = annotate_regions_and_gris(
        rbh,
        genes_fwd_tsv, genes_rev_tsv,
        gris_fwd_tsv, gris_rev_tsv,
        sp1=sp1, sp2=sp2
    )

    # How often are forward vs reverse region GRIS different?
    tmp = rbh.dropna(subset=["region_gris_fwd", "region_gris_rev"]).copy()
    if not tmp.empty:
        tmp["gris_delta"] = (tmp["region_gris_fwd"] - tmp["region_gris_rev"]).abs()
        n_diff = int((tmp["gris_delta"] > 1e-9).sum())
        print(f"{sp1}-{sp2}: pairs with different fwd vs rev region GRIS = {n_diff:,} / {len(tmp):,}")

    # Sort by total reciprocal GRIS (forward + reverse) before 1:many resolution
    rbh_sorted = rbh.sort_values(
        by=["region_gris_total", "Query_Gene", "Target_Gene"],
        ascending=[False, True, True]
    )

    # Resolve 1:many by greedy selection
    one2one, inparalogs = pick_one_to_one(rbh_sorted)

    # Select and order columns
    cols = [
        "Query_Gene", "Target_Gene",
        "Query_Protein", "Target_Protein",
        "Region_Fwd", "Region_Rev",
        "region_gris_fwd", "region_gris_rev", "region_gris_total"
    ]
    for c in cols:
        if c not in one2one.columns:
            one2one[c] = np.nan
        if c not in inparalogs.columns:
            inparalogs[c] = np.nan

    one2one = one2one[cols]
    inparalogs = inparalogs[cols]

    # Save results
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    one2one.to_csv(out_tsv, sep="\t", index=False)
    inparalogs_out = Path(out_tsv).with_name(f"{Path(out_tsv).stem.replace('one2one', 'inparalogs')}.tsv")
    inparalogs.to_csv(inparalogs_out, sep="\t", index=False)

    print(f"{sp1}-{sp2}: 1:1 pairs = {len(one2one):,}  -> {out_tsv}")
    print(f"{sp1}-{sp2}: in-paralogs = {len(inparalogs):,}  -> {inparalogs_out}")




def main():
    
    # Input files

    # Uncomment path needed 

    # diamond and map

    # For mammals
    #DIAMOND_GEOMEAN = "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/geometric_mean.txt"
    #PROT_GENE_MAP = "/storage/EasyVectorOmics/synteny_algorithm/material/mammalian/gene_to_protein_map.tsv"

    # For plants
    #DIAMOND_GEOMEAN = "/storage/EasyVectorOmics/synteny_algorithm/material/cardamine/geometric_mean.txt"
    #PROT_GENE_MAP = "/storage/EasyVectorOmics/synteny_algorithm/material/cardamine/gene_to_protein_map.tsv"

    # For drosophila
    DIAMOND_GEOMEAN = "/storage/EasyVectorOmics/synteny_algorithm/material/drosophila/geometric_mean.txt"
    PROT_GENE_MAP = "/storage/EasyVectorOmics/synteny_algorithm/material/drosophila/gene_to_protein_map.tsv"

    # species pairs

    # For mammals
    '''
    PAIRS = [
        {
            "sp1": "can", "sp2": "mac",
            "genes_fwd_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/can_mac_synteny_tables_genes_fix.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/mac_can_synteny_tables_genes_fix.tsv",
            "gris_fwd_tsv":  "results/region_gris/can_mac_region_gris_mean_fix.tsv",
            "gris_rev_tsv":  "results/region_gris/mac_can_region_gris_mean_fix.tsv",
            "out_tsv":       "results/reciprocal_pairs/can_mac_one2one_fix.tsv",
        },
        {
            "sp1": "can", "sp2": "rat",
            "genes_fwd_tsv": "results/genes_in_blocks/mammalia/region_gene_tables//can_rat_synteny_tables_genes_fix.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/mammalia/region_gene_tables//rat_can_synteny_tables_genes_fix.tsv",
            "gris_fwd_tsv":  "results/region_gris/can_rat_region_gris_mean_fix.tsv",
            "gris_rev_tsv":  "results/region_gris/rat_can_region_gris_mean_fix.tsv",
            "out_tsv":       "results/reciprocal_pairs/can_rat_one2one_fix.tsv",
        },
        {
            "sp1": "can", "sp2": "mus",
            "genes_fwd_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/can_mus_synteny_tables_genes_fix.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/mus_can_synteny_tables_genes_fix.tsv",
            "gris_fwd_tsv":  "results/region_gris/can_mus_region_gris_mean_fix.tsv",
            "gris_rev_tsv":  "results/region_gris/mus_can_region_gris_mean_fix.tsv",
            "out_tsv":       "results/reciprocal_pairs/can_mus_one2one_fix.tsv",
        },
        {
            "sp1": "mac", "sp2": "rat",
            "genes_fwd_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/mac_rat_synteny_tables_genes_fix.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/rat_mac_synteny_tables_genes_fix.tsv",
            "gris_fwd_tsv":  "results/region_gris/mac_rat_region_gris_mean_fix.tsv",
            "gris_rev_tsv":  "results/region_gris/rat_mac_region_gris_mean_fix.tsv",
            "out_tsv":       "results/reciprocal_pairs/mac_rat_one2one_fix.tsv",
        },
        {
            "sp1": "mac", "sp2": "mus",
            "genes_fwd_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/mac_mus_synteny_tables_genes_fix.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/mus_mac_synteny_tables_genes_fix.tsv",
            "gris_fwd_tsv":  "results/region_gris/mac_mus_region_gris_mean_fix.tsv",
            "gris_rev_tsv":  "results/region_gris/mus_mac_region_gris_mean_fix.tsv",
            "out_tsv":       "results/reciprocal_pairs/mac_mus_one2one_fix.tsv",
        },
        {
            "sp1": "rat", "sp2": "mus",
            "genes_fwd_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/rat_mus_synteny_tables_genes_fix.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/mammalia/region_gene_tables/mus_rat_synteny_tables_genes_fix.tsv",
            "gris_fwd_tsv":  "results/region_gris/rat_mus_region_gris_mean_fix.tsv",
            "gris_rev_tsv":  "results/region_gris/mus_rat_region_gris_mean_fix.tsv",
            "out_tsv":       "results/reciprocal_pairs/rat_mus_one2one_fix.tsv",
        }
    ]
    '''

    # For plants
    '''
    PAIRS = [
        {
            "sp1": "ar", "sp2": "card",
            "genes_fwd_tsv": "results/genes_in_blocks/brassicaceae/region_gene_tables/ar_card_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/brassicaceae/region_gene_tables/card_ar_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/ar_card_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/card_ar_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/ar_card_one2one.tsv",
        },
        {
            "sp1": "card", "sp2": "ar",
            "genes_fwd_tsv": "results/genes_in_blocks/brassicaceae/region_gene_tables/card_ar_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/brassicaceae/region_gene_tables/ar_card_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/card_ar_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/ar_card_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/card_ar_one2one.tsv",
        }
    ]
    '''

    # For drosophila

    PAIRS = [
        {
            "sp1": "dmel", "sp2": "dsec",
            "genes_fwd_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dmel_dsec_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsec_dmel_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/dmel_dsec_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/dsec_dmel_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/dmel_dsec_one2one.tsv",
        },
        {
            "sp1": "dsec", "sp2": "dmel",
            "genes_fwd_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsec_dmel_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dmel_dsec_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/dsec_dmel_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/dmel_dsec_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/dsec_dmel_one2one.tsv",
        },
         {
            "sp1": "dmel", "sp2": "dsim",
            "genes_fwd_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dmel_dsim_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsim_dmel_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/dmel_dsim_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/dsim_dmel_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/dmel_dsim_one2one.tsv",
        },
        {
            "sp1": "dsim", "sp2": "dmel",
            "genes_fwd_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsim_dmel_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dmel_dsim_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/dsim_dmel_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/dmel_dsim_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/dsim_dmel_one2one.tsv",
        },
        {
            "sp1": "dsec", "sp2": "dsim",
            "genes_fwd_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsec_dsim_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsim_dsec_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/dsec_dsim_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/dsim_dsec_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/dsec_dsim_one2one.tsv",
        },
        {
            "sp1": "dsim", "sp2": "dsec",
            "genes_fwd_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsim_dsec_synteny_tables_genes.tsv",
            "genes_rev_tsv": "results/genes_in_blocks/drosophila/region_gene_tables/dsec_dsim_synteny_tables_genes.tsv",
            "gris_fwd_tsv":  "results/region_gris/dsim_dsec_region_gris_mean.tsv",
            "gris_rev_tsv":  "results/region_gris/dsec_dsim_region_gris_mean.tsv",
            "out_tsv":       "results/reciprocal_pairs/dsim_dsec_one2one.tsv",
        }
    ]
    



    print("Loading DIAMOND…")
    diamond_df = extract_from_blast(DIAMOND_GEOMEAN)
    print(f"DIAMOND rows: {len(diamond_df):,}")

    print("Loading protein->gene map…")
    protein2gene = extract_from_map(PROT_GENE_MAP)
    print(f"Map entries:  {len(protein2gene):,}")

    print("Mapping proteins to genes…")
    gene_pairs = map_proteins_to_genes(diamond_df, protein2gene)
    print(f"Gene-level pairs: {len(gene_pairs):,}")

    

    for cfg in PAIRS:
        build_one_to_one_for_pair(
            cfg["sp1"], cfg["sp2"],
            cfg["genes_fwd_tsv"], cfg["genes_rev_tsv"],
            cfg.get("gris_fwd_tsv"), cfg.get("gris_rev_tsv"),
            gene_pairs,
            cfg["out_tsv"]
        )

if __name__ == "__main__":
    main()
