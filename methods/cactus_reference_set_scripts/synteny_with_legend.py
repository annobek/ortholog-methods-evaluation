#!/usr/bin/env python3
"""
Enhanced Synteny Plot Generator with Legend
Shows which numbers correspond to which scaffold/chromosome names
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np
from collections import defaultdict
import seaborn as sns

def create_bezier_ribbon(x1_start, x1_end, y1, x2_start, x2_end, y2):
    """Create a smooth ribbon connecting two regions"""
    ctrl = abs(y2 - y1) * 0.4
    
    verts = [
        (x1_start, y1), (x1_start, y1 - ctrl), (x2_start, y2 + ctrl), (x2_start, y2),
        (x2_end, y2), (x2_end, y2 + ctrl), (x1_end, y1 - ctrl), (x1_end, y1),
        (x1_start, y1)
    ]
    codes = [Path.MOVETO] + [Path.CURVE4]*3 + [Path.LINETO] + [Path.CURVE4]*3 + [Path.CLOSEPOLY]
    
    return Path(verts, codes)

def plot_synteny_with_legend(data_files, species_pairs, output='synteny_with_legend.png', 
                             top_scaffolds=15, min_length=1000000,
                             show_scaffold_labels=True, color_scaffolds=True,
                             label_fontsize=8, show_legend=True,
                             legend_file='scaffold_legend.tsv'):
    """
    Create synteny plot with a legend file showing scaffold names
    
    Args:
        data_files: List of file paths to synteny data
        species_pairs: List of tuples (species1, species2) for each file
        output: Output filename for plot
        top_scaffolds: Number of top scaffolds to show per species
        min_length: Minimum synteny region length to display
        show_scaffold_labels: Whether to show scaffold numbers on scaffolds
        color_scaffolds: Whether to color scaffolds individually
        label_fontsize: Font size for scaffold labels
        show_legend: Whether to create a legend file
        legend_file: Output filename for legend text file
    """
    
    # Load and combine data
    all_dfs = []
    for fpath, (sp1, sp2) in zip(data_files, species_pairs):
        df = pd.read_csv(fpath, sep='\t')
        df['Species1'] = sp1
        df['Species2'] = sp2
        df['Region_Length'] = (df['Query_End'] - df['Query_Start'] + 
                               df['Target_End'] - df['Target_Start']) / 2
        all_dfs.append(df)
    
    combined = pd.concat(all_dfs, ignore_index=True)
    combined = combined[combined['Region_Length'] >= min_length].copy()
    
    print(f"Total regions after filtering (min length {min_length}): {len(combined)}")
    
    # Get all unique species
    all_species = sorted(set(combined['Species1'].unique()) | set(combined['Species2'].unique()))
    n_species = len(all_species)
    
    print(f"Species found: {all_species}")
    
    # Calculate scaffold lengths
    scaffold_lengths = defaultdict(dict)
    for _, row in combined.iterrows():
        sp1, sc1 = row['Species1'], row['Query_Scaffold']
        scaffold_lengths[sp1][sc1] = max(scaffold_lengths[sp1].get(sc1, 0), row['Query_End'])
        
        sp2, sc2 = row['Species2'], row['Target_Scaffold']
        scaffold_lengths[sp2][sc2] = max(scaffold_lengths[sp2].get(sc2, 0), row['Target_End'])
    
    # Select top scaffolds per species and create mapping
    top_scaffolds_dict = {}
    scaffold_order = {}
    scaffold_name_map = {}  # Maps number to actual scaffold name
    
    for species in all_species:
        scaffolds = scaffold_lengths[species]
        sorted_scaffolds = sorted(scaffolds.items(), key=lambda x: x[1], reverse=True)
        top_scaffolds_dict[species] = {s[0]: s[1] for s in sorted_scaffolds[:top_scaffolds]}
        
        # Create mapping: number -> scaffold name
        scaffold_order[species] = {}
        scaffold_name_map[species] = {}
        for idx, (scaffold_name, length) in enumerate(sorted_scaffolds[:top_scaffolds]):
            number = idx + 1
            scaffold_order[species][scaffold_name] = number
            scaffold_name_map[species][number] = {
                'name': scaffold_name,
                'length': length
            }
        
        print(f"{species}: {len(sorted_scaffolds)} scaffolds, showing top {min(top_scaffolds, len(sorted_scaffolds))}")
    
    # Create legend file if requested
    if show_legend:
        # Create TSV format legend
        legend_data = []
        for species in all_species:
            for num in sorted(scaffold_name_map[species].keys()):
                info = scaffold_name_map[species][num]
                legend_data.append({
                    'Label': num,
                    'Species': species,
                    'Scaffold_Name': info['name'],
                    'Length_bp': info['length']
                })
        
        # Save as TSV
        legend_df = pd.DataFrame(legend_data)
        legend_df.to_csv(legend_file, sep='\t', index=False)
        
        print(f"\n✓ Legend saved to: {legend_file}")
        print(f"  Format: TSV with columns [Label, Species, Scaffold_Name, Length_bp]")
        print(f"  Note: Labels are numbered within each species (1=longest scaffold)")
    
    # Filter data to only include top scaffolds
    def is_in_top(row):
        return (row['Query_Scaffold'] in top_scaffolds_dict[row['Species1']] and
                row['Target_Scaffold'] in top_scaffolds_dict[row['Species2']])
    
    combined = combined[combined.apply(is_in_top, axis=1)].copy()
    print(f"Regions after scaffold filtering: {len(combined)}")
    
    # Generate color palette for scaffolds
    max_scaffolds = max(len(scaffold_order[sp]) for sp in all_species)
    
    if color_scaffolds:
        scaffold_colors = {}
        base_colors = sns.color_palette("husl", max_scaffolds)
        
        for species in all_species:
            scaffold_colors[species] = {}
            for scaffold, idx in scaffold_order[species].items():
                scaffold_colors[species][scaffold] = base_colors[(idx-1) % len(base_colors)]
    else:
        scaffold_colors = {sp: {sc: 'lightgray' for sc in scaffold_order[sp]} 
                          for sp in all_species}
    
    # Calculate positions for scaffolds
    scaffold_positions = {}
    for species in all_species:
        scaffolds = top_scaffolds_dict[species]
        total_len = sum(scaffolds.values())
        gap = total_len * 0.02
        
        pos = 0
        positions = {}
        for scaffold in sorted(scaffolds.keys(), key=lambda x: scaffolds[x], reverse=True):
            positions[scaffold] = {
                'start': pos,
                'end': pos + scaffolds[scaffold],
                'length': scaffolds[scaffold],
                'index': scaffold_order[species][scaffold],
                'name': scaffold
            }
            pos += scaffolds[scaffold] + gap
        
        scaffold_positions[species] = {
            'scaffolds': positions,
            'total': pos - gap
        }
    
    # Create plot
    fig, ax = plt.subplots(figsize=(26, 8))
    
    # Y positions for species
    y_positions = {sp: n_species - i for i, sp in enumerate(all_species)}
    scaffold_height = 0.2
    
    # Draw scaffolds
    for species, y_pos in y_positions.items():
        total = scaffold_positions[species]['total']
        
        for scaffold, pos in scaffold_positions[species]['scaffolds'].items():
            x_start = pos['start'] / total
            x_end = pos['end'] / total
            width = x_end - x_start
            
            color = scaffold_colors[species][scaffold]
            
            rect = patches.Rectangle(
                (x_start, y_pos - scaffold_height/2), width, scaffold_height,
                linewidth=1.5, edgecolor='black', facecolor=color, alpha=0.7
            )
            ax.add_patch(rect)
            
            # Add scaffold label
            if show_scaffold_labels:
                label_x = x_start + width / 2
                label_y = y_pos
                label_text = str(pos['index'])
                
                ax.text(label_x, label_y, label_text,
                       ha='center', va='center',
                       fontsize=label_fontsize, fontweight='bold',
                       color='black',
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor='white', 
                                edgecolor='none',
                                alpha=0.7))
    
    # Draw ribbons
    for idx, row in combined.iterrows():
        sp1, sp2 = row['Species1'], row['Species2']
        y1, y2 = y_positions[sp1], y_positions[sp2]
        
        q_scaffold = row['Query_Scaffold']
        t_scaffold = row['Target_Scaffold']
        
        if (sp1 not in scaffold_positions or 
            q_scaffold not in scaffold_positions[sp1]['scaffolds'] or
            sp2 not in scaffold_positions or
            t_scaffold not in scaffold_positions[sp2]['scaffolds']):
            continue
        
        total1 = scaffold_positions[sp1]['total']
        total2 = scaffold_positions[sp2]['total']
        
        q_offset = scaffold_positions[sp1]['scaffolds'][q_scaffold]['start']
        t_offset = scaffold_positions[sp2]['scaffolds'][t_scaffold]['start']
        
        x1_start = (q_offset + row['Query_Start']) / total1
        x1_end = (q_offset + row['Query_End']) / total1
        x2_start = (t_offset + row['Target_Start']) / total2
        x2_end = (t_offset + row['Target_End']) / total2
        
        ribbon_path = create_bezier_ribbon(
            x1_start, x1_end, y1 - scaffold_height/2,
            x2_start, x2_end, y2 + scaffold_height/2
        )
        
        ribbon_color = scaffold_colors[sp1][q_scaffold]
        
        ribbon = patches.PathPatch(
            ribbon_path,
            facecolor=ribbon_color,
            alpha=0.4,
            edgecolor='none'
        )
        ax.add_patch(ribbon)
    
    # Format plot
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.5, n_species + 0.5)
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels([s.replace('_', ' ').title() for s in all_species], 
                       fontsize=13, fontweight='bold')
    ax.set_xlabel('')
    ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    
    title = f'Synteny Plot - {len(combined)} Regions across {n_species} Species'
    #if show_legend:
    #    title += f'\n(See {legend_file} for scaffold names)'
    
    plt.title(title, fontsize=16, pad=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved to: {output}")
    
    return fig, ax, scaffold_name_map


if __name__ == "__main__":
    print("Enhanced Synteny Plot with Legend")
    print("=" * 70)
    print("\nThis version creates:")
    print("  1. The synteny plot (PNG)")
    print("  2. A legend file (TXT) showing what each number represents")
    print("\nUsage:")
    print("""
from synteny_with_legend import plot_synteny_with_legend

fig, ax, legend = plot_synteny_with_legend(
    data_files=['file1.tsv', 'file2.tsv', ...],
    species_pairs=[('sp1', 'sp2'), ...],
    output='synteny.png',
    legend_file='scaffold_legend.txt',  # Shows number -> name mapping
    top_scaffolds=15,
    min_length=1000000
)
""")
