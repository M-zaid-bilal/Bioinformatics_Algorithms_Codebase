import pandas as pd
import matplotlib.pyplot as plt
import logomaker

# Import functions directly from your motif_profile module
from motif_profile import profile_matrix_df, column_entropy


def information_content_df(motifs):
    """
    Computes Information Content matrix for logomaker:
    Height(nuc, col) = Frequency(nuc, col) * (2 - Column_Entropy)
    """
    df_profile = profile_matrix_df(motifs)
    
    # Calculate (2 - H) for each column position
    column_entropies = df_profile.apply(column_entropy, axis=0)
    info_content_scales = 2.0 - column_entropies
    
    # Scale nucleotide frequencies by (2 - H) and transpose so rows are positions (0..k-1)
    df_info = (df_profile * info_content_scales).T
    return df_info


if __name__ == "__main__":
    # 10 NF-κB motifs transcribed from motif_profile setup
    nf_kb_motifs = [
        "TCGGGGgTTTtt",
        "cCGGtGAcTTaC",
        "aCGGGGATTTtC",
        "TtGGGGAcTTtt",
        "aaGGGGAcTTCC",
        "TtGGGGAcTTCC",
        "TCGGGGATTcat",
        "TCGGGGATTcCt",
        "TaGGGGAacTaC",
        "TCGGGtATaaCC"
    ]

    # Compute Information Content DataFrame
    df_info = information_content_df(nf_kb_motifs)

    # Plot Sequence Logo using logomaker
    fig, ax = plt.subplots(figsize=(10, 3.5))
    
    logo = logomaker.Logo(
        df_info,
        ax=ax,
        color_scheme='classic',   # Classic DNA scheme (A: green, C: blue, G: orange/yellow, T: red)
        vpad=0.05
    )

    # Clean styling
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    
    ax.set_ylabel('Information Content (bits)', fontsize=12)
    ax.set_xlabel('Column Position', fontsize=12)
    ax.set_title('NF-κB Motif Logo', fontsize=14, fontweight='bold')
    
    # Set 1-indexed string tick labels (1 to k)
    k = len(nf_kb_motifs[0])
    ax.set_xticks(range(k))
    ax.set_xticklabels([str(i) for i in range(1, k + 1)])
    
    plt.tight_layout()
    plt.savefig("nf_kb_motif_logo.png", dpi=300)
    print("Motif logo saved successfully to 'nf_kb_motif_logo.png'")
    plt.show()