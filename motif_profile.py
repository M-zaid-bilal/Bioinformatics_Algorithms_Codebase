import math
import pandas as pd


def count_matrix_df(motifs):
    """
    Computes the 4 x k Count matrix for a list of motif strings.
    Returns a pandas DataFrame with nucleotides 'A', 'C', 'G', 'T' as rows.
    """
    k = len(motifs[0])
    counts = {nuc: [0] * k for nuc in "ACGT"}

    for seq in motifs:
        for col, nuc in enumerate(seq.upper()):
            counts[nuc][col] += 1

    return pd.DataFrame(counts, index=range(k)).T


def profile_matrix_df(motifs):
    """
    Computes the 4 x k Profile matrix (frequencies) for a list of motif strings.
    Divides the Count DataFrame by t (number of motif rows).
    """
    t = len(motifs)
    df_count = count_matrix_df(motifs)
    return df_count / t


def consensus_string_df(motifs):
    """
    Forms Consensus(Motifs) using pandas idxmax to extract 
    the nucleotide with the maximum count in each column.
    """
    df_count = count_matrix_df(motifs)
    return "".join(df_count.idxmax(axis=0).values)


def column_entropy(col_series):
    """
    Computes Shannon Entropy for a pandas Series of nucleotide probabilities.
    Handles 0 * log2(0) = 0 safely.
    """
    H = 0.0
    for p in col_series:
        if p > 0:
            H -= p * math.log2(p)
    return H


def entropy_motifs_df(motifs):
    """
    Computes total entropy by applying the entropy function across 
    each column of the Profile DataFrame.
    """
    df_profile = profile_matrix_df(motifs)
    column_entropies = df_profile.apply(column_entropy, axis=0)
    return column_entropies.sum()


if __name__ == "__main__":
    # Exact 10 NF-κB motifs transcribed directly from the image grid (all 12 chars)
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

    print("--- Count Matrix (pandas) ---")
    df_count = count_matrix_df(nf_kb_motifs)
    print(df_count)

    print("\n--- Profile Matrix (pandas) ---")
    df_profile = profile_matrix_df(nf_kb_motifs)
    print(df_profile)

    print("\n--- Consensus String ---")
    consensus = consensus_string_df(nf_kb_motifs)
    print(consensus)

    print("\n--- Total Matrix Entropy ---")
    total_entropy = entropy_motifs_df(nf_kb_motifs)
    print(f"Entropy: {total_entropy:.8f}")