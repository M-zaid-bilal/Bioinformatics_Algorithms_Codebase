from Hamming_Distance import hamming_distance

def neighbors(pattern: str, d: int) -> set:
    """Generates all k-mers within Hamming distance d from Pattern (d-neighborhood)."""
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}

    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)

    for suffix in suffix_neighbors:
        if hamming_distance(pattern[1:], suffix) < d:
            for nuc in ["A", "C", "G", "T"]:
                neighborhood.add(nuc + suffix)
        else:
            neighborhood.add(pattern[0] + suffix)

    return neighborhood


if __name__ == "__main__":
    pattern = "GGGGGGCGT"
    d=2
    result = neighbors(pattern, d)
    print(" ".join(result))