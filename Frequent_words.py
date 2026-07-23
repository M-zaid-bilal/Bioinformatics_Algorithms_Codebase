import PatternCount_Chapter_1_Task1 as pc

def FrequentWords(Text, k):
    # FrequentPatterns ← an empty set
    FrequentPatterns = set()
    
    # Boundary check to prevent errors with large k or empty text
    if k <= 0 or k > len(Text):
        return FrequentPatterns

    n = len(Text) - k + 1
    
    # count(i) ← PatternCount(Text, Pattern)
    counts = [pc.pattern_count(Text, Text[i:i+k]) for i in range(n)]
    
    # maxCount ← maximum value in array Count
    maxCount = max(counts) if counts else 0
    
    # for i ← 0 to |Text| − k
    for i in range(n):
        # if Count(i) = maxCount
        if counts[i] == maxCount:
            # add Text(i, k) to FrequentPatterns
            FrequentPatterns.add(Text[i:i+k])
            
    # Sets automatically remove duplicates, so no explicit removal is needed
    return FrequentPatterns