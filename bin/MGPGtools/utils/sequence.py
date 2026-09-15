def nucl_complement(sequence):
    trantab = str.maketrans("ACGTacgt", "TGCAtgca")
    transtr = sequence.translate(trantab)
    return transtr[::-1]
