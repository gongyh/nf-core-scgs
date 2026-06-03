# 在vambtools中会用到

def myoverwrite_matrix(matrix, mask):
    """Given a float32 matrix and Uint8 mask, does the same as setting the first
    rows of matrix to matrix[mask], but in-place.
    This is only important to save on memory.
    """
    i = 0
    j = 0
    matrixindex = 0
    length = matrix.shape[1]
    masklength = len(mask)

    # First skip to the first zero in the mask, since the matrix at smaller
    # indices than this should remain untouched.
    for i in range(masklength):
        if mask[i] == 0:
            break

    # If the mask is all true, don't touch array.
    if i == masklength:
        return masklength

    matrixindex = i

    for i in range(matrixindex, masklength):
        if mask[i] == 1:
            for j in range(length):
                matrix[matrixindex][j] = matrix[i][j]
            matrixindex += 1

    return matrixindex

def mykmercounts(bytesarray, k, counts):
    """Count tetranucleotides of contig and put them in counts vector.

    The bytearray is expected to be np.uint8 of bytevalues of the contig.
    Only values 64, 67, 71, 84 are accepted, all others are skipped.
    The counts is expected to be an array of 4^k 32-bit integers with value 0.
    """
    kmer = 0
    countdown = k-1
    contiglength = len(bytesarray)
    mask = (1 << (2 * k)) - 1
    lut = [
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    ]

    for i in range(contiglength):
        character = bytesarray[i]
        charvalue = lut[character]

        if charvalue == 4:
            countdown = k

        kmer = ((kmer << 2) | charvalue) & mask

        if countdown == 0:
            counts[kmer] += 1
        else:
            countdown -= 1
