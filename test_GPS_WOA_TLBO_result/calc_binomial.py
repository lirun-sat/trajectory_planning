

def calc_binomial(n, k):
    LUT = [[1],
           [1, 1],
           [1, 2, 1],
           [1, 3, 3, 1],
           [1, 4, 6, 4, 1],
           [1, 5, 10, 10, 5, 1],
           [1, 6, 15, 20, 15, 6, 1]]
    while n >= len(LUT):
        s = len(LUT)
        prev = s - 1
        nextRow = [1]
        for i in range(1, s):
            nextRow.append(LUT[prev][i-1] + LUT[prev][i])
        nextRow.append(1)
        LUT.append(nextRow)

    return LUT[n][k]

