def bm_approx_mat(pat, tx, maxmis):
    patlen = pat.getSize()
    txlen = tx.getSize()

    if patlen == 0 or txlen == 0 or patlen > txlen:
        return []

    pos = []

    i = 0
    while i <= txlen - patlen:
        mismc = 0

        # Separated loop for mismc calculation
        for j in range(patlen):
            if i + j < txlen and pat.getItr(j) != tx.getItr(i + j):
                mismc += 1

        if mismc <= maxmis:
            pos.append(i)

        if i + patlen < txlen:
            lastoc = pat.find(tx.getItr(i + patlen))
            i += max(1, patlen - lastoc)
        else:
            break

    return pos

