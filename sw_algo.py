def sw_algo(seq1, seq2, st):
    r = seq1.getSize()
    c = seq2.getSize()

    # First loop: Creating an empty matrix with zeros
    sm = []
    for _ in range(r):
        sm.append([0] * c)

    maxs = 0
    maxsp = (0, 0)

    # The rest of the code remains unchanged
    for i in range(r):
        for j in range(c):
            curr_1 = seq1.getItr(i)
            curr_2 = seq2.getItr(j)

            ms = st[curr_1][curr_2]
            gappen = st["gap"]["gap"]

            diags = sm[i - 1][j - 1] + ms
            ups = sm[i - 1][j] + gappen
            les = sm[i][j - 1] + gappen

            sm[i][j] = max(0, diags, ups, les)

            if sm[i][j] > maxs:
                maxs = sm[i][j]
                maxsp = (i, j)

    alseq1 = ""
    alseq2 = ""
    i, j = maxsp

    while i >= 0 and j >= 0 and sm[i][j] > 0:
        curr_1 = seq1.getItr(i)
        curr_2 = seq2.getItr(j)

        ms = st[curr_1][curr_2]
        gappen = st["gap"]["gap"]

        if sm[i][j] == sm[i - 1][j - 1] + ms:
            alseq1 = curr_1 + alseq1
            alseq2 = curr_2+ alseq2
            i -= 1
            j -= 1
        elif sm[i][j] == sm[i - 1][j] + gappen:
            alseq1 = curr_1 + alseq1
            alseq2 = '-' + alseq2
            i -= 1
        else:
            alseq1 = '-' + alseq1
            alseq2 = curr_2 + alseq2
            j -= 1

    return alseq1, alseq2