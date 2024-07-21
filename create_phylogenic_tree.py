def create_phylogenic_tree(species):
    def check_valid_sequences(species, keys):
        length = species[keys[0]].getSize()
        for k in keys:
            if length > species[k].getSize():
                length = species[k].getSize()
        return length
        
    def minmatch(seq1, seq2, length):
        return sum(1 for i in range(length) if seq1.getItr(i) != seq2.getItr(i))
    
    def generate_initial_scores(species, keys, num_species, length):
        mismatch_score = {}
        for k in keys:
            mismatch_score[k] = {}
       
        for i in range(num_species):
            for j in range(i+1, num_species):
                s1 = keys[i]
                s2 = keys[j]
                mismatch_score[s1][s2] = minmatch(species[s1], species[s2], length)

        return mismatch_score
    
    def find_min_pos(score, num_species, seq_len):
        min_value = seq_len
        min_pos = (0, 0)

        for i in range(num_species):
            for j in range(i+1, num_species):
                s1 = keys[i]
                s2 = keys[j]
                if score[s1][s2] < min_value:
                    min_value = score[s1][s2]
                    min_pos = (i, j)
        
        return min_pos
    
    def update_score(score, keys, pos1, pos2):
        m = "[" + keys[pos1] + ", " + keys[pos2] + "]";

        mismatch_score = score.copy()
        mismatch_score[m] = {}

        for i in range(num_species):
            if i < pos1:
                mismatch_score[keys[i]][m] = (mismatch_score[keys[i]][keys[pos1]] + mismatch_score[keys[i]][keys[pos2]]) / 2
                del mismatch_score[keys[i]][keys[pos1]]
                del mismatch_score[keys[i]][keys[pos2]]
            elif i > pos1 and i < pos2:
                mismatch_score[keys[i]][m] = (mismatch_score[keys[pos1]][keys[i]] + mismatch_score[keys[i]][keys[pos2]]) / 2
                del mismatch_score[keys[i]][keys[pos2]]
            elif pos2 < i:
                mismatch_score[keys[i]][m] = (mismatch_score[keys[pos1]][keys[i]] + mismatch_score[keys[pos2]][keys[i]]) / 2

        del mismatch_score[keys[pos1]]
        del mismatch_score[keys[pos2]]
        return mismatch_score

    keys = list(species)
    num_species = len(keys)
    seq_length = check_valid_sequences(species, keys)

    if seq_length == 0:
        print("INVALID: All Sequences should be of a length greater than 0.")
        return None

    if num_species == 1:
        return keys[0]
    else:
        score = generate_initial_scores(species, keys, num_species, seq_length)

        while(num_species > 1):
            (i, j) = find_min_pos(score, num_species, seq_length)
            score = update_score(score, keys, i, j)
            keys = list(score)
            num_species = len(keys)
            
        return keys[0]