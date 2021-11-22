import numpy as np
import sys



def get_score(aa1, aa2, amino_acids_list, scoring_matrix):
    return int(scoring_matrix[amino_acids_list.index(aa1)][amino_acids_list.index(aa2)])


def fill_tables(row_index, col_index, partial_score_table, alignment_algorithm):
    #diagonal (d)
    matching_score = get_score(seq1[col_index-1], seq2[row_index-1], amino_acids_list, scoring_matrix)
    diagonal_score = partial_score_table[row_index-1][col_index-1] + matching_score
    
    #horizontal (h)
    if (direction_table[row_index][col_index-1] != 0) and ('h' in direction_table[row_index][col_index-1]):  #gap extension
        horizontal_score = partial_score_table[row_index][col_index-1] + gap_extension_penalty
    else: #gap opening: 'v' or 'd' or empty(first row)
        horizontal_score = partial_score_table[row_index][col_index-1] + gap_opening_penalty
        
    #vertical (v)
    if (direction_table[row_index-1][col_index] != 0) and ('v' in direction_table[row_index-1][col_index]):  #gap extension
        vertical_score = partial_score_table[row_index-1][col_index] + gap_extension_penalty
    else: #gap opening: 'h' or 'd' or empty(first column)
        vertical_score = partial_score_table[row_index-1][col_index] + gap_opening_penalty
    
    possible_dirs = {'d': diagonal_score, 'h': horizontal_score, 'v': vertical_score}
    max_score_dir = ','.join([k for k, v in possible_dirs.items() if v == max(possible_dirs.values())])
    
    if alignment_algorithm == 'global':
        partial_score_table[row_index, col_index] = max(possible_dirs.values())
        direction_table[row_index, col_index] = max_score_dir
    elif alignment_algorithm == 'local':
        if max(possible_dirs.values()) > 0:
            partial_score_table[row_index, col_index] = max(possible_dirs.values())
            direction_table[row_index, col_index] = max_score_dir
        else:
            partial_score_table[row_index, col_index] = 0
            direction_table[row_index, col_index] = 0
            
            
def traceback(alignment_algorithm, row, col, table, path):
    
    if alignment_algorithm == 'global':
        if row == 0 and col == 0:
            return path

        elif row == 0 and col != 0:
            path = path + ['h']*col
            return path

        elif row != 0 and col == 0:
            path = path + ['v']*row
            return path
    
    elif alignment_algorithm == 'local':
        if table[row][col] == 0:
            return path
    
    directions = table[row][col].split(',')
    direction = directions[0] 
    
    if direction == 'd':
        path.append('d')
        return traceback(alignment_algorithm, row-1, col-1, table, path)
    elif direction == 'h':
        path.append('h')
        return traceback(alignment_algorithm, row, col-1, table, path)
    elif direction == 'v':
        path.append('v')
        return traceback(alignment_algorithm, row-1, col, table, path)
    
    

    
# ARG0 : FILE NAME OF SEQUENCES
protein_sequences_file_name = sys.argv[1]  
with open(protein_sequences_file_name) as file:
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]
    
seq1 = lines[0]
seq2 = lines[1]


# ARG1 : ALIGNMENT ALGORITHM
alignment_algorithm = sys.argv[2] 
if alignment_algorithm not in ['local', 'global']:
    raise Exception("Please provide a valid alignment algorithm: either 'local' or 'global' ")
    
    
# ARG2 : FILE NAME OF SCORING MATRIX 
scoring_matrix_file_name = sys.argv[3] 
with open(scoring_matrix_file_name) as file:
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]

amino_acids_list = lines[0].split()
scoring_matrix = []
for i, line in enumerate(lines):
    if i == 0:
        pass
    else:
        scoring_matrix.append(line.split()[1:])

        
# ARG3 : GAP OPENING PENALTY
try:
    gap_opening_penalty = int(sys.argv[4])
    if gap_opening_penalty >= 0:
        raise Exception("Please provide a negative integer value for the gap opening penalty.")
except:
    raise Exception("Please provide a negative integer value for the gap opening penalty.")

    
# ARG3 : GAP EXTENSION PENALTY    
try:
    gap_extension_penalty = sys.argv[5] 
    if gap_extension_penalty == '-':
        gap_extension_penalty = gap_opening_penalty
    else:
        gap_extension_penalty = int(gap_extension_penalty)
        if gap_extension_penalty >= 0:
            raise Exception("Please provide a negative integer value for the gap extension penalty.")
except:
    raise Exception("Please provide a negative integer value for the gap extension penalty.")

    
    
    
# CREATE PARTIAL-SCORE AND DIRECTION TABLES
partial_score_table = np.zeros([len(seq2)+1, len(seq1)+1])
direction_table = np.zeros([len(seq2)+1, len(seq1)+1], dtype='object')

if alignment_algorithm == 'global':
    for i in range(partial_score_table.shape[0]):
        if (i == 0): partial_score_table[i][0] = 0
        else:partial_score_table[i][0] = gap_opening_penalty + (i-1) * gap_extension_penalty
    for i in range(partial_score_table.shape[1]):
        if (i == 0): partial_score_table[0][i] = 0
        else: partial_score_table[0][i] = gap_opening_penalty + (i-1) * gap_extension_penalty
            


# FILL PARTIAL-SCORE AND DIRECTION TABLES
for row_index in range(1, partial_score_table.shape[0]):
    for col_index in range(1, partial_score_table.shape[1]):
        fill_tables(row_index, col_index, partial_score_table, alignment_algorithm)
        
        
        
# TRACEBACK TO GET ALIGNING PATH
if alignment_algorithm == 'local':    
    max_raw_alignment_score = partial_score_table.max()
    max_raw_alignment_score_index = np.argwhere(partial_score_table == max_raw_alignment_score)[0]
    a_path = traceback(alignment_algorithm, max_raw_alignment_score_index[0], max_raw_alignment_score_index[1], direction_table, path=[])
elif alignment_algorithm == 'global': 
    a_path = traceback(alignment_algorithm, direction_table.shape[0]-1, direction_table.shape[1]-1, direction_table, path=[])
    
    
    
# OBTAIN ALIGNED SEQUENCES
if alignment_algorithm == 'global':
    seq1_rev = seq1[::-1]
    seq2_rev = seq2[::-1]
elif alignment_algorithm == 'local':
    seq1_rev = seq1[-(len(seq1) - max_raw_alignment_score_index[1] + 1)::-1]
    seq2_rev = seq2[-(len(seq2) - max_raw_alignment_score_index[0] + 1)::-1]

aligned_seq1 = []
aligned_seq2 = []

counter_seq1 = 0
for direction in a_path:
    if (direction == 'd') or (direction == 'h'):
        aligned_seq1.append(seq1_rev[counter_seq1])
        counter_seq1 += 1
    elif direction == 'v':
        aligned_seq1.append('-')
        
counter_seq2 = 0
for direction in a_path:
    if (direction == 'd') or (direction == 'v'):
        aligned_seq2.append(seq2_rev[counter_seq2])
        counter_seq2 += 1
    elif direction == 'h':
        aligned_seq2.append('-')
            
aligned_seq1 = aligned_seq1[::-1]
aligned_seq2 = aligned_seq2[::-1]



# CREATE INTERMEDIATE LINE
matching_arr = []
total_matched = 0
for as1, as2 in zip(aligned_seq1, aligned_seq2):
    if as1 == as2:
        matching_arr.append('|')
        total_matched += 1
    else:
        matching_arr.append(' ')
        
        
        
# WRITE TO OUTPUT.txt (aligned sequences, raw alignment score, percent identity)
with open("output.txt", "w") as text_file:
    if alignment_algorithm == 'global':
        text_file.write(f'raw alignment score: {partial_score_table[partial_score_table.shape[0]-1][partial_score_table.shape[1]-1]}')
        text_file.write('\n')
        text_file.write(f'percent identity between the two aligned sequences: {(total_matched*100)/len(aligned_seq1)}')
    elif alignment_algorithm == 'local':
        text_file.write(f'raw alignment score: {max_raw_alignment_score}')
        text_file.write('\n')
        text_file.write(f'percent identity between the two aligned sequences: {(total_matched*100)/len(aligned_seq1)}')        
    
    text_file.write('\n')
    text_file.write('\n')
    text_file.write(''.join(aligned_seq1))
    text_file.write('\n')
    text_file.write(''.join(matching_arr))
    text_file.write('\n')
    text_file.write(''.join(aligned_seq2))

print("output.txt has been created successfully.")