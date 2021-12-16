**Run command**<br>
`python main.py sequences.txt global scoring_matrices\BLOSUM62.txt -10 -5`

<br>

**Command Line Arguments** <br>
- *Arg1 (sequences.txt) :* text file including 2 sequences of any length, written in 2 different lines. 
- *Arg2 (global) :* alignment algorithm: local or global. 
- *Arg3 (scoring_matrices\BLOSUM62.txt) :* text file including scoring matrix* 
- *Arg4 (-10) :* gap opening penalty (a negative integer). 
- *Arg5 (-5) :* gap extension penalty (a negative integer).

<br>

**Output** <br>
An output.txt file will be created, which includes the alignments of the sequences in 3-line format, corresponding raw alignment score and percent identity.

**Requirements** <br>
numpy==1.21.4
