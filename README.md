Updated to work with python3.9+ and output a file of the user's choosing

# BPP - branch point prediction 
BPP is a method that can identify branch point based on only the intron sequence.  

### How BPP works?

BPP predicts the branch point sequence by integrating the degenerative motif of BPS and PPT characteristics. Specifically, BPP uses a mixture model to infer the BPS motif and a set of weighted octanucleotides to estimate the contribution of the 65kDa subunit of U2AF (U2AF65). 

### Detailed information on BPSP and citation
A paper describing BPP is under review.  


### Dependent libraries or software
- Python3.9

### Quickstart

Usage: 
```bash
python3 BP_PPT.py -b -p -i -r -o -h
```

#### Required:
```bash
-b, --pwm file    STR     The file including PWM of BPS
-p, --ppt file    STR     The file including the PPT score
-i, --FASTA file  STR     The file including the fasta sequence
-r, --report nu   INT     The reported sites; default=1; 0: print all positions
-o, --output      STR     Output file to write the results
-h, --help
```

#### Example:
```bash
python3 BP_PPT_updated.py \
-b demo/pwmBP_human.txt \
-p demo/scPPT_human.txt \
-i demo/example.fa \
-o demo/example_output.txt
```

#### Format of the output file:
- `id`:      ID of the intron
- `bps`:     the branch point sequence
- `bp_pos`:  the position of the branch point relative to 3'SS upstream
- `sc_bps`:  the score of the BPS
- `sc_ppt`:  the score of the PPT
- `sc`:      the score of the BPS and PPT
- `zsc_bps`: the z-score of the BPS
- `zsc_ppt`: the z-score of the PPT
- `zsc`:     the z-score of the BPS and PPT
