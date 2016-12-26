# bpp 
BPP is a method that can identify branch point based on only the intron sequence.  

### How BPP works?

BPP predicts the branch point sequence by integrating the degenerative motif of BPS and PPT characteristics. Specifically, BPP uses a mixture model to infer the BPS motif and a set of weighted octanucleotides to estimate the contribution of the 65kDa subunit of U2AF (U2AF65). 

### Detailed information on BPSP and citation

A paper describing BPP is under review.  


### Dependent libraries or software

- Python3.5

### Quickstart

Usage: BP_PPT.py -b -p -i -r -h
#### Required:
- `-b, --pwm file    STR     The file including PWM of BPS`
- `-p, --ppt file    STR     The file including the PPT score`
- `-i, --FASTA file  STR     The file including the fasta sequence`
- `-r, --report nu   INT     The reported sites; default=1; 0: print all positions`
- `-h, --help`

