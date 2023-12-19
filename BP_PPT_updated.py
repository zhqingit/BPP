#!/usr/bin/python

import sys
import getopt
import math

# Function to read and process PWM (Position Weight Matrix) from a file
def bppt_get_pwm(pwmf):
    index = 0
    PWMBP = {}
    with open(pwmf, 'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            PWMBP[index] = {}
            line = tmp.split('\t')
            PWMBP[index]['A'] = float(line[1])
            PWMBP[index]['C'] = float(line[2])
            PWMBP[index]['G'] = float(line[3])
            PWMBP[index]['T'] = float(line[4])
            index += 1
    return PWMBP

# Function to read and process PPT (Polypyrimidine Tract) scores from a file
def bppt_get_ppt(pptf):
    PPTS = {}
    with open(pptf, 'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            line = tmp.split('\t')
            PPTS[line[0]] = float(line[4])
    return PPTS

# Function to calculate the score of a BPS (Branch Point Sequence)
def bppt_bpscore(cbps):
    cbpsc = 1
    cbpsN = list(cbps)
    for i in range(len(cbpsN)):
        if i in PWMBP:
            cbpsc *= PWMBP[i][cbpsN[i]]
        else:
            print(f"Error: the input bps is longer than the PWM")
            break
    return cbpsc

# Function to generate and score all possible BPS combinations
def bppt_get_bpscore(l, basebp):
    nn = ['A', 'C', 'G', 'T']  # Nucleotides
    NNS = ['A', 'C', 'G', 'T']  # Nucleotide combinations with fixed length
    NN = []

    for i in range(2, l + 1):
        for nns in NNS:
            for n in nn:
                newN = nns + n
                NN.append(newN)
        NNS = NN
        NN = []

    basebpsc = bppt_bpscore(basebp)
    cBPSC = {}
    for nn in NNS:
        cBPSC[nn] = bppt_bpscore(nn) / basebpsc

    return cBPSC

# Function to get the AGEZ (AG Exclusion Zone)
def bppt_get_AGEZ(seq, offset=12, maxL=-1):
    sL = len(seq)
    ss = seq[0:-offset].split("AG")
    if maxL > 0:
        pAG = sL - maxL
    else:
        pAG = sL - (offset + len(ss[-1]) + 14)
    if pAG < 0:
        pAG = 0
    return pAG

# Function to calculate PPT scores
def bppt_get_pptsc(pptS, lppt, l_max_ppt, baseppt):
    end = len(pptS) - lppt
    if end > l_max_ppt:  # If the ppt sequence is larger than l_max_ppt, only check the l_max_ppt
        end = l_max_ppt
    pptbasesc = (l_max_ppt - lppt + 1) * PPTS[baseppt]
    pptSC = 0
    for i in range(end):
        cppts = pptS[i:i + lppt]
        pptSC += PPTS[cppts]
    if pptbasesc == 0:
        return pptSC
    else:
        return (pptSC / pptbasesc)


# Function to calculate the distribution probability
def bppt_dis_pro(pAG, offset=22, interval=4):
 # 4 is best
    if abs(pAG - offset) > 700:
        return 0
    else:
        return 1 / (math.exp((abs(pAG - offset)) / interval) + 1)  # Best formula


# Main function to get BPPT (Branch Point and Polypyrimidine Tract) scores
def bppt_get_BPPTsc(seq, maxL, baseppt):
    lmotif = 7
    lppt = 8
    l_max_ppt = 20
    pstart = bppt_get_AGEZ(seq=seq, offset=12, maxL=maxL)
    sL = len(seq)

    totbpsc = 0
    totpptsc = 0
    totsc = 0
    npos = 0

    orinp = []
    orsc = []
    # Main loop to calculate scores
    for ipos in range(pstart, sL - 3 - lmotif):  # The BPS could be close to the 3' end with only one nucleotide distance
        pAG = sL - ipos - 5
        bpS = seq[ipos:ipos + lmotif]  # BPS sequence
        bpSC = cBPSC[bpS]  # BPS sequence score

        pptSC = 0
        dis3 = pAG - 1  # Distance of BPS last nucleotide to the 3' end
        if dis3 > lppt + 3:  # 3 means the AG + one nucleotide distance, no U2AF can bind the downstream sequence of BPS if this is true
            pptS = seq[ipos + lmotif:sL - 3]  # BPS sequence
            pptSC = bppt_get_pptsc(pptS, lppt, l_max_ppt, baseppt)

        SC = (bpSC * pptSC)  # Combined score
        inp = f"{bpS}\t{pAG}\t{bpSC}\t{pptSC}\t{SC}"

        # Sorting the scores
        if not orinp:
            orinp.append(inp)
        else:
            flag_in = 0
            for i in range(len(orinp)):
                line = orinp[i].split("\t")
                scold = float(line[4])
                if scold < SC:
                    orinp.insert(i, inp)
                    flag_in = 1
                    break
            if flag_in == 0:
                orinp.append(inp)

        totsc += SC
        totbpsc += bpSC
        totpptsc += pptSC
        npos += 1

    msc = totsc / npos
    mbpsc = totbpsc / npos
    mpptsc = totpptsc / npos

    # Calculating standard deviation and z-scores
    dsc = []
    dbpsc = []
    dpptsc = []
    sdsc = 0
    sdbpsc = 0
    sdpptsc = 0
    for i in range(len(orinp)):
        line = orinp[i].split("\t")
        sc = float(line[4])
        bpsc = float(line[2])
        pptsc = float(line[3])
        dd = sc - msc
        dsc.append(dd)
        sdsc += dd * dd
        dd = bpsc - mbpsc
        dbpsc.append(dd)
        sdbpsc += dd * dd
        dd = pptsc - mpptsc
        dpptsc.append(dd)
        sdpptsc += dd * dd

    sdsc = math.sqrt(sdsc / npos)
    sdbpsc = math.sqrt(sdbpsc / npos)
    sdpptsc = math.sqrt(sdpptsc / npos)

    zsc = [d / sdsc for d in dsc]
    zbps = [d / sdbpsc for d in dbpsc]
    zppt = [d / sdpptsc for d in dpptsc]

    return orinp, zbps, zppt, zsc

# Function to print results to the console or write to a file
def bppt_print(idd, orinp, zbps, zppt, zsc, REPORTN, outfile=None):
    output = []
    output.append("#id\tbps\tbp_pos\tsc_bps\tsc_ppt\tsc\tzsc_bps\tzsc_ppt\tzsc")
    for i in range(min(len(orinp), REPORTN)):
        output.append(f"{idd}\t{orinp[i]}\t{zbps[i]}\t{zppt[i]}\t{zsc[i]}")

    if outfile:
        with open(outfile, 'a') as out:  # Open in append mode
            out.write('\n'.join(output) + '\n')
    else:
        for line in output:
            print(line)

# Function to display help message
def bppt_print_help():
    print('Usage: BP_PPT.py -b -p -i -r -o -h')
    print(' -b, --pwm file    STR     The file including PWM of BPS')
    print(' -p, --ppt file    STR     The file including the PPT score')
    print(' -i, --FASTA file  STR     The file including the fasta sequence')
    print(' -r, --report nu   INT     The reported sites; default=1; 0: print all positions')
    print(' -o, --output      STR     Output file to write the results')
    print(' -h, --help')
    sys.exit()

# Function to parse command line arguments
def bppt_get_options(argv):
    try:
        opts, args = getopt.getopt(argv, "b:r:p:i:o:h", ["bppt"])
    except getopt.GetoptError as err:
        sys.exit(2)

    pwmf = pptf = fastaF = outputFile = ""
    repn = 1
    for opt, arg in opts:
        if opt == '-h':
            bppt_print_help()
        elif opt == '-b':
            pwmf = arg
        elif opt == '-p':
            pptf = arg
        elif opt == '-r':
            repn = arg
        elif opt == '-i':
            fastaF = arg
        elif opt == '-o':
            outputFile = arg
    if not pwmf or not pptf or not fastaF:
        bppt_print_help()

    return pwmf, pptf, fastaF, int(repn), outputFile

# Main function
def main(argv):
    mL = 7
    basebp = "TACTAAC"
    baseppt = "TTTTTTTT"

    global REPORTN
    pwmf, pptf, fastaF, REPORTN, outputFile = bppt_get_options(argv)
    REPORTN = int(REPORTN)
    global PWMBP, PPTS, cBPSC
    PWMBP = bppt_get_pwm(pwmf)
    PPTS = bppt_get_ppt(pptf)
    cBPSC = bppt_get_bpscore(mL, basebp)

    idd = ""
    seq = ""
    with open(fastaF, 'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith(">"):
                if not idd:
                    idd = tmp
                else:
                    orinp, zbps, zppt, zsc = bppt_get_BPPTsc(seq, maxL=-1, baseppt=baseppt)
                    bppt_print(idd, orinp, zbps, zppt, zsc, REPORTN, outputFile)
                    idd = tmp
                    seq = ""
            else:
                seq += tmp
        orinp, zbps, zppt, zsc = bppt_get_BPPTsc(seq, maxL=-1, baseppt=baseppt)
        bppt_print(idd, orinp, zbps, zppt, zsc, REPORTN, outputFile)


if __name__ == "__main__":
    main(sys.argv[1:])