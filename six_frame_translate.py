#!/usr/bin/env python
"""
This is pretty basic, just takes a fasta files and does a six frame translation
for it. Currently it's set up for genetic codes 1, 4, and 11 (plus a version of
each that separately labels start codons).
"""

import sys
from Bio.Seq import Seq

if len(sys.argv) != 3:
    print "Usage: six_frame_translate.py <fasta file> <genetic code (either 1,\
           4 or 11 currently)>\n"
    quit()

fasta = open(sys.argv[1], "r")
gcode = sys.argv[2]

if gcode == "1":
    out = open(sys.argv[1] + ".six.frame", "w")
    aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
              "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
              "TGT":"C","TGC":"C","TGA":"*","TGG":"W", "CTT":"L","CTC":"L",
              "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
              "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
              "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
              "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
              "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
              "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A",
              "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
              "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

elif gcode == "1_start":
    out = open(sys.argv[1] + ".six.frame.startcodons", "w")
    aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"-","TCT":"S","TCC":"S",
              "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
              "TGT":"C","TGC":"C","TGA":"*","TGG":"W", "CTT":"L","CTC":"L",
              "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
              "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
              "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"-",
              "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
              "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
              "GTT":"V","GTC":"V","GTA":"V","GTG":"-","GCT":"A","GCC":"A",
              "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
              "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

elif gcode == "4":
    out = open(sys.argv[1] + ".six.frame", "w")
    aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
                  "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
                  "TGT":"C","TGC":"C","TGA":"W","TGG":"W", "CTT":"L","CTC":"L",
                  "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
                  "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
                  "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
                  "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
                  "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
                  "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A",
                  "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
                  "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

elif gcode == "4_start":
    out = open(sys.argv[1] + ".six.frame.startcodons", "w")
    aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"-","TCT":"S","TCC":"S",
                  "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
                  "TGT":"C","TGC":"C","TGA":"W","TGG":"W", "CTT":"L","CTC":"L",
                  "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
                  "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
                  "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"-",
                  "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
                  "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
                  "GTT":"V","GTC":"V","GTA":"V","GTG":"-","GCT":"A","GCC":"A",
                  "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
                  "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

elif gcode == "11":
    out = open(sys.argv[1] + ".six.frame", "w")
    aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
                  "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
                  "TGT":"C","TGC":"C","TGA":"*","TGG":"W", "CTT":"L","CTC":"L",
                  "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
                  "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
                  "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
                  "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
                  "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
                  "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A",
                  "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
                  "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

elif gcode == "11_start":
    out = open(sys.argv[1] + ".six.frame.startcodons", "w")
    aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"-","TCT":"S","TCC":"S",
                  "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
                  "TGT":"C","TGC":"C","TGA":"*","TGG":"W", "CTT":"L","CTC":"L",
                  "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
                  "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
                  "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"-",
                  "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
                  "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
                  "GTT":"V","GTC":"V","GTA":"V","GTG":"-","GCT":"A","GCC":"A",
                  "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
                  "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}


# My usual bit to read the fasta file into a dictionary
seq_dict = {}
for line in fasta:
    if line.startswith(">"):
        header = line.strip()
        seq_dict[header] = ""
    else:
        seq_dict[header] += line.strip()


def translate(offset, seq):
    """Does the actual translating"""

    for x in range(int(offset), len(seq), 3):
        codon = seq[x:x + 3]
        if len(codon) < 3:
            break
        elif "N" in codon:
            out.write("X")
        else:
            out.write("%s" % aa_dict[codon])
    return


# Goes through each sequence and calls 'translate' for all reading frames
for item in seq_dict:
    seq = seq_dict[item].upper()
    seq_rc = Seq(seq_dict[item].upper())
    seq_rc = seq_rc.reverse_complement()
    seq_rc = str(seq_rc)
    out.write("%s_frame1\n" % item)
    translate(0, seq)
    out.write("\n")
    out.write("%s_frame2\n" % item)
    translate(1, seq)
    out.write("\n")
    out.write("%s_frame3\n" % item)
    translate(2, seq)
    out.write("\n")
    out.write("%s_frame4\n" % item)
    translate(0, seq_rc)
    out.write("\n")
    out.write("%s_frame5\n" % item)
    translate(1, seq_rc)
    out.write("\n")
    out.write("%s_frame6\n" % item)
    translate(2, seq_rc)
    out.write("\n")
