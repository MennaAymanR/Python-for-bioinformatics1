N=["A","C","G","T"]
def DNAtool(DNA_Seq):
    for i in DNA_Seq:
        DNA=DNA_Seq.upper()

        if i not in DNA:
            return False
    return DNA


def count_nuc(Seq):
    nuc_dict={"A":0, "C":0, "G":0, "T":0}
    for i in Seq:
        nuc_dict[i]+=1
    return nuc_dict

def transcriptionDNA(DNASEQ):
    dna=""

    for i in DNASEQ:
        if i =="A":
            dna +="T"
        elif i=="T":
            dna+="A"
        elif i =="G":
            dna+="C"
        elif i == "C":
            dna+="G"
        else:
            print("invalid")
    print(dna) 

def GCpercentage(DNA_SEQ):
    total=len(DNA_SEQ)
    C=DNA_SEQ.count("C")
    G=DNA_SEQ.count("G")
    GC=((C+G)/total)*100
    print(GC)


def AminoAcidSearch(AminoAcids):
    AminoAcids=AminoAcids.upper()
    aminoacids=["TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA"
    ,"TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT"
    ,"CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA"
    ,"ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT"
    ,"GCC","GCA","GCG","GAT","GAC","GGT","GGC","GGA","GGG"]
    aminoacidsname=["Phe(F)","Phe(F)","Leu(L)","Leu(L)","Ser(S)","Ser(S)","Ser(S)","Ser(S)"
    ,"Tyr(Y)","Tyr(Y)","Stop","Stop","Cys(C)","Cys(C)", "Stop", "Trp(W)" ,"leu(L)","leu(L)"
    ,"leu(L)","leu(L)","Pro(P)","Pro(P)","Pro(P)","Pro(P)","His(H),His(H)","Gln(G)","Gln(G)"
    ,"Arg(R)","Arg(R)","Arg(R)","Arg(R)","Ile(I)","Ile(I)","Ile(I)","Met(M)","Thr(T)","Thr(T)"
    ,"Thr(T)","Thr(T)","Asn(N)","Asn(N)","Lys(K)","Lsy(K)","Ser(S)","Ser(S)","Ser(S)","Arg(R)"
    ,"Arg(R)","Val(v)","Val(v)","Val(v)","Val(v)","Asp(D)","Asp(D)","Glu(E)","Glu(E)","Gly(G)"
    ,"Gly(G)","Gly(G)","Gly(G)"]
    i=0
    find=False
    for X in aminoacids:
        if(X==AminoAcids):
            print("The amino acid you entered is:"+ aminoacidsname[i])    
            find=True
        i+=1        
    if (find==False):
        print("Wrong Amino acid")    
user_input=input("Please enter DNA seq")
def translate(user_input):
    dictCodon = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M','ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T','AAC': 'N',
        'AAT': 'N', 'AAA': 'K', 'AAG': 'K','AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L','CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q','CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V','GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E','GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S','TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_','TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = " "
    if len(user_input) % 3 == 0:
        for i in range(0, len(user_input), 3):
            codon = user_input[i:i + 3]
            protein += dictCodon[codon]
        return protein