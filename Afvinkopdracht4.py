import re

def readfile(bestands_naam):
    """
    This function reads a fasta file and returns two variables with the headers and sequence of the file as strings
    input: protein fasta file
    output: Headers and sequence
    """
    try:
        bestand = open(bestands_naam)
        headers = []

        seqs = ""
        seq = ""
        for line in bestand:
            line = line.strip()
            if ">" in line:
                if seq != "":
                    seqs+=seq
                    seq = ""
                headers.append(line)
                if headers == '':
                    print('Please insert a fasta file')
            else:
                seq += line.strip()
        return headers, seqs

    except FileNotFoundError:
        print('file not found')
        exit()

def itereer_en_controleer(sequenties,nieteiwit):
    for i in sequenties:
        if i == "":
            print("The sequence is empty")
            exit()
    else:
        consensus(sequenties)



def consensus(seq):
    """
    This function checks if the consensus pattern from the p53 protein is in the a sequence
    if the pattern is found, a boolean with True is returned else a False is returned
    input: sequence
    output: Boolean
    """

    p53pattern = "MCNSSCMGGMNRR|MCNSSCVGGMNRR"
    x = re.search(str(p53pattern), seq)
    print(x)
    if (x):
        print("YES! We have a match!")
    else:
        print("No match")
    return x


def is_eiwitsequentie(nieteiwit,line):
    """
    This function checks if the sequence is a protein by checking if it contains and characters that do not represent a protein
    in standard fasta format encoding, it does this by using the search function form regular expressions
    input: string with characters that are not proteins and a string
    output: print stating whetever the sequence is a protein or not

    """
    if re.search(nieteiwit,line,flags=0) == True:
        print("It is not a protein")
    else:
        print("This is a valid protein")


def main():
    bestandsnaam = "Mus_musculus.GRCm38.pep.all.fa"
    headers, seqs = readfile(bestandsnaam)
    nieteiwit = 'B'
    itereer_en_controleer(seqs,nieteiwit)

main()