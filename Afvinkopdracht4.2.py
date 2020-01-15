import re
import time
 # op dit moment worden van het fasta bestand de sequentie en de headers geopend
 # de sequentie staat in een string en de headers zijn onderdelen van een list


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
                    print("Please insert a fasta file")
            else:
                seq += line.strip()
        return headers, seqs

    except FileNotFoundError:
        print('file not found')
        exit()


def is_dnasequentie(nietdna,line):
    """
    This function checks if the sequence is a protein by checking if it contains and characters that do not represent a protein
    in standard fasta format encoding, it does this by using the search function form regular expressions
    input: string with characters that are not proteins and a string
    output: print stating whetever the sequence is a protein or not

    """
    seconds = time.time()
    if re.search(nietdna,line,flags=0) == True:
        print("Regular expression does not confirm this to be DNA")
    else:
        print("Regular expression confirms this to be DNA")
    return seconds
def is_dnasequentie_iterate(nietdna,line):
    dnacheck = True
    second = time.time()
    for i in line:
        if nietdna in line:
            dnacheck = False
    return dnacheck, second


def main():
    startfunctie = time.time()
    bestandsnaam = "Mus_Musculus_GRCm38_Chromosome1.fa"
    headers, seqs = readfile(bestandsnaam)
    nietdna = 'B'
    startrexp = is_dnasequentie(nietdna,seqs)
    dna, startiterate = is_dnasequentie_iterate(nietdna,seqs)
    iteratetime = startrexp - startfunctie
    rexptime = startiterate - startrexp
    if dna == True:
        print("Iteration confirms this as DNA")
    else:
        print("Iteration does not confirm this to be DNA")
    print("The time it took to Iterate = ",iteratetime,"\n","The time it took to use regular experss = ",rexptime)

main()