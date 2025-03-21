from Bio import SeqIO
import sys

def Quote_Remover(input_string):
    Out = ''
    for char in input_string:
        if char != '"':
            Out += char
    return Out

def All_Core_Finder(gene_presence_absence):
    f = open(gene_presence_absence, 'r')
    String1 = f.readline()
    List1 = String1.split(',')
    Total = len(List1)
    All = len(List1) - 14
    Core_IDs = []
    for line in f:
        List1 = line.split('","')
        Count = Quote_Remover(List1[3])
        if int(Count) == All:
            for entry in range(14, Total):
                Allele = Quote_Remover(List1[entry])
                Core_IDs.append(Allele)
    f.close()
    return Core_IDs

def Gene_List_Extractor(input_list, input_genome, output_genome):
    Genome = list(SeqIO.parse(input_genome, 'fasta'))
    Out = open(output_genome, 'w')
    for gene in Genome:
        if (gene.id in input_list):
            gene.id = gene.id + '_1'
            SeqIO.write(gene, Out, 'fasta')
    Out.close()

Genes = All_Core_Finder(sys.argv[1])
Gene_List_Extractor(Genes, sys.argv[2], 'Core_Genome.fasta')
