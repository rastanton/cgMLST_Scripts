import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from decimal import *
getcontext().prec = 3
from multiprocessing import Pool
import dendropy
import shutil
import subprocess
import math

##Written by Richard Stanton (rstanton@cdc.gov; stanton.rich@gmail.com)
##Requires Python, biopython, dendropy, and blat
##Usage: $ python General_cgMLST_Exe.py Path/To/cgMLST_Alleles.fasta Path/To/cgMLST_Individual_Alleles_Fasta/ Path/to/new/assemblies/
##https://github.com/rastanton/PA_cgMLST

def Folder_Name_Remover(input_string):
    Temp_string =''
    for characters in input_string:
        if characters == '/':
            Temp_string = ''
        else:
            Temp_string = Temp_string + characters
    return Temp_string

def cgMLST_Similarity(match_file1, match_file2):
    f = open(match_file1, 'r')
    g = open(match_file2, 'r')
    string1 = f.readline()
    string2 = g.readline()
    Matches = 0
    while string1 != '':
        if string1 == string2:
            Matches = Matches + 1
            string2 = g.readline()
            string1 = f.readline()
        else:
            string2 = g.readline()
            string1 = f.readline()
    return Matches

##def cgMLST_Distance(match_file1, match_file2):
##    f = open(match_file1, 'r')
##    g = open(match_file2, 'r')
##    string1 = f.readline()
##    string2 = g.readline()
##    Total = 0
##    Matches = 0
##    while string1 != '':
##        if match_file1 == match_file2:
##            break
##        elif string1 == string2 and string1[-3:] != '_0\n' and string1[-7:] != '_Novel\n':
##            Total = Total + 1
##            Matches = Matches + 1
##            string2 = g.readline()
##            string1 = f.readline()
##        elif string1[-3:] == '_0\n' or string2[-3:] == '_0\n':
##            string2 = g.readline()
##            string1 = f.readline()
##        else:
##            Total = Total + 1
##            string2 = g.readline()
##            string1 = f.readline()
##    return int(Total - Matches)

def cgMLST_Lister(input_match_file):
    Out = []
    f = open(input_match_file, 'r')
    for line in f:
        Info = line.strip()
        Info = Info.split(',')
        Out.append(Info)
    f.close()
    return Out

def cgMLST_Distance(match_file1, match_file2):
    List1 = cgMLST_Lister(match_file1)
    List2 = cgMLST_Lister(match_file2)
    Difference = 0
    for entry in range(len(List1)):
        Diff = 1
        for data in List1[entry]:
            if (data in List2[entry]):
                Diff = 0
        Difference += Diff
    return Difference

def Jaccard_Difference(match_file1, match_file2):
    f = open(match_file1, 'r')
    g = open(match_file2, 'r')
    string1 = f.readline()
    string2 = g.readline()
    Union = 0
    Intersection = 0
    while string1 != '':
        if string1 == string2 and string1[-3:] != '_0\n' and string1[-7:] != '_Novel\n':
            Union = Union + 1
            Intersection = Intersection + 1
            string2 = g.readline()
            string1 = f.readline()
        elif string1[-3:] == '_0\n' or string2[-3:] == '_0\n':
            string2 = g.readline()
            string1 = f.readline()
        else:
            Union = Union + 2
            string2 = g.readline()
            string1 = f.readline()
    Jaccard = Decimal(Intersection) / Decimal(Union)
    Difference = 1 - Jaccard
    return Difference

def cgMLST_Jaccard_Difference_Matrix_Maker(input_index, output_file):
    File_List = glob.glob(input_index)
    File_List.sort()
    f = open(output_file, 'w')
    String1 = ''
    f.write('\t')
    for file_names in File_List:
        String1 = String1 + Folder_Name_Remover(file_names) + '\t'
    f.write(String1[0:-1] + '\n')
    for files in File_List:
        String1 = Folder_Name_Remover(files) + '\t'
        for match_files in File_List:
            String1 = String1 + str(Jaccard_Difference(files, match_files)) + '\t'
        f.write(String1[0:-1] + '\n')
    f.close()
    
def List_Overlap(List1, List2):
    for entry in List1:
        if (entry in List2):
            return True
    return False

def Allele_Compare(Allele1, Allele2, Multifasta):
    Genes = SeqIO.to_dict(SeqIO.parse(Multifasta, 'fasta'))
    Gene1 = str(Genes[Allele1].seq)
    Gene2 = str(Genes[Allele2].seq)
    Gene2_Rev = str(Genes[Allele2].seq.reverse_complement())
    Diff = 1
    if (Gene1 in Gene2) or (Gene2 in Gene1) or (Gene2_Rev in Gene1) or (Gene1 in Gene2_Rev):
        Diff = 0
    return Diff

def Allele_Differences(file1, file2, fasta_folder):
    Difference = 0
    if file1 == file2:
        return Difference
    f = open(file1, 'r')
    g = open(file2, 'r')
    for line1 in f:
        line1 = line1.strip()
        Diff = 0
        line2 = g.readline()
        line2 = line2.strip()
        if ('_0' in line1) or ('_0' in line2):
            continue
        elif line1 != line2:
            List1 = line1.split(',')
            List2 = line2.split(',')
            Allele = List1[0].split('_')[0:-1]
            Allele = '_'.join(Allele)
            if List_Overlap(List1, List2) == True:
                break
            else:
                All_Diffs = []
                for entry in List1:
                    for entry2 in List2:
                        Diff = Allele_Compare(entry, entry2, fasta_folder + Allele + '_alleles.fasta')
                        All_Diffs.append(Diff)
                if (0 in All_Diffs):
                    Diff = 0
                else:
                    Diff = 1
        Difference += Diff
    f.close()
    g.close()
    return Difference

def Allele_Difference_List(input_list):
##    print(str(input_list))
    Out = Allele_Differences(input_list[0], input_list[1], input_list[2])
    return Out

def Folder_List_Maker(input_list, folder):
    New = []
    for entry in range(len(input_list)):
        for entry2 in range(entry, len(input_list)):
            New.append([input_list[entry], input_list[entry2], folder])
    return New

def List_Matrix(length):
    Out = []
    Length = int(math.sqrt(length * 2))
    New = []
    for entry in range(Length):
        New.append(entry)
    Out.append(New)
    for entry in range(1, Length):
        New = []
        for entry2 in range(Length):
##            if (entry2) == 0:
##                New.append(entry)
            if entry2 < entry:
                New.append(Out[entry2][entry])
            elif entry == entry2:
                New.append(Out[entry - 1][-1] + 1)
            else:
                New.append(New[-1] + 1)
        Out.append(New)
    return Out

def Matrix_Maker(input_list):
    Matches = List_Matrix(len(input_list))
    Out = []
    for entry in Matches:
        New = []
        for entry2 in entry:
            New.append(input_list[entry2])
        Out.append(New)
    return Out

def Matrix_Writer(input_matrix, input_index, output_file):
    List1 = glob.glob(input_index)
    Names = []
    for files in List1:
        ID = Folder_Name_Remover(files)
        Names.append(ID)
    Out = open(output_file, 'w')
    Line1 = '\t'.join(Names)
    Out.write('\t' + Line1 + '\n')
    for entry in range(len(Names)):
        Out.write(Names[entry] + '\t')
        Data_Line = []
        for entry2 in input_matrix[entry]:
            Data_Line.append(str(entry2))
        New = '\t'.join(Data_Line)
        Out.write(New + '\n')
    Out.close()
    
def cgMLST_Distance_List_Writer(input_index, Allele_Folder, Output_File):
    List1 = glob.glob(input_index)
    Matches = Folder_List_Maker(List1, Allele_Folder)
    with Pool() as pool:
        Data = pool.map(Allele_Difference_List, Matches)
    Matrix = Matrix_Maker(Data)
    Matrix_Writer(Matrix, input_index, Output_File)

def cgMLST_Distance_Matrix_Maker(input_index, output_file, Allele_Folder):
    File_List = glob.glob(input_index)
    File_List.sort()
    f = open(output_file, 'w')
    String1 = ''
    f.write('\t')
    for file_names in File_List:
        String1 = String1 + Folder_Name_Remover(file_names) + '\t'
    f.write(String1[0:-1] + '\n')
    for files in File_List:
        String1 = Folder_Name_Remover(files) + '\t'
        for match_files in File_List:
            String1 = String1 + str(Allele_Differences(files, match_files, Allele_Folder)) + '\t'
        f.write(String1[0:-1] + '\n')
    f.close()
    
def cgMLST_Distance_Matrix_Maker2(input_index, output_file):
    File_List = glob.glob(input_index)
    File_List.sort()
    f = open(output_file, 'w')
    String1 = ''
    f.write('\t')
    for file_names in File_List:
        String1 = String1 + Folder_Name_Remover(file_names) + '\t'
    f.write(String1[0:-1] + '\n')
    Matches = []
    Distances = []
    for files in File_List:
        String1 = Folder_Name_Remover(files) + '\t'
        File1 = Folder_Name_Remover(files)
        for match_files in File_List:
            File2 = Folder_Name_Remover(match_files)
            if ([File1, File2] in Matches):
                Pos = Matches.index([File1, File2])
                distance = Distances[Pos]
            elif ([File2, File1] in Matches):
                Pos = Matches.index([File2, File1])
                distance = Distances[Pos]
            else:
                distance = cgMLST_Distance(files, match_files)
                Distances.append(distance)
                Matches.append([File1, File2])
            String1 = String1 + str(distance) + '\t'
        f.write(String1[0:-1] + '\n')
    f.close()


def Replace_All(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

def Replace_All_File(text_file, dic):
    import fileinput
    for line in fileinput.input(text_file, inplace=True):
        line = Replace_All(line, dic)
        print(line, end="")

def Tree_Maker(input_distance_matrix, NJ_file, UPGMA_file):
    f = open(NJ_file, 'w')
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open(input_distance_matrix),
        delimiter="\t")
    nj_tree = pdm.nj_tree()
    f.write(nj_tree.as_string("newick"))
    f.close()
    f = open(UPGMA_file, 'w')
    upgma_tree = pdm.upgma_tree()
    f.write(upgma_tree.as_string("newick"))
    f.close()
    mydic = {"'":""}
    Replace_All_File(NJ_file, mydic)
    Replace_All_File(UPGMA_file, mydic)

def String_Converter(Input_String):
    Counter = 0
    Character = Input_String[0]
    Current_String = ''
    Out_List = []
    while Counter < len(Input_String):
        if Character == '\n':
            Out_List.append(Current_String)
            return Out_List
        elif Character == '\t':
            Out_List.append(Current_String)
            Current_String = ''
            Counter = Counter + 1
            Character = Input_String[Counter]
        else:
            Current_String = Current_String + Character
            Counter = Counter + 1
            Character = Input_String[Counter]

def Folder_Name_Remover(input_string):
    Temp_string =''
    for characters in input_string:
        if characters == '/':
            Temp_string = ''
        else:
            Temp_string = Temp_string + characters
    return Temp_string

def Single_Allele_Matches(input_fasta, allele_fasta):
    """Finds Match to single allele"""
    Allele_List = list(SeqIO.parse(allele_fasta, 'fasta'))
    Genome = list(SeqIO.parse(input_fasta, 'fasta'))
    for genes in Genome:
        for alleles in Allele_List:
            if str(alleles.seq) == str(genes.seq) or str(alleles.seq) == str(genes.seq.reverse_complement()):
                return alleles.id
                break
            else:
                continue
            
def Single_Gene_Match(input_gene, allele_fasta):
    """Finds Match to single gene"""
    Allele_List = list(SeqIO.parse(allele_fasta, 'fasta'))
    Allele = Folder_Name_Remover(allele_fasta)
    Allele_Match = Allele[0:-14] + '_Novel'
    for alleles in Allele_List:
        if str(alleles.seq) == str(input_gene.seq) or str(alleles.seq) == str(input_gene.seq.reverse_complement()):
            Allele_Match = alleles.id
            break
        else:
            continue
    return Allele_Match

def Single_Gene_Match_2(input_gene, allele_fasta):
    """Finds Match to single gene"""
    Allele_List = list(SeqIO.parse(allele_fasta, 'fasta'))
    Allele = Folder_Name_Remover(allele_fasta)
    Allele_Match = Allele[0:-14] + '_Novel'
    for alleles in Allele_List:
        if str(alleles.seq) == str(input_gene.seq) or str(alleles.seq) == str(input_gene.seq.reverse_complement()):
            Allele_Match = alleles.id
            break
        else:
            continue
    if Allele_Match[-5:] == 'Novel':
        Novel_Allele_Adder(input_gene, allele_fasta)
        Allele_List = list(SeqIO.parse(allele_fasta, 'fasta'))
        for alleles in Allele_List:
            if str(alleles.seq) == str(input_gene.seq) or str(alleles.seq) == str(input_gene.seq.reverse_complement()):
                Allele_Match = alleles.id
                break
            else:
                continue
    return Allele_Match

def Novel_Allele_Adder(novel_gene, allele_fasta):
    Gene_List = list(SeqIO.parse(allele_fasta, 'fasta'))
    output_handle = open(allele_fasta, 'w')
    for genes in Gene_List:
        SeqIO.write(genes, output_handle, 'fasta')
    Allele_Number = len(Gene_List) + 1
    New_Gene = novel_gene
    New_Gene.id = Gene_List[0].id[0:-1] + str(Allele_Number)
    SeqIO.write(New_Gene, output_handle, 'fasta')
    output_handle.close()

def Zero_Finder(allele_list_file):
    f = open(items, 'r')
    string1 = f.readline()
    while string1 != '':
	    if string1[-3:-1] == '_0':
		    print(items[77:] + ' ' + string1)
		    string1 = f.readline()
	    else:
		    string1 = f.readline()
    f.close()

def Allele_Name(input_name):
    output_name = ''
    for letters in input_name:
        if letters == '_':
            return output_name
        else:
            output_name = output_name + letters

##def Match_List_Output(allele_list, match_list, output_file):
##    out = open(output_file, 'w')
##    allele_list.sort()
##    for items in allele_list:
##        Item_Matches = []
##        for item_matches in match_list:
##            if Allele_Name(item_matches) == items:
##                Item_Matches.append(item_matches)
##            else:
##                continue
##        if len(Item_Matches) == 0:
##            Item_Matches.append(items + '_0')
##        out.write(List_Converter(Item_Matches) + '\n')
##    out.close()

##def Match_List_Output(allele_list, match_list, output_file):
##    out = open(output_file, 'w')
##    allele_list.sort()
##    for items in allele_list:
##        Item_Matches = []
##        for item_matches in match_list:
##            if ('_'.join(item_matches.split('_')[0:2])) == items:
##                Item_Matches.append(item_matches)
##            else:
##                continue
##        if len(Item_Matches) == 0:
##            Item_Matches.append(items + '_0')
##        Info = ','.join(Item_Matches)
##        out.write(Info + '\n')
##    out.close()

def Match_List_Output(allele_list, match_list, output_file):
    out = open(output_file, 'w')
    allele_list.sort()
    for items in allele_list:
        Item_Matches = []
        for matches in match_list:
            match = matches.split('_')
            match = '_'.join(match[0:-1])
            if (items == match):
                Item_Matches.append(matches)
            else:
                continue
        if len(Item_Matches) == 0:
            Item_Matches.append(items + '_0')
        Item_Matches = Repeat_Remover(Item_Matches)
        Info = ','.join(Item_Matches)
        out.write(Info + '\n')
    out.close()


def List_Converter(input_list):
    if len(input_list) == 1:
        output_string = input_list[0]
    else:
        output_string = input_list[0]
        for items in range(1, len(input_list)):
            new_item = input_list[items]
            output_string = output_string + ', ' + new_item
    return output_string  

def Repeat_Remover(any_list):
    """Removes repeats for any list"""
    new_list = []
    for items in any_list:
        if (items in new_list) == False:
            new_list.append(items)
    return new_list

def Match_List_Test(allele_list, match_list):
    allele_list.sort()
    for items in allele_list:
        Item_Matches = []
        for item_matches in match_list:
            if Allele_Name(item_matches) == items:
                Item_Matches.append(item_matches)
            else:
                continue
        if len(Item_Matches) == 0:
            Item_Matches.append(items + '_0')
        print(items + '\t' + List_Converter(Item_Matches))

def PSL_Allele_Match_List(input_PSL, input_fasta, allele_folder):
    Match_List = []
    f = open(input_PSL, 'r')
    String1 = f.readline()
    Isolate_Genome = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    while String1 != '':
        List1 = String_Converter(String1)
        PA_Gene_Name = List1[13][0:-2]
        Isolate_Gene_Name = List1[9]
        if int(List1[0]) > (int(List1[14]) / 2) and List1[8] == '+':
            Gene_Match = Isolate_Genome[Isolate_Gene_Name][int(List1[11]):int(List1[12])]
            Match_List.append(Single_Gene_Match_2(Gene_Match, allele_folder
                                                + PA_Gene_Name + '_alleles.fasta'))
            String1 = f.readline()
        elif int(List1[0]) > (int(List1[14]) / 2) and List1[8] == '-':
            Sequenced_Gene = Isolate_Genome[Isolate_Gene_Name][int(List1[11]):int(List1[12])]
            Rev_Gene = Sequenced_Gene.seq.reverse_complement()
            Gene_Match = SeqRecord(Seq(str(Rev_Gene)))
            Gene_Match.id = Sequenced_Gene.id
            Gene_Match.description = Sequenced_Gene.description
            Match_List.append(Single_Gene_Match_2(Gene_Match, allele_folder
                                                + PA_Gene_Name + '_alleles.fasta'))
            String1 = f.readline()
        else:
            String1 = f.readline()
    Match_List.sort()
    return Match_List

def PSL_Allele_Writer(intput_PSL, input_fasta, out_file, allele_folder):
    Allele_List = []
    List1 = glob.glob(allele_folder + '*.fasta')
    for items in List1:
        Allele_List.append(Folder_Name_Remover(items)[0:-14])
    Matches = PSL_Allele_Match_List(intput_PSL, input_fasta, allele_folder)
    Match_List_Output(Allele_List, Matches, out_file)

def cgMLST_Blat(input_list):
    subprocess.call('blat ' + input_list[1] + ' ' + input_list[0] +' -noHead cgMLST_PSL/' + input_list[0][0:-6] + '_cgMLST.psl', shell=True)

Fasta = sys.argv[1]
Folder = sys.argv[2]
FFN_Folder = sys.argv[3]
#if FFN_Folder[-1] != '/':
#    FFN_Folder = FFN_Folder + '/'
if Folder[-1] != '/':
    Folder = Folder + '/'

subprocess.call('cd ' + FFN_Folder, shell=True)
subprocess.call('mkdir cgMLST_PSL', shell=True)
subprocess.call('mkdir cgMLST_Allele_Matches', shell=True)
subprocess.call('mkdir cgMLST_Trees', shell=True)   

List1 = glob.glob('*.fasta')
##for items in List1:
##    print(items)
Lines = []
for items in List1:
    Lines.append([items, Fasta])

with Pool() as pool:
    pool.map(cgMLST_Blat, Lines)

for files in List1:
    PSL_Allele_Writer('cgMLST_PSL/' + files[0:-6] + '_cgMLST.psl', files, 'cgMLST_Allele_Matches/' + files[0:-6] + '_Allele_Matches.txt', Folder)

cgMLST_Distance_List_Writer('cgMLST_Allele_Matches/*.txt', Folder, 'cgMLST_Trees/Distance_Matrix.txt')
Tree_Maker('cgMLST_Trees/Distance_Matrix.txt', 'cgMLST_Trees/Distance_NJ.tre', 'cgMLST_Trees/Distance_UPGMA.tre')
