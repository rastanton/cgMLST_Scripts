import glob
import sys

##Usage: python Matrix_Converter_Exe.py input_Matrix_file output_list_file

def In_List(item, list1):
    """determines if item is in list1"""
    for x in list1:
        if x == item:
            return True
    else: return False

def String_Converter(Input_String):
    Counter = 0
    Character = Input_String[0]
    Current_String = ''
    Out_List = []
    while Counter < len(Input_String):
        Character = Input_String[Counter]
        if Character == '\t':
            Out_List.append(Current_String)
            Current_String = ''
            Counter = Counter + 1
        else:
            Current_String = Current_String + Character
            Counter = Counter + 1
    Out_List.append(Current_String)
    if Out_List[-1][-1] == '\n':
        Out_List[-1] = Out_List[-1][0:-1]
    return Out_List

def Underscore_Converter(Input_String):
    Counter = 0
    Character = Input_String[0]
    Current_String = ''
    Out_List = []
    while Counter < len(Input_String):
        Character = Input_String[Counter]
        if Character == '_':
            Out_List.append(Current_String)
            Current_String = ''
            Counter = Counter + 1
        else:
            Current_String = Current_String + Character
            Counter = Counter + 1
    Out_List.append(Current_String)
    if Out_List[-1][-1] == '\n':
        Out_List[-1] = Out_List[-1][0:-1]
    return Out_List

def Matrix_Converter(input_matrix, output_file):
    """Converts a distance matrix to a list of distances"""
    f = open(input_matrix, 'r')
    g = open(output_file, 'w')
    String1 = f.readline()
    Match_List = String_Converter(String1)
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        for values in range(1, len(Match_List)):
            g.write(List1[0] + '\t' + Match_List[values] + '\t' + List1[values] + '\n')
##            print(List1[0] + '\t' + Match_List[values] + '\t' + List1[values])
        String1 = f.readline()
    f.close()
    g.close()

def Matrix_Converter_Non_Self(input_matrix, output_file):
    """Converts a distance matrix to a list of distances"""
    f = open(input_matrix, 'r')
    g = open(output_file, 'w')
    String1 = f.readline()
    Match_List = String_Converter(String1)
    String1 = f.readline()
    while String1 != '':
        List1 = String_Converter(String1)
        for values in range(1, len(Match_List)):
            if List1[0] != Match_List[values]:
                g.write(List1[0] + '\t' + Match_List[values] + '\t' + List1[values] + '\n')
##            print(List1[0] + '\t' + Match_List[values] + '\t' + List1[values])
        String1 = f.readline()
    f.close()
    g.close()

def Matrix_Converter_Short(input_matrix, output_file):
    """Converts a distance matrix to a list of distances"""
    f = open(input_matrix, 'r')
    g = open(output_file, 'w')
    String1 = f.readline()
    Match_List = String_Converter(String1)
    String1 = f.readline()
    Counter = 2
    while String1 != '':
        List1 = String_Converter(String1)
        for values in range(Counter, len(Match_List)):
            g.write(List1[0] + '\t' + Match_List[values] + '\t' + List1[values] + '\n')
##            print(List1[0] + '\t' + Match_List[values] + '\t' + List1[values])
        Counter = Counter + 1
        String1 = f.readline()
    f.close()
    g.close()

Matrix_Converter_Short(sys.argv[1], sys.argv[2])
