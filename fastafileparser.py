from os.path import exists as file_exists
import fileinput

class FASTAFileParser:
    def __init__(self, filePath):
        if file_exists(filePath) != True:
            FileNotFoundError

        self.file = open(filePath) 
        self.sequence = ""

    def get_sequence(self):
        sequenceFound = False

        for line in self.file:
            line = line.strip()   

            if line.startswith('>') == False and sequenceFound == False:
                #print("No Sequence Found!")
                return

            if line.startswith('>') == True and sequenceFound == True:
                #print("Found another Sequence in file ... Stopping!")
                return self.sequence

            if line.startswith('>') == True and sequenceFound == False:
                #print("Found a Sequence in file!")
                sequenceFound = True
            else:
                self.sequence += line

        return self.sequence

