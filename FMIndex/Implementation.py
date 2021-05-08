from FMIndex.fmindex import *

class Implementation:
    def __init__(self, location): 
        self.fileloc = location
        with open(location) as file:
            data = file.read().replace('\n', '')
            self.fmindex = create_fm_index(data)
            self.n = len(data)
    def search(self, query):
        ''' 
        returns list of all occurences of query/substring/gene
        '''
        return self.fmindex.search(query)

    def correction(self, string1, string2, approxpos=None):
        ''' 
        If there are no syntax and formatting errors then it searches for string1 in the data
        If found, it replaces string1 with string2 and overwrites the file
        If there are multiple occurrences and approximate position is given by user
        then it replaces string that is closest to position. 
        '''
        if string2.isupper() == False:
            return 'Your replacement DNA is not in uppercase!'
        for letter in string2:
            if letter not in 'ATGC':
                return 'Your replacement seems to be invalid! You can only replace the gene with another gene (so input DNA letters)'
        Rs1 = self.fmindex.search(string1)
        originalf = open(self.fileloc, "r")
        lines = originalf.read()
        lines = str(lines)
        if approxpos==None:
            if len(Rs1)==1:
                lines = lines[:Rs1[-1]] + str(string2) + lines[Rs1[-1]+len(string2):]
                originalf = open(self.fileloc, "w")
                originalf.writelines(lines)
                originalf.close()
                return string1 + " Now changed to " + string2
            else:
                return "Multiple Genes found, please specify what gene to change or pass in additional approximate positional argument"
        else:
            if Rs1!=[]:
                d = math.inf
                pos = -1
                for i in Rs1:
                    if abs(approxpos - i) < d:
                        d = abs(approxpos - i)
                        pos = i
                lines = lines[:pos] + str(string2) + lines[pos+len(string2):]
                originalf = open(self.fileloc, "w")
                originalf.writelines(lines)
                originalf.close()
                return string1 + " Now changed to " + string2 + " at position " + str(pos)
            else:
                return "No matching found"

    def gene_analysis(self, string1, string2, approxlength=None):
        '''
        If there are no syntax and formatting errors then it searches for the gene between string1 and string2
        '''
        if string2.isupper() == False or string1.isupper() == False:
            return 'Your DNA letters are not in uppercase!'
        for letter in string2:
            if letter not in 'ATGC':
                return 'Your input seems to be invalid! Please make sure that you have input DNA letters'
        Rs1 = self.fmindex.search(string1)
        Rs2 = self.fmindex.search(string2)
        if approxlength!=None:
            d = math.inf
            lst = []
            for i in range(len(Rs1)):
                for j in range(len(Rs2)):
                    if Rs1[i] < Rs2[j] and Rs2[j]-Rs1[i]-len(string1)-approxlength < d and Rs2[j]-Rs1[i]-len(string1)-approxlength > 0:
                        d = Rs2[j]-len(string1) -approxlength - Rs1[i]
            for i in range(len(Rs1)):
                for j in range(len(Rs2)):
                    if Rs1[i] < Rs2[j] and Rs2[j]-Rs1[i]-len(string1)-approxlength == d:
                        lst.append((Rs1[i]+len(string1), Rs2[j]))
            if len(lst)>1:
                return "Multiple instances found please specify gene more correctly"
            else:
                if lst==[]:
                    return "No instance found"
                else:
                    file = open(self.fileloc, "r")
                    s = file.read()
                    file.close()
                    return s[ lst[-1][0] : lst[-1][1] ]
        else:
            if len(Rs1) > 1 or len(Rs2) > 1:
                return "Multiple instances found please specify gene more correctly"
            elif Rs1==[] or Rs2==[]:
                return "No instance found"
            else:
                file = open(self.fileloc, "r")
                s = file.read()
                file.close()
                return s[ Rs1[-1]+len(string1) : Rs2[-1] ]
