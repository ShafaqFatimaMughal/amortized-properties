from FMIndex.fmindex import *
from SuffixTree.suffixtree import *
from matplotlib import pyplot as plt
import numpy as np
import time, math, sys

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
                
def get_FM_build_time(filename):
    '''
    Returns the time taken to build the FM-Index for a given file with a certain number of letters of DNA
    '''
    start_time = time.time()
    I = Implementation(filename) # creating an implementation object that builds an fm-index for a given file
    end_time = time.time()
    return  (end_time-start_time)

def get_suffix_build_time(filename):
    '''
    Returns the time taken to build the Suffix Tree for a given file with a certain number of letters of DNA
    '''
    start_time = time.time()
    with open(filename, 'r') as file:
        data = file.read().replace('\n', '')
        tree = suffix_tree(data) # building suffix tree for the given data set
    end_time = time.time()
    return  (end_time-start_time)    

# def search(txt, query):
#     res = [i for i in range(len(txt)) if txt.startswith(query, i)]
#     return res
def linear_substring_search(string, substring):
    '''
    Finds each occurence of the substring in the string and then returns the indexes for each
    '''
    indexes = []
    for i in range(len(string)):
        if string[i] == substring[0]: # if the current element is the same as the first element of the substring 
            if string[i:i+len(substring)] == substring: # then check the rest of the string. if it matches add the index of first letter to tthe list
                indexes.append(i)
    return indexes

def linear_find(filename, substring, time_bool=False):
    '''
    Reads a file and does linearSubStringSearch on it
    Either time taken for this substring search is returned or the list of indexes
    This depends on third argument
    '''
    with open(filename, 'r') as file:
        data = file.read().replace('\n', '')
        start_time = time.time()
        result =  linear_substring_search(data, substring)
        print(result)
        end_time = time.time()
        if time_bool == False:
            return result
        else: # if timebool argument is given then even return the time taken to search 
            return (end_time - start_time)

def suffix_tree_find(filename, substring, time_bool=False):
    ''' 
    Reads a file and builds a suffix tree for it. Then the substring search is done on the file
    Either time taken for this substring search is returned or the list of indexes
    This depends on third argument
    '''
    with open(filename, 'r') as file:
        data = file.read().replace('\n', '')
        tree = suffix_tree(data) # building the suffix tree
        start_time = time.time()
        result = tree.find_all(substring) # searching for first index of all occurrences of substring
        print(result)
        end_time = time.time()
        if time_bool == False:
            return result
        else: # if timebool argument is given then even return the time taken to search
            return (end_time - start_time)

def fm_index_find(filename, substring, time_bool=False):
    ''' 
    Reads a file and builds a FM-Index for it. Then the substring search is done on the file
    Either time taken for this substring search is returned or the list of indexes
    This depends on third argument
    '''
    I = Implementation(filename) # creating an implementation object that builds an fm-index for a given file
    start_time = time.time()
    result = I.search(substring) # searching for first index of all occurrences of substring 
    print(result)
    end_time = time.time()
    if time_bool == False:
        return result
    else: # if timebool argument is given then even return the time taken to search
        return (end_time - start_time)

def plot(expected, title, n, color):
    ''' 
    Plots a histogram for a list of expected values
    Resulting graph will give the most expected value at mean
    '''
    bins = 20 
    binWidth = (max(expected) - min(expected)) / bins 
    plt.hist(expected, bins=bins , weights=np.ones(len(expected))/(len(expected)*binWidth), color=color) # produce histogram 
    plt.xlabel("Expected Time when n="+str(n))
    plt.ylabel("Frequency")
    plt.title(title)
    plt.show() 

def expected_linear_time(numExp, sample, filename, substring, title, color):
    ''' 
    Reads a file and does linear substring search on it multiple times
    Time taken for this substring search is stored each time to get a sample and then the mean/expected value is calculated
    This is repeated several times to get several expected values
    The Expected values are plotted using a histogram using helper function
    '''
    expected = np.zeros(numExp) # numpy array of all expected values
    with open(filename, 'r') as file:
        data = file.read().replace('\n', '')
        for i in range(numExp): # To find several expectations. numExp is the number of experiments
            timeTaken = [] 
            for j in range(sample): # To find expectation using several simulations of samples
                start_time = time.time()
                print( linear_substring_search(data, substring) )
                end_time = time.time()
                timeTaken.append(end_time-start_time)
                # print(i,j)
            expected[i] = sum(timeTaken)/sample
        plot(expected, title, len(data), color)

def expected_suffix_time(numExp, sample, filename, substring, title, color):
    ''' 
    Reads a file, builds its suffix tree and does a substring search on it multiple times
    Time taken for this substring search is stored each time to get a sample and then the mean/expected value is calculated
    This is repeated several times to get several expected values
    The Expected values are plotted using a histogram using helper function
    '''
    expected = np.zeros(numExp)  # numpy array of all expected values
    with open(filename, 'r') as file:
        data = file.read().replace('\n', '')
        tree = suffix_tree(data)
        for i in range(numExp): # To find several expectations. numExp is the number of experiments
            timeTaken = []
            for j in range(sample): # To find expectation using several simulations of samples
                start_time = time.time()
                print(tree.find_all(substring))
                end_time = time.time()
                timeTaken.append(end_time-start_time)
                # print(i,j)
            expected[i] = sum(timeTaken)/sample
        plot(expected, title, len(data), color)

def expected_fm_time(numExp, sample, filename, substring, title, color):
    ''' 
    Reads a file, builds its FM-Index and does a substring search on it multiple times
    Time taken for this substring search is stored each time to get a sample and then the mean/expected value is calculated
    This is repeated several times to get several expected values
    The Expected values are plotted using a histogram using helper function
    '''
    expected = np.zeros(numExp)  # numpy array of all expected values
    I = Implementation(filename) # creating an implementation object that builds an fm-index for a given file
    for i in range(numExp): # To find several expectations. numExp is the number of experiments
        timeTaken = []
        for j in range(sample): # To find expectation using several simulations of samples
            start_time = time.time()
            print( I.search(substring) )
            end_time = time.time()
            timeTaken.append(end_time-start_time)
            # print(i,j)
        expected[i] = sum(timeTaken)/sample
    plot(expected, title, I.n, color)

def line_plot(x_axis, suffix_y, fm_y, line_y=False):
    '''
    Plots a line graph for the expected time taken for substring search using suffix and fm-index
    Linear search plot is optional
    '''
    if type(line_y) == type(fm_y):
        plt.plot(x_axis, line_y, label='Linear Search', color='orange')
    plt.plot(x_axis, suffix_x, label='Suffix Tree', color='purple')
    plt.plot(x_axis, fm_y, label='FM-Index', color=None)
    plt.legend()
    plt.xlabel('log( Number of DNA letters )')
    plt.ylabel('log( Time (s) )')
    plt.title('Expected Time w.r.t file size (number of DNA letters)')
    plt.show()

def build_time_analysis(x_axis, lst, n):
    '''
    Plots a line graph for the expected time taken for fm index and sufflix tree to be build based on size of data set
    '''
    expected_fm_build = []
    expected_suffix_build = []
    for name in lst: # for each file
        fm_build = []
        suffix_build = []
        for i in range(n): # getting expected build time for fm index and suffix tree, for each file
            fm_build.append(get_FM_build_time(name))
            suffix_build.append(get_suffix_build_time(name))
        # print(fm_build, suffix_build)
        expected_fm_build.append(sum(fm_build)/n)
        expected_suffix_build.append(sum(suffix_build)/n)
    plt.plot(x_axis, np.array(expected_suffix_build), label='Suffix Tree', color='purple')
    plt.plot(x_axis, np.array(expected_fm_build), label='FM-Index', color=None)
    plt.legend()
    plt.xlabel('Number of DNA letters (in millions)')
    plt.ylabel('Time (s)')
    plt.title('Expected build Time w.r.t file size (number of DNA letters)')
    plt.show()


######################################## Testing for 1 search #########################################################################################################
# print(get_FM_build_time('DataSets/100000.txt'))

# linear_find('DataSets/10000.txt', 'AAT')
# suffix_tree_find('DataSets/10000.txt', 'AAT')
# fm_index_find('DataSets/10000.txt', 'AAT')

# print('Linear Search Time', linear_find('DataSets/1000.txt', 'AAT', True))
# print('Suffix Tree Search Time', suffix_tree_find('DataSets/1000.txt', 'AAT', True))
# print('FM Index Search Time', fm_index_find('DataSets/100000.txt', 'AAGTCT', True))

################################### Testing for expected time ####################################################################################################
######################################## LINEAR SEARCH ##########################################
# expected_linear_time(250, 50, 'DataSets/1000.txt', 'GGAATT', 'Linear Search', 'orange')
# expected_linear_time(250, 50, 'DataSets/10000.txt', 'CTCGTGA', 'Linear Search', 'orange')
# expected_linear_time(250, 50, 'DataSets/100000.txt', 'TATGCAC', 'Linear Search', 'orange')
# expected_linear_time(250, 20, 'DataSets/1000000.txt', 'AGTACAGC', 'Linear Search', 'orange')
# expected_linear_time(250, 20, 'DataSets/2500000.txt', 'CACATTT', 'Linear Search', 'orange')

####################################### SUFFIX-TREE ############################################# 
# expected_suffix_time(250, 50, 'DataSets/1000.txt', 'GGAATT', 'Suffix-Tree', 'purple')
# expected_suffix_time(250, 50, 'DataSets/10000.txt', 'CTCGTGA', 'Suffix-Tree', 'purple')
# expected_suffix_time(250, 50, 'DataSets/100000.txt', 'TATGCAC', 'Suffix-Tree', 'purple')
# expected_suffix_time(250, 25, 'DataSets/1000000.txt', 'AGTACAGC', 'Suffix-Tree', 'purple')
# expected_suffix_time(250, 25, 'DataSets/2500000.txt', 'CACATTT', 'Suffix Tree', 'purple')
# expected_suffix_time(250, 50, 'DataSets/5000000.txt', 'GGACTACT', 'Suffix Tree', 'purple')

######################################## FM-INDEX ###############################################
# expected_fm_time(250, 50, 'DataSets/1000.txt', 'GGAATT', 'FM/-Index', None)
# expected_fm_time(250, 50, 'DataSets/10000.txt', 'CTCGTGA', 'FM-Index', None)
# expected_fm_time(250, 50, 'DataSets/100000.txt', 'TATGCAC', 'FM-Index', None)
# expected_fm_time(250, 25, 'DataSets/1000000.txt', 'AGTACAGC', 'FM-Index', None)
# expected_fm_time(250, 25, 'DataSets/2500000.txt', 'CACATTT', 'FM-Index', None)
# expected_fm_time(250, 50, 'DataSets/5000000.txt', 'GGACTACT', 'FM-Index', None)

############################################ Line Plot Analysis of time to Search ######################################################################################################
x_axis = np.array([math.log(1000, 10), math.log(10000, 10), math.log(100000, 10), math.log(1000000, 10), math.log(2500000, 10)])
line_y = np.array([math.log(0.0013, 10), math.log(0.0056, 10), math.log(0.040, 10), math.log(0.225, 10), math.log(0.6, 10)])
suffix_x = np.array([math.log(0.00072, 10), math.log(0.0008, 10), math.log(0.00116, 10), math.log(0.00120, 10), math.log(0.0077, 10)])
fm_y = np.array([math.log(0.0006, 10), math.log(0.00065, 10), math.log(0.0009, 10), math.log(0.0015, 10), math.log(0.0050, 10)])
line_plot(x_axis, suffix_x, fm_y, line_y)    
# line_plot(x_axis, suffix_x, fm_y)

############################################ Line Plot Analysis of time to build #######################################################################################
# x = np.array([1000, 10000, 100000, 1000000, 2500000])
# file_lst = ['DataSets/1000.txt', 'DataSets/10000.txt', 'DataSets/100000.txt', 'DataSets/1000000.txt', 'DataSets/2500000.txt']
# build_time_analysis(x, file_lst, 5)

################################ Testing Gene Correction And Analysis #####################################################################################
# I = Implementation('DataSets/1000000.txt')

# print(I.correction('TGCT', 'genomic data'))
# print(I.correction('TGCT', 'BIOINFORMATICS', 50))
# print(I.correction('TGCT', 'ATGC', 50))

# print(I.gene_analysis('CCCC', 'Cmkwx'))
# print(I.gene_analysis('CCCC', 'CTMK'))
# print(I.gene_analysis('CCCC', 'AGTCA'))
# print(I.gene_analysis('CCCC', 'AGTCA', 150), '150')