#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
A somewhat memory-efficient blastword generator.
Optimized for simplicity.

Uses neither https://xkcd.com/208/ nor Perl. Nor Java.
"""
import itertools
#What a pity we can't use numpy in Jython...

__author__  = "Uli Koehler"
__license__ = "Apache License v2.0"
__version__ = "3.14"

def newMatNxN(n):
    """Generate a new nxn matrix. Jython workaround for lack of Numpy"""
    return [[None] * n for i in range(n)]

def readMatblasAlignmentMatrix(filename):
    """
    Read a substitution matrix in matblas format.

    Keyword arguments:
        filename: The filename to read the matrix from

    Returns a tuple (column/row list, substitution matrix)
    """
    infile = open(filename)
    currentRow = 0
    for line in infile:
        if line.startswith("#"): continue
        elif line.startswith(" "): #Column indicator
            columns = line.split()
            matrix = newMatNxN(len(columns))
        else: #Matrix row
            parts = line.split()
            assert(len(parts) == len(columns) + 1)
            #Assume rows are in the same order as column
            assert(columns[currentRow] == parts[0])
            for x in range(len(columns)):
                matrix[x][currentRow] = int(parts[x + 1])
            currentRow += 1
    #Postcondition: Matrix filled completely (numpy would be *quite* nice)
    for x in range(len(columns)):
        for y in range(len(columns)):
            assert matrix[x][y] is not None
            #Assert symmetry
            assert matrix[x][y] == matrix[y][x]
    infile.close()
    return (columns, matrix)

def wordDistance(a, b, mat, alphabet):
    assert len(a) == len(b)
    totalScore = 0
    for i in range(len(a)):
        aIdx = alphabet.index(a[i])
        bIdx = alphabet.index(b[i])
        substScore = mat[aIdx][bIdx]
        totalScore += substScore
    return totalScore

def genWords(alphabet, length):
    """
    Iterate over cartesian product of n times the alphabet
    
    >>> sorted(list(genWords("abc", 2)))
    ['aa', 'ab', 'ac', 'ba', 'bb', 'bc', 'ca', 'cb', 'cc']
    """
    for word in itertools.product(alphabet, repeat=length):
        if "*" in word: continue #Skip wildcard character
        yield "".join(word)

def getSubsequences(sequence, n):
    """
    Generate all n-long subsequence of the given sequence
    
    >>> list(getSubsequences("ab", 3))
    []
    >>> list(getSubsequences('abc', 3))
    ['abc']
    >>> list(getSubsequences('abcdef', 3))
    ['abc', 'bcd', 'cde', 'def']
    >>> list(getSubsequences('abcdefabc', 3))
    ['abc', 'bcd', 'cde', 'def', 'efa', 'fab', 'abc']
    """
    #Zip with itself n times
    return ("".join(ss) for ss in itertools.izip(*[sequence[i:] for i in range(n)]))

def getBLASTWords(sequence, length, alphabet, mat, threshold):
    """
    Given a sequence, a word length, an alphabet and a substitution matrix, calculates all BLAST words below the given threshold
    """
    blastWords = set()
    for word in genWords(alphabet, length):
        for subsequence in getSubsequences(sequence, length):
            distance = wordDistance(word, subsequence, mat, alphabet)
            # Similarity metric = operator >=
            # --> Distance metric = operator <
            # --> But we have a distance matrix with negative values!
            if distance >= threshold:
                if not word in blastWords:
                    blastWords.add(word)
                    yield word

def writeBLASTWords(outfilename, sequence, length, alphabet, mat, threshold):
    outfile =  open(outfilename, "w");
    for word in getBLASTWords(sequence, length, alphabet, mat, threshold):
        print >>outfile, word
    outfile.close()

def readFASTASimple(filename):
    """A simple FASTA reader that concats all sequences and ignores other lines"""
    f = open(filename)
    sequence = ""
    for line in f:
        #We want to be ultra-compatible
        if line.startswith(">") or line.startswith(";"):
            continue
        sequence += line.strip()
    f.close()
    return sequence

def main():
    import sys
    #The QnD way
    #Our Jython wrapper makes arguments begin at index 0 instead of (cython) 1
    wordLength = int(sys.argv[1])
    threshold = float(sys.argv[2])
    sequenceFile = sys.argv[3]
    scoringMatrixFile = sys.argv[4]
    outputFile = sys.argv[5]
    sequence = readFASTASimple(sequenceFile)
    alphabet, mat = readMatblasAlignmentMatrix(scoringMatrixFile)
    #Do it now.
    writeBLASTWords(outputFile, sequence, wordLength, alphabet, mat, threshold)

#Run main both in Jython and in standalone mode
main()

def run_tests():
    #Unit tests. Don't work in Jython
    import doctest
    doctest.testmod()
    #Some random testing.
    alphabet, mat = readMatblasAlignmentMatrix("blosum62.txt")
    sequence = "ARNDCQEGHILKMFPSTWYVBZX"
    for blastword in getBLASTWords(sequence, 4, alphabet, mat, 20.0):
        print blastword