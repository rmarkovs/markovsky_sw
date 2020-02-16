import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s1', dest ='s1', help = 'input the first sequence to be used')
    parser.add_argument('-s2', dest = 's2', help = 'input the second sequence to be used')
    parser.add_argument('-m', dest='m', help='input the matrix to be used')
    parser.add_argument('-g', dest = 'g', help='input the gap penalty to be used')
    parser.add_argument('-o', dest = 'o',help='input the output file')
    args = parser.parse_args()

    def matrixRead():
        matrix = open(args.m, 'r')
        matrix = matrix.readlines()
        return matrix

    matrixRead()

    def matrixValues():
        matrixData = matrixRead()
        newSeries = pd.Series([])
        columnAndRowNames = []
        lineCount =0
        columnCount = 0
        data = []
        df1 = pd.DataFrame(data, columns = [])
        frames = [df1]
        for line in matrixData:
            if lineCount == 15:
                i = 0
                temp = ""
                temp = line
                temp = temp.rstrip(' ')
                for a in temp:
                    if (a != ' ') & (a != '*') & (a != '\\') & (a != '\\n') :
                        columnAndRowNames.append(a)
            if lineCount == 16:
                lineArray = []
                i = 1
                temp = ""
                temp = line
                tempPositionString = ''
                while i<len(temp):
                    if (temp[i:i+1] != ' ') & (temp[i:i+1] != '*') & (temp[i:i+1] != '\\') & (temp[i:i+1] != '\n'):
                        if temp[i:i+1] == '-':
                            tempPositionString = temp[i:i+2]
                            lineArray.append(tempPositionString)
                            i = i +1
                        else:
                            lineArray.append(temp[i:i+1])
                    i = i + 1
                tempString = ""
                tempString = columnAndRowNames[columnCount]
                df3 = pd.DataFrame({tempString: lineArray})
                frames.append(df3)
                columnCount +=1
            if lineCount > 16:
                lineArray = []
                i = 1
                temp = ""
                temp = line
                tempPositionString = ""
                while i<len(temp):
                     if (temp[i:i + 1] != ' ') & (temp[i:i + 1] != '*') & (temp[i:i + 1] != '\\') & (temp[i:i + 1] != '\n'):
                        if temp[i:i + 1] == '-':
                            tempPositionString = temp[i:i + 2]
                            lineArray.append(tempPositionString)
                            i = i + 1
                        else:
                            lineArray.append(temp[i:i + 1])
                     i = i + 1
                tempString = ""
                tempString = columnAndRowNames[columnCount]
                frames.append(pd.DataFrame({tempString: lineArray}))
                columnCount +=1
            lineCount +=1
        result = pd.concat(frames, axis = 1)
        columnAndRowNames.append('ZERO_ROW')
        result.index = columnAndRowNames
        return result
    def matrixToScore():
        zeroData = []
        #makes a zero row
        differentMatrices = []
        #set of zero row matrices
        seq1 = open(args.s1,'r')
        seq1 = seq1.readlines()
        #column
        columnElementsArray = []
        rowElementsArray = []
        colNum = 0
        rowNum = 0
        for line in seq1:
            temp = ""
            temp = line
            if(temp[0:1] != ">"):
                for a in temp:
                    if(a != ' ') & (a != '\n'):
                        columnElementsArray.append(a)
                        colNum +=1
        #creates a list of what the columns are
        seq2 = open(args.s2,'r')
        seq2 = seq2.readlines()
        for line in seq2:
            temp2 = ""
            temp2 = line
            if (temp2[0:1] != ">"):
                for a in temp2:
                    if (a != ' ') & (a != '\n'):
                        rowElementsArray.append(a)
                        rowNum +=1
        #creates a list of what the rows are
        columnSize = colNum
        #establishes number of columns
        rowSize = rowNum
        #establishes number of columns
        i = 0
        while (i < colNum):
            zeroData.append(0)
            i += 1
        #appens a number of zero columns that is equal to the number of columns
        n = 0
        while (n < rowNum):
            temp = ""
            temp = rowElementsArray[n]
            differentMatrices.append(pd.DataFrame({temp: zeroData}))
            n += 1
        #names the rows and creates an array of the elements of the rows
        dataMatrix = pd.concat(differentMatrices, axis = 1)
        dataMatrix.index = columnElementsArray


        sourceMatrix = matrixValues()
        #calls the source matrix
        alignValues = []
        #creates an array for the values for the diagnoal path
        values = {}
        #creates the dictionary of pairs and their alignment values
        if(len(rowElementsArray) > len(columnElementsArray)):
            #checks if there are more rows than columns
            i = 0
            #sets counter
            #was columnElements
            while i<=len(columnElementsArray):


                #sets the number of the row
                column = columnElementsArray[len(columnElementsArray) -1 - i]
                #sets the column element
                columnCount = len(columnElementsArray) -1 - i
                rowCount = columnCount
                row = rowElementsArray[rowCount]

                #sets the column number
                neededValue = sourceMatrix.at[row, column]
                #gets the alignment value for the row and column
                alignValues.append(neededValue)
                #alignValues keeps track of the diagonal alignments
                rowCol = ""
                rowCol += row + column
                #creates a string combination of the row and the column
                if rowCol not in values.keys():
                    values[rowCol] == neededValue
                #if such a combination doesnt exist, puts it into the dictionary with the alighment value
                i +=1
        else:
            i = 0
            while i <= len(rowElementsArray):
                #if there are not more rows than columns
                row = rowElementsArray[len(rowElementsArray) - 1 - i]
                #the row is the element of the array at that spot
                rowCount = len(rowElementsArray) -1 -i
                #gets the number of the row


                columnCount = rowCount
                column = columnElementsArray[columnCount]
                #gets the element of the column



                #gets the count of the column
                neededValue = sourceMatrix.at[row,column]
                #gets the alignment value for those amino acids
                alignValues.append(neededValue)
                #puts the alignment value into the list
                rowCol = ""
                rowCol += row + column
                #creates a string of the row and column
                if rowCol not in values.keys():
                    values[rowCol] = neededValue
                    #appends the pair if it isnt in the dictionary
                    #print(rowCol)
                i += 1
        rowAlignments = []
        rowAlignString = ""
        for i in columnElementsArray:
            column = i
            for x in rowElementsArray:
                row = x
                alignVal = sourceMatrix.at[row,column]
                rowAlignString = row + column
                if rowAlignString not in values:
                    values[rowAlignString] = alignVal
        #print(values)
        columnAlignString = ""
        for i in rowElementsArray:
            row = i
            for x in columnElementsArray:
                column = x
                alignVal = sourceMatrix.at[row, column]
                columnAlignString = row + column
                if columnAlignString not in values:
                    values[columnAlignString] = alignVal
        #print(values)
        bestAlignmentValue = 0
        for i in values.keys():
            if int(values[i]) <0:
                values[i] = 0
        #print(values)
        bestAlignmentTrace = []
        currentRow = ""
        currentColumn = ""
        maxSize = 0
        gapPen = int(args.g)
        if(len(columnElementsArray)<len(rowElementsArray)):
            currentRow = rowElementsArray[len(columnElementsArray) -1]
            currentColumn = columnElementsArray[len(columnElementsArray) -1]
            bestAlignmentValue = -1 * (gapPen * (len(rowElementsArray) - len(columnElementsArray)))
            maxSize = len(columnElementsArray)
        else:
            currentRow = rowElementsArray[len(rowElementsArray) -1]
            currentColumn = columnElementsArray[len(rowElementsArray) -1]
            bestAlignmentValue = -1 * (gapPen * (len(columnElementsArray) - len(rowElementsArray)))
            maxSize = len(rowElementsArray)
        currentPair = ""
        currentPair = currentRow + currentColumn
        bestAlignmentValue += int(values[currentPair])
        pairLeft = ""
        pairUp = ""
        pairDi = ""
        rowUp = ""
        columnUp = ""
        z = 1
        while z <= maxSize - 1:
            rowUp = rowElementsArray[maxSize - 1 - z]
            columnUp = columnElementsArray[maxSize  - 1 -z]
            pairLeft = currentRow + columnUp
            pairUp = rowUp + currentColumn
            pairDi = rowUp + columnUp
            leftVal = int(values[pairLeft])
            upVal = int(values[pairUp])
            diVal = int(values[pairDi])
            if(leftVal >= upVal) & (leftVal >= diVal):
                currentPair = pairLeft
                currentColumn = columnUp
                bestAlignmentValue += leftVal
                bestAlignmentTrace.append('left')
            if(upVal>=leftVal) & (upVal >= diVal):
                currentPair = pairUp
                currentRow = rowUp
                bestAlignmentValue += upVal
                bestAlignmentTrace.append('up')
            if(diVal>=upVal) & (diVal >= leftVal):
                currentPair = pairDi
                currentRow = rowUp
                currentColumn = columnUp
                bestAlignmentValue += diVal
                bestAlignmentTrace.append('diagonal')
            z +=1
            # print(bestAlignmentValue)
            # print(bestAlignmentTrace)



        outFile = open(args.o, 'w')
        outFile.write(seq1)
        outFile.write("\n")
        outFile.write(seq2)
        outFile.write("\n")
        outFile.write(bestAlignmentValue)
        outFile.write("\n")
        for i in bestAlignmentTrace:
            outFile.write(i)
        # print(seq1)
        # print(seq2)
        # print(bestAlignmentValue)
        # print(bestAlignmentTrace)
        #print(values)

        #return dataMatrix
    matrixToScore()


    matrixValues()



main()
