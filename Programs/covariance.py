import xlrd
import xlwt
import csv
from itertools import repeat
import scipy
import numpy
from scipy import stats
import pickle

#function to read the data from an excel spread sheet. 
#returns a list of all of the comparisons in a probe with the P value and fold change for every probe. I adjusted this so it just pickled the data (python module to save an object to a file) and so this function is no longer being called. 
def getData(myData):
	#I added a few print lines so I know which functions are currently running
	print "get data"

        #initialize empty 2D matrix. In python, a 2D matrix is just a list of lists.
	#Here I have a list of 15131 lists. Each sublist corresponds with one probe.
	cells=[[] for i in repeat(None, 15131)]
	assert len(cells)==15131
	probeNum=0

	
	#read in data from excel file and get column 0
	data=xlrd.open_workbook(myData)
	sheet=data.sheets()[0]
	c=sheet.col_values(0)

	#get data for each probe
	for probe in range(0,len(c)):
		#all probes end with 'at' so I used this to distinguish cells with probes
		if str(c[probe]).endswith("at"):

                        #add the corresponding P and C values to the location in the matrix corresponding with the probe

			#if sheet.col(3)[probe].value==0:
				#cells[probeNum].append([float('nan'),0])
			#else:
			cells[probeNum].append([sheet.col(1)[probe].value,sheet.col(2)[probe].value])
			
			#I use this to put info about each probe in the correct spot in the matrix
			probeNum=probeNum+1
			if probeNum==15131:
				probeNum=0
	
    #pickle the python list			
	pickle.dump(cells,open('/home/smparker/Desktop/HeadRawNoCutoffs.p','wb'))
	return cells	

#create a covariance matrix of all of the data
def makeMatrix(dataList):
	print "make matrix"

	cutoffs=[1.0]#, 0.5, 0.05, 0.01]
	count=0

	#for each cutoff value, get all of the comparisons to use in matrix and construct a matrix
	for pVal in cutoffs:

        #initialize an empty covariance matrix. scipy.zeros creates an array in a shape specified by a tuple. 
		matrix = scipy.zeros((len(dataList),len(dataList)))
		assert scipy.shape(matrix)==(13948,13948)
		for i in range(0,len(dataList)):
			for j in range(0,len(dataList)):
				
				toUseI=[]
				toUseJ=[]
				dictionary={}

				#fill dictionary with appropriate fold changes

				for gene in range(0,len(dataList[i])):
                                        #add fold changes to dicitonary if they are less than p value cutoff
					if float(dataList[i][gene][0])<=pVal or float(dataList[i][gene][0])!=pVal:
						dictionary[gene]=numpy.float64(dataList[i][gene][1])

					#else:
						#dictionary[gene]=None
					else:
						print dataList[i][gene][0]
				#use dictionary to add in overlap of fold changes in j to the list toUse

				for gene in range(0,len(dataList[j])):

					if float(dataList[j][gene][0])<=pVal or float(dataList[j][gene][0])!=pVal:
                       
                        #if the corresponding spot in the dictionary is full, there must be overlap

						if gene in dictionary:
							toUseI.append(dictionary[gene])
							toUseJ.append(numpy.float64(dataList[j][gene][1]))
						#else:
							#print dataList[j][gene][1]
					
                #We cannot compute a covariance if the lists have only one element, so this is a special case
				#that must be considered

				#if len(toUseJ)>1:
				assert len(toUseJ)==len(toUseI)


				covariance=numpy.cov(toUseI,toUseJ)[0][1]

					#count+=1
					
				#else:
					#covariance=None
					
				#round to 5 decimal places to save space

				matrix[i][j]=round(covariance,5)
		
        #pickle the final matrix
		pickle.dump(matrix,open('/home/smparker/Desktop/HeadMatrixNoNAN.p','wb'))
			




def main():

    #data='/home/smparker/Desktop/Cooper Lab/Data/HeadRevised.xlsx'
    #getData(data)
	data=pickle.load(open('/home/smparker/Desktop/HeadRawNoNAN.p','rb'))
	makeMatrix(data)
	

if __name__ == '__main__':
	main()
