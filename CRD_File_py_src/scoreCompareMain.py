'''
compares the local scores of two proteins
'''
import matplotlib.pyplot as plt
import math
import time
import InternalEnergyMethods

'''
@brief creates a sub-map of an atomMap from [startIndex, endIndex)
'''
def submap(atomMap, startIndex, endIndex):
    subMap = {}
    i=0
    for j in range(startIndex, endIndex):
        if j in atomMap:
            subMap[i] = atomMap[j]
        i+=1
    return subMap

def main():
    rangeN = 10
    filePath1 = "../data/model1.crd"
    filePath2 = "../data/model2.crd"
    outputPath = "proteinScoreDifference" + str(rangeN) + ".csv"
    # process protein #1
    runtimeProfile1 = InternalEnergyMethods.RuntimeProfile()
    atomMap1 = InternalEnergyMethods.atomMapFromCRD(filePath1, runtimeProfile1)
    print("Protein #1 file processing time:",
          runtimeProfile1.getDurationMS(), "milliseconds")
    # Uncomment next line to debug hash table structure
    # toFile(atomMap1, "protein1.txt")
    # process protein #2
    runtimeProfile2 = InternalEnergyMethods.RuntimeProfile()
    atomMap2 = InternalEnergyMethods.atomMapFromCRD(filePath2, runtimeProfile2)
    print("Protein #2 file processing time:",
          runtimeProfile2.getDurationMS(), "milliseconds")
    # Uncomment next line to debug hash table structure
    # toFile(atomMap2, "protein2.txt")
    differenceMap = {}
    n = len(atomMap1) + 2 -rangeN
    for i in range(0, n):
        print("Beginning Iteration", i, "to", i+rangeN)
        # find submap
        subMap1 = submap(atomMap1, i, i+rangeN)
        subMap2 = submap(atomMap2, i, i+rangeN)
        #print(subMap1)
        #print(subMap2)
        # calculate protein #1
        runtimeProfile3 = InternalEnergyMethods.RuntimeProfile()
        internalEnergy1 = InternalEnergyMethods.calculateInternalEnergy(subMap1, runtimeProfile3)
        print("Protein #1 calculate internal energy processing time:",
              runtimeProfile3.getDurationMS(), "milliseconds")
        print("The internal energy of the protein #1 is",
              int(round(internalEnergy1)), "kcal/mol")
        # calculate protein #2
        runtimeProfile4 = InternalEnergyMethods.RuntimeProfile()
        internalEnergy2 = InternalEnergyMethods.calculateInternalEnergy(subMap2, runtimeProfile4)
        print("Protein #2 calculate internal energy processing time:",
              runtimeProfile4.getDurationMS(), "milliseconds")
        print("The internal energy of the protein #2 is",
              int(round(internalEnergy2)), "kcal/mol")
        # difference map data
        differenceMap[i] = {}
        differenceMap[i]["startIndex"] = i
        differenceMap[i]["endIndex"] = i+rangeN
        differenceMap[i]["protein1"] = internalEnergy1
        differenceMap[i]["protein2"] = internalEnergy2
        differenceMap[i]["difference"] = internalEnergy1-internalEnergy2
    # write to file
    n = len(differenceMap)
    with open(outputPath, "w") as outputFile:
        outputFile.write("startIndex,endIndex,protein1,protein2,difference\n")
        for i in range(0, n):
            if i in differenceMap:
                outputFile.write(str(differenceMap[i]["startIndex"]) + ",")
                outputFile.write(str(differenceMap[i]["endIndex"]) + ",")
                outputFile.write(str(int(round(differenceMap[i]["protein1"]))) + ",")
                outputFile.write(str(int(round(differenceMap[i]["protein2"]))) + ",")
                outputFile.write(str(int(round(differenceMap[i]["difference"]))) + "\n")

if __name__ == "__main__":
    main()


