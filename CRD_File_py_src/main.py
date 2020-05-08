import math
import time
import InternalEnergyMethods
import matplotlib

def main():
    filePath1 = "../data/model1.crd"
    filePath2 = "../data/model2.crd"
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
    # calculate protein #1
    runtimeProfile3 = InternalEnergyMethods.RuntimeProfile()
    internalEnergy1 = InternalEnergyMethods.calculateInternalEnergy(atomMap1, runtimeProfile3)
    print("Protein #1 calculate internal energy processing time:",
          runtimeProfile3.getDurationMS(), "milliseconds")
    print("The internal energy of the protein #1 is",
          int(round(internalEnergy1)), "kcal/mol")
    # calculate protein #2
    runtimeProfile4 = InternalEnergyMethods.RuntimeProfile()
    internalEnergy2 = InternalEnergyMethods.calculateInternalEnergy(atomMap2, runtimeProfile4)
    print("Protein #2 calculate internal energy processing time:",
          runtimeProfile4.getDurationMS(), "milliseconds")
    print("The internal energy of the protein #2 is",
          int(round(internalEnergy2)), "kcal/mol")
    # Protein comparison
    if(internalEnergy1 < internalEnergy2):
        print(
            "The internal energy of the protein #1 is less, so it is more likely to occur")
    else:
        print(
            "The internal energy of the protein #2 is less, so it is more likely to occur")

if __name__ == "__main__":
    main()
