import math
import time

"""
@name Atom
@brief atomic information
PYTHON: For some reason, changing an Atom from a dictionary to an object made the code
slightly slower.
__slots__ slightly sped up the code, but dictionary method still faster. 
representing Atom as an onject will still be done for semantics 
and consistency with C++ implementation.
"""


class Atom:
    __slots__ = ("num", "x", "y", "z", "R", "Epsilon", "Sigma",
                 "Charge", "ASP", "Atm_name", "Res_name", "Res_num", "exclude")

    def __init__(self):
        self.num = int(0)
        self.x = self.y = self.z = self.R = self.Epsilon = self.Sigma = self.Charge = self.ASP = float(
            0)
        self.Atm_name = self.Res_name = self.Res_num = ""
        self.exclude = []

"""
@name RuntimeProfile
@brief A utility class passed along for function profiling for 
process time and operations performed
"""


class RuntimeProfile:
    processDuration = False

    def setDuration(self, start, stop):
        self.processDuration = stop - start

    def getDurationMS(self):
        return int(round(self.processDuration * 1000))

# Opens the .crd file and stores atom data into a python dictionary


def atomMapFromCRD(filePath, runtimeProfile):
    processStart = time.process_time()
    # Initialize a python dictionary, which is implemented as a hash table
    # internally
    atomMap = {}
    with open(filePath) as fp:
        line = fp.readline()
        i = 0
        atomNumberMax = int(line)
        line = fp.readline()
        lastAtomFound = False
        # iterate through each atom
        while ((i < atomNumberMax) and not lastAtomFound):
            # Ignore commented out lines
            if(line.strip()[:1] == "#"):
                line = fp.readline()
                continue
            # Let"s map each column in the file to a variable
            dataArray = line.split()
            atomNum = int(dataArray[0])
            atomMap[atomNum] = Atom()
            atomMap[atomNum].num = int(dataArray[0])
            atomMap[atomNum].x = float(dataArray[1])  # Angstroms
            atomMap[atomNum].y = float(dataArray[2])  # Angstroms
            atomMap[atomNum].z = float(dataArray[3])  # Angstroms
            # Van der Waals radius in Angstroms
            atomMap[atomNum].R = float(dataArray[4])
            atomMap[atomNum].Epsilon = float(dataArray[5])
            atomMap[atomNum].Sigma = float(dataArray[6])
            atomMap[atomNum].Charge = float(dataArray[7])
            atomMap[atomNum].ASP = float(dataArray[8])
            atomMap[atomNum].Atm_name = dataArray[9]
            atomMap[atomNum].Res_name = dataArray[10]
            atomMap[atomNum].Res_num = dataArray[11]
            atomMap[atomNum].exclude = []
            # print (atomMap[i]["#"], atomMap[i]["x"], atomMap[i]["y"],
            # atomMap[i]["z"])
            if(atomNum >= atomNumberMax):
                # Last atom found prematurely so we will stop looking for atoms
                lastAtomFound = True
            else:
                line = fp.readline()
            i += 1

        # Now we get the exlude list
        line = fp.readline()
        while (line):
            # Ignore commented out lines
            if(line.strip()[:1] == "#"):
                line = fp.readline()
                continue
            # Let"s map each column in the file to a variable
            dataArray = line.split()
            i = int(dataArray[0])
            exludeNumber = int(dataArray[1])
            # iterate through the exclude list
            for j in range(0, exludeNumber):
                atomMap[i].exclude.append(int(dataArray[j + 2]))
                # print (i, j, exludeNumber, int(dataArray[j+2]))
            line = fp.readline()
    processEnd = time.process_time()
    runtimeProfile.setDuration(processStart, processEnd)
    return atomMap


def toFile(atomMap, outputPath):
    n = len(atomMap) + 2
    with open(outputPath, "w") as outputFile:
        for i in range(0, n):
            if i in atomMap:
                outputFile.write(str(atomMap[i].num) + " ")
                outputFile.write("{:.4f} ".format(atomMap[i].x))
                outputFile.write("{:.4f} ".format(atomMap[i].y))
                outputFile.write("{:.4f} ".format(atomMap[i].z))
                outputFile.write("{:.4f} ".format(atomMap[i].R))
                outputFile.write("{:.4f} ".format(atomMap[i].Epsilon))
                outputFile.write("{:.4f} ".format(atomMap[i].Sigma))
                outputFile.write("{:.4f} ".format(atomMap[i].Charge))
                outputFile.write("{:.4f} ".format(atomMap[i].ASP))
                outputFile.write(str(atomMap[i].Atm_name) + " ")
                outputFile.write(str(atomMap[i].Res_name) + " ")
                outputFile.write(str(atomMap[i].Res_num) + " ")
                outputFile.write("\n")
                for atomNum in atomMap[i].exclude:
                    outputFile.write(str(atomNum) + " ")
                outputFile.write("\n")


# calculates the internal energy of a protein
# atomMap = a map of atoms
def calculateInternalEnergy(atomMap, runtimeProfile):
    processStart = time.process_time()
    """ Brad:
    I"m just implementing the equation that was in the PDF here.
    I"m not exactly sure of the physical meaning though.
    """
    internalEnergy = 0
    # radius of a water molecule in angstroms
    radius_h2o = 1.4
    coulomb_constant = 83
    # Van der waal + Electrostatic Force
    # Go through each combination of two atoms once, order doesn"t matter
    n = len(atomMap) + 2
    for i in range(0, n - 1):
        if i in atomMap:
            for j in range(i + 1, n):
                # only non-bonded atoms are calculated
                if (j in atomMap and j not in atomMap[i].exclude):
                    epiilon_ij = math.sqrt(
                        atomMap[i].Epsilon * atomMap[j].Epsilon)
                    sigma_ij = 0.5 * (atomMap[i].Sigma + atomMap[j].Sigma)
                    r_ij = math.sqrt(math.pow(atomMap[i].x - atomMap[j].x, 2) + math.pow(
                        atomMap[i].y - atomMap[j].y, 2) + math.pow(atomMap[i].z - atomMap[j].z, 2))
                    van_der_waal_energy = epiilon_ij * \
                        (math.pow(sigma_ij / r_ij, 12) -
                         2 * math.pow(sigma_ij / r_ij, 6))
                    electrostatic_force = coulomb_constant * \
                        atomMap[i].Charge * atomMap[j].Charge / r_ij
                    # energy units are kcal/mol
                    internalEnergy += van_der_waal_energy
                    internalEnergy += electrostatic_force

    # Solvation Energy
    for i in range(0, n):
        if i in atomMap:
            ASA = 0.2 * 4 * 3.1415 * math.pow(atomMap[i].R + radius_h2o, 2)
            internalEnergy += atomMap[i].ASP * ASA
    processEnd = time.process_time()
    runtimeProfile.setDuration(processStart, processEnd)
    return internalEnergy