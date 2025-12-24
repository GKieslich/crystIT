import ase
from ase.spacegroup import crystal
from ase.units import kB,mol,kJ

import spglib

import pyxtal
from pyxtal.symmetry import Group

import numpy                # arrays
import math                 # log
import os.path              # isfile, isdir
import copy                 # copy dictionary
import glob                 # iterate through dir
import time                 # for batch processing
import io                   # creating file from string
import multiprocessing      # for batch mode
import warnings             # control warning output
import traceback            # detailed error messages
warningCache = ''           # 

# Default Settings
symmetryTolerance = 5e-3                  # distance tolerance in cartesian coordinates to find crystal symmetry
occupancy = False                         # show menu to correct occupancy values
maxThreads = multiprocessing.cpu_count()  # maximum no of parallel threads
decimalSeparator = '.'
entropyOptions = False                    # calculation of entropy values from Krivovichev (2016)
recursive = False                         # subdirectory scanning in batch mode
# except for userMenu() these settings are usually forwarded through function parameters, as nested functions sometimes do not realize that global variables have been changed

# Program Information
programName = 'crystIT'
paper = 'Kaußler, Kieslich (2021). J. Appl. Cryst. 54, DOI: 10.1107/S1600576720016386'
versionNumber = '0.3'
releaseDate = '2023-04-24'
authors = 'Clemens Kaußler and Gregor Kieslich'
institution = 'Technical University of Munich'


def getComplexity(structure, pathToCif, verbose, entropy, sym):
    """
    calculates complexity of crystal structure based on an ASE Atoms object (including tags, storing CIF data-block)

    Parameters:
    arg1 (Atoms):   ASE Atoms object, including CIF data tags (store_tags = True)
    arg2 (string):  path to CIF
    arg3 (bool):    output results to console (True) or suppress console output and return result array (False)
    arg4 (bool):    entropy options
    arg5 (float):   symmetry tolerance value in cartesian coordinates

    Returns:        if (arg3 == False): array will be returned; most important variables given below:
                    if (arg4 == True):  values in {brackets} are returned additionally
    array:
        warningCache,                                                                                       errors and warnings
        chemical_formula,                                                                                   chemical formula composed from CIF-entry, ignoring dummy entries
        aSG,                                                                                                spacegroup assumed by spglib
        SG,                                                                                                 spacegroup given in CIF
        atomsPerUnitCell,                                                                                   number of atoms per crystallographic unit cell (vacancies do not count as atoms)
        atomsPerPrimitiveUnitCell,                                                                          number of atoms per primitive unit cell (vacancies do not count as atoms)
        positionsPerPrimitiveUnitCell,                                                                      amount of positions per primitive unit cell, corresponding to the sum over the crystallographic orbits' multiplicities
        uniqueSpecies,                                                                                      number of unique species, defined by combination of element (vacancies count as elements too) and crystallographic orbit
        aritySum,                                                                                           number of coordinational degrees of freedom (arities) per reduced unit cell
        I_comb, I_comb_mix, I_comb_max, I_comb_norm, I_comb_tot, I_comb_density                             combinatorial information, as defined by S. Krivovichev in 2014 and 2022 (corresponding to I_G, I_mix, I_G,max, I_G,norm, I_G,total, rho_inf), but extended by partial occupancies
        I_coor, I_coor_max, I_coor_norm, I_coor_tot, I_coor_density                                         coordinational information, as defined by W. Hornfeck in 2020
        I_conf, I_conf_max, I_conf_norm, I_conf_tot, I_conf_density                                         configurational information, as defined by W. Hornfeck in 2020, but extended by partial occupancies
        S_cfg, S_mix, S_cryst                                                                               configurational entropy (sensu strictu), mixing entropy, total entropy within the crystal
    """

    if not verbose:
        global warningCache

    # direct input of ASE Atoms object into spglib is deprecated!
    cell = (
        structure.get_cell(),
        structure.get_scaled_positions(),
        structure.get_atomic_numbers()
    )
    # find reduced unit cell
    primitiveCell = spglib.find_primitive(cell, symprec = sym)
    # get symmetry from reduced unit cell
    primitiveDataset = spglib.get_symmetry_dataset(primitiveCell, symprec = sym)
    primitiveCrystallographicOrbits = primitiveDataset['crystallographic_orbits']
    primitiveWyckoff = primitiveDataset['wyckoffs']

    # compare spacegroup set in CIF (SG) with assumed spacegroup (aSG)
    cifTags = structure.info.copy()
    try:
        iSG = cifTags['_symmetry_space_group_name_h-m']
    except:
        try:
            iSG = cifTags['_space_group_name_h-m_alt']
        except:
            iSG = 'not set'
    try:
        iSGNo = str(cifTags['_symmetry_int_tables_number'])
    except:
        try:
            iSGNo = str(cifTags['_space_group_it_number'])
        except:
            iSGNo = 'not set'
    SG = iSG + ' (' + iSGNo + ')'
    aSG = spglib.get_spacegroup(cell, symprec = sym)
    groupnumber = aSG[aSG.index('(')+1:aSG.index(')')]

    if not iSGNo == 'not set' and not iSGNo == groupnumber:
        if verbose:
            print(f'Wrong space group detected by spglib: {groupnumber} vs. {iSGNo} given in CIF. Try to alter the symmetry tolerance value. Continuing with fingers crossed.')
        else:
            warningCache += f'Wrong space group detected by spglib: {groupnumber} vs. {iSGNo} given in CIF. Try to alter the symmetry tolerance value. Continuing with fingers crossed. '

    # gather some more info about publication (batch documentation)
    try:
        journal = str(cifTags['_journal_name_full']).replace('\n', ' ').replace(';', ',')
    except:
        journal = ''
    try:
        year = str(cifTags['_journal_year'])
    except:
        year = ''
    try:
        doi = str(cifTags['_journal_paper_doi']).replace(';', '')
    except:
        doi = ''

    # compose matrix of wyckoff letters, multiplicities and arities for all crystallographic orbits
    g = Group(int(groupnumber))
    iCrystallographicOrbits = {}
    crystallographicOrbitCount = 0
    for x in numpy.unique(primitiveCrystallographicOrbits):
        iCrystallographicOrbits[crystallographicOrbitCount, 0] = numpy.count_nonzero(primitiveCrystallographicOrbits == x)          # 0 - multiplicity (in context of red uc)
        wyckoffLetter = primitiveWyckoff[list(primitiveCrystallographicOrbits).index(x)]
        iCrystallographicOrbits[crystallographicOrbitCount, 1] = wyckoffLetter                                                      #1 - wyckoff letter
        iCrystallographicOrbits[crystallographicOrbitCount, 2] = getArity(g[wyckoffLetter])                                         #2 - arity
        crystallographicOrbitCount += 1

    # identify duplicate atoms (same x,y,z coordinates = same cryst orbit) from structure in order to condense occupancyDict for all entries with identical coordinates!
    try: 
        atomSiteTypeSymbol = []
        for entry in cifTags['_atom_site_type_symbol']:
            if len(entry) > 1 and entry[1].islower():
                atomSiteTypeSymbol.append(entry[0:2])
            else:
                atomSiteTypeSymbol.append(entry[0])
    except:
        # sometimes _atom_site_type_symbol isn't set, usually when there are no fractional occupancies to consider -> extract atom species from _atom_site_label
        atomSiteTypeSymbol = []
        for entry in cifTags['_atom_site_label']:
            if len(entry) > 1 and entry[1].islower():
                atomSiteTypeSymbol.append(entry[0:2])
            else:
                atomSiteTypeSymbol.append(entry[0])

    duplicateArray = []
    identPos = []
    for x in range(0, len(atomSiteTypeSymbol)):
        XYZInfo = [
            cifTags['_atom_site_fract_x'][x],
            cifTags['_atom_site_fract_y'][x],
            cifTags['_atom_site_fract_z'][x]
        ]
        # check whether coordinates of current atom are already contained in identPos
        for y in range(0, len(identPos)):
            if numpy.allclose(XYZInfo, identPos[y], atol = sym):
                duplicateArray.append([x, y])
                break
        identPos.append(XYZInfo)

    discrepancy = len(atomSiteTypeSymbol) - crystallographicOrbitCount - len(duplicateArray)
    if discrepancy > 0:
    # same crystallographic orbit has probably been reached with different coordinates (e.g. GITWIQ) 
    # ==> construct all symmetrically equivalent positions & compare with priors. Requires significantly more computing power, therefore only executed in second step...
        duplicateArray = []
        symEquivPos = []
        for x in range(0, len(atomSiteTypeSymbol)):
            duplicate = False

            XYZInfo = [
                cifTags['_atom_site_fract_x'][x],
                cifTags['_atom_site_fract_y'][x],
                cifTags['_atom_site_fract_z'][x]
            ]

            # check whether coordinates of current atom are already contained in symEquivPos
            for y in range(0, len(symEquivPos)):
                for pos in symEquivPos[y]:
                    if numpy.allclose(XYZInfo, pos, atol = sym):
                        duplicateArray.append([x, y])
                        duplicate = True
                        break
                if duplicate:
                    break
            if not duplicate:
                # generate all symmetrically equivalent positions
                offset = len(duplicateArray) # if duplicates were identified, x has to be reduced
                wyckoffLetter = iCrystallographicOrbits[x-offset, 1]
                arity = iCrystallographicOrbits[x-offset, 2]
                # using partially parametrized positions ==> find out which wyckoff instance is present and isolate actual (x,y,z)
                if arity > 0:
                    lineNo = -1
                    for line in str(g[wyckoffLetter]).split('\n'):
                        if lineNo == -1:
                            lineNo += 1
                            continue
                        elements = line.split(',')
                        matches = 0
                        for y in range(0, 3):
                            if(
                                'x' not in elements[y]
                                and 'y' not in elements[y]
                                and 'z' not in elements[y]
                                and XYZInfo[y] == eval(elements[y])
                            ):
                                matches += 1
                        if matches == (3 - arity):
                            correctedXYZInfo = [0, 0, 0]
                            for z in range (0, 3):
                                if 'x' in elements[z]:
                                    correctedXYZInfo[0] = correctCoordinates(elements[z], 'x', XYZInfo[z])
                                elif 'y' in elements[z]:
                                    correctedXYZInfo[1] = correctCoordinates(elements[z], 'y', XYZInfo[z])
                                elif 'z' in elements[z]:
                                    correctedXYZInfo[2] = correctCoordinates(elements[z], 'z', XYZInfo[z])
                            XYZInfo = correctedXYZInfo
                            break
                        lineNo += 1
                    symEquivPos.append(
                        pyxtal.operations.filtered_coords(
                            pyxtal.operations.apply_ops(XYZInfo, g[wyckoffLetter])
                        )
                    )
            else:
                symEquivPos.append([])
        discrepancy = len(atomSiteTypeSymbol) - crystallographicOrbitCount - len(duplicateArray)

    if discrepancy == 0:
        # compose own occupancyDict, as too many errors may occur while correcting the one given by ASE (structure.info['occupancy'])
        try:
            siteOccupancy = cifTags['_atom_site_occupancy']
        except:
            siteOccupancy = []
            for i in range(0, len(atomSiteTypeSymbol)):
                siteOccupancy.append(1)

        occupancyDict = {}
        offset = 0
        for i in range(0, crystallographicOrbitCount):
            # ignore duplicates
            for entry in duplicateArray:
                if entry[0] == (i + offset):
                    offset += 1
            # add value
            occupancyDict[i] = {}
            occupancyDict[i][atomSiteTypeSymbol[i + offset]] = siteOccupancy[i + offset]
            # add all duplicates
            for entry in duplicateArray:
                if entry[1] == (i + offset):
                    try:
                        occupancyDict[i][atomSiteTypeSymbol[entry[0]]] += siteOccupancy[entry[0]]
                    except:
                        occupancyDict[i][atomSiteTypeSymbol[entry[0]]] = siteOccupancy[entry[0]]
            # double check for too high occupancy value at current crystallographic orbit
            occupancySum = 0
            for element in occupancyDict[i]:
                occupancySum += occupancyDict[i][element]
            if occupancySum > 1:
                if verbose:
                    print(f'Warning: Occupancy sum {occupancySum} at Wyckoff {iCrystallographicOrbits[i, 0]}{iCrystallographicOrbits[i, 1]}, crystallographic orbit #{i}: {occupancyDict[i]}.')
                else:
                    warningCache += f'Warning: Occupancy sum {occupancySum} at Wyckoff {iCrystallographicOrbits[i, 0]}{iCrystallographicOrbits[i, 1]}, crystallographic orbit #{i}: {occupancyDict[i]}. '
    elif verbose:
        print(f'Error: discrepancy of {discrepancy} positions between crystallographic orbits calculated by spglib and given CIF-entries. Wrong space group detected? Try to adjust symmetry tolerance!')
        return
    else:
        warningCache += f'Error: discrepancy of {discrepancy} positions between crystallographic orbits calculated by spglib and given CIF-entries. Wrong space group detected? Try to adjust symmetry tolerance! '
        return [warningCache, pathToCif]

    # allow corrections if occupancy options are enabled
    if occupancy:
        if '[' in pathToCif or verbose == False:
            print('\n\n'+pathToCif)
        occupancyDict = correctOccupancy(occupancyDict, iCrystallographicOrbits)

    # determine number of atoms in primitive unit cell and thereby compose sum formula
    # w/ occupancy (find gcd of crystal orbit muliplicities, consider occupancy)
    wyckoffSum = 0.0
    chemicalFormulaDict = {}
    numbers = []
    for i in range(0, crystallographicOrbitCount):
        numbers.append(iCrystallographicOrbits[i, 0])

    divisor = gcd(numbers)
    if divisor < 0:
        divisor = 1

    counter = 0
    for x in occupancyDict:
        multiplicity = iCrystallographicOrbits[counter, 0]
        for element in occupancyDict[x]:
            try:
                chemicalFormulaDict[element] += occupancyDict[x][element] * multiplicity / divisor
            except:
                chemicalFormulaDict[element] = occupancyDict[x][element] * multiplicity / divisor
            wyckoffSum += occupancyDict[x][element] * multiplicity
        counter += 1    

    # sometimes gcd of multiplicities does not yield empirical formula (e.g. Cu2P6O18Li2 / MnN10C18H28)
    # better safe than sorry: try to reduce formula a second time
    # (multiplicity approach still implemented bc fractional occupancies often complicate computation of gcd)
    numbers = []
    for element in chemicalFormulaDict:
        # suppose: a) lacking precision
        if abs(chemicalFormulaDict[element] - round(chemicalFormulaDict[element])) < 0.1:
            numbers.append(round(chemicalFormulaDict[element]))
        # or b) more severe defects
        else:
            numbers.append(math.ceil(chemicalFormulaDict[element]))
    if not numbers:
        divisor = 1
    else:
        divisor = gcd(numbers)
        if divisor < 0:
            divisor = 1

    # compose assumed chemical formula
    chemical_formula = ''
    for element in sorted(chemicalFormulaDict):
        stoichiometry = chemicalFormulaDict[element] / divisor
        if stoichiometry == 1:
            stoichiometry = ''
        elif stoichiometry % 1 == 0:
            stoichiometry = str(int(stoichiometry))
        else:
            stoichiometry = str(stoichiometry)
        chemical_formula = chemical_formula + element + stoichiometry

    atomsPerPrimitiveUnitCell = wyckoffSum
    atomsPerUnitCell = wyckoffSum * len(structure) / len(primitiveCrystallographicOrbits) 

    # sum over multiplicities / arities of all crystallographic orbits
    positionsPerPrimitiveUnitCell = aritySum = 0
    for x in range(0, crystallographicOrbitCount):
        positionsPerPrimitiveUnitCell += iCrystallographicOrbits[x, 0]
        aritySum += iCrystallographicOrbits[x, 2]

    # calculate information contents
    I_comb = I_coor = I_conf = I_comb_str = I_comb_mix = 0.0
    uniqueSpecies = 0

    for x in range(0, crystallographicOrbitCount):
        # the coordinational sum is formed over all crystallographic orbits
        if aritySum > 0:
            probability = iCrystallographicOrbits[x, 2] / aritySum
            if probability > 0:
                I_coor -= probability * math.log(probability, 2)

            # coordinational part of configurational sum
            probability = iCrystallographicOrbits[x, 2] / (aritySum + positionsPerPrimitiveUnitCell)
            if probability > 0:
                I_conf -= probability * math.log(probability, 2)

        # the combinatorial sum is formed over each element in a crystallographic orbit individually (in other words: over unique species)
        # mixing complexity is calculated according to Krivovichev 2022
        # vacancies count as elements too -> probability according to positionsPerPrimitiveUnitCell 
        occupancySum = 0
        multiplicity = iCrystallographicOrbits[x, 0]
        
        probability = multiplicity / positionsPerPrimitiveUnitCell
        if probability > 0: 
            I_comb_str -= probability * math.log(probability, 2)
            uniqueSpecies += 1
            
            for element in occupancyDict[x]:
                occupancyValue = occupancyDict[x][element]
                occupancySum += occupancyDict[x][element]
                if occupancyValue == 0:
                    break
                if occupancyValue < 1:
                    I_comb_mix -= probability*occupancyValue*math.log(occupancyValue, 2)   
            # add vacancy species for every crystallographic orbit that is only partially occupied
            if occupancySum < 1:
                I_comb_mix -= probability*(1 - occupancySum)*math.log(1 - occupancySum, 2)
                uniqueSpecies += 1
            
            #combinatorial part of configurational sum
            probability = multiplicity / (aritySum + positionsPerPrimitiveUnitCell)
            I_conf -= probability * math.log(probability,2)
            for element in occupancyDict[x]:
                occupancyValue = occupancyDict[x][element]
                if occupancyValue == 0:
                    break
                if occupancyValue < 1:
                    I_conf += probability*occupancyValue*math.log(occupancyValue, 2)
            if occupancySum < 1:
                I_conf += probability*(1 - occupancySum)*math.log(1 - occupancySum, 2)
            
        elif verbose:
            print(f'Probability <= 0 was skipped: {element} at pos. {x}')
        else:
            warningCache += f'Probability <= 0 was skipped: {element} at pos. {x} '


    I_comb = I_comb_str - I_comb_mix
    I_comb_tot = positionsPerPrimitiveUnitCell * I_comb
    I_coor_tot = aritySum * I_coor
    I_conf_tot = (aritySum + positionsPerPrimitiveUnitCell) * I_conf
    

    # maximum combinatorial information content based on number of unique species which are defined by a combination of crystallographic orbit and element (vacancies obviously count too).
    # otherwise: I_comb > I_comb_max for alloys (in general: cases w/ all occupancies < 1)
    # maximum coordinational information content based on number of different crystallographic orbits (in general: cardinality)
    I_comb_max = math.log(uniqueSpecies, 2)
    I_coor_max = math.log(crystallographicOrbitCount, 2)
    I_conf_max = math.log(uniqueSpecies + crystallographicOrbitCount, 2)

    if I_comb_max != 0:
        I_comb_norm = I_comb / I_comb_max
    else:
        I_comb_norm = 0
    if I_coor_max != 0:
        I_coor_norm = I_coor / I_coor_max
    else:
        I_coor_norm = 0
    if I_conf_max != 0:
        I_conf_norm = I_conf / I_conf_max
    else:
        I_conf_norm = 0

    # correct cell volume to primitive cell volume
    perVolume = atomsPerUnitCell / (atomsPerPrimitiveUnitCell * structure.cell.volume)
    I_comb_density = perVolume * I_comb_tot
    I_coor_density = perVolume * I_coor_tot
    I_conf_density = perVolume * I_conf_tot

    if entropy:
        gasConstantR = mol * kB / (kJ / 1000)
        conversionFactor = math.log(2, math.e)

        S_mix = gasConstantR * I_comb_mix * positionsPerPrimitiveUnitCell * conversionFactor #per reduced unit cell
        
        # exact configurational entropy using Boltzmann equation on crystallographic orbits
        S_cfg = 0.0
        for x in range(0, crystallographicOrbitCount):
            if iCrystallographicOrbits[x, 0] <= 90:
                S_cfg += math.log(math.factorial(iCrystallographicOrbits[x, 0]), math.e)
            else:
                S_cfg += iCrystallographicOrbits[x, 0] * (math.log(iCrystallographicOrbits[x, 0], math.e) - 1)
        S_cfg = gasConstantR * S_cfg
        S_cryst = S_cfg + S_mix


    if verbose:
        print(f'\n\n------------ {pathToCif} ------------')
        print(f'assumed formula\t {chemical_formula}')
        print(f'assumed SG\t {aSG}')
        print(f'SG from CIF\t {SG}')
        print(
            'lattice [A] \t a: {:.2f}, b: {:.2f}, c: {:.2f}'.format(
                structure.get_cell_lengths_and_angles()[0],
                structure.get_cell_lengths_and_angles()[1],
                structure.get_cell_lengths_and_angles()[2]
            ).replace('.', decimalSeparator)
        )
        print(
            'angles [°] \t b,c: {:.2f}, a,c: {:.2f}, a,b: {:.2f}'.format(
                structure.get_cell_lengths_and_angles()[3],
                structure.get_cell_lengths_and_angles()[4],
                structure.get_cell_lengths_and_angles()[5]
            ).replace('.', decimalSeparator)
        )

        print('---')
        print('{:.6f} \t atoms / unit cell'.format(atomsPerUnitCell).replace('.', decimalSeparator))
        print('{:.6f} \t atoms / reduced unit cell'.format(atomsPerPrimitiveUnitCell).replace('.', decimalSeparator))
        print('{:.6f} \t positions / reduced unit cell'.format(positionsPerPrimitiveUnitCell).replace('.', decimalSeparator))
        print('{:.6f} \t crystallographic orbits'.format(crystallographicOrbitCount).replace('.', decimalSeparator))
        print('{:.6f} \t unique species'.format(uniqueSpecies).replace('.', decimalSeparator))
        print('{:.6f} \t coordinational degrees of freedom (arities)'.format(aritySum).replace('.', decimalSeparator))

        print('--- combinatorial (extended Krivovichev) ---')
        print('{:.6f} \t I_comb \t\t [bit / position]'.format(I_comb).replace('.', decimalSeparator))
        print('{:.6f} \t I_comb_mix \t\t [bit / position]'.format(I_comb_mix).replace('.', decimalSeparator))
        print('{:.6f} \t I_comb_max \t\t [bit / position]'.format(I_comb_max).replace('.', decimalSeparator))
        print('{:.6f} \t I_comb_norm \t\t [-]'.format(I_comb_norm).replace('.', decimalSeparator))
        print('{:.6f} \t I_comb_tot \t\t [bit / reduced unit cell]'.format(I_comb_tot).replace('.', decimalSeparator))
        print('{:.6f} \t I_comb_dens \t\t [bit / A^3]'.format(I_comb_density).replace('.', decimalSeparator))

        print('--- coordinational (Hornfeck) ---')
        print('{:.6f} \t I_coor \t\t [bit / freedom]'.format(I_coor).replace('.', decimalSeparator))
        print('{:.6f} \t I_coor_max \t\t [bit / freedom]'.format(I_coor_max).replace('.', decimalSeparator))
        print('{:.6f} \t I_coor_norm \t\t [-]'.format(I_coor_norm).replace('.', decimalSeparator))
        print('{:.6f} \t I_coor_tot \t\t [bit / reduced unit cell]'.format(I_coor_tot).replace('.', decimalSeparator))
        print('{:.6f} \t I_coor_dens \t\t [bit / A^3]'.format(I_coor_density).replace('.', decimalSeparator))

        print('--- configurational (extended Hornfeck) ---')
        print('{:.6f} \t I_conf \t\t [bit / (position + freedom)]'.format(I_conf).replace('.', decimalSeparator))
        print('{:.6f} \t I_conf_max \t\t [bit / (position + freedom)]'.format(I_conf_max).replace('.', decimalSeparator))
        print('{:.6f} \t I_conf_norm \t\t [-]'.format(I_conf_norm).replace('.', decimalSeparator))
        print('{:.6f} \t I_conf_tot \t\t [bit / reduced unit cell]'.format(I_conf_tot).replace('.', decimalSeparator))
        print('{:.6f} \t I_conf_dens \t\t [bit / A^3]'.format(I_conf_density).replace('.', decimalSeparator))
        
        if entropy:
            print('--- entropy ---')
            print('{:.6f} \t S_cfg \t\t [J / (mol (reduced unit cell) * K)]'.format(S_cfg).replace('.', decimalSeparator))
            print('{:.6f} \t S_mix \t\t [J / (mol (reduced unit cell) * K)]'.format(S_mix).replace('.', decimalSeparator))
            print('{:.6f} \t S_cryst \t\t [J / (mol (reduced unit cell) * K)]'.format(S_cryst).replace('.', decimalSeparator))
        
        return

    elif entropy:
        returnArray = [
            warningCache,
            pathToCif,
            doi, journal, year,
            chemical_formula,
            aSG,
            SG,
            structure.get_cell_lengths_and_angles()[0],
            structure.get_cell_lengths_and_angles()[1],
            structure.get_cell_lengths_and_angles()[2],
            structure.get_cell_lengths_and_angles()[3],
            structure.get_cell_lengths_and_angles()[4],
            structure.get_cell_lengths_and_angles()[5],
            atomsPerUnitCell,
            atomsPerPrimitiveUnitCell,
            positionsPerPrimitiveUnitCell,
            uniqueSpecies,
            aritySum,
            I_comb, I_comb_mix, I_comb_max, I_comb_norm, I_comb_tot, I_comb_density, 
            I_coor, I_coor_max, I_coor_norm, I_coor_tot, I_coor_density,
            I_conf, I_conf_max, I_conf_norm, I_conf_tot, I_conf_density,
            S_cfg, S_mix, S_cryst
        ]
    else:
        returnArray = [
            warningCache,
            pathToCif,
            doi, journal, year,
            chemical_formula,
            aSG,
            SG,
            structure.get_cell_lengths_and_angles()[0],
            structure.get_cell_lengths_and_angles()[1],
            structure.get_cell_lengths_and_angles()[2],
            structure.get_cell_lengths_and_angles()[3],
            structure.get_cell_lengths_and_angles()[4],
            structure.get_cell_lengths_and_angles()[5],
            atomsPerUnitCell,
            atomsPerPrimitiveUnitCell,
            positionsPerPrimitiveUnitCell,
            uniqueSpecies,
            aritySum,
            I_comb, I_comb_mix, I_comb_max, I_comb_norm, I_comb_tot, I_comb_density,
            I_coor, I_coor_max, I_coor_norm, I_coor_tot, I_coor_density,
            I_conf, I_conf_max, I_conf_norm, I_conf_tot, I_conf_density
        ]
    return returnArray


def correctCoordinates(coordinateDescription, parameter, coordinate):
    """
    extracts x/y/z parameter of a wyckoff position's individual coordinates. e.g. the z-coordinate of a wyckoff position 4c in SG 24 might be defined as (-z+1/2) = 0.3 --> returns (z) = 0.2

    Parameters
    arg1 (string)       parametrized description of the coordinate                                      e.g. '-z+1/2'
    arg2 (string)       'x', 'y' or 'z' as parameter to isolate from arg1 (coordinateDescription)       e.g. 'z'
    arg3 (float)        fractional coordinate on x/y/z axis                                             e.g. 0.3

    Returns
    float               fractional coordinate, corresponding to the isolated parameter (x, y or z)      e.g. 0.2
    """
    if coordinateDescription.split(parameter)[0] == '-':
        factor = -1
    else:
        factor = +1
    if coordinateDescription.split(parameter)[1] != '':
        summand = eval(coordinateDescription.split(parameter)[1])
    else:
        summand = 0
    return (factor * (coordinate - summand)) % 1


def getArity(pyxtalWyckoff):
    """
    calculates the arity of a given wyckoff position

    Parameters
    arg1 (Wyckoff_position)     pyxtal Wyckoff_position class object

    Returns
    int                         arity
    """
    firstSymmOp = str(pyxtalWyckoff).splitlines()[1] # line 0 contains general description: 'wyckoff pos nA in SG xx with site symmetry xx'
    arity = 0
    if 'x' in firstSymmOp:
        arity += 1
    if 'y' in firstSymmOp:
        arity += 1
    if 'z' in firstSymmOp:
        arity += 1
    return arity


def correctOccupancy(occupancyDict, iCrystallographicOrbits):
    """
    a menu that allows for on-the-fly editing of occupancy values
    
    Parameters
    arg1 (dictionary)       dictionary, containing {Element1 : occupancy1, Element2 : occupancy2} for every crystallographic orbit
    arg2 (array)            array, containing the multiplicities [x, 0], wyckoff letters [x, 1] and arities [x, 2] of every crystallographic orbit
    
    Returns
    dictionary              updated occupancyDict
    """

    corrOccupancyDict = copy.deepcopy(occupancyDict)

    while True:
        print('\n\nEnter a number on the left to correct the species\' occupancy. \'c\' to continue with current values. \'d\' to discard changes.')
        print('#\t Element \t Wyckoff \t arity \t original \t current')

        positions = []
        for x in corrOccupancyDict:
            for element in corrOccupancyDict[x]:
                positions.append([x,element])
                print(f'{len(positions) - 1} \t {element} \t\t {iCrystallographicOrbits[x, 0]}{iCrystallographicOrbits[x, 1]} \t\t {iCrystallographicOrbits[x, 2]} \t {occupancyDict[x][element]} \t\t {corrOccupancyDict[x][element]}')
            print('')

        userInput = input()
        if userInput == 'c':
            return corrOccupancyDict
        elif userInput == 'd':
            return occupancyDict
        elif RepresentsInt(userInput) and 0 <= int(userInput) < len(positions):
            x = positions[int(userInput)][0]
            element = positions[int(userInput)][1]
            print(f'\n\nInput the new stoichiometry for {element} at Wyckoff {iCrystallographicOrbits[x, 0]}{iCrystallographicOrbits[x, 1]} with \'.\' as decimal separator. Currently: {corrOccupancyDict[x][element]}')

            userInput2 = input()
            if RepresentsFloat(userInput2) and 0 < float(userInput2) <= 1:
                corrOccupancyDict[x][element] = float(userInput2)
            else:
                print(f'\n\nPlease only insert occupancy values  0 < x <= 1')
            continue
        else:
            print(f'\n\nPlease only enter integer numbers in the range of 0 to {len(positions) - 1}')
            continue


def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


def RepresentsFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False


def gcd(numbers):
    """
    calculates the greatest common divisor of a given array of integers
    """
    divisor = numbers[0]
    while True:
        try:
            for number in numbers:
                rest = number % divisor
                if not rest:
                    pass
                else:
                    raise
            break
        except:
            pass
        divisor -= 1
    return divisor


def customWarnings(message, category, filename, lineno, file, random):
    """
    redirects warnings into the global variable warningCache (batch mode)
    """
    global warningCache
    warningCache += str(message) + ', in: ' + str(filename) + ', line: ' + str(lineno) + ' '


def processFile(pathToCif, verbose, entropy, symprec):
    """
    open CIF from given path, perform corrections that enhance ASE-compatibility and facilitate calculations in getComplexity()
    let ASE parse the file and forward the data blocks in form of Atoms objects to getComplexity()
    
    Parameters:
    arg1 (string)   path to valid CIF
    arg2 (Boolean)  verbosity: (True) --> output to console <=> (False) --> output to .csv-file in respective folder
    arg3 (Boolean)  entropy options
    arg4 (float)    symmetry tolerance in cartesian coordinates
    
    Returns:
    returns return valuess of getComplexity() as an array
    """

    # redirect warnings for batch mode
    if not verbose:
        resultArray = []
        global warningCache
        warnings.showwarning = customWarnings

    # get contents from CIF-file and thereby correct spacegroups that are written with brackets (ASE will throw errors)
    # crystal water is often denominated as "Wat", ASE hates that, replace "Wat" with "O" as hydrogen atoms are missing anyway
    # ignore dummy atoms completely as they will cause problems and should not contribute to any information content
    # filter fractional coordinates with modulo operator (should be between 0 and 1!), thereby discard of uncertainty values
    input = open(pathToCif)
    output = ''
    xPos = yPos = zPos = counter = -1
    for line in input:
        low = line.lower()
        if line[0] == '#':
            continue
        elif '_' in line:
            if (
                '_symmetry_space_group_name_h-m' in low 
                or '_space_group_name_h-m_alt' in low
            ):
                output += line.replace('(', '').replace(')', '')
            elif 'loop_' in low:
                output += line
                xPos = yPos = zPos = counter = -1
            elif '_atom_site_fract_x' in low:
                output += line
                xPos = counter
            elif '_atom_site_fract_y' in low:
                output += line
                yPos = counter
            elif '_atom_site_fract_z' in low:
                output += line
                zPos = counter
            else:
                output += line
            counter += 1
        elif xPos >= 0 and yPos >=0 and zPos >= 0:
            if 'dum' in low:
                continue

            segments = line.split()
            if len(segments) > max([xPos, yPos, zPos]):
                if '(' in segments[xPos]:
                    segments[xPos] = segments[xPos][0:segments[xPos].find('(')]
                if '(' in segments[yPos]:
                    segments[yPos] = segments[yPos][0:segments[yPos].find('(')]
                if '(' in segments[zPos]:
                    segments[zPos] = segments[zPos][0:segments[zPos].find('(')]
                if RepresentsFloat(segments[xPos]):
                    segments[xPos] = str(float(segments[xPos]) % 1)
                if RepresentsFloat(segments[yPos]):
                    segments[yPos] = str(float(segments[yPos]) % 1)
                if RepresentsFloat(segments[zPos]):
                    segments[zPos] = str(float(segments[zPos]) % 1)

                for segment in segments:
                    output += '    '
                    output += segment.replace('Wat', 'O')
                output += '\n'
            else:
                output += line.replace('Wat', 'O')
        else:
            output += line

    cifFile = io.StringIO(output)

    #let ase read adjusted CIF-file
    try:
        structureList = ase.io.read(cifFile, format = 'cif', index = ':', store_tags = True, reader = 'ase') #, fractional_occupancies = True
    except Exception as e:
        errorMessage = 'File is either empty or corrupt. ' + traceback.format_exc().replace('\n', ' ')
        if verbose:
            print(errorMessage)
            return
        else:
            errorMessage += warningCache
            warningCache = ''
            resultArray.append([errorMessage, pathToCif])
            return resultArray

    # iterate through entries in CIF-file
    index = 0
    for structure in structureList:
        outputPath = pathToCif
        if len(structureList) > 1:
            outputPath = outputPath + ' [' + str(index) + ']'
        try:
            if verbose:
                getComplexity(structure, outputPath, verbose, entropy, symprec)
            else:
                resultArray.append(getComplexity(structure, outputPath, verbose, entropy, symprec))
        except Exception as e:
            errorMessage = 'Error: ' + traceback.format_exc().replace('\n', ' ')
            if verbose:
                print(errorMessage)
            else:
                warningCache += errorMessage
                resultArray.append([warningCache, outputPath])
        warningCache = ''
        index += 1

    if not verbose:
        return resultArray


def processDirectory(dir, recursive, entropy, symprec):
    """
    iterates through all .cif-files in a given directory with multithreading and compiles results into .csv-file

    Parameters:
    arg1 (string):  path to directory 
    arg2 (Boolean): iterate through subdirs as well?
    arg3 (Boolean): entropy options
    arg4 (float):   symmetry tolerance in cartesian coordinates

    Returns: results as .csv-file into dir
    """

    start = time.time()
    if not dir[-1] == '/' and not dir[-1] == '\\':
        dir += '\\'

    if recursive:
        extension = '**/*.cif'
    else:
        extension = '*.cif'

    resultArray = []
    fileList = glob.glob(dir + extension, recursive = recursive)
    numFiles = len(fileList)
    if numFiles == 0:
        print(f'\n\n{dir} does not contain .cif-files. Try to activate recursive subdirectory scanning in the settings.')
        return

    if numFiles > maxThreads:
        numProcesses = maxThreads
    else:
        numProcesses = numFiles
    pool = multiprocessing.Pool(processes = numProcesses)
    for file in fileList:
        resultArray.append(pool.apply_async(processFile, args = (file, False, entropy, symprec)))

    output = ''
    numEntries = 0
    for fileResult in resultArray:
        for cifResult in fileResult.get():
            counter = 0
            numEntries += 1
            for string in cifResult:
                if counter > 7:
                    if decimalSeparator == ',':
                        output += '{:.6f}; '.format(string).replace('.', ',')
                    else:
                        output += '{:.6f}; '.format(string)
                else:
                    output += string + '; '
                counter += 1
            output += '\n '

    if entropyOptions:
        header = 'Errors; Path; DOI; Journal; Year; Assumed Formula; assumed SG; SG from CIF; a [A]; b [A]; c [A]; b,c [°]; a,c [°]; a,b [°]; atoms / uc; atoms / reduc; pos / reduc; unique species; coor freedom (aritySum); I_comb; I_comb_mix; I_comb_max; I_comb_norm; I_comb_tot; I_comb_density; I_coor; I_coor_max; I_coor_norm; I_coor_tot; I_coor_density; I_conf; I_conf_max; I_conf_norm; I_conf_tot; I_conf_density;  S_cfg; S_mix; S_cryst; \n '
    else:
        header = 'Errors; Path; DOI; Journal; Year; Assumed Formula; assumed SG; SG from CIF; a [A]; b [A]; c [A]; b,c [°]; a,c [°]; a,b [°]; atoms / uc; atoms / reduc; pos / reduc; unique species; coor freedom (aritySum); I_comb; I_comb_mix; I_comb_max; I_comb_norm; I_comb_tot; I_comb_density; I_coor; I_coor_max; I_coor_norm; I_coor_tot; I_coor_density; I_conf; I_conf_max; I_conf_norm; I_conf_tot; I_conf_density; \n '

    finish = time.time()
    outputFile = dir + f'batch_{int(finish)}.csv'
    f = open(outputFile, 'w', encoding = 'utf-8')
    f.write(header + output)
    f.close()

    timer = '{:.3f}'.format(finish - start)
    print(f'\n\nProcessed {numFiles} files ({numEntries} entries) in {timer} s. Results written into {outputFile}')


def userMenu():
    global symmetryTolerance
    global occupancy
    global maxThreads
    global decimalSeparator
    global entropyOptions
    global recursive

    print(
        f'Welcome to {programName} -- A Crystal Structure Complexity Analyzer Based on Information Theory\n'
        + f'Version {versionNumber}, release date: {releaseDate}\n'
        + f'Written by {authors} ({institution})\n'
        + f'Please cite the following paper if {programName} is utilized in your work:\n'
        + f'\t {paper}'
    )

    while True:
        print(f'\n\nInput path of .cif file or directory for complexity analysis. \'s\' for settings. \'e\' to exit.')
        userInput = input().replace('\"', '')
        if userInput == 'exit' or userInput == 'e':
            break
        elif userInput == 's':
            while True:
                print(
                    f'\n\nInput float as symmetry tolerance 0 < x < 1\t (currently {symmetryTolerance}).'
					+ f'\nInput int as maximum number of threads\t\t (currently {maxThreads})'
					+ f'\n\'d\' to toggle between decimal separators\t (currently \'{decimalSeparator}\').'
					+ f'\n\'o\' to toggle occupancy editing options\t\t (currently {occupancy}).'
					+ f'\n\'r\' to toggle recursive subdir scan\t\t (currently {recursive}). '
					+ f'\n\'s\' to toggle entropy calculation\t\t (currently {entropyOptions}).'
					+ '\n\'e\' exit to main menu:'
                )
                userInput = input()
                if userInput == 'o':
                    occupancy = not occupancy
                elif userInput == 'r':
                    recursive = not recursive
                elif userInput == 's':
                    entropyOptions = not entropyOptions
                elif userInput == 'd':
                    if decimalSeparator == '.':
                        decimalSeparator = ','
                    else:
                        decimalSeparator = '.'
                elif userInput == 'e' or userInput == 'exit':
                    break
                elif RepresentsFloat(userInput) and 0 < float(userInput) < 1:
                    symmetryTolerance = float(userInput)
                elif RepresentsInt(userInput) and int(userInput) > 0:
                    maxThreads = int(userInput)
                else:
                    print('\n\nInvalid input')
                continue
            continue
        elif os.path.isdir(userInput):
            processDirectory(userInput, recursive, entropyOptions, symmetryTolerance)
            continue
        elif '.' in userInput:
            extension = userInput.split('.')[-1]
            if extension != 'cif':
                userInput += '.cif'
        else:
            userInput += '.cif'
        if os.path.isfile(userInput):
            processFile(userInput, True, entropyOptions, symmetryTolerance)
        else:
            print('\n\nInvalid path')
        continue


if __name__ == '__main__':
    userMenu()