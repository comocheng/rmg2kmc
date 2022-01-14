# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo'],
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006_adjusted', False), 'BurkeH2O2inArHe'],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies =['surface', 'default'],
    kineticsEstimator = 'rate rules',
)
catalystProperties(
    bindingEnergies={
                        'H': (-2.479, 'eV/molecule'),  # katrin
                        'C': (-6.750, 'eV/molecule'),  # katrin
                        'N': (-4.352, 'eV/molecule'),  # katrin
                        'O': (-3.586, 'eV/molecule'),  # katrin
    },
    # I have no idea where these come from, so I'll use Katrin's values
    # bindingEnergies={   # default values for Pt(111)
    #                    'H': (-2.479, 'eV/molecule'),  # katrin
    #                    'C': (-6.568, 'eV/molecule'),
    #                    'N': (-4.352, 'eV/molecule'),
    #                    'O': (-4.610, 'eV/molecule'),
    # },
    # bindingEnergies = {
    #                   'C':(-2.000000, 'eV/molecule'),
    #                   'O':(-6.500000, 'eV/molecule'),
    #                   'N':(-4.352, 'eV/molecule'),
    #                   'H':(-2.479, 'eV/molecule'),
    #                   },
    surfaceSiteDensity=(2.72e-9, 'mol/cm^2'),
)

# List of species
species(
    label='X',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
   label='O2',
   reactive=True,
   structure=adjacencyList(
       """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)

species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]"),
)

species(
    label='CO2',
    reactive=True,
    structure=SMILES("O=C=O"),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)


#----------
# Reaction systems
surfaceReactor(
    temperature=(600,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.041866,
        "O2": 0.03488,
        "Ar": 0.131246,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=(600,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.108574,
        "O2": 0.02088,
        "Ar": 0.78547,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=(2000,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.041866,
        "O2": 0.03488,
        "Ar": 0.131246,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=(2000,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.108574,
        "O2": 0.02088,
        "Ar": 0.78547,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.95,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=500000,
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)

#generatedSpeciesConstraints(
#    allowed=['input species','reaction libraries'],
#)
