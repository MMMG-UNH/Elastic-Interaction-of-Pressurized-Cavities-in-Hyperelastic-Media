#-*- coding: utf-8 -*-
# =============================================================================
# Two_pressurized_holes_ABAQUS_script_2D.py
#
# Full pipeline for simulating two equal, pressurized cylindrical cavities in
# an infinite 2D elastic medium using ABAQUS/Standard.
#
# The script is organized into three sequential sections:
#   1. Model creation and input file generation
#   2. Job submission and execution
#   3. Post-processing: strain energy and deformed cavity coordinates
#
# Problem setup:
#   - Two circular holes of radius R, centered symmetrically on the x-axis
#     at x = ±(d/2), embedded in a large circular domain of radius R_m.
#   - Uniform internal pressure P is applied to both hole boundaries.
#   - The outer boundary is traction-free.
#   - Plane-strain conditions are assumed (2D).
#   - The separation η (center-to-center distance) is swept from d_min to
#     d_max in steps of d_increment.
#
# Supported material models (set via `material_model`):
#   'NH'  - Neo-Hookean hyperelastic (C10, D1)
#   'LE'  - Linear elastic (E, nu)
#   'MR'  - Mooney-Rivlin hyperelastic (C10, C01, D1)
#   'AB'  - Arruda-Boyce hyperelastic (C2, lambda_m, D)
#
# Outputs:
#   - ABAQUS .inp files (one per separation distance)
#   - StrainEnergy.rpt: total strain energy at each separation (appended)
#   - deformed_coordinates-{index}.csv: deformed (x,y) coordinates of left-
#     hole boundary nodes (used to compute cavity area change ΔA for PE)
#
# Dependencies: ABAQUS 2022 or compatible, NumPy (bundled with ABAQUS Python)
#
# Usage: Run inside ABAQUS/CAE via File > Run Script, or:
#   abaqus cae noGUI=Two_pressurized_holes_ABAQUS_script_2D.py
#
# Reference:
#   Saeedi & Kothari (2025). J. Appl. Mech., 92(5), 051008.
#   DOI: 10.1115/1.4067855
# =============================================================================

# --- ABAQUS module imports (must run within ABAQUS Python environment) ---
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import os
from abaqus import *
from abaqusConstants import *
import regionToolset
from odbAccess import *
import csv
import numpy as np

# Set working directory for all file I/O (input files, ODB files, reports)
base_dir = r"YOUR PATH\files"
os.chdir(base_dir)


# =============================================================================
# SECTION 1: MODEL CREATION AND INPUT FILE GENERATION
#
# Creates one ABAQUS model per separation distance, applies material,
# boundary conditions, mesh, and writes a .inp file for each.
# =============================================================================

# --- Simulation parameters ---
P           = 30    # Internal pressure applied to both hole surfaces
d_min       = 0.22  # Minimum center-to-center distance (must be > 2R = 0.2)
d_max       = 2.0 + 0.0001  # Maximum center-to-center distance (inclusive bound)
d_increment = 0.02  # Step size for sweeping separation distances
R_m, R      = 5.0, 0.1  # Domain radius (truncation of infinite domain), hole radius

# Material model selector: 'NH', 'LE', 'MR', or 'AB'
material_model = 'NH'

# Mesh seed sizes (element edge lengths at different regions of the domain)
# Finer near the holes, coarser toward the outer boundary
HoleSeed, MiddleSeed1, MiddleSeed2, OuterSeed = 0.002, 0.06, 0.12, 0.2

# --- Material property definitions ---
# Each model requires different constitutive parameters.
# D1 = 0 enforces incompressibility for hyperelastic models.
if material_model == 'NH':
    # Neo-Hookean: W = C10*(I1 - 3), shear modulus mu = 2*C10
    C10, D1 = 10.0, 0.0
    material_type = NEO_HOOKE
    table = ((C10, D1),)
elif material_model == 'LE':
    # Linear elastic: standard Hooke's law
    E, nu = 52.0, 0.3
    material_type = ELASTIC
    table = ((E, nu),)
elif material_model == 'MR':
    # Mooney-Rivlin: W = C10*(I1-3) + C01*(I2-3), mu = 2*(C10+C01)
    # alpha = C10/(C10+C01) controls strain-stiffening (here alpha=0.5)
    C10, C01, D1 = 5.0, 5.0, 0.0
    material_type = MOONEY_RIVLIN
    table = ((C10, C01, D1),)
elif material_model == 'AB':
    # Arruda-Boyce (8-chain model): C2 and lambda_m are material constants.
    # lambda_m is the limiting chain stretch (controls strain-stiffening onset).
    # Shear modulus: mu = C2*(1 + 3/(5*lambda_m^2) + ...)  (consistency cond.)
    C2, lambda_m, D = 16.7091, 2.0, 0.0
    material_type = ARRUDA_BOYCE
    table = ((C2, lambda_m, D),)
else:
    raise ValueError("Unknown material model specified.")

# --- Generate the list of hole-center x-positions (right hole only; left is mirrored) ---
# circle_values is a list of (x_center, 0) for the right hole at each separation
# The right hole center is at x = d/2, left hole at x = -d/2
circle_values = [(d_min/2 + d_increment/2 * i, 0.0)
                 for i in range(int((d_max/2 - d_min/2) / (d_increment/2) + 1))]

# --- Loop: create one ABAQUS model per separation distance ---
for index, circle in enumerate(circle_values):
    model_name = 'Model-' + str(index + 1)
    part_name  = 'Part-'  + str(index + 1)
    step_name  = 'Step-'  + str(index + 1)

    # dR = eta/R = center-to-center distance / hole radius (non-dimensional separation)
    dR       = 2 * circle[0] / R
    job_name = "P{0}_dR{1:.1f}".format(int(P), dR)
    job_name = job_name.replace('.', '_')  # replace dot for valid filename

    # --- Create a new ABAQUS model ---
    mdb.Model(name=model_name, modelType=STANDARD_EXPLICIT)

    # --- Sketch the 2D geometry ---
    # Outer circle: large domain approximating an infinite medium
    mdb.models[model_name].ConstrainedSketch(name='__profile__', sheetSize=20.0)
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(
        center=(0.0, 0.0), point1=(-R_m, 0.0))

    # Right hole centered at (circle[0], 0), radius R
    center1_x = circle[0]
    center2_x = -center1_x  # left hole (symmetric about x=0)
    point1_x  = center1_x + R
    point2_x  = center2_x - R
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(
        center=(center1_x, 0.0), point1=(point1_x, 0.0))
    # Left hole (mirrored)
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(
        center=(center2_x, 0.0), point1=(point2_x, 0.0))

    # Create the 2D planar deformable part from the sketch
    mdb.models[model_name].Part(
        dimensionality=TWO_D_PLANAR, name=part_name, type=DEFORMABLE_BODY)
    mdb.models[model_name].parts[part_name].BaseShell(
        sketch=mdb.models[model_name].sketches['__profile__'])
    del mdb.models[model_name].sketches['__profile__']

    # --- Assign material and section properties ---
    mdb.models[model_name].Material(name='Material-1')
    if material_model == 'NH':
        mdb.models[model_name].materials['Material-1'].Hyperelastic(
            materialType=ISOTROPIC, table=table, testData=OFF,
            type=material_type, volumetricResponse=VOLUMETRIC_DATA)
    elif material_model == 'LE':
        mdb.models[model_name].materials['Material-1'].Elastic(table=table)
    elif material_model in ('MR', 'AB'):
        mdb.models[model_name].materials['Material-1'].Hyperelastic(
            materialType=ISOTROPIC, table=table, testData=OFF,
            type=material_type, volumetricResponse=VOLUMETRIC_DATA)

    # Plane-strain section with unit thickness (2D problem)
    mdb.models[model_name].HomogeneousSolidSection(
        material='Material-1', name='Section-1', thickness=1.0)
    mdb.models[model_name].parts[part_name].SectionAssignment(
        offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
        region=Region(faces=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#1 ]',))),
        sectionName='Section-1', thicknessAssignment=FROM_SECTION)

    # --- Assembly ---
    mdb.models[model_name].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models[model_name].rootAssembly.Instance(
        dependent=ON, name=part_name,
        part=mdb.models[model_name].parts[part_name])

    # --- Analysis step ---
    # nlgeom=ON for nonlinear materials (large-deformation formulation)
    # nlgeom=OFF for linear elastic (small-deformation)
    if material_model == 'LE':
        mdb.models[model_name].StaticStep(
            adaptiveDampingRatio=0.05, continueDampingFactors=False,
            initialInc=0.1, name=step_name, nlgeom=OFF,
            previous='Initial', stabilizationMagnitude=0.0002,
            stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
    else:
        mdb.models[model_name].StaticStep(
            adaptiveDampingRatio=0.05, continueDampingFactors=False,
            initialInc=0.1, name=step_name, nlgeom=ON,
            previous='Initial', stabilizationMagnitude=0.0002,
            stabilizationMethod=DISSIPATED_ENERGY_FRACTION)

    # Request energy and nodal force output (used for strain energy extraction)
    mdb.models[model_name].FieldOutputRequest(
        name='Energy_output', createStepName=step_name,
        variables=('ENER', 'ELEN', 'ELEDEN'))
    mdb.models[model_name].FieldOutputRequest(
        name='NF_output', createStepName=step_name, variables=('NFORC',))
    mdb.models[model_name].steps[step_name].setValues(
        maxNumInc=500, initialInc=0.1, maxInc=0.1)

    # --- Boundary conditions: uniform pressure on both hole edges ---
    # mask '[#3 ]' selects both hole boundary edges in the assembly
    mdb.models[model_name].Pressure(
        amplitude=UNSET, createStepName=step_name,
        distributionType=UNIFORM, field='', magnitude=P,
        name='InternalPressure',
        region=Region(side1Edges=mdb.models[model_name].rootAssembly.instances[part_name].edges.getSequenceFromMask(mask=('[#3 ]',))))

    # --- Mesh generation ---

    # Partition the domain into three radial zones to allow graded meshing:
    #   Inner zone (radius ~2.4): fine mesh near the holes
    #   Middle zone (radius ~3.7): transition region
    #   Outer zone: coarse mesh toward the far field
    mdb.models[model_name].ConstrainedSketch(
        gridSpacing=0.7, name='__profile__', sheetSize=28.28,
        transform=mdb.models[model_name].parts[part_name].MakeSketchTransform(
            sketchPlane=mdb.models[model_name].parts[part_name].faces[0],
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models[model_name].parts[part_name].projectReferencesOntoSketch(
        filter=COPLANAR_EDGES, sketch=mdb.models[model_name].sketches['__profile__'])
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(
        center=(0.0, 0.0), point1=(2.4, 0.0))   # inner partition circle
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(
        center=(0.0, 0.0), point1=(3.7, 0.0))   # outer partition circle
    mdb.models[model_name].parts[part_name].PartitionFaceBySketch(
        faces=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#1 ]',)),
        sketch=mdb.models[model_name].sketches['__profile__'])
    del mdb.models[model_name].sketches['__profile__']

    # Assign seed sizes to each region (finer near holes, coarser far away)
    mdb.models[model_name].parts[part_name].seedEdgeBySize(
        constraint=FINER, deviationFactor=0.1,
        edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#18 ]',)),
        size=HoleSeed)      # hole boundary edges (finest)
    mdb.models[model_name].parts[part_name].seedEdgeBySize(
        constraint=FINER, deviationFactor=0.1,
        edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#1 ]',)),
        size=MiddleSeed1)   # inner partition circle
    mdb.models[model_name].parts[part_name].seedEdgeBySize(
        constraint=FINER, deviationFactor=0.1,
        edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#2 ]',)),
        size=MiddleSeed2)   # outer partition circle
    mdb.models[model_name].parts[part_name].seedEdgeBySize(
        constraint=FINER, deviationFactor=0.1,
        edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#4 ]',)),
        size=OuterSeed)     # outer domain boundary (coarsest)

    # Use quad elements with free meshing; CPE4RH = 4-node reduced-integration
    # hybrid (incompatible modes) plane-strain element, suited for
    # incompressible/nearly-incompressible hyperelastic materials
    mdb.models[model_name].parts[part_name].setMeshControls(
        elemShape=QUAD,
        regions=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#7 ]',)),
        technique=FREE)
    mdb.models[model_name].parts[part_name].setElementType(
        regions=(mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#7 ]',)),),
        elemTypes=(ElemType(elemCode=CPE4RH, elemLibrary=STANDARD),))
    mdb.models[model_name].parts[part_name].setMeshControls(
        regions=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#7 ]',)),
        elemShape=QUAD)
    mdb.models[model_name].parts[part_name].generateMesh()

    # --- Identify nodes and elements on the left hole boundary ---
    # These are used to extract deformed coordinates for the ΔA calculation.
    # A bounding sphere is used to locate nodes/elements near the left hole center.
    center_LH          = (center2_x, 0.0, 0.0)
    radius_LH_elements = R + 5.0 * HoleSeed   # slightly larger than R to catch adjacent elements
    radius_LH_nodes    = R + 0.1 * HoleSeed   # very tight to capture only boundary nodes

    assembly  = mdb.models[model_name].rootAssembly
    part      = mdb.models[model_name].parts[part_name]
    instance  = assembly.instances[part_name]

    # Get nodes exactly on the left hole boundary
    hole_nodes = part.nodes.getByBoundingSphere(center=center_LH, radius=radius_LH_nodes)
    part.Set(nodes=hole_nodes, name='H-Nodes')
    node_set = part.sets['H-Nodes']

    # Get elements in the vicinity of the left hole boundary
    elements_hole_edge = part.elements.getByBoundingSphere(
        center=center_LH, radius=radius_LH_elements)
    part.Set(elements=elements_hole_edge, name='H-Elements')
    element_set = part.sets['H-Elements']

    # Keep only elements that share at least one node with the hole boundary node set
    # This filters out elements that are nearby but not actually adjacent to the hole
    filtered_element_ids = []
    for element in element_set.elements:
        element_nodes = element.connectivity
        if any(part.nodes[node] in node_set.nodes for node in element_nodes):
            filtered_element_ids.append(element.label)
    part.SetFromElementLabels(elementLabels=filtered_element_ids, name='Hole')

    # --- Write ABAQUS input file ---
    mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS,
            atTime=None, waitMinutes=0, waitHours=0, queue=None,
            memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
            echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
            userSubroutine='', scratch='', resultsFormat=ODB,
            numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT,
            numCpus=1, numGPUs=0)
    mdb.jobs[job_name].writeInput(consistencyChecking=OFF)


# =============================================================================
# SECTION 2: JOB SUBMISSION
#
# Reads the .inp files written above and submits them to ABAQUS/Standard.
# Jobs are submitted sequentially and the script waits for each to finish
# before proceeding. Adjust numCpus and numDomains for parallel execution.
# =============================================================================

dR_min       = d_min / R
dR_max       = d_max / R + 0.0001
dR_increment = d_increment / R

# Reconstruct the list of dR values (must match Section 1)
dR_values = [dR_min + i * dR_increment
             for i in range(int((dR_max - dR_min) / dR_increment) + 1)]

for dR in dR_values:
    dR_str     = "{:.1f}".format(dR).replace(".", "_")
    job_name   = "P{}_dR{}".format(int(P), dR_str)
    input_file = os.path.join(base_dir, 'P{0}_dR{1}.inp'.format(int(P), dR_str))

    mdb.JobFromInputFile(
        name=job_name, inputFileName=input_file,
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
        userSubroutine='', scratch='', resultsFormat=ODB,
        numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT,
        numCpus=2, numDomains=2, numGPUs=0)
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()  # block until this job is done


# =============================================================================
# SECTION 3: POST-PROCESSING
#
# For each completed simulation:
#   (a) Extract total strain energy history → StrainEnergy.rpt
#   (b) Extract deformed coordinates of left-hole boundary nodes at the
#       final load step → deformed_coordinates-{index}.csv
#
# The deformed coordinates are later used by Deformed_area_calculator_2D.m
# to compute the cavity area change ΔA needed for the potential energy:
#   PE = SE - P * ΔA  (Eq. 9 of paper)
# =============================================================================

counter = 1  # file index for deformed_coordinates CSV files

for dR in dR_values:
    dR_str   = "{:.1f}".format(dR).replace(".", "_")
    odb_path = os.path.join(base_dir, 'P{0}_dR{1}.odb'.format(int(P), dR_str))

    # Open the ODB in the ABAQUS Viewer session
    o3 = session.openOdb(name=odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()
    odb       = session.odbs[odb_path]
    step_name = 'Step-'  + str(counter)
    PART_name = 'PART-'  + str(counter)

    # --- (a) Extract total strain energy time history ---
    # ALLSE is the ABAQUS variable for total strain energy of the whole model.
    # We use the final value (at load factor = 1.0).
    session.XYDataFromHistory(
        name='ALLSE', odb=odb,
        outputVariableName='Strain energy: ALLSE for Whole Model',
        steps=(step_name,), __linkedVpName__='Viewport: 1')
    x0 = session.xyDataObjects['ALLSE']
    session.xyReportOptions.setValues(numDigits=7)
    report_strain_energy = os.path.join(base_dir, 'StrainEnergy.rpt')
    session.writeXYReport(fileName=report_strain_energy, xyData=(x0,))

    # --- (b) Extract deformed nodal coordinates of left hole boundary ---
    # Node IDs are specific to the mesh generated in Section 1.
    # Node 4 is the node at (center2_x, 0) (leftmost point of left hole).
    # Nodes 607–919 form the circumference of the left hole boundary.
    # If the mesh changes, these IDs must be updated accordingly.
    node_ids = [4] + list(range(607, 920))

    # Re-open ODB directly (odbAccess API, distinct from session.odbs)
    odb        = openOdb(odb_path)
    last_step  = odb.steps.keys()[-1]
    last_frame = odb.steps[last_step].frames[-1]  # final load increment

    node_coordinates = {}
    for node_id in node_ids:
        # Get displacement field at this node
        node         = last_frame.fieldOutputs['U'].getSubset(
            region=odb.rootAssembly.instances[PART_name].nodes[node_id - 1])
        displacement = node.values[0].data

        # Get reference (undeformed) coordinates
        initial_coordinates = odb.rootAssembly.instances[PART_name].nodes[node_id - 1].coordinates

        # Compute deformed position: x_deformed = x_ref + u
        displacement = displacement[:2]              # keep only X and Y
        displacement = np.append(displacement, 0.0) # pad Z component (2D problem)
        deformed_coordinates = (initial_coordinates + displacement)[:2]
        node_coordinates[node_id] = deformed_coordinates

    # Write deformed coordinates to CSV for use by Deformed_area_calculator_2D.m
    csv_file = 'deformed_coordinates-' + str(counter) + '.csv'
    with open(csv_file, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['Node ID', 'X Coordinate', 'Y Coordinate'])
        for node_id, coordinates in node_coordinates.items():
            writer.writerow([node_id, coordinates[0], coordinates[1]])

    odb.close()
    counter += 1
