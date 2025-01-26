# Import necessary Abaqus modules
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

# Change directory to the desired path
base_dir = r"D:\ElasticInteraction\files"
os.chdir(base_dir)


############################## 1. Creating the conputational model and write input files ##############################

# Initialize parameters
P = 30 # pressure of holes
d_min = 0.22 # starting value of center-to-center holes' distance
d_max = 2.0 + 0.0001 # ending value of center-to-center holes' distance, ensuring the upper bound is inclusive
d_increment = 0.02 # distance increment
R_m, R = 5.0, 0.1 # domain radius, hole radius
material_model = 'NH' # Define the material model to be used # Change this to 'LE', 'MR', or 'AB' for other models
HoleSeed, MiddleSeed1, MiddleSeed2, OuterSeed = 0.002, 0.06, 0.12, 0.2 # mesh seeding sizes

# Define material properties based on the selected model
if material_model == 'NH':
    C10, D1 = 10.0, 0.0
    material_type = NEO_HOOKE
    table = ((C10, D1),)
elif material_model == 'LE':
    E, nu = 52.0, 0.3
    material_type = ELASTIC
    table = ((E, nu),)
elif material_model == 'MR':
    C10, C01, D1 = 5.0, 5.0, 0.0
    material_type = MOONEY_RIVLIN
    table = ((C10, C01, D1),)
elif material_model == 'AB':
    C2, lambda_m, D = 16.7091, 2.0, 0.0
    material_type = ARRUDA_BOYCE
    table = ((C2, lambda_m, D),)
else:
    raise ValueError("Unknown material model specified.")

# Generate circle values as (d/2, 0.0) pairs
circle_values = [(d_min / 2 + d_increment / 2 * i, 0.0) for i in range(int((d_max / 2 - d_min / 2) / (d_increment / 2) + 1))]

# Loop over circle values to create models
for index, circle in enumerate(circle_values):
    model_name = 'Model-' + str(index + 1)
    part_name = 'Part-' + str(index + 1)
    step_name = 'Step-' + str(index + 1)
    # Calculate dR value
    dR = 2 * circle[0] / R
    job_name = "P{0}_dR{1:.1f}".format(int(P), dR)  # Format job_name as "P10_dR3.0"
    # Ensure job_name does not contain invalid characters
    job_name = job_name.replace('.', '_')  # Replace dot with 'p'
    
    # Model creation
    mdb.Model(name=model_name, modelType=STANDARD_EXPLICIT)

    # Part creation
    mdb.models[model_name].ConstrainedSketch(name='__profile__', sheetSize=20.0)
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(-R_m, 0.0))
    center1_x = circle[0]
    center2_x = -center1_x
    point1_x = center1_x + R
    point2_x = center2_x - R
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(center=(center1_x, 0.0), point1=(point1_x, 0.0))
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(center=(center2_x, 0.0), point1=(point2_x, 0.0))
    mdb.models[model_name].Part(dimensionality=TWO_D_PLANAR, name=part_name, type=DEFORMABLE_BODY)
    mdb.models[model_name].parts[part_name].BaseShell(sketch=mdb.models[model_name].sketches['__profile__'])
    del mdb.models[model_name].sketches['__profile__']

    # Material and section properties
    mdb.models[model_name].Material(name='Material-1')
    if material_model == 'NH':
        mdb.models[model_name].materials['Material-1'].Hyperelastic(materialType=ISOTROPIC, table=table, testData=OFF, type=material_type, volumetricResponse=VOLUMETRIC_DATA)
    elif material_model == 'LE':
        mdb.models[model_name].materials['Material-1'].Elastic(table=table)
    elif material_model == 'MR':
        mdb.models[model_name].materials['Material-1'].Hyperelastic(materialType=ISOTROPIC, table=table, testData=OFF, type=material_type, volumetricResponse=VOLUMETRIC_DATA)
    elif material_model == 'AB':
        mdb.models[model_name].materials['Material-1'].Hyperelastic(materialType=ISOTROPIC, table=table, testData=OFF, type=material_type, volumetricResponse=VOLUMETRIC_DATA)

    mdb.models[model_name].HomogeneousSolidSection(material='Material-1', name='Section-1', thickness=1.0)
    mdb.models[model_name].parts[part_name].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                                                               region=Region(faces=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#1 ]',))),
                                                               sectionName='Section-1', thicknessAssignment=FROM_SECTION)
    
    # Assembly
    mdb.models[model_name].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models[model_name].rootAssembly.Instance(dependent=ON, name=part_name,
                                                 part=mdb.models[model_name].parts[part_name])

    # Step
    if material_model == 'LE':
        mdb.models[model_name].StaticStep(adaptiveDampingRatio=0.05, continueDampingFactors=False, initialInc=0.1,name=step_name, nlgeom=OFF,
                                      previous='Initial', stabilizationMagnitude=0.0002, stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
    else:
        mdb.models[model_name].StaticStep(adaptiveDampingRatio=0.05, continueDampingFactors=False, initialInc=0.1, name=step_name, nlgeom=ON,
                                      previous='Initial', stabilizationMagnitude=0.0002, stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
    mdb.models[model_name].FieldOutputRequest(name='Energy_output', createStepName=step_name, variables=('ENER', 'ELEN', 'ELEDEN'))
    mdb.models[model_name].FieldOutputRequest(name='NF_output', createStepName=step_name, variables=('NFORC', ))
    mdb.models[model_name].steps[step_name].setValues(maxNumInc=500, initialInc=0.1, maxInc=0.1)

    # Loading
    mdb.models[model_name].Pressure(amplitude=UNSET, createStepName=step_name, distributionType=UNIFORM, field='', magnitude=P,
                                    name='InternalPressure', region=Region(side1Edges=mdb.models[model_name].rootAssembly.instances[part_name].edges.getSequenceFromMask(mask=('[#3 ]',))))
    
    # Mesh
    # Partitioning
    mdb.models[model_name].ConstrainedSketch(gridSpacing=0.7, name='__profile__', sheetSize=28.28,
                                             transform=mdb.models[model_name].parts[part_name].MakeSketchTransform(sketchPlane=mdb.models[model_name].parts[part_name].faces[0],
                                                                                                                   sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models[model_name].parts[part_name].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models[model_name].sketches['__profile__'])
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(2.4, 0.0))
    mdb.models[model_name].sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(3.7, 0.0))
    mdb.models[model_name].parts[part_name].PartitionFaceBySketch(faces=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#1 ]',)),
                                                                  sketch=mdb.models[model_name].sketches['__profile__'])
    del mdb.models[model_name].sketches['__profile__']
    
    # Seeding and generating mesh
    mdb.models[model_name].parts[part_name].seedEdgeBySize(constraint=FINER, deviationFactor=0.1,
                                                           edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#18 ]',)),
                                                           size=HoleSeed)
    mdb.models[model_name].parts[part_name].seedEdgeBySize(constraint=FINER, deviationFactor=0.1,
                                                           edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#1 ]',)),
                                                           size=MiddleSeed1)
    mdb.models[model_name].parts[part_name].seedEdgeBySize(constraint=FINER, deviationFactor=0.1,
                                                           edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#2 ]',)),
                                                           size=MiddleSeed2)
    mdb.models[model_name].parts[part_name].seedEdgeBySize(constraint=FINER, deviationFactor=0.1,
                                                           edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(mask=('[#4 ]',)),
                                                           size=OuterSeed)
    mdb.models[model_name].parts[part_name].setMeshControls(elemShape=QUAD, regions=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#7 ]',)), technique=FREE)
    mdb.models[model_name].parts[part_name].setElementType(regions=(mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#7 ]',)),),
        elemTypes=(ElemType(elemCode=CPE4RH, elemLibrary=STANDARD),))
    mdb.models[model_name].parts[part_name].setMeshControls(regions=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(mask=('[#7 ]',)),elemShape=QUAD)
    mdb.models[model_name].parts[part_name].generateMesh()

    # Define the center and radius of the bounding sphere for the hole
    center_LH = (center2_x, 0.0, 0.0)
    radius_LH_elements = R + 5.0 * HoleSeed
    radius_LH_nodes = R + 0.1 * HoleSeed
    
    # Get elements within the bounding sphere from the part
    assembly = mdb.models[model_name].rootAssembly
    part = mdb.models[model_name].parts[part_name]
    instance = assembly.instances[part_name]
    hole_nodes = part.nodes.getByBoundingSphere(center=center_LH, radius=radius_LH_nodes)
    part.Set(nodes=hole_nodes, name='H-Nodes')
    node_set = part.sets['H-Nodes']
    elements_hole_edge = part.elements.getByBoundingSphere(center=center_LH, radius=radius_LH_elements)
    part.Set(elements=elements_hole_edge, name='H-Elements')
    element_set = part.sets['H-Elements']

    # Initialize an empty list to store element IDs to keep
    filtered_element_ids = []
    # Iterate over each element in the element set
    for element in element_set.elements:
        # Get the nodes of the element
        element_nodes = element.connectivity
        # Check if any of the element's nodes are in the node set
        if any(part.nodes[node] in node_set.nodes for node in element_nodes):
            filtered_element_ids.append(element.label)
    # Create a new element set with the filtered element IDs
    part.SetFromElementLabels(elementLabels=filtered_element_ids, name='Hole')
    
    # Create the Job and write the input file
    mdb.Job(name=job_name, model=model_name, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
            queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
            userSubroutine='', scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    mdb.jobs[job_name].writeInput(consistencyChecking=OFF)




############################################## 2. Submitting input files ##############################################

dR_min = d_min/R  # starting value of d/R
dR_max = d_max/R + 0.0001  # ending value of d/R, ensuring the upper bound is inclusive
dR_increment = d_increment/R  # interval

# Loop through the range of dR values
dR_values = [dR_min + i * dR_increment for i in range(int((dR_max - dR_min) / dR_increment) + 1)]

for dR in dR_values:
    dR_str = "{:.1f}".format(dR).replace(".", "_")  # Replace period with hyphen for dR string
    job_name = "P{}_dR{}".format(int(P), dR_str)  # Format job name
    input_file = os.path.join(base_dir, 'P{0}_dR{1}.inp'.format(int(P), dR_str))

    mdb.JobFromInputFile(name=job_name,
                         inputFileName=input_file,
                         type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None,
                         memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                         explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
                         userSubroutine='', scratch='', resultsFormat=ODB,
                         numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=2,
                         numDomains=2, numGPUs=0)
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()




############## 3. Extracting Strain Energy of the whole system, and Deformed Coordinates of the Left Hole ##############

# Loop through the range of dR values
counter = 1

for dR in dR_values:
    dR_str = "{:.1f}".format(dR).replace(".", "_")  # Replace period with hyphen for dR string

    # Reading odb files
    odb_path = os.path.join(base_dir, 'P{0}_dR{1}.odb'.format(int(P), dR_str))
    o3 = session.openOdb(name=odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()
    odb = session.odbs[odb_path]
    step_name = 'Step-' + str(counter)
    PART_name = 'PART-' + str(counter)

    # Strain Energy Extraction
    session.XYDataFromHistory(name='ALLSE', odb=odb,
                              outputVariableName='Strain energy: ALLSE for Whole Model', steps=(step_name,),__linkedVpName__='Viewport: 1')
    x0 = session.xyDataObjects['ALLSE']
    session.xyReportOptions.setValues(numDigits=7)
    report_strain_energy = os.path.join(base_dir, 'StrainEnergy.rpt')
    session.writeXYReport(fileName=report_strain_energy, xyData=(x0,))

    # Deformed Coordinates Extraction
    # Node IDs for which you want to extract deformed coordinates
    node_ids = [4] + list(range(607, 920))  # The range(start, stop) function generates numbers from start to stop-1
    # Open the output database
    odb = openOdb(odb_path)
    # Get the last frame (i.e., the final step) in the output database
    last_step = odb.steps.keys()[-1]
    last_frame = odb.steps[last_step].frames[-1]
    # Get the deformed coordinates for the specified nodes
    node_coordinates = {}
    for node_id in node_ids:
        node = last_frame.fieldOutputs['U'].getSubset(region=odb.rootAssembly.instances[PART_name].nodes[node_id - 1])
        displacement = node.values[0].data
        initial_coordinates = odb.rootAssembly.instances[PART_name].nodes[node_id - 1].coordinates
        # Reshape the displacement array to match the shape of the initial coordinates array
        displacement = displacement[:2]  # Extract only X and Y components for 2D models
        displacement = np.append(displacement, 0.0)  # Add a zero value for the Z component
        # Calculate deformed coordinates by adding displacement to initial coordinates
        deformed_coordinates = initial_coordinates + displacement
        deformed_coordinates = deformed_coordinates[:2]  # Extract only X and Y coordinates for 2D models
        node_coordinates[node_id] = deformed_coordinates
    # Write the deformed coordinates to a CSV file
    csv_file = 'deformed_coordinates-' + str(counter) + '.csv'  # Update with your desired CSV file name
    with open(csv_file, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['Node ID', 'X Coordinate', 'Y Coordinate'])

        for node_id, coordinates in node_coordinates.items():
            writer.writerow([node_id, coordinates[0], coordinates[1]])
    # Close the output database
    odb.close()

    counter += 1  # Increment the counter for the next iteration








