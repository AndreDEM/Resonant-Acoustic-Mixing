"""This package is used to write macros for the design files in StarCCM+. This will be used to determine the dimensionalisation of the problem at hand.""" 
import numpy as np
import os
import pandas as pd
from pandas.tools.plotting import table
import matplotlib.pyplot as plt

def write_header(name,DIR):
    """Write header files for the macro"""
    os.chdir(DIR)
    file = open(name+'.java', 'w')
    file.write('// STAR-CCM+ macro: '+name+'\n')
    file.write('// Written by STAR-CCM+ 11.06.010\n')
    file.write('\n')
    file.write('package macro;\n')
    file.write('import java.util.*;\n')
    file.write('\n')
    file.write('import star.cadmodeler.*;\n')
    file.write('import star.common.*;\n')
    file.write('import star.base.neo.*;\n')
    file.write('import star.meshing.*;\n')
    file.write('import star.lagrangian.*;\n')
    file.write('import star.material.*;\n')
    file.write('import star.lagrangian.dem.*;\n')
    file.write('import star.flow.*;\n')
    file.write('import star.vis.*;\n')
    file.write('import star.motion.*;\n')
    file.write('\n')
    file.write('public class '+name+' extends StarMacro {\n')
    file.write('\n')
    file.write('  public void execute() {\n')
    file.write('    execute0();\n')
    file.write('  }')
    file.write('\n')
    file.close()
    
def volume_set(name,H,R,mD,DIR):
    """This section changes the volume of the current state:
        H: Cup height [mm]
        R: Cup base diameter [mm]
        mD: Minimum mesh base size [m]"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    file.write('  private void execute0() {\n')
#Gets Active simulation to edit
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//-----------------------------EDIT VOLUME--------------------------------------\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('//Gets Active simulation to edit\n')
    file.write('    Simulation simulation_0 =\n')
    file.write('      getActiveSimulation();\n')
#Gets Active CAD model
    file.write('\n//Gets Active CAD model\n')
    file.write('    CadModel cadModel_0 = \n')
    file.write('      ((CadModel) simulation_0.get(SolidModelManager.class).getObject("3D-CAD Model 1"));\n')

#Cup Height change
    file.write('\n//Cup Height change\n')
    file.write('    ScalarQuantityDesignParameter scalarQuantityDesignParameter_1 = \n')
    file.write('      ((ScalarQuantityDesignParameter) cadModel_0.getDesignParameterManager().getObject("Height"));\n')
    file.write('    scalarQuantityDesignParameter_1.getQuantity().setValue('+str(H*1000)+');\n')
    
#Cup Base radius change
    file.write('\n//Cup Base radius change\n')
    file.write('    ScalarQuantityDesignParameter scalarQuantityDesignParameter_0 = \n')
    file.write('      ((ScalarQuantityDesignParameter) cadModel_0.getDesignParameterManager().getObject("Radius"));\n')
    file.write('    scalarQuantityDesignParameter_0.getQuantity().setValue('+str(R*1000)+');\n')
    
#Remesh operation
    file.write('\n//Remesh operation\n')
    file.write('    AutoMeshOperation autoMeshOperation_0 = \n')
    file.write('      ((AutoMeshOperation) simulation_0.get(MeshOperationManager.class).getObject("Automated Mesh"));\n')
    file.write('    autoMeshOperation_0.getDefaultValues().get(BaseSize.class).setValue('+str(mD)+');\n')
    file.write('    MeshPipelineController meshPipelineController_0 = \n')
    file.write('      simulation_0.get(MeshPipelineController.class);\n')
    file.write('    meshPipelineController_0.generateSurfaceMesh();\n')
    file.write('    meshPipelineController_0.generateVolumeMesh();\n')
    file.close()

def add_particles(name,pName,k,rho,E,Ew,DIR):
    """This section changes the volume of the current state:
        name: File/Class name
        pName: Phase Name
        k: Phase counter start at 0
        rho: Particle Density [kg/m^3]
        E: Particle Youngs Modulus [Pa]
        Ew: Wall Youngs Modulus [Pa]"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    if k is 0:
        file.write('\n//------------------------------------------------------------------------------\n')
        file.write('\n//-----------------------------ADD PARTICLE-------------------------------------\n')
        file.write('\n//------------------------------------------------------------------------------\n')
        file.write('//Get current physics model\n')
        file.write('    PhysicsContinuum physicsContinuum_0 = \n')
        file.write('      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));\n')
        
        file.write('\n//Get current Lagrangian Phase model\n')
        file.write('    LagrangianMultiphaseModel lagrangianMultiphaseModel_0 = \n')
        file.write('      physicsContinuum_0.getModelManager().getModel(LagrangianMultiphaseModel.class);\n')
    
#Create the respective Lagrangian Phase
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//Create Lagrangian Phase '+str(k)+'\n')
    file.write('    LagrangianPhase lagrangianPhase_'+str(k)+' = \n')
    file.write('      lagrangianMultiphaseModel_0.createPhase();\n')
    file.write('    lagrangianPhase_'+str(k)+'.setPresentationName("'+pName+'");\n')
    file.write('    lagrangianPhase_'+str(k)+'.enable(DemModel.class);\n')
    file.write('    lagrangianPhase_'+str(k)+'.enable(SphericalParticleModel.class);\n')
    file.write('    lagrangianPhase_'+str(k)+'.enable(SingleComponentParticleModel.class);\n')
    file.write('    lagrangianPhase_'+str(k)+'.enable(ConstantDensityModel.class);\n')
    
#Set Material Properties
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//Material Properties: Particle '+str(k)+'\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('    SingleComponentParticleModel singleComponentParticleModel_'+str(k)+' = \n')
    file.write('      lagrangianPhase_'+str(k)+'.getModelManager().getModel(SingleComponentParticleModel.class);\n')
    file.write('    SingleComponentParticleMaterial singleComponentParticleMaterial_'+str(k)+' = \n')
    file.write('      ((SingleComponentParticleMaterial) singleComponentParticleModel_'+str(k)+'.getMaterial());\n')
    file.write('    singleComponentParticleMaterial_'+str(k)+'.setPresentationName("Material'+str(k)+'");\n')
    
#Particle Density
    file.write('\n//Particle Density\n')
    file.write('    ConstantMaterialPropertyMethod constantMaterialPropertyMethodDens_'+str(k)+' = \n')
    file.write('      ((ConstantMaterialPropertyMethod) singleComponentParticleMaterial_'+str(k)+'.getMaterialProperties().getMaterialProperty(ConstantDensityProperty.class).getMethod());\n')
    file.write('    constantMaterialPropertyMethodDens_'+str(k)+'.getQuantity().setValue('+str(rho)+');\n')
    
#Particle Youngs Modulus
    file.write('\n//Particle Youngs Modulus\n')
    file.write('    ConstantMaterialPropertyMethod constantMaterialPropertyMethodE_'+str(k)+' = \n')
    file.write('      ((ConstantMaterialPropertyMethod) singleComponentParticleMaterial_'+str(k)+'.getMaterialProperties().getMaterialProperty(YoungsModulusProperty.class).getMethod());\n')
    file.write('    constantMaterialPropertyMethodE_'+str(k)+'.getQuantity().setValue('+str(E)+');\n')
    
#Wall
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//Material Properties: Particle '+str(k)+' wall\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('    DemMultiphaseModel demMultiphaseModel_'+str(k)+' = \n')
    file.write('      physicsContinuum_0.getModelManager().getModel(DemMultiphaseModel.class);\n')
    file.write('    DemBoundaryPhase demBoundaryPhase_'+str(k)+' = \n')
    file.write('      ((DemBoundaryPhase) demMultiphaseModel_'+str(k)+'.getDemBoundaryPhaseManager().getPhase("Wall"));\n')
    file.write('\n')
    file.write('    SolidModel solidModel_'+str(k)+' = \n')
    file.write('      demBoundaryPhase_'+str(k)+'.getModelManager().getModel(SolidModel.class);\n')
    file.write('    Solid solid_'+str(k)+' = \n')
    file.write('      ((Solid) solidModel_'+str(k)+'.getMaterial());\n')
    file.write('    solid_'+str(k)+'.setPresentationName("wall");\n')
    
    file.write('\n//Wall Youngs Modulus\n')
    file.write('    ConstantMaterialPropertyMethod constantMaterialPropertyMethod_'+str(k)+' = \n')
    file.write('      ((ConstantMaterialPropertyMethod) solid_'+str(k)+'.getMaterialProperties().getMaterialProperty(YoungsModulusProperty.class).getMethod());\n')
    file.write('    constantMaterialPropertyMethod_'+str(k)+'.getQuantity().setValue('+str(Ew)+');\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.close()
    
def particle_interaction(name,I1,I2,mu,Rn,Rt,W,k,DIR):
    """This section changes the volume of the current state:
    name: File/Class name
    I1: Phase Name 1 [same as particle name, if wall set I2 to "Wall"]
    I1: Phase Name 2 [same as particle name, if wall: "Wall"]
    mu: Static Friction Coef.
    Rn: Restitution Coef. (Normal)
    Rt: Restitution Coef. (Tangential)
    W: Linear Cohesion Coef. (Cohesion: Similar particles, Adhesion: Dissimilar particles)
    k: Phase counter start at 0"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    if k is 0:
        file.write('\n//------------------------------------------------------------------------------\n')
        file.write('\n//-------------------------PARTICLE INTERACTION---------------------------------\n')
        file.write('\n//------------------------------------------------------------------------------\n')
        file.write('//Get Multiphase Interaction Tree\n')
        file.write('    MultiPhaseInteractionModel multiPhaseInteractionModel_0 = \n')
        file.write('      physicsContinuum_0.getModelManager().getModel(MultiPhaseInteractionModel.class);\n')
    
    file.write('\n//------------------------------------------------------------------------------\n')
#Get Phase Interaction
    file.write('//Get Phase Interaction '+str(k)+'\n')
    file.write('    PhaseInteraction phaseInteraction_'+str(k)+' = \n')
    file.write('      multiPhaseInteractionModel_0.createPhaseInteraction();\n')
    file.write('    phaseInteraction_'+str(k)+'.setPresentationName("'+I1+I2+'");\n')
    file.write('\n')
    
    file.write('\n//Set Hertz Mendlin Model\n')
    file.write('    phaseInteraction_'+str(k)+'.enable(DemPhaseInteractionModel.class);\n')
    file.write('    phaseInteraction_'+str(k)+'.enable(HertzMindlinNoSlipContactModel.class);\n')
    file.write('    phaseInteraction_'+str(k)+'.enable(LinearCohesionModel.class);\n')
    file.write('    DemPhaseInteractionModel demPhaseInteractionModel_'+str(k)+' = \n')
    file.write('      phaseInteraction_'+str(k)+'.getModelManager().getModel(DemPhaseInteractionModel.class);\n')
    file.write('\n')

    file.write('\n//Link the respective interactions\n')
#Phase 1
    file.write('    LagrangianPhase lagrangianPhase_1_'+str(k)+' = \n')
    file.write('      ((LagrangianPhase) lagrangianMultiphaseModel_0.getPhaseManager().getPhase("'+I1+'"));\n')
    file.write('    demPhaseInteractionModel_'+str(k)+'.setPhase0(lagrangianPhase_1_'+str(k)+');\n')
    
    if I2 is 'Wall':
    #Wall
        file.write('   demPhaseInteractionModel_'+str(k)+'.setPhase1(demBoundaryPhase_0);\n')    
    else:
    #Phase 2
        file.write('    LagrangianPhase lagrangianPhase_2_'+str(k)+' = \n')
        file.write('      ((LagrangianPhase) lagrangianMultiphaseModel_0.getPhaseManager().getPhase("'+I2+'"));\n')
        file.write('    demPhaseInteractionModel_'+str(k)+'.setPhase1(lagrangianPhase_2_'+str(k)+');\n')
        
#HM
    file.write('\n//Set HM Model parameters\n')
    file.write('    HertzMindlinNoSlipContactModel hertzMindlinNoSlipContactModel_'+str(k)+' = \n')
    file.write('      phaseInteraction_'+str(k)+'.getModelManager().getModel(HertzMindlinNoSlipContactModel.class);\n')
    
#Friction
    file.write('\n//Static Friction Coef.\n')
    file.write('    hertzMindlinNoSlipContactModel_'+str(k)+'.getCoeffStaticFriction().setValue('+str(mu)+');\n')
#N Rest. Coef
    file.write('\n//Normal Restitution Coef.\n')
    file.write('    hertzMindlinNoSlipContactModel_'+str(k)+'.getCoeffRestitutionNormal().setValue('+str(Rn)+');\n')
#T Rest. Coef
    file.write('\n//Tangential Restitution Coef.\n')
    file.write('    hertzMindlinNoSlipContactModel_'+str(k)+'.getCoeffRestitutionTangential().setValue('+str(Rt)+');\n')
    
#Set Linear Cohesion Model
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//Set Linear Cohesion Model\n')
    file.write('    LinearCohesionModel linearCohesionModel_'+str(k)+' = \n')
    file.write('      phaseInteraction_'+str(k)+'.getModelManager().getModel(LinearCohesionModel.class);\n')
    
    #Cohesion Energy
    file.write('\n//Cohesion Energy\n')
    file.write('    linearCohesionModel_'+str(k)+'.getWorkOfCohesion().setValue('+str(W)+');\n')
    file.close()

def injector(name,H,R,pD,I,k,DIR):
    """This section changes the volume of the current state:
    name: File/Class name
    H: Cup Height [m]
    R: Cup Base radius [m]
    pD: Particle Diameter
    I: Particle Name
    k: Phase counter start at 0"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    file.write('    Injector injector_'+str(k)+' = \n')
    file.write('      simulation_0.get(InjectorManager.class).createInjector();\n')
    if k is 0:
        file.write('\n//------------------------------------------------------------------------------\n')
        file.write('\n//-----------------------------ADD INJECTOR-------------------------------------\n')
        file.write('\n//------------------------------------------------------------------------------\n')
        file.write('\n//Set Probe plane\n')
        file.write('    PlaneProbePart planeProbePart_0 = \n')
        file.write('      ((PlaneProbePart) simulation_0.getPartManager().getObject("Presentation Grid"));\n')
        file.write('    Coordinate coordinate_origin = \n')
        file.write('      planeProbePart_0.getOriginCoordinate();\n')
        file.write('    Units units_0 = \n')
        file.write('      ((Units) simulation_0.getUnitsManager().getObject("m"));\n')
        
        #Plane Origin 1
        file.write('\n//Set Plane origin coordinates\n')
        file.write('    coordinate_origin.setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {'+str(-R)+', '+str(-R)+', '+str(H-0.005)+'}));\n')
        
        #Plane Point 1
        file.write('\n//Set Plane point 1\n')
        file.write('    Coordinate coordinate_point1 = \n')
        file.write('      planeProbePart_0.getPoint1Coordinate();\n')
        file.write('    coordinate_point1.setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {'+str(R)+', '+str(-R)+', '+str(H-0.005)+'}));\n')
        
        #Plane Point 2
        file.write('\n//Set Plane point 2\n')
        file.write('    Coordinate coordinate_point2 = \n')
        file.write('      planeProbePart_0.getPoint2Coordinate();\n')
        file.write('    coordinate_point2.setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {'+str(-R)+', '+str(R)+', '+str(H-0.005)+'}));\n')
        
        
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//Injector '+str(k)+'\n')
    file.write('    injector_'+str(k)+'.setPresentationName("'+I+'");\n')
    file.write('    injector_'+str(k)+'.setLagrangianPhase(lagrangianPhase_'+str(k)+');\n')
    file.write('    PartInjector partInjector_'+str(k)+' = \n')
    file.write('      ((PartInjector) simulation_0.get(ConditionTypeManager.class).get(PartInjector.class));\n')
    file.write('    injector_'+str(k)+'.setInjectorType(partInjector_'+str(k)+');\n')
    file.write('    injector_'+str(k)+'.getPartGroup().setObjects(planeProbePart_0);\n')
    
    file.write('\n//Flow Rate specification\n')
    file.write('    injector_'+str(k)+'.getInjectorConditions().get(InjectorFlowRateOption.class).setSelected(InjectorFlowRateOption.Type.PARTICLE_FLOW_RATE);\n')
    file.write('    InjectorParticleFlowRate injectorParticleFlowRate_'+str(k)+' = \n')
    file.write('      injector_'+str(k)+'.getInjectorValues().get(InjectorParticleFlowRate.class);')
    file.write('    injectorParticleFlowRate_'+str(k)+'.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${inject_'+str(k)+'}");\n')
    
    file.write('\n//Particle Diameter\n')
    file.write('    InjectorParticleDiameter injectorParticleDiameter_'+str(k)+' = \n')
    file.write('      injector_'+str(k)+'.getInjectorValues().get(InjectorParticleDiameter.class);\n')
    file.write('    injectorParticleDiameter_'+str(k)+'.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue('+str(pD)+');\n')
    file.close()
    
def flowrate(name,FR,k,DIR):
    """This section changes the flowrate of the current state:
    name: File/Class name
    FR: Flowrate particles per injector per second so (num of particles needed/4) = TOTAL PARTICLES ADDED as each injection timestep is set to 0.25s
    k: Phase counter start at 0"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//-------------------------------FLOW RATE--------------------------------------\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//Define injector '+str(k)+'\n')
    file.write('    UserFieldFunction userFieldFunctionFR_'+str(k)+' = \n')
    file.write('      ((UserFieldFunction) simulation_0.getFieldFunctionManager().createFieldFunction());\n')
    file.write('    userFieldFunctionFR_'+str(k)+'.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);\n')
    file.write('    userFieldFunctionFR_'+str(k)+'.setPresentationName("inject_'+str(k)+'");\n')
    file.write('    userFieldFunctionFR_'+str(k)+'.setFunctionName("inject_'+str(k)+'");\n')
    file.write('\n')
    
    if k is 0:
        file.write('    userFieldFunctionFR_'+str(k)+'.setDefinition("${Time} > 0.25? 0:'+str(FR)+'");\n')

    else:
        file.write('    userFieldFunctionFR_'+str(k)+'.setDefinition("${Time} > '+str(k*0.25)+'? (${Time} < '+str((k+1)*0.25)+'? '+str(FR)+':0):0");\n')
    
    file.close()
    
def speed(name,A,f,start,TS,DIR):
    """This section changes the flowrate of the current state:
    name: File/Class name
    A: Cup amplitude [m]
    f: Cup frequency [Hz]
    start: Excitation Start Time [s]
    TS: Total Solution Time [s] (Settling Time to be assumed at 1 s)
    k: Phase counter start at 0"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')

    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//---------------------------------SPEED----------------------------------------\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('    UserFieldFunction userFieldFunction_S = \n')
    file.write('      ((UserFieldFunction) simulation_0.getFieldFunctionManager().createFieldFunction());\n')
    file.write('    userFieldFunction_S.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);\n')
    file.write('    userFieldFunction_S.setPresentationName("shake");\n')
    file.write('    userFieldFunction_S.setFunctionName("shake");\n')
    file.write('\n')
    
    w = 2*np.pi*f
    
    file.write('    userFieldFunction_S.setDefinition("${Time} > '+str(TS+1)+'? 0 :(${Time} < '+str(start)+'? 0 : ('+str(A)+'*'+str(w)+'*sin((${Time}-'+str(start)+')*'+str(w)+')))");\n')
    
    #Add Motion
    file.write('    TranslatingMotion translatingMotion_0 = \n')
    file.write('      simulation_0.get(MotionManager.class).createMotion(TranslatingMotion.class, "Translation");\n')
    file.write('    translatingMotion_0.getTranslationVelocity().setDefinition("[0.0, 0.0, ${shake}]");\n')
    file.write('    Region region_0 = \n')
    file.write('      simulation_0.getRegionManager().getRegion("inner");\n')
    file.write('    MotionSpecification motionSpecification_0 = \n')
    file.write('      region_0.getValues().get(MotionSpecification.class);\n')
    file.write('    motionSpecification_0.setMotion(translatingMotion_0);\n')
    file.write('\n')
    
    file.write('//Edit Inner Iterations (FREEZE FLOW) [iter = 1]\n')
    file.write('    InnerIterationStoppingCriterion innerIterationStoppingCriterion_0 = \n')
    file.write('      ((InnerIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Inner Iterations"));\n')
    file.write('    innerIterationStoppingCriterion_0.setMaximumNumberInnerIterations(1);\n')
    file.write('\n')
    
    file.write('//Edit total solution time\n')
    file.write('    PhysicalTimeStoppingCriterion physicalTimeStoppingCriterion_0 = \n')
    file.write('      ((PhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));\n')
    file.write('    physicalTimeStoppingCriterion_0.getMaximumTime().setValue('+str(TS+1)+');\n')
    
    file.write('//Stopping Criteria\n')
    file.write('    StepStoppingCriterion stepStoppingCriterion_0 = \n')
    file.write('      ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));\n')
    file.write('    stepStoppingCriterion_0.setIsUsed(false);\n')
    file.close()

def add_table(name,I,start,k,DIR):
    """This section adds a export table to be used in post processing:
    name: File/Class name
    sDIR: Simulation directory [str]
    I: Iteration export frequency
    start: Excitation Start Time [s]
    k: Phase range (only max value now)"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//-------------------------------ADD TABLE--------------------------------------\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    
    file.write('    XyzInternalTable xyzInternalTable_0 = \n')
    file.write('      simulation_0.getTableManager().createTable(XyzInternalTable.class);')
    
    file.write('\n//Export PARCEL ID\n')
    file.write('    PrimitiveFieldFunction primitiveFieldFunction_0 = \n')
    file.write('      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("ParcelId"));\n')
    
    file.write('\n//Export PARCEL MASS\n')
    file.write('    PrimitiveFieldFunction primitiveFieldFunction_1 = \n')
    file.write('      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("ParcelMass"));\n')
    
    file.write('\n//Export INJECTOR INDEX\n')
    file.write('    PrimitiveFieldFunction primitiveFieldFunction_2 = \n')
    file.write('      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("InjectorIndex"));\n')
    
    file.write('\n//Export PARTICLE VELOCITY\n')
    file.write('    PrimitiveFieldFunction primitiveFieldFunction_3 = \n')
    file.write('      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("ParticleVelocity"));\n')
    file.write('    VectorComponentFieldFunction vectorComponentFieldFunction_0 = \n')
    file.write('      ((VectorComponentFieldFunction) primitiveFieldFunction_3.getComponentFunction(0));\n')
    file.write('    VectorComponentFieldFunction vectorComponentFieldFunction_1 = \n')
    file.write('      ((VectorComponentFieldFunction) primitiveFieldFunction_3.getComponentFunction(1));\n')
    file.write('    VectorComponentFieldFunction vectorComponentFieldFunction_2 = \n')
    file.write('      ((VectorComponentFieldFunction) primitiveFieldFunction_3.getComponentFunction(2));\n')
    
    file.write('\n//Export TIME\n')
    file.write('    PrimitiveFieldFunction primitiveFieldFunction_4 = \n')
    file.write('      ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("Time"));\n')
    
    file.write('//Write Table\n')
    file.write('    xyzInternalTable_0.setFieldFunctions(new NeoObjectVector(new Object[] {primitiveFieldFunction_0, primitiveFieldFunction_2, vectorComponentFieldFunction_0, vectorComponentFieldFunction_1, vectorComponentFieldFunction_2, primitiveFieldFunction_1,primitiveFieldFunction_4}));\n')    
        
    file.write('\n//Table Update\n')
    file.write('    TableUpdate tableUpdate_0 = \n')
    file.write('      xyzInternalTable_0.getTableUpdate();\n')
    
    file.write('    tableUpdate_0.setSaveToFile(true);\n')
    
    #Check wheter simulation directory exists and if not create
    if not os.path.exists(DIR+'/Table'):
        #Windows
        #os.makedirs(DIR+'\\Table')
        #Linux or Mac
        os.makedirs(DIR+'/Table')
    
    file.write('    tableUpdate_0.setFilePath("'+DIR+'/Table");\n')
    file.write('    tableUpdate_0.getUpdateModeOption().setSelected(StarUpdateModeOption.Type.ITERATION);\n')
    file.write('    IterationUpdateFrequency iterationUpdateFrequency_0 = \n')
    file.write('      tableUpdate_0.getIterationUpdateFrequency();\n')
    
    file.write('    iterationUpdateFrequency_0.setIterations('+str(I)+');\n')
    file.write('    iterationUpdateFrequency_0.setStart('+str(int(np.round((start/0.00025)-1, decimals = 0)))+');\n')

    S = 'lagrangianPhase_0';
    if k > 0:
        for s in range(0,k):
                S = S+', lagrangianPhase_'+str(s+1)
    
    file.write('    xyzInternalTable_0.getParts().setObjects('+S+');\n')
    file.close()
    
def save_and_run(name,DIR):
    """This section saves the simulation file to the simulation directory and runs the simulation:
    name: File/Class name
    DIR: Simulation directory [str]"""
    
    os.chdir(DIR)
    file = open(name+'.java', 'a+')
    file.write('\n//------------------------------------------------------------------------------\n')
    file.write('\n//-----------------------------SAVE AND RUN-------------------------------------\n')
    file.write('\n//------------------------------------------------------------------------------\n')
    
    file.write('//Save File as: '+name+'.sim\n')
    file.write('    simulation_0.saveState(resolvePath("'+DIR+'/'+name+'.sim"));\n')
    
    file.write('\n//Get NEW simulation and RUN\n')
    file.write('    Simulation simulation_1 =\n')
    file.write('      getActiveSimulation();\n')
    file.write('    simulation_1.getSimulationIterator().run();\n')
    file.write('  }\n')
    file.write('}')
    file.close()