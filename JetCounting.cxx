/**
	JetCounting
	JetCounting.cxx 
	Purpose: Defines the entry point for the application. 

	@author Gordon Stevenson
	@version 0.1 16/01/2017
*/
#pragma once

#include <stdio.h>
#include <vtkActor.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkImageData.h>
#include <vtkImageThreshold.h>
#include <vtkMarchingCubes.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkProbeFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkTriangleFilter.h>

#include "JetCountingMeshToolbox.h"

using namespace std;

int main(int argc, char* argv[])
	/**
		Runs the JetCounting main Application. Using 4 images that are pre-processed the application will extract the poly data meshes at defined distance in millimetrs from the utero-placenta interface
		and then output the 3D image files and meshes to an output directory.

		Choice of 3D axis to extract and orientation can be user defined by the optional arguments.

		**this program assumes all volumes are in the same resolution and spacing, in MetaImage format. NB Resampling may be required. **

		@param inputBModeVolume path to the 3D bmode volume
		@param inputPowerDopplerVolume path to the 3D PD volume
		@param inputSegVolume path to the 3D binary segmentation volume
		@param inputDTVolume path to the 3D Signed Euclidean Distance Transform volume
		@param dist floating point argument of distance from the UPI to measure
		@param outputDir output directory to store volumes and meshes for analysis

		@param axis optional integer parameter of axis (0,1,2) in which extraction should occur
		@param interface optional integer parameter of negative or positive normal (0,1) extraction of the UPI from the total mesh so occur.

	*/
	{
	if ((argc != 7) && (argc!=8) && (argc!=9)){
		std::cerr << "Usage: ./" << argv[0] << " inputBModeVolume inputPowerDopplerVolume inputSegVolume inputDTVolume dist outputDir [axis] [interface]" << std::endl;			
		return EXIT_FAILURE;
	}
	
	//load in data
	vtkSmartPointer<vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
	vtkSmartPointer<vtkMetaImageReader> pdreader = vtkSmartPointer<vtkMetaImageReader>::New();
	vtkSmartPointer<vtkMetaImageReader> dtreader = vtkSmartPointer<vtkMetaImageReader>::New();
	vtkSmartPointer<vtkMetaImageReader> segreader = vtkSmartPointer<vtkMetaImageReader>::New();
	vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();

	//convert to mesh format.
	vtkSmartPointer<vtkImageData> batchinput = vtkSmartPointer<vtkImageData>::New();
	vtkSmartPointer<vtkImageData> seeding    = vtkSmartPointer<vtkImageData>::New();
	vtkSmartPointer<vtkImageData> seg = vtkSmartPointer<vtkImageData>::New();
	vtkSmartPointer<vtkImageData> result     = vtkSmartPointer<vtkImageData>::New();
		
	cout << "Greyscale File: " << argv[1] <<"\n";
	cout << "Doppler File: " << argv[2] <<"\n";
	cout << "Seg File : " << argv[3] << "\n";
	cout << "EDT File : " << argv[4] << "\n";
	cout << "Distance (mm) : " << argv[5] << "\n";

	reader->SetFileName(argv[1]);
	reader->Update();

        //pdreader->SetFileName(argv[2]);
	//pdreader->Update();

	segreader->SetFileName(argv[3]);
	segreader->Update();

	dtreader->SetFileName(argv[4]);
	dtreader->Update();

	std::string out_path = argv[6];
	std::string pd_path = out_path;
	std::string edt_path = out_path;
	std::string seg_path = out_path;
	std::string mesh_path = out_path;
	std::string mesh_path1 = out_path;
	std::string full_mesh_path = out_path;
	std::string bmode_path = out_path;

	pd_path.append("/rot_pd.mha");
	edt_path.append("/rot_edt.mha");
	seg_path.append("/rot_seg.mha");
	bmode_path.append("/rot_bmode.mha");
	mesh_path.append("/mesh.vtk");
	mesh_path1.append("/mesh_otherside.vtk");
	full_mesh_path.append("/mesh_nocut.vtk");


	JetCountingMeshToolbox *toolbox = new JetCountingMeshToolbox();
	
	toolbox->VolumeData = reader->GetOutput();
	cout << "Test1" << "\n";
	toolbox->SegmentationData = segreader->GetOutput();
	cout << "Test2" << "\n";
	//toolbox->DopplerData = pdreader->GetOutput();
	toolbox->EDTData = dtreader->GetOutput();

	if ((argc == 8) && (atoi(argv[7]) > 2) && (atoi(argv[7])<0)){				
		{toolbox->RotatePlacentaBasedOnMinorAxis(atoi(argv[7]));}
	}
	else{
		if ((argc == 9)  && (atoi(argv[8]) <= 2) && (atoi(argv[8])>=0)){
			if((atoi(argv[7]) <= 2) && (atoi(argv[7])>=0)){
				toolbox->RotatePlacentaBasedOnMinorAxis(atoi(argv[7]),atoi(argv[8]));
			}
			else{
				toolbox->RotatePlacentaBasedOnMinorAxis(3,atoi(argv[8]));
			}
		}
	}
	if (argc == 7){
		toolbox->RotatePlacentaBasedOnMinorAxis();
	}
	cout << "Test3" << "\n";
	std::cout << pd_path<< "\n";
	//toolbox->WriteImage(toolbox->Rot_DopplerData,pd_path);
	toolbox->WriteImage(toolbox->Rot_EDTData,edt_path);
	toolbox->WriteImage(toolbox->Rot_SegData,seg_path);
	toolbox->WriteImage(toolbox->Rot_BModeData,bmode_path);

	vtkSmartPointer<vtkPolyData> spurtPD = vtkSmartPointer<vtkPolyData>::New();
	spurtPD = toolbox->CreateMesh(toolbox->Rot_EDTData,atof(argv[5]));
	spurtPD = toolbox->FinalClean(spurtPD);
	
	vtkSmartPointer<vtkPolyData> fullmesh = vtkSmartPointer<vtkPolyData>::New();
	fullmesh->DeepCopy(spurtPD);

	vtkSmartPointer<vtkPolyData> spurtPD_other = vtkSmartPointer<vtkPolyData>::New();
	spurtPD_other->DeepCopy(fullmesh);
	
	if (argc == 9){
		toolbox->Orientation = atoi(argv[8]);
	}

	int other_Ori =1;
	if (toolbox->Orientation){
		other_Ori = 0;
	}
	std::cout << other_Ori << " Ori \t" << toolbox->Orientation << std::endl;
	if ((argc == 9)  && (atoi(argv[7]) <= 2) && (atoi(argv[7])>=0)){
		spurtPD_other = toolbox->ExtractInterface(spurtPD_other,atoi(argv[7]),other_Ori);
		spurtPD = toolbox->ExtractInterface(spurtPD,atoi(argv[7]),toolbox->Orientation);
	}
	else{
		spurtPD_other = toolbox->ExtractInterface(spurtPD_other,toolbox->Axis,other_Ori);
		spurtPD = toolbox->ExtractInterface(spurtPD,toolbox->Axis,toolbox->Orientation);
	}
	//spurtPD = toolbox->PaintMesh(spurtPD,toolbox->Rot_DopplerData);
	//spurtPD_other = toolbox->PaintMesh(spurtPD_other,toolbox->Rot_DopplerData);		
	//fullmesh = toolbox->PaintMesh(fullmesh,toolbox->Rot_DopplerData);

	toolbox->WriteMesh(toolbox->RotateMeshBack(spurtPD,toolbox->transform),mesh_path);
	toolbox->WriteMesh(toolbox->RotateMeshBack(spurtPD_other,toolbox->transform), mesh_path1);
	toolbox->WriteMesh(toolbox->RotateMeshBack(fullmesh,toolbox->transform),full_mesh_path);

	std::cout << "Complete!!" << std::endl;
}
