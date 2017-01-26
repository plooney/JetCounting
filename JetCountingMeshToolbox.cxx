/**
	JetCounting
	JetCountingMeshToolbox.cxx
	Purpose: Class containing swiss army-knife of functions and tools for manipulating meshes and 3D volume data to count and analyse 3D PD data.

	@author Gordon Stevenson
	@version 0.1 16/01/2017

*/

#include "JetCountingMeshToolbox.h"

JetCountingMeshToolbox::JetCountingMeshToolbox(void)
{
	Rot_SegData = vtkSmartPointer<vtkImageData>::New();
	Rot_EDTData = vtkSmartPointer<vtkImageData>::New();
	Rot_DopplerData = vtkSmartPointer<vtkImageData>::New();
	Rot_BModeData = vtkSmartPointer<vtkImageData>::New();

	VolumeData = vtkSmartPointer<vtkImageData>::New();
	DopplerData = vtkSmartPointer<vtkImageData>::New();
	SegmentationData = vtkSmartPointer<vtkImageData>::New();
	EDTData = vtkSmartPointer<vtkImageData>::New();
	
	transform = NULL;
}

JetCountingMeshToolbox::~JetCountingMeshToolbox(void){}

 vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::ReadMesh(std::string fname){
	 /**
	 Read mesh file from string fname use PolyDataReader. Assumes PolyData format (.vtk)

	 @param fname The filename to the mesh
	 @return SmartPointer to vtkPolyData object, poly.
	 */
	  
	vtkSmartPointer<vtkPolyDataReader> reader;
	vtkSmartPointer<vtkPolyData> poly;
	
	size_t dotpos;
	std::string format;

	dotpos  = fname.find_last_of(".");

	if (dotpos == -1) {
		return NULL;
	}

	format = fname.substr(dotpos,4);
	//got the format string.

	if (format.compare(".vtk") == 0){
		reader = vtkSmartPointer<vtkPolyDataReader>::New();	
		reader->SetFileName(fname.c_str());
		reader->Update();
		poly = reader->GetOutput();
		return poly;
	}
	return NULL;
}

 vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::RotateMeshBack(vtkSmartPointer<vtkPolyData> polydata, TransformType::Pointer transform){
	 
	vtkSmartPointer<vtkMatrix4x4> vtk_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	TransformType::MatrixType matrix = transform->GetMatrix(); 
	TransformType::OffsetType offset = transform->GetOffset(); 
	TransformType::InputPointType center = transform->GetCenter();

	for(unsigned int i=0; i<3; i++ ){ 
		for(unsigned int j=0; j<3; j++ ){ 
			vtk_matrix->SetElement(i,j, 
			matrix.GetVnlMatrix().get(i,j));   
		} 
		vtk_matrix->SetElement( i, 3, offset[i]); 
	} 

	vtkSmartPointer<vtkTransform> poly_transform = vtkSmartPointer<vtkTransform>::New();
	poly_transform->SetMatrix(vtk_matrix);
	poly_transform->Update();
	
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetInputData(polydata);
	transformFilter->SetTransform(poly_transform);
	transformFilter->Update();
	return transformFilter->GetOutput();
}

vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::SmoothMesh(vtkSmartPointer<vtkPolyData> polydata, vtkSmartPointer<vtkImageData> scalars){
	 /** 
		 Take a MC Mesh and add greyscale data to it by probe, then decimate, smooth and produce the
		 mesh to then be cut based on curvature. 

		 @param polydata 3D vtkPolyData that needs to be coloured and smoothed
		 @param scalars 3D volumetric power Doppler data used to colour the vtkPolyData input, polydata.
		 
		 @return 3D polydata
	 */

	vtkSmartPointer<vtkTriangleFilter> triang = vtkSmartPointer<vtkTriangleFilter>::New();
	triang->SetInputData(polydata);
	triang->SetPassLines(0);
	triang->SetPassVerts(0);
	triang->Update();

	vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
	probe->SetInputConnection(triang->GetOutputPort());
	probe->SetSourceData(scalars);
	probe->Update();

	vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
	deci->SetInputConnection(triang->GetOutputPort());
	deci->SetTargetReduction(0.7);
	deci->PreserveTopologyOn();
	deci->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> polyconn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	polyconn->SetInputConnection(deci->GetOutputPort());
	polyconn->SetExtractionModeToLargestRegion();
	polyconn->Update();

	vtkSmartPointer<vtkPolyData> pd = polyconn->GetOutput();
	return pd;
 }

vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::ExtractInterface(vtkSmartPointer<vtkPolyData> poly, int axis, int orientation){

	/**
		Extracts from a polydata a cut-surface based on the surface normals chosen by axis and surface normal orientation.

		@param poly The 3D polydata - assumes we need to cut this.
		@param axis integer describing x,y or z. (0 to 2)
		@param orientation negative or positive selection based on bool (0 or 1) respectively.

		@return The cut 3D polydata.

	*/

	vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputData(poly);
	normals->ComputeCellNormalsOn();
	normals->ComputePointNormalsOn();
	normals->Update();

	vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);

	vtkSmartPointer<vtkFloatArray> normArray;
	normArray = vtkFloatArray::SafeDownCast(normals->GetOutput()->GetCellData()->GetNormals());
	int i;

	float norm;
	if (orientation){
		for (i = 0; i < normArray->GetNumberOfTuples(); i++){		
		norm = normArray->GetTuple(i)[axis];  // x,y,z axis based on input
			if (norm > 0.000000){
				ids->InsertNextValue(i);		 
			}	
		}
	}
	else{
		for (i = 0; i < normArray->GetNumberOfTuples(); i++){		
		norm = normArray->GetTuple(i)[axis];  // x,y,z axis based on input
			if (norm < 0.000000){
				ids->InsertNextValue(i);		 
			}	
		}
	}
			
	vtkSmartPointer<vtkSelectionNode> selectionnode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionnode->SetFieldType(vtkSelectionNode::CELL);
	selectionnode->SetContentType(vtkSelectionNode::INDICES);
	selectionnode->SetSelectionList(ids);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionnode);	

	vtkSmartPointer<vtkExtractSelectedPolyDataIds> extractSelection = vtkSmartPointer<vtkExtractSelectedPolyDataIds>::New();
	extractSelection->SetInputData(0, poly);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();  
	
	vtkSmartPointer<vtkPolyDataConnectivityFilter> polyfilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	polyfilter->SetInputData(extractSelection->GetOutput());
	polyfilter->SetExtractionModeToLargestRegion();  
	polyfilter->Update();	

	vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
	tri->SetInputConnection(polyfilter->GetOutputPort());
	tri->SetPassLines(1);
	tri->SetPassVerts(1);
	tri->Update();

	vtkSmartPointer<vtkFillHolesFilter> filler = vtkSmartPointer<vtkFillHolesFilter>::New();
	filler->SetInputData(tri->GetOutput());  
	filler->SetHoleSize(20.0);
	filler->Update();

	vtkSmartPointer<vtkPolyData> outpoly = vtkSmartPointer<vtkPolyData>::New();
	outpoly = filler->GetOutput();
	return outpoly;
}

//write to a fixed file "C:\\spurt_info.csv" the information on path, axis and orientation for this particular run.
void JetCountingMeshToolbox::WriteInfo(std::string pd_path, int axis, int orientation){	
	ofstream csvfile;
	csvfile.open("C:\\spurt_info.csv", ios::app);
	csvfile << pd_path << "," << axis << "," << orientation << "\r";
	csvfile.close();
}

//Take in a polydata pd, and remove any cells (triangles) that are
//either non-manifold (no neighbours) or are hanging. (one point only in the triangle) - so returned "clean")
vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::CleanMesh(vtkSmartPointer<vtkPolyData> pd){
	//list for cells to keep.
	vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();	
	
	//GetEulerCharacteristic(pd);

	//for each cell in the mesh
	unsigned short ncells;
	int lonevertex,badcellcount;
	//variable to store number of cells that are deleted.
	badcellcount = 0.0;
	int a = 0;
	bool cornersfound = true;
	while(cornersfound){
		//create connectivity list for cells. --required for GetPointCells calls.
		pd->BuildLinks();
		cornersfound = false;
		
		int b;
		//for each cell.
		for(vtkIdType i = 0; i<pd->GetNumberOfCells(); i++){
			//generate the points each that cell.
			vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
			pointIds = pd->GetCell(i)->GetPointIds();			
			if (pointIds->GetNumberOfIds() > 3) {
				std::cout << "Polygon got more than 3 vertices!" << std::endl;
			}

			vtkIdType *dummyIDT;
			lonevertex = 0;
			ncells = 0.0;
			//so for each point in the triangle, check that the vertex has more than one cell
			//that uses it.
				
				//loop over each vertex in the cell. i.e 3 times.
				for(vtkIdType j = 0; j < pointIds->GetNumberOfIds(); j++){				
					//get the number of cells that use this vertex.
					b = pointIds->GetId(j);
					pd->GetPointCells(pointIds->GetId(j),ncells,dummyIDT);
					if(ncells == 2.0){
					//may be a duplicate cell lurking underneath.
						if( pd->IsPointUsedByCell(pointIds->GetId(0),dummyIDT[0]) &&
							pd->IsPointUsedByCell(pointIds->GetId(1),dummyIDT[0]) &&
							pd->IsPointUsedByCell(pointIds->GetId(2),dummyIDT[0]) &&
							pd->IsPointUsedByCell(pointIds->GetId(0),dummyIDT[1]) &&
							pd->IsPointUsedByCell(pointIds->GetId(1),dummyIDT[1]) &&
							pd->IsPointUsedByCell(pointIds->GetId(2),dummyIDT[1]) ){
							lonevertex++;
						}
					}
					//if vertex is not used by any other cells then -- including itself.
					if (ncells <= 1.0){				
						//vertex is a corner, and therefore a candidate for deletion.
						lonevertex++;
					}
				}
					
				//or, check that for an edge there is no neighbour.
				// but we don't need to check any cells that have already been targetted
				// for deletion
				if(lonevertex == 0){
					int noEdgeNeighbors = 0;
					pd->GetCellEdgeNeighbors(i,pointIds->GetId(0),pointIds->GetId(1),neighbors);
						if (neighbors->GetNumberOfIds() == 0){
							noEdgeNeighbors += 1;
						}
					neighbors->Reset();
					pd->GetCellEdgeNeighbors(i,pointIds->GetId(1),pointIds->GetId(2),neighbors);
						if (neighbors->GetNumberOfIds() == 0){
							noEdgeNeighbors += 1;
						}
					neighbors->Reset();
					pd->GetCellEdgeNeighbors(i,pointIds->GetId(2),pointIds->GetId(0),neighbors);
						if (neighbors->GetNumberOfIds() == 0){
							noEdgeNeighbors += 1;
						}

					if(noEdgeNeighbors >= 2){
						lonevertex++;
					}
				}

				//have we got any lone vertices?
				if(lonevertex!=0){
					badcellcount++;
					//put up cell as candidate for deletion.
					pd->DeleteCell(i);					
					//need to loop back over the cells again checking for more corners
					//that may have been revealed.
					cornersfound = true;
				}
		}
		pd->RemoveDeletedCells();
		pd->Modified();
	}

	vtkSmartPointer<vtkPolyDataConnectivityFilter> polyconnect = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	polyconnect->SetInputData(pd);
	polyconnect->SetExtractionModeToLargestRegion();
	polyconnect->Update();

	//GetEulerCharacteristic(polyconnect->GetOutput());
	return polyconnect->GetOutput();
}

//need to try and clean, fill holes and remove connected components. 
vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::FinalClean(vtkSmartPointer<vtkPolyData> pd){

	vtkSmartPointer<vtkPolyDataConnectivityFilter> polyconnect = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	polyconnect->SetInputData(pd);
	polyconnect->SetExtractionModeToLargestRegion();
	polyconnect->Update();

	vtkSmartPointer<vtkFillHolesFilter> fholes = vtkSmartPointer<vtkFillHolesFilter>::New();
	fholes->SetInputConnection(polyconnect->GetOutputPort());
	fholes->SetHoleSize(2000.0);
	fholes->Update();

	return fholes->GetOutput();
}

//helper function to extract Euler Characteristic - used to check for handles (holes)
int JetCountingMeshToolbox::GetEulerCharacteristic(vtkSmartPointer<vtkPolyData> pd){

	vtkSmartPointer<vtkExtractEdges> extractEdges = 
	  vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputData(pd);
	extractEdges->Update();
 
	vtkSmartPointer<vtkCellArray> lines = extractEdges->GetOutput()->GetLines();
	int linesnum = extractEdges->GetOutput()->GetLines()->GetNumberOfCells();
	
	vtkSmartPointer<vtkPoints> points = extractEdges->GetOutput()->GetPoints();
	int pointnum = (int) extractEdges->GetOutput()->GetPoints()->GetNumberOfPoints();
 
	int EC = pointnum - linesnum + pd->GetNumberOfCells();
	std::printf("Euler Characteristic V - E + F = 2 (for a closed a surface)\n");
	std::printf("%d -%d + %d = %d \n",pointnum,linesnum,pd->GetNumberOfCells(),EC);
	
	return EC;
}


void JetCountingMeshToolbox::ColorScalarCopy(std::string infname, std::string outfname){
//takes the scalar information stored in the original surface file, and copies it to the another (parameterised) file 
//ASSUMES that only one ``POINT DATA'' portion in the vtk name is given.
	
	//open the file
	FILE * sFile;
	sFile = fopen(infname.c_str(),"r");

	char buff[100];

	std::string cscalars = "COLOR_SCALARS";
	std::string crtscalars = "SCALARS";
	std::string pdata = "POINT_DATA";
	std::string sfileStore;

//	long pos;
	long lSize;

	// obtain file size:
	fseek (sFile , 0 , SEEK_END);
	lSize = ftell (sFile);
	rewind (sFile);

	//std::stringstream sstream(std::stringstream::in | std::stringstream::out);
	std::ofstream outfile(outfname.c_str(),ios::ate | ios::app);
	//input and output stream.

	//start reading..
	if (sFile==NULL) perror ("Error opening file");
		else
		{
			while (!feof(sFile)) {
				//read the line in...
				fgets (buff,100,sFile);
				sfileStore = buff;				

				//in the line, we found an occurance of POINT_DATA.
				if(sfileStore.find(pdata) != -1){
					//now need to pull all the data from this point to the end of the file into a buffer.
					//move the pointer to the previous line.
					fseek(sFile,(ftell(sFile)-sfileStore.size()-1),ios::beg);
					outfile << "\n";

					while (!feof(sFile)) {
					//we then loop through the file from this point and pass every line into a string stream.
						fgets (buff,100,sFile);
						//replace the colors line.
						sfileStore = buff;
						//if we find a COLOR_SCALARS in the line....
						if(sfileStore.find(cscalars) != -1){
							//sfileStore.replace(sfileStore.find(cscalars),cscalars.length(),crtscalars);
							outfile << "SCALARS scalars float 1 \n";
							outfile << "LOOKUP_TABLE default \n";
						}
						else{
						//put to file. 
							if(!feof(sFile)){
								outfile << sfileStore.c_str();						
							}
						}
					}

					//when done, break.
					break;
				};
				sfileStore.clear();
			}
	}


	//close the files.
	fclose(sFile);
	outfile.close();
	//paramfile.close();
}

//simple counter
int JetCountingMeshToolbox::CountSpurts(vtkSmartPointer<vtkPolyData> pd, double thresh){
	vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
	contourFilter->SetInputData(pd);
	contourFilter->SetValue(0,thresh);
	contourFilter->SetGenerateTriangles(1);
	contourFilter->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> polyconnectFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	polyconnectFilter->SetInputData(contourFilter->GetOutput());
	polyconnectFilter->SetExtractionModeToAllRegions();
	//polyconnectFilter->SetScalarConnectivity(1);
	polyconnectFilter->ColorRegionsOn();
	polyconnectFilter->Update();
	cout << "RegionIDs # " << polyconnectFilter->GetNumberOfExtractedRegions() << "\n";
	if  (polyconnectFilter->GetNumberOfExtractedRegions() == 0){
		return 0;
	}
	else{
	double *range = polyconnectFilter->GetOutput()->GetPointData()->GetArray("RegionId")->GetRange();
	//std::string o_str = "D:\\mesh_test.vtk";
	//WriteMesh(polyconnectFilter->GetOutput(),o_str);
	return static_cast<int>(range[1]-range[2]);
	}
}

//display out polyData
void JetCountingMeshToolbox::Display(vtkSmartPointer<vtkPolyData> polyData){

  vtkSmartPointer<vtkPolyDataNormals> norm = vtkSmartPointer<vtkPolyDataNormals>::New();
  norm->SetInputData(polyData);
  norm->SetFeatureAngle(45);

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(norm->GetOutput());

  vtkSmartPointer<vtkLookupTable> lut1 = vtkSmartPointer<vtkLookupTable>::New();
  lut1->SetHueRange(0.0,0.0);
  lut1->SetSaturationRange(0.0,0.0);
  lut1->SetValueRange(0,1.0);
  lut1->Build();

  mapper->SetLookupTable(lut1);
  mapper->SetScalarVisibility(1);
  mapper->SetScalarRange(0,255);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
  ren->AddActor(actor);
  ren->SetBackground(1,1,1);

  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren);
  renWin->SetSize(300,300);

  vtkSmartPointer<vtkRenderWindowInteractor> inter = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  inter->SetRenderWindow(renWin);
  renWin->Render();
  inter->Start();
}

void JetCountingMeshToolbox::WriteMesh(vtkSmartPointer<vtkPolyData> poly){
	//DEBUG write mesh file
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName("C:/test.vtk");
	writer->Write();
}

void JetCountingMeshToolbox::WriteMesh(vtkSmartPointer<vtkPolyData> poly, std::string fname){
	//write mesh file, poly to location fname.
	//outputs the mesh file based on the file extension....
	vtkSmartPointer<vtkPolyDataWriter> writer;
	
	size_t dotpos;
	std::string format;

	dotpos  = fname.find_last_of(".");

	if (dotpos == -1) {
		return;
	}

	format = fname.substr(dotpos,4);
	//got the format string.

	if (format.compare(".vtk") == 0){
	writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fname.c_str());
	writer->Write();
	}
}

void JetCountingMeshToolbox::GetAxes(ImageMomentsCalculatorType::MatrixType m, int axes[], double angles[]){
	/*given a matrix of the principal axes in physical co-ordinates
	//which are sorted from smallest to largest
	//return an int array which contains the membership of the rows of the matrix
	//to each axis in order
	//for each row find the highest abs value in i,0 i,1 i,2*/
	for (int i=0; i<m.GetVnlMatrix().rows(); i++){
		double max_val = std::abs(m(i,0));
		double val = m(i,0);
		int axis = 0;
		//for each value in the row - 
		for(int j=1; j<m.GetVnlMatrix().cols(); j++){
			//if the next value is greater than the one just now.
			//store the values in their respective variables.
			if (max_val < std::abs(m(i,j))){
				max_val = std::abs(m(i,j));
				val = m(i,j);
				axis = j;
			}
		}
		//now moves these values into 
		axes[i] = axis;
		angles[i] = val;
	}
}

int JetCountingMeshToolbox::BasalOrientation(ImageType::Pointer dp_image, ImageType::Pointer seg_image, int direction){
	CastImageFilterType::Pointer seg_castFilter = CastImageFilterType::New();
	seg_castFilter->SetInput(seg_image);
	seg_castFilter->Update();

  //simple gaussian filter and threshold.
  GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
  gaussianFilter->SetInput(seg_castFilter->GetOutput());
  gaussianFilter->SetUseImageSpacingOn();
  gaussianFilter->SetVariance(3.0);
  gaussianFilter->Update();
  //threshold the gaussian filter.
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput(gaussianFilter->GetOutput());
  thresholdFilter->SetLowerThreshold(0.5);
  thresholdFilter->SetUpperThreshold(1.0);
  thresholdFilter->SetInsideValue(1.0);
  thresholdFilter->SetOutsideValue(0.0);
  thresholdFilter->Update();

  //now, extract the gradient of the image.
  GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
  gradientFilter->SetInput(thresholdFilter->GetOutput());
  //gradientFilter->SetUseImageSpacingOn();
  gradientFilter->SetDirection(direction);
  gradientFilter->Update();

	//returns 1 if positive normals contain more vascularity.
	//0 if it is negative.
	ThresholdFilterType::Pointer neg_thresholdFilter = ThresholdFilterType::New();
	neg_thresholdFilter->SetInput(gradientFilter->GetOutput());
	neg_thresholdFilter->SetUpperThreshold(0.0);
	neg_thresholdFilter->SetLowerThreshold(-1.0);
	neg_thresholdFilter->SetInsideValue(1.0);
	neg_thresholdFilter->SetOutsideValue(0.0);
	neg_thresholdFilter->Update();

	FloatingImageType::Pointer negativeNormals = neg_thresholdFilter->GetOutput();
	
	ThresholdFilterType::Pointer pos_thresholdFilter = ThresholdFilterType::New();
	pos_thresholdFilter->SetInput(gradientFilter->GetOutput());
	pos_thresholdFilter->SetUpperThreshold(1.0);
	pos_thresholdFilter->SetLowerThreshold(0.0);
	pos_thresholdFilter->SetOutsideValue(0.0);
	pos_thresholdFilter->SetInsideValue(1.0);
	pos_thresholdFilter->Update();
	//pos_thresholdFilter->SetInsideValue;
	FloatingImageType::Pointer positiveNormals = pos_thresholdFilter->GetOutput();

	//cast the normals to integral type - AND filter requirment.
	CastImageToCharFilterType::Pointer pos_castFilter = CastImageToCharFilterType::New();
	pos_castFilter->SetInput(positiveNormals);
	pos_castFilter->Update();

	//cast the normals to integral type - AND filter requirment.
	CastImageToCharFilterType::Pointer neg_castFilter = CastImageToCharFilterType::New();
	neg_castFilter->SetInput(negativeNormals);
	neg_castFilter->Update();

	//CastImageToCharFilterType::Pointer doppler_castFilter = CastImageToCharFilterType::New();
	//doppler_castFilter->SetInput(dp_image);
	//doppler_castFilter->Update();
		
	//std::cout << "Casted\n";
	MultiplyImageFilter::Pointer pos_multiFilter = MultiplyImageFilter::New();
	pos_multiFilter->SetInput1(pos_castFilter->GetOutput());
	pos_multiFilter->SetInput2(dp_image);
	pos_multiFilter->Update();

	MultiplyImageFilter::Pointer neg_multiFilter = MultiplyImageFilter::New();
	neg_multiFilter->SetInput1(neg_castFilter->GetOutput());
	neg_multiFilter->SetInput2(dp_image);
	neg_multiFilter->Update();
	//std::cout << "Multiply \n";
	
	ImageMomentsCalculatorType::Pointer pos_MomentCalc = ImageMomentsCalculatorType::New();
	pos_MomentCalc->SetImage(pos_multiFilter->GetOutput());
	pos_MomentCalc->Compute();

	ImageMomentsCalculatorType::Pointer neg_MomentCalc = ImageMomentsCalculatorType::New();
	neg_MomentCalc->SetImage(neg_multiFilter->GetOutput());
	neg_MomentCalc->Compute();
	
	float posMass = pos_MomentCalc->GetTotalMass();
	float negMass = neg_MomentCalc->GetTotalMass();

	if (posMass > negMass){
		std::cout << "Found Orientation 1\n";
		return 1;
	}

	if (posMass < negMass){
		std::cout << "Found Orientation 0\n";
		return 0;
	}

	if (posMass == negMass){
		std::cout << "blrgh!!\n";
		return EXIT_FAILURE; 
	}
	std::cout << "blrgh!!\n";
	return EXIT_FAILURE;
}

void JetCountingMeshToolbox::WriteImage(ImageType::Pointer img, std::string path){
	  ImageWriterType::Pointer writer = ImageWriterType::New();
	  writer->SetFileName(path.c_str());
	  writer->SetInput(img);
	  writer->Update();
}

void JetCountingMeshToolbox::WriteImage(FloatingImageType::Pointer img, std::string path){
	  FloatImageWriterType::Pointer writer = FloatImageWriterType::New();
	  writer->SetFileName(path.c_str());
	  writer->SetInput(img);
	  writer->Update();
}

void JetCountingMeshToolbox::RotatePlacentaBasedOnMinorAxis(int axis,int interface_axis){
  //convert to ITK
  InputVTKType::Pointer vtkToITKConverter = InputVTKType::New();
  vtkToITKConverter->SetInput(SegmentationData);
  vtkToITKConverter->Update();

  InputVTKType::Pointer vtkToITKConverter_BMode = InputVTKType::New();
  vtkToITKConverter_BMode->SetInput(VolumeData);
  vtkToITKConverter_BMode->Update();

  InputVTKType::Pointer vtkToITKConv_Doppler = InputVTKType::New();
  vtkToITKConv_Doppler->SetInput(DopplerData);
  vtkToITKConv_Doppler->Update();

  InputFloatVTKType::Pointer vtkToITKConv_EDT = InputFloatVTKType::New();
  vtkToITKConv_EDT->SetInput(EDTData);
  vtkToITKConv_EDT->Update();

  //get the prinicpal axes of the segmentation image
  ImageMomentsCalculatorType::Pointer calculatorFilter = ImageMomentsCalculatorType::New();
  calculatorFilter->SetImage(vtkToITKConverter->GetOutput());
  calculatorFilter->Compute();
  
  //print out principal axes,
  std::cout <<	"Principal Axes "<< std::endl;
  std::cout << calculatorFilter->GetPrincipalAxes()<< std::endl;

  std::cout << "Principal Moments"<< std::endl;
  std::cout << calculatorFilter->GetPrincipalMoments() << std::endl;
  std::cout << calculatorFilter->GetSecondMoments() << std::endl;

  //take the smallest principal axis, from the calculator -  this is the axis where the two sides of the
  //placenta are at either end.
  VectorType principal_axis;
  principal_axis.SetVnlVector(calculatorFilter->GetPrincipalAxes().GetVnlMatrix().get_row(0));
  
  //print out the angle of the largest vector.
  double min_axis = std::max(std::abs(principal_axis.Get_vnl_vector()[0]) , std::max(std::abs(principal_axis.Get_vnl_vector()[1]),std::abs(principal_axis.Get_vnl_vector()[2])));
  std::cout << "Minimum Axis Computed" << min_axis << std::endl;
  //centered or typical affine transform?
  transform = TransformType::New();

  ImageType::PointType pt_center;
  pt_center[0]  = calculatorFilter->GetCenterOfGravity()[0];
  pt_center[1]  = calculatorFilter->GetCenterOfGravity()[1];
  pt_center[2]  = calculatorFilter->GetCenterOfGravity()[2];

  //center the transform relative to the center of mass.
  transform->SetCenter(pt_center);
  std::cout << "center -" << pt_center[0] << ", " << pt_center[1] << ", " << pt_center[2] << "\n";
  //transform->ComputeOffset();
  //TransformType::OutputVectorType axis;
  
  char ax[3] = {'x','y','z'};
	//get the largest vector component for each angle and rotate by the smallest principal axis
	//unless manually overridden
	int axes[3];
	double angles[3];
	GetAxes(calculatorFilter->GetPrincipalAxes(), axes,angles);
	std::cout << "axes - " << axes[0] << " " << axes[1] << " " << axes[2] << "\n";
	std::cout << "angles - " << angles[0] << " " << angles[1] << " " << angles[2] << "\n";;  

	if (axis == 3){
	  Axis = axes[0];
	  std::cout << ax[axes[0]] << "-axis selected \n";
	  std::cout << "Angle of minimal principal axis = "  << (angles[0]*180.0)/3.14 << std::endl;
	  transform->Rotate3D(angles[0],axes[0],false);
	  //transform->Rotate3D(angles[1],axes[1],false);
	  //transform->Rotate3D(angles[2],axes[2],false);
	  if (interface_axis == 3){
	   Orientation = BasalOrientation(vtkToITKConv_Doppler->GetOutput(),vtkToITKConverter->GetOutput(),axes[0]);
	  }
	  else{
	   Orientation = BasalOrientation(vtkToITKConv_Doppler->GetOutput(),vtkToITKConverter->GetOutput(),interface_axis);
	  }
  }
  else{
		for(int a_x = 0; a_x <3; a_x++){
			if (axes[a_x] == axis){
				Axis = axes[a_x];
				std::cout << ax[Axis] << "-axis selected \n";
				std::cout << "Angle of axis = "  << (angles[a_x]*180.0)/3.14 << std::endl;
				transform->Rotate3D(angles[a_x],axes[a_x],false);

				if (interface_axis == 3){
					Orientation = BasalOrientation(vtkToITKConv_Doppler->GetOutput(),vtkToITKConverter->GetOutput(),axes[a_x]);
				}
				else{
					Orientation = BasalOrientation(vtkToITKConv_Doppler->GetOutput(),vtkToITKConverter->GetOutput(),interface_axis);
				}
			}
		}
  }

  //std::cout << "Orientation Based on Vascularity " <<



  // Without this, the program crashes == need to adjust the output to match the input images parameters.  
  //need to pad the size in order to avoid cropping out the images.
  ImageType::SizeType size;
  size = vtkToITKConverter->GetOutput()->GetLargestPossibleRegion().GetSize();
  std::cout << "Input Origin of Volume " << vtkToITKConverter->GetOutput()->GetOrigin() << "\n";

  ImageType::SpacingType output_spacings;
  output_spacings = vtkToITKConverter->GetOutput()->GetSpacing();
  //output_spacings[0] = 0.33;
  //output_spacings[1] = 0.33;
  //output_spacings[2] = 0.33;

  //calculate diagonal
  double diag =  std::sqrt((double) ((size[0]*size[0]) + (size[1]*size[1]) + (size[2]*size[2])) );
  int xdiag = (int) diag/10;
  int ydiag = (int) diag/10;
  int zdiag = (int) diag/10;

  //size[0] = size[0] + xdiag;
  //size[1] = size[1] + ydiag; 
  //size[2] = size[2] + zdiag;

  ImageType::PointType origin = vtkToITKConverter->GetOutput()->GetOrigin();
  //origin[0] = origin[0] - (xdiag/2);
  //origin[1] = origin[1] - (ydiag/2);
  //origin[2] = origin[2] - (zdiag/2);
  std::cout << "Output Origin of Volume " << origin << "\n";

  //resample again for image, need to set the same transform and the correct output.
  ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  resampleFilter->SetTransform(transform);
  resampleFilter->SetInput(vtkToITKConverter->GetOutput());
  resampleFilter->SetOutputSpacing(output_spacings);
  resampleFilter->SetSize ( size );
  resampleFilter->SetOutputOrigin(origin);
  resampleFilter->Update();

  //resample again for image, need to set the same transform and the correct output.
  ResampleFilterType::Pointer resampleBModeFilter = ResampleFilterType::New();
  resampleBModeFilter->SetTransform(transform);
  resampleBModeFilter->SetInput(vtkToITKConverter_BMode->GetOutput());
  resampleBModeFilter->SetOutputSpacing(output_spacings);
  resampleBModeFilter->SetSize ( size );
  resampleBModeFilter->SetOutputOrigin(origin);
  resampleBModeFilter->Update();


  ResampleFilterType::Pointer resampleDopplerFilter = ResampleFilterType::New();
  resampleDopplerFilter->SetTransform(transform);
  resampleDopplerFilter->SetInput(vtkToITKConv_Doppler->GetOutput());
  resampleDopplerFilter->SetSize ( size );
  resampleDopplerFilter->SetOutputSpacing(output_spacings);
  resampleDopplerFilter->SetOutputOrigin(origin);
  resampleDopplerFilter->Update();

  FloatResampleFilterType::Pointer resampleEDTFilter = FloatResampleFilterType::New();
  resampleEDTFilter->SetTransform(transform);
  resampleEDTFilter->SetInput(vtkToITKConv_EDT->GetOutput());
  resampleEDTFilter->SetSize ( size );
  resampleEDTFilter->SetOutputSpacing(output_spacings);
  resampleEDTFilter->SetOutputOrigin(origin);
  resampleEDTFilter->Update();
  
  OutputVTKType::Pointer rotPD_toVTK = OutputVTKType::New();
  OutputVTKType::Pointer rotEDT_toVTK = OutputVTKType::New();
  OutputVTKType::Pointer rotSeg_toVTK = OutputVTKType::New();
  OutputVTKType::Pointer rotBMode_toVTK = OutputVTKType::New();

  rotPD_toVTK->SetInput(resampleDopplerFilter->GetOutput());
  rotPD_toVTK->Update();
  
  rotEDT_toVTK->SetInput(resampleEDTFilter->GetOutput());
  rotEDT_toVTK->Update();
  
  rotSeg_toVTK->SetInput(resampleFilter->GetOutput());
  rotSeg_toVTK->Update();
  
  rotBMode_toVTK->SetInput(resampleBModeFilter->GetOutput());
  rotBMode_toVTK->Update();

  
  Rot_DopplerData->DeepCopy(rotPD_toVTK->GetOutput());
  Rot_SegData->DeepCopy(rotSeg_toVTK->GetOutput());
  Rot_EDTData->DeepCopy(rotEDT_toVTK->GetOutput());
  Rot_BModeData->DeepCopy(rotBMode_toVTK->GetOutput());
}


vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::CreateMesh(vtkSmartPointer<vtkImageData> edt, double distance){
	vtkSmartPointer<vtkMarchingCubes> MCMesh = vtkSmartPointer<vtkMarchingCubes>::New();
	MCMesh->SetInputData(edt);
	MCMesh->SetValue(0,distance);
	MCMesh->SetComputeNormals(1);
	MCMesh->Update();

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoother->SetNumberOfIterations(100);
	smoother->SetInputConnection(MCMesh->GetOutputPort());
	smoother->Update();

	vtkSmartPointer<vtkFillHolesFilter> fholes = vtkSmartPointer<vtkFillHolesFilter>::New();
	fholes->SetInputConnection(smoother->GetOutputPort());
	fholes->SetHoleSize(1e6);	
	fholes->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> polyFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	polyFilter->SetInputConnection(fholes->GetOutputPort());
	polyFilter->SetExtractionModeToLargestRegion();
	polyFilter->Update();

	/*vtkSmartPointer<vtkPolyDataWriter> dbgWriter  = vtkSmartPointer<vtkPolyDataWriter>::New();
	dbgWriter->SetInputData(MCMesh->GetOutput());
	dbgWriter->SetFileName("D:\\test.vtk");
	dbgWriter->Write();*/

	return polyFilter->GetOutput();

}

vtkSmartPointer<vtkPolyData> JetCountingMeshToolbox::PaintMesh(vtkSmartPointer<vtkPolyData> pd, vtkSmartPointer<vtkImageData> scalars){
	 ///Take a MC Mesh and add greyscale data to it by probe, then decimate, smooth and produce the
	 //mesh to then be cut based on curvature. 

	scalars->AllocateScalars(VTK_DOUBLE,1);

	vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
	probe->SetInputData(pd);
	probe->SetSourceData(scalars);
	probe->Update();	

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoother->SetNumberOfIterations(250);
	smoother->SetInputConnection(probe->GetOutputPort());
	smoother->Update();

	vtkSmartPointer<vtkTriangleFilter> trifilter = vtkSmartPointer<vtkTriangleFilter>::New();
	trifilter->SetInputConnection(smoother->GetOutputPort());
	trifilter->SetPassLines(0);
	trifilter->SetPassVerts(0);
	trifilter->Update();

	vtkSmartPointer<vtkPolyData> polyd = trifilter->GetOutput();
	return polyd;
 }

void JetCountingMeshToolbox::WriteImage(vtkSmartPointer<vtkImageData> img, std::string path){
	vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
	writer->SetInputData(img);
	writer->SetFileName(path.c_str());
	writer->Write();	
}