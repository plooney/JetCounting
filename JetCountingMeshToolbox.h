#pragma once

#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkCellLinks.h>
#include <vtkCleanPolyData.h>
#include <vtkContourFilter.h>
#include <vtkDecimatePro.h>
#include <vtkDelaunay2D.h>
#include <vtkDiscretizableColorTransferFunction.h>
#include <vtkExtractEdges.h>
#include <vtkExtractSelection.h>
#include <vtkExtractSelectedPolyDataIds.h>
#include <vtkFeatureEdges.h>
#include <vtkFillHolesFilter.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkLookupTable.h>
#include <vtkMatrix4x4.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkImageMask.h>
#include <vtkImageConstantPad.h>
#include <vtkMarchingCubes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkProbeFilter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkStripper.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>

#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCenteredAffineTransform.h"
#include "itkDerivativeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageMomentsCalculator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"

typedef itk::Image<double,3> FloatingImageType;
typedef itk::Image<unsigned char,3> ImageType;

//caster from uchar to float
typedef itk::CastImageFilter<ImageType,FloatingImageType> CastImageFilterType;
//caster from float to uchar
typedef itk::CastImageFilter<FloatingImageType,ImageType> CastImageToCharFilterType;

//get the first moment, center of gravity and principal axes.
typedef itk::ImageMomentsCalculator<FloatingImageType> FloatingImageMomentsCalculatorType;
typedef itk::ImageMomentsCalculator<ImageType> ImageMomentsCalculatorType;
typedef itk::Vector<double, 3> VectorType;

//ITK image filters
typedef itk::DiscreteGaussianImageFilter<FloatingImageType,FloatingImageType> GaussianFilterType;
typedef itk::BinaryThresholdImageFilter<FloatingImageType,FloatingImageType> ThresholdFilterType;
typedef itk::DerivativeImageFilter<FloatingImageType,FloatingImageType> GradientFilterType;
typedef itk::MultiplyImageFilter<ImageType,ImageType> MultiplyImageFilter;

//resampling types for the image
typedef itk::ResampleImageFilter<ImageType,FloatingImageType>  ResampleFilterType;  
typedef itk::ResampleImageFilter<FloatingImageType,FloatingImageType>  FloatResampleFilterType;  
typedef itk::CenteredAffineTransform< double, 3 > TransformType;

//ITK to VTK
typedef itk::ImageToVTKImageFilter<FloatingImageType> OutputVTKType;
typedef itk::VTKImageToImageFilter<ImageType> InputVTKType;
typedef itk::VTKImageToImageFilter<FloatingImageType> InputFloatVTKType;

//writers
typedef itk::ImageFileWriter<ImageType> ImageWriterType;
typedef itk::ImageFileWriter<FloatingImageType> FloatImageWriterType;

class JetCountingMeshToolbox
{
	public:
		JetCountingMeshToolbox(void);
		~JetCountingMeshToolbox(void);

		TransformType::Pointer transform;
	
		vtkSmartPointer<vtkImageData> VolumeData;
		vtkSmartPointer<vtkImageData> DopplerData;
		vtkSmartPointer<vtkImageData> SegmentationData;
		vtkSmartPointer<vtkImageData> EDTData;

		vtkSmartPointer<vtkImageData> Rot_SegData;
		vtkSmartPointer<vtkImageData> Rot_EDTData;
		vtkSmartPointer<vtkImageData> Rot_DopplerData;
		vtkSmartPointer<vtkImageData> Rot_BModeData;
		int Orientation;
		int Axis;

		void DebugOn();

		//readers and writers.
		vtkSmartPointer<vtkPolyData> ReadMesh(std::string fname);
		void WriteMesh(vtkSmartPointer<vtkPolyData>);
		void WriteMesh(vtkSmartPointer<vtkPolyData> poly, std::string fname);

		//classes to display a vtkPolyData in standalongRenderWindow.
		void Display(vtkSmartPointer<vtkPolyData>);
	
		//classes used to perform surface smoothing, extracting and cleaning/DT before being able to parameterise.
		vtkSmartPointer<vtkPolyData> SmoothMesh(vtkSmartPointer<vtkPolyData> polydata, vtkSmartPointer<vtkImageData> scalars);		
		vtkSmartPointer<vtkPolyData> ExtractInterface(vtkSmartPointer<vtkPolyData> poly, int axis, int orientation);
		vtkSmartPointer<vtkPolyData> CleanMesh(vtkSmartPointer<vtkPolyData> pd);
		vtkSmartPointer<vtkPolyData> FinalClean(vtkSmartPointer<vtkPolyData> pd);

	
		//copies scalar information from one vtk mesh into one which doesn't.
		void ColorScalarCopy(std::string infname, std::string outfname);

		//given some scalars paint the scalars onto the triangular mesh.
		vtkSmartPointer<vtkPolyData> PaintMesh(vtkSmartPointer<vtkPolyData> pd, vtkSmartPointer<vtkImageData> scalars);
	
		//metric for calculating genus/Euler characterisc of mesh
		int GetEulerCharacteristic(vtkSmartPointer<vtkPolyData> pd);
	
		void GetMinorAxis(vtkSmartPointer<vtkImageData> imgData);

		void GetAxes(ImageMomentsCalculatorType::MatrixType m, int axes[], double angles[]);

		//Takes the existing binary volumes of segmentation and rotate the B-Mode and PD volumes 
		//appropriately to provide a surface to then take a surface based on normals.
		void RotatePlacentaBasedOnMinorAxis(int axis = 3, int interface_axis=0);

		vtkSmartPointer<vtkPolyData> RotateMeshBack(vtkSmartPointer<vtkPolyData> polydata, TransformType::Pointer transform);
	
		int BasalOrientation(ImageType::Pointer dp_image, ImageType::Pointer seg_image, int axis);

		int CountSpurts(vtkSmartPointer<vtkPolyData> pd, double thresh=128.0);

		void WriteInfo(std::string pd_path, int axis, int orientation);

		vtkSmartPointer<vtkPolyData> CreateMesh(vtkSmartPointer<vtkImageData> edt, double distance);
	
		//image writers
		void WriteImage(vtkSmartPointer<vtkImageData> img, std::string path);
		void WriteImage(ImageType::Pointer img, std::string path);
		void WriteImage(FloatingImageType::Pointer img, std::string path);
};
