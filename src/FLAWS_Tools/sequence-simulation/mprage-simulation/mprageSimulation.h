#ifndef __MPRAGESimulationImageFilter_h
#define __MPRAGESimulationImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkImage.h"

namespace flaws
{

/** \class MPRAGESimulationImageFilter
 *  \brief Simulate MPRAGE images from given tissue segmentations and sequence parameters.
 *
 *  \par Description
 *  This filter simulates MPRAGE images from partial volume maps
 *  and sequence parameters using Bloch equations. 
 * 
 *  \par Template parameters
 *  TInputImage: the input image type. TOutputImage: the output image type.
 *
 *  \author Jeremy Beaumont
 *  \version 0.2
 *  \date Thu 31 May 2018 16:11:00
 *  \ingroup sequence-simulation
 */
    template <class TInputImage, class TOutputImage>
	class MPRAGESimulationImageFilter:
	public itk::ImageToImageFilter<TInputImage, TOutputImage>
    {
    public:
	/** Standard class typedefs */
	typedef MPRAGESimulationImageFilter                      Self;
	typedef itk::ImageToImageFilter<TInputImage,TOutputImage> Superclass;
	typedef itk::SmartPointer<Self>                           Pointer;
	typedef itk::SmartPointer<const Self>                     ConstPointer;

	/** Method for creation through the object factory */
	itkNewMacro(Self);

	/** Run-time type information */
	itkTypeMacro(MPRAGESimulationImageFilter, itk::ImageToImageFilter);

	/** WMMap */
	void SetWMMap(const TInputImage *WMMap)
	{
	    this->SetInput(WMMap);
	}

	const TInputImage* GetWMMap()
	{
	    return this->GetInput();
	}

	/** GMMap */
	void SetGMMap(const TInputImage *GMMap)
	{
	    this->SetNthInput(1, const_cast<TInputImage *>(GMMap) );
	}

	const TInputImage* GetGMMap() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(1));
	}

	/** CSFMap */
	void SetCSFMap(const TInputImage *CSFMap)
	{
	    this->SetNthInput(2, const_cast<TInputImage *>(CSFMap) );
	}

	const TInputImage* GetCSFMap() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(2));
	}

	/** B1PlusMap */
	void SetB1PlusMap(const TInputImage *B1PlusMap)
	{
	    this->SetNthInput(3, const_cast<TInputImage *>(B1PlusMap) );
	}

	const TInputImage* GetB1PlusMap() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(3));
	}

	/** B1MinusMap */
	void SetB1MinusMap(const TInputImage *B1MinusMap)
	{
	    this->SetNthInput(4, const_cast<TInputImage *>(B1MinusMap) );
	}

	const TInputImage* GetB1MinusMap() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(4));
	}

	
	/** T1Map */
	void SetT1Map(const TInputImage *T1Map)
	{
	    this->SetNthInput(5, const_cast<TInputImage *>(T1Map) );
	}

	const TInputImage* GetT1Map() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(5));
	}
	
	/* Add B1+ inhomogeneity in the images */
	itkSetMacro(UseB1PlusMap,bool);
	itkGetConstMacro(UseB1PlusMap,bool);

	itkBooleanMacro(UseB1PlusMap);

	/* Add B1- inhomogeneity in the images */
	itkSetMacro(UseB1MinusMap,bool);
	itkGetConstMacro(UseB1MinusMap,bool);

	itkBooleanMacro(UseB1MinusMap);
	
	/** Level of additive gaussian noise (%)*/
	itkSetMacro(NoiseLevel, unsigned short);
	itkGetConstMacro(NoiseLevel, unsigned short);
		
	/** Field to simulate the signal */
	itkSetMacro(Field, float);
	itkGetConstMacro(Field, float);
	
	/** MPRAGE sequence information: Alpha (Degrees) */
	itkSetMacro(Alpha, unsigned char);
	itkGetConstMacro(Alpha, const unsigned char);

	/** MPRAGE sequence information: Inversion time (s) */
	itkSetMacro(InversionTime, float);
	itkGetConstMacro(InversionTime, const float);

	/** MPRAGE sequence information: Echo time (s) */
	itkSetMacro(EchoTime, float);
	itkGetConstMacro(EchoTime, const float);

	/** MPRAGE sequence information: Echo spacing (s) */
	itkSetMacro(EchoSpacing, float);
	itkGetConstMacro(EchoSpacing, const float);
  
	/** MPRAGE sequence information: Sequence repetition time (s) */
	itkSetMacro(RepetitionTime, float);
	itkGetConstMacro(RepetitionTime, const float);

	/** MPRAGE sequence information: Number of slices */
	itkSetMacro(NumberOfSlices, unsigned short);
	itkGetConstMacro(NumberOfSlices, const unsigned short);
      
	/** MPRAGE sequence information: Partial Fourier (6,7,8) */
	itkSetMacro(PartialFourier, unsigned char);
	itkGetConstMacro(PartialFourier, const unsigned char);

    protected:
	/** Constructor */
	MPRAGESimulationImageFilter();
  
	/** Destructor */
	virtual ~MPRAGESimulationImageFilter() {};

	/** Check the coherence in sequence parameters */
	void CheckSequenceParameters();

	/** Simulation of the MPRAGE sequence to generate the lookup table (Bloch equations) */
	float signalSimulation(float,float,float,float);

	/** Run the simulation from tissue segmentations */
	void SimulationFromSegmentations(void);

	/** Run the simulation from t1 map */
	void SimulationFromT1Map(void);
		  
	/** Run the filter */
	virtual void GenerateData();
	
    private:
	MPRAGESimulationImageFilter(const Self&);  // purposely not implemented
	void operator=(const Self&);  // purposely not implemented

	/* Sequence parameters  */
	unsigned char m_PartialFourier,m_Alpha;
	unsigned short m_NumberOfSlices,m_nbEx,m_nbExBefore,m_nbExAfter,m_NoiseLevel;
	float m_Field,m_InversionTime,m_RepetitionTime,m_EchoTime,m_EchoSpacing,m_AlphaRad,m_timeGre,m_timeGreBefore,m_timeGreAfter,m_invEff,m_t1WM,m_t1GM,m_t1CSF,m_t2StarWM,m_t2StarGM,m_t2StarCSF,m_protonDensityWM,m_protonDensityGM,m_protonDensityCSF;

	/* Mode to simulate the combination image */
	bool m_UseB1PlusMap,m_UseB1MinusMap,m_UseVar;

    };  // end class

}  // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "mprageSimulation.hxx"
#endif

#endif
