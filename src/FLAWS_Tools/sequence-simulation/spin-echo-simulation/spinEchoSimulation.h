#ifndef __SpinEchoSimulationImageFilter_h
#define __SpinEchoSimulationImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkImage.h"

namespace flaws
{

/** \class SpinEchoSimulationImageFilter
 *  \brief Simulate Spin echo images from given tissue partial volume maps and sequence parameters.
 *
 *  \par Description
 *  This filter simulates Spin echo images from given partial volume maps
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
	class SpinEchoSimulationImageFilter:
	public itk::ImageToImageFilter<TInputImage, TOutputImage>
    {
    public:
	/** Standard class typedefs */
	typedef SpinEchoSimulationImageFilter                     Self;
	typedef itk::ImageToImageFilter<TInputImage,TOutputImage> Superclass;
	typedef itk::SmartPointer<Self>                           Pointer;
	typedef itk::SmartPointer<const Self>                     ConstPointer;

	/** Method for creation through the object factory */
	itkNewMacro(Self);

	/** Run-time type information */
	itkTypeMacro(SpinEchoSimulationImageFilter, itk::ImageToImageFilter);

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

      	/** FatMap */
	void SetFatMap(const TInputImage *FatMap)
	{
	    this->SetNthInput(3, const_cast<TInputImage *>(FatMap) );
	}

	const TInputImage* GetFatMap() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(3));
	}

	/** Field to simulate the signal */
	itkSetMacro(Field, float);
	itkGetConstMacro(Field, const float);
      
        /** Simulation with fat */
	itkSetMacro(SimuFat, bool);
	itkGetConstMacro(SimuFat, const bool);
	
	/** SpinEcho sequence information: Echo time (s) */
	itkSetMacro(EchoTime, float);
	itkGetConstMacro(EchoTime, const float);

	/** SpinEcho sequence information: Sequence repetition time (s) */
	itkSetMacro(RepetitionTime, float);
	itkGetConstMacro(RepetitionTime, const float);

    protected:
	/** Constructor */
	SpinEchoSimulationImageFilter();
  
	/** Destructor */
	virtual ~SpinEchoSimulationImageFilter() {};

	/** Check the coherence in sequence parameters */
	void CheckSequenceParameters();

	/** Simulation of the SpinEcho sequence to generate the lookup table (Bloch equations) */
	float signalSimulation(float,float,float);

	/** Run the simulation from tissue segmentations */
	void SimulationFromSegmentations(void);

        /** Run the simulation from tissue segmentations */
	void SimulationFromSegmentationsWithFat(void);

	/** Run the filter */
	virtual void GenerateData();
	
    private:
	SpinEchoSimulationImageFilter(const Self&);  // purposely not implemented
	void operator=(const Self&);  // purposely not implemented

	/* Sequence parameters  */
        float m_Field,m_SimuFat,m_RepetitionTime,m_EchoTime,m_t1WM,m_t1GM,m_t1CSF,m_t1Fat,m_t2WM,m_t2GM,m_t2CSF,m_t2Fat,m_protonDensityWM,m_protonDensityGM,m_protonDensityCSF,m_protonDensityFat;
    };  // end class

}  // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "spinEchoSimulation.hxx"
#endif

#endif
