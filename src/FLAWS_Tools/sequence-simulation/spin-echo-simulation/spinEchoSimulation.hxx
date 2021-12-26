#ifndef __SpinEchoSimulationImageFilter_hxx
#define __SpinEchoSimulationImageFilter_hxx

#include <math.h>

#include "spinEchoSimulation.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace flaws
{

    /* Constructor */
    template <class TInputImage, class TOutputImage>
    SpinEchoSimulationImageFilter<TInputImage, TOutputImage>
    ::SpinEchoSimulationImageFilter():
	m_EchoTime(0),
	m_RepetitionTime(0),
	m_Field(0),
	m_SimuFat(false),
	m_protonDensityWM(0.7),
	m_protonDensityGM(0.8),	
	m_protonDensityCSF(1),
	m_protonDensityFat(1)
    {	
	this->SetNumberOfRequiredInputs(0);
	this->SetNumberOfRequiredOutputs(1);	

	/* The output image type allow decimal numbers */
	if(!std::is_floating_point<typename TOutputImage::PixelType>::value)
	{
	    std::cerr << "Error: only floating points image types are supported for simulation." << std::endl;
	    exit(EXIT_FAILURE);
	}
    }

    /* Check the coherence in sequence parameters */
    template <class TInputImage, class TOutputImage>
    void
    SpinEchoSimulationImageFilter<TInputImage,TOutputImage>
    ::CheckSequenceParameters()
    {	
	if(m_EchoTime==0 || m_RepetitionTime==0)
	{
	    std::cerr << "Error: one of the following parameters has not been set: echo time, repetition time" << std::endl;
	    exit(EXIT_FAILURE);
	}

	if(m_EchoTime>=m_RepetitionTime)
	{
	    std::cerr << "Error: the echo time shouldn't be higher than the repetition time." << std::endl;
	    exit(EXIT_FAILURE);
	}

	/** T1 values from "Water proton T1 measurements in brain tissue at 7,3 and 1.5T using IR-EPI, IR-TSE, and MPRAGE: results and optimization", Wright et al. 2008
	 *  SpinEcho, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field", Marques et al. 2010
	 *  T2 values from "What are normal relaxation times of  tissues at 3T?", Bojorquez et al. 2017 **/
	if(m_Field==1.5)
	{
	    m_t1WM=0.650;
	    m_t1GM=1.2;
	    m_t1CSF=4;

	    std::cerr << "Currently not available" << std::endl;
	    exit(EXIT_FAILURE);
	}
	else if(m_Field==3)
	{
	    m_t1WM=0.850;
	    m_t1GM=1.35;
	    m_t1CSF=4;
	    m_t1Fat=0.385;

	    m_t2WM=0.075;
	    m_t2GM=0.083;
	    m_t2CSF=0.250;
	    m_t2Fat=0.103;
	}
	else if(m_Field==7)
	{
	    m_t1WM=1.050;
	    m_t1GM=1.85;
	    m_t1CSF=4;

	    std::cerr << "Currently not available" << std::endl;
	    exit(EXIT_FAILURE);
	}
	else
	{
	    std::cerr << "Error: invalid field: should be 1.5, 3 or 7T." << std::endl;
	    exit(EXIT_FAILURE);
	}
    }
        
    /** Simulation of the Spin echo sequence to generate the T1 lookup table */
    template <class TInputImage, class TOutputImage>
    float
    SpinEchoSimulationImageFilter<TInputImage, TOutputImage>
    ::signalSimulation(float _t1, float _t2, float _protonDensity)
    {	
	float E1 = exp(-1*m_RepetitionTime/_t1);
	float E2 = exp(-1*m_EchoTime/_t2);

	float seIntensity = _protonDensity * E2 * (1 - E1);

	return seIntensity;
    }    

    /* Run the simulation from tissue segmentations */
    template <class TInputImage, class TOutputImage>
    void
    SpinEchoSimulationImageFilter<TInputImage, TOutputImage>
    ::SimulationFromSegmentations()
    {
	itk::ImageRegionConstIterator<TInputImage> wmIt(this->GetWMMap(), this->GetWMMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> gmIt(this->GetGMMap(), this->GetGMMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> csfIt(this->GetCSFMap(), this->GetCSFMap()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> outputIt(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

	float wmSignal = this->signalSimulation(m_t1WM,m_t2WM,m_protonDensityWM);
	float gmSignal = this->signalSimulation(m_t1GM,m_t2GM,m_protonDensityGM);
	float csfSignal = this->signalSimulation(m_t1CSF,m_t2CSF,m_protonDensityCSF);
	
	/* Simulations */
	for(wmIt.GoToBegin(),gmIt.GoToBegin(),csfIt.GoToBegin(),outputIt.GoToBegin() ;
	    !wmIt.IsAtEnd() ;
	    ++wmIt,++gmIt,++csfIt,++outputIt)
	{
	    typename TOutputImage::PixelType wmM0 = (typename TInputImage::PixelType) wmIt.Get();
	    typename TOutputImage::PixelType gmM0 = (typename TInputImage::PixelType) gmIt.Get();
	    typename TOutputImage::PixelType csfM0 = (typename TInputImage::PixelType) csfIt.Get();

	    float seSignal=0;
	    
	    if(wmM0>0)
	    {
		float wmSignal=this->signalSimulation(m_t1WM,m_t2WM,m_protonDensityWM); 
		seSignal+=wmM0*wmSignal;
	    }

	    if(gmM0>0)
	    {
		float gmSignal=this->signalSimulation(m_t1GM,m_t2GM,m_protonDensityGM); 
		seSignal+=gmM0*gmSignal;
	    }

	    if(csfM0>0)
	    {
		float csfSignal=this->signalSimulation(m_t1CSF,m_t2CSF,m_protonDensityCSF); 
		seSignal+=csfM0*csfSignal;
	    }

	    
	    outputIt.Set(fabs(seSignal));
	}
    }

      /* Run the simulation from tissue segmentations */
    template <class TInputImage, class TOutputImage>
    void
    SpinEchoSimulationImageFilter<TInputImage, TOutputImage>
    ::SimulationFromSegmentationsWithFat()
    {
	itk::ImageRegionConstIterator<TInputImage> wmIt(this->GetWMMap(), this->GetWMMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> gmIt(this->GetGMMap(), this->GetGMMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> csfIt(this->GetCSFMap(), this->GetCSFMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> fatIt(this->GetFatMap(), this->GetFatMap()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> outputIt(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

	float wmSignal = this->signalSimulation(m_t1WM,m_t2WM,m_protonDensityWM);
	float gmSignal = this->signalSimulation(m_t1GM,m_t2GM,m_protonDensityGM);
	float csfSignal = this->signalSimulation(m_t1CSF,m_t2CSF,m_protonDensityCSF);
	float fatSignal = this->signalSimulation(m_t1Fat,m_t2Fat,m_protonDensityFat);
	
	/* Simulations */
	for(wmIt.GoToBegin(),gmIt.GoToBegin(),csfIt.GoToBegin(),fatIt.GoToBegin(),outputIt.GoToBegin() ;
	    !wmIt.IsAtEnd() ;
	    ++wmIt,++gmIt,++csfIt,++fatIt,++outputIt)
	{
	    typename TOutputImage::PixelType wmM0 = (typename TInputImage::PixelType) wmIt.Get();
	    typename TOutputImage::PixelType gmM0 = (typename TInputImage::PixelType) gmIt.Get();
	    typename TOutputImage::PixelType csfM0 = (typename TInputImage::PixelType) csfIt.Get();
	    typename TOutputImage::PixelType fatM0 = (typename TInputImage::PixelType) fatIt.Get();

	    float seSignal=0;
	    
	    if(wmM0>0)
	    {
		float wmSignal=this->signalSimulation(m_t1WM,m_t2WM,m_protonDensityWM); 
		seSignal+=wmM0*wmSignal;
	    }

	    if(gmM0>0)
	    {
		float gmSignal=this->signalSimulation(m_t1GM,m_t2GM,m_protonDensityGM); 
		seSignal+=gmM0*gmSignal;
	    }

	    if(csfM0>0)
	    {
		float csfSignal=this->signalSimulation(m_t1CSF,m_t2CSF,m_protonDensityCSF); 
		seSignal+=csfM0*csfSignal;
	    }

	    if(fatM0>0)
	    {
		float fatSignal=this->signalSimulation(m_t1Fat,m_t2Fat,m_protonDensityFat); 
		seSignal+=fatM0*fatSignal;
	    }

	    outputIt.Set(fabs(seSignal));
	}
    }
  
    /* Run the filter */
    template <class TInputImage, class TOutputImage>
    void
    SpinEchoSimulationImageFilter<TInputImage, TOutputImage>
    ::GenerateData()
    {
	/* Check the coherence in sequence parameters */
	this->CheckSequenceParameters();
	
	/* Memory allocation of the filters output and image iterator declaration */
	this->AllocateOutputs();

	if(m_SimuFat)
	{
	  this->SimulationFromSegmentationsWithFat();
	}
	else
	{
	  this->SimulationFromSegmentations();
	}
    }
}  // end namespace

#endif
