#ifndef __MPRAGESimulationImageFilter_hxx
#define __MPRAGESimulationImageFilter_hxx

#define G_PI 3.14159265358979323846 /* value of pi */

#include <math.h>
#include <random>

#include "mprageSimulation.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkAdditiveGaussianNoiseImageFilter.h"

namespace flaws
{

    /* Constructor */
    template <class TInputImage, class TOutputImage>
    MPRAGESimulationImageFilter<TInputImage, TOutputImage>
    ::MPRAGESimulationImageFilter():
	m_Alpha(0),
	m_InversionTime(0),
	m_EchoTime(0),
	m_EchoSpacing(0),
	m_RepetitionTime(0),
	m_NumberOfSlices(0),
	m_PartialFourier(0),
	m_timeGre(0),
	m_timeGreBefore(0),
	m_timeGreAfter(0),
	m_nbEx(0),
	m_nbExBefore(0),
	m_nbExAfter(0),
	m_AlphaRad(0),
	m_invEff(0.96),
	m_Field(0),
	m_NoiseLevel(0),
	m_protonDensityWM(0.7),
	m_protonDensityGM(0.8),	
	m_protonDensityCSF(1),
	m_UseB1PlusMap(false),
	m_UseB1MinusMap(false)
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
    MPRAGESimulationImageFilter<TInputImage,TOutputImage>
    ::CheckSequenceParameters()
    {
	if((this->GetWMMap()==nullptr || this->GetGMMap()==nullptr || this->GetCSFMap()==nullptr) && this->GetT1Map()==nullptr)
	{
	    std::cerr << "Error: If no T1 Map provided, wm, gm and csf segmentations needed." << std::endl;
	    exit(EXIT_FAILURE);
	}
	
	if(m_Alpha==0 || m_InversionTime==0 || m_EchoTime==0 || m_RepetitionTime==0 || m_NumberOfSlices==0 || m_PartialFourier==0)
	{
	    std::cerr << "Error: one of the following parameters has not been set: alpha, inversion time, echo time, repetition time, number of slices, partial Fourier" << std::endl;
	    exit(EXIT_FAILURE);
	}

	if(m_EchoTime>=m_EchoSpacing)
	{
	    std::cerr << "Error: the echo time shouldn't be higher than the echo spacing." << std::endl;
	    exit(EXIT_FAILURE);
	}
	
	if(m_Alpha > 15)
	{
	    std::cout << "Warning: uncommon flip angle value: alpha = " << m_Alpha << std::endl; 
	}

	if( !(m_PartialFourier==6 || m_PartialFourier==7 || m_PartialFourier==8) )
	{
	    std::cerr << "Error: partial Fourier must be set to {6,7,8}." << std::endl;
	    exit(EXIT_FAILURE);
	}

	if( ((float) m_NumberOfSlices)*m_PartialFourier/8 != std::floor(((float) m_NumberOfSlices)*m_PartialFourier/8) )
	{
	    std::cerr << "Error: number of slices * partial fourier bust be an integer." << std::endl;
	    exit(EXIT_FAILURE);
	}

	m_nbEx = ((float) m_PartialFourier)/8 * m_NumberOfSlices;
	m_nbExBefore = ((float) m_PartialFourier-4)/8 * m_NumberOfSlices;
	m_nbExAfter = ((float) m_NumberOfSlices)/2;
	
	m_timeGre = m_EchoSpacing*m_nbEx;
	m_timeGreBefore = m_EchoSpacing*m_nbExBefore;
	m_timeGreAfter = m_EchoSpacing*m_nbExAfter;

	if(m_InversionTime+m_timeGreAfter >= m_RepetitionTime)
	{
	    std::cerr << "Error: invalid parameter combination: inversion time + GRE acquisition time/2 >= repetition time." << std::endl;
	    exit(EXIT_FAILURE);
	}

	m_AlphaRad = ((float) m_Alpha) * G_PI/180;

	/** T1 values from "Water proton T1 measurements in brain tissue at 7,3 and 1.5T using IR-EPI, IR-TSE, and MPRAGE: results and optimization", Wright et al. 2008
	 *  MPRAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field", Marques et al. 2010
	 *  T2 star values from "T2 star measurements in human brain at 1.5, 3 and 7T", Peters et al. 2007 **/
	if(m_Field==1.5)
	{
	    m_t1WM=0.650;
	    m_t1GM=1.2;
	    m_t1CSF=4;

	    m_t2StarWM=0.0662;
	    m_t2StarGM=0.084;
	    m_t2StarCSF=0.150;
	}
	else if(m_Field==3)
	{
	    m_t1WM=0.850;
	    m_t1GM=1.35;
	    m_t1CSF=4;

	    m_t2StarWM=0.0532;
	    m_t2StarGM=0.066;
	    m_t2StarCSF=0.150;
	}
	else if(m_Field==7)
	{
	    m_t1WM=1.050;
	    m_t1GM=1.85;
	    m_t1CSF=4;
	    m_t1CSF=4.425;

	    m_t2StarWM=0.0268;
	    m_t2StarGM=0.0332;
	    m_t2StarCSF=0.150;
	}
	else if(this->GetT1Map()==nullptr)
	{
	    std::cerr << "Error: invalid field: should be 1.5, 3 or 7T." << std::endl;
	    exit(EXIT_FAILURE);
	}
    }
        
    /** Simulation of the MPRAGE sequence to generate the T1 lookup table */
    template <class TInputImage, class TOutputImage>
    float
    MPRAGESimulationImageFilter<TInputImage, TOutputImage>
    ::signalSimulation(float _t1, float _t2Star, float _protonDensity, float _b1Plus)
    {
	float E1 = exp(-m_EchoSpacing/_t1);
	float E2 = exp(-m_EchoTime/_t2Star);
	
	float TA = m_InversionTime - m_timeGreBefore;
	float TB = m_RepetitionTime - m_timeGre - TA;
       
	float EAB = exp((TA+TB)/_t1);
	float EB = exp(TB/_t1);

	float cosAlphaE1 = cos(m_AlphaRad*_b1Plus)*E1;
	
	float sinAlpha = sin(m_AlphaRad*_b1Plus);

	float signal = EAB*(-1+E1);
	signal += (-EAB*E1-m_invEff*E1+(1+m_invEff)*EB)*pow(cosAlphaE1,m_nbExBefore);
	signal += m_invEff*(-1+E1)*pow(cosAlphaE1,m_nbEx);
	signal += (EAB+m_invEff-(1+m_invEff)*EB)*pow(cosAlphaE1,m_nbExBefore+1);
	signal /= (-1+cosAlphaE1)*(EAB+m_invEff*pow(cosAlphaE1,m_nbEx));
	signal *= sinAlpha*E2*_protonDensity;

	return fabs(signal);
    }    

    /* Run the simulation from tissue segmentations */
    template <class TInputImage, class TOutputImage>
    void
    MPRAGESimulationImageFilter<TInputImage, TOutputImage>
    ::SimulationFromSegmentations()
    {
	itk::ImageRegionConstIterator<TInputImage> wmIt(this->GetWMMap(), this->GetWMMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> gmIt(this->GetGMMap(), this->GetGMMap()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> csfIt(this->GetCSFMap(), this->GetCSFMap()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> outputIt(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());
	
	/* Simulate the effect of B1+ inhomogeneity on the signal */
	if(m_UseB1PlusMap)
	{
	    itk::ImageRegionConstIterator<TInputImage> b1PlusIt(this->GetB1PlusMap(), this->GetB1PlusMap()->GetLargestPossibleRegion());
	    
	    if(m_UseB1MinusMap)
	    {
		itk::ImageRegionConstIterator<TInputImage> b1MinusIt(this->GetB1MinusMap(), this->GetB1MinusMap()->GetLargestPossibleRegion());

		/* Simulations */
		for(wmIt.GoToBegin(),gmIt.GoToBegin(),csfIt.GoToBegin(),b1PlusIt.GoToBegin(),b1MinusIt.GoToBegin(),outputIt.GoToBegin() ;
		    !wmIt.IsAtEnd() ;
		    ++wmIt,++gmIt,++csfIt,++b1PlusIt,++b1MinusIt,++outputIt)
		{
		    typename TOutputImage::PixelType outputSignal=0;

		    typename TOutputImage::PixelType wmM0 = (typename TInputImage::PixelType) wmIt.Get();
		    typename TOutputImage::PixelType gmM0 = (typename TInputImage::PixelType) gmIt.Get();
		    typename TOutputImage::PixelType csfM0 = (typename TInputImage::PixelType) csfIt.Get();
		    
		    unsigned int b1PlusValue = ((float) b1PlusIt.Get())/1000;
		    unsigned int b1MinusValue = ((float) b1MinusIt.Get())/1000;

		    
		    if(wmM0>0)
		    {
			outputSignal += wmM0*this->signalSimulation(m_t1WM,m_t2StarWM,m_protonDensityWM, b1PlusValue);
		    }

		    if(gmM0>0)
		    {
			outputSignal += gmM0*this->signalSimulation(m_t1GM,m_t2StarGM,m_protonDensityGM, b1PlusValue);
		    }

		    if(csfM0>0)
		    {
			outputSignal += csfM0*this->signalSimulation(m_t1CSF,m_t2StarCSF,m_protonDensityCSF, b1PlusValue);
		    }

		    outputSignal *= b1MinusValue;

		    outputIt.Set(outputSignal);
		}
	    }
	    else
	    {
                /* Simulations */
		for(wmIt.GoToBegin(),gmIt.GoToBegin(),csfIt.GoToBegin(),b1PlusIt.GoToBegin(),outputIt.GoToBegin() ;
		    !wmIt.IsAtEnd() ;
		    ++wmIt,++gmIt,++csfIt,++b1PlusIt,++outputIt)
		{
		    typename TOutputImage::PixelType outputSignal=0;

		    typename TOutputImage::PixelType wmM0 = (typename TInputImage::PixelType) wmIt.Get();
		    typename TOutputImage::PixelType gmM0 = (typename TInputImage::PixelType) gmIt.Get();
		    typename TOutputImage::PixelType csfM0 = (typename TInputImage::PixelType) csfIt.Get();
		    
		    unsigned int b1PlusValue = ((float) b1PlusIt.Get())/1000;
		    
		    if(wmM0>0)
		    {
			outputSignal += wmM0*this->signalSimulation(m_t1WM,m_t2StarWM,m_protonDensityWM, b1PlusValue);
		    }

		    if(gmM0>0)
		    {
			outputSignal += gmM0*this->signalSimulation(m_t1GM,m_t2StarGM,m_protonDensityGM, b1PlusValue);
		    }

		    if(csfM0>0)
		    {
			outputSignal += csfM0*this->signalSimulation(m_t1CSF,m_t2StarCSF,m_protonDensityCSF, b1PlusValue);
		    }

		    outputIt.Set(outputSignal);
		}
	    }
	}
	else if(m_UseB1MinusMap)
	{
	    /* Pure signal simulation */
	    float wmSignal = this->signalSimulation(m_t1WM,m_t2StarWM,m_protonDensityWM, 1);
	    float gmSignal = this->signalSimulation(m_t1GM,m_t2StarGM,m_protonDensityGM, 1);
	    float csfSignal = this->signalSimulation(m_t1CSF,m_t2StarCSF,m_protonDensityCSF, 1);

	    itk::ImageRegionConstIterator<TInputImage> b1MinusIt(this->GetB1MinusMap(), this->GetB1MinusMap()->GetLargestPossibleRegion());
	    
	    /* Simulations */
	    for(wmIt.GoToBegin(),gmIt.GoToBegin(),csfIt.GoToBegin(),b1MinusIt.GoToBegin(),outputIt.GoToBegin() ;
		!wmIt.IsAtEnd() ;
		++wmIt,++gmIt,++csfIt,++b1MinusIt,++outputIt)
	    {
		typename TOutputImage::PixelType outputSignal=0;

		typename TOutputImage::PixelType wmM0 = (typename TInputImage::PixelType) wmIt.Get();
		typename TOutputImage::PixelType gmM0 = (typename TInputImage::PixelType) gmIt.Get();
		typename TOutputImage::PixelType csfM0 = (typename TInputImage::PixelType) csfIt.Get();
		    
		unsigned int b1MinusValue = ((float) b1MinusIt.Get())/1000;

		if(wmM0>0)
		{
		    outputSignal += wmSignal*wmM0;
		}

		if(gmM0>0)
		{
		    outputSignal += gmSignal*gmM0;
		}

		if(csfM0>0)
		{
		    outputSignal += csfSignal*csfM0;
		}

		outputSignal *= b1MinusValue;

		outputIt.Set(outputSignal);
	    }
	}
	/* The effect of B1+ inhomogeneity is not taken into account to simulate the signal */
	else
	{
	    /* Pure signal simulation */
	    float wmSignal = this->signalSimulation(m_t1WM,m_t2StarWM,m_protonDensityWM, 1);
	    float gmSignal = this->signalSimulation(m_t1GM,m_t2StarGM,m_protonDensityGM, 1);
	    float csfSignal = this->signalSimulation(m_t1CSF,m_t2StarCSF,m_protonDensityCSF, 1);
    
	    /* Simulations */
	    for(wmIt.GoToBegin(),gmIt.GoToBegin(),csfIt.GoToBegin(),outputIt.GoToBegin() ;
		!wmIt.IsAtEnd() ;
		++wmIt,++gmIt,++csfIt,++outputIt)
	    {
		typename TOutputImage::PixelType outputSignal=0;

		typename TOutputImage::PixelType wmM0 = (typename TInputImage::PixelType) wmIt.Get();
		typename TOutputImage::PixelType gmM0 = (typename TInputImage::PixelType) gmIt.Get();
		typename TOutputImage::PixelType csfM0 = (typename TInputImage::PixelType) csfIt.Get();		    

		// if(wmM0>0)
		// {
		//     outputSignal += wmSignal*wmM0;
		// }

		// if(gmM0>0)
		// {
		//     outputSignal += gmSignal*gmM0;
		// }

		// if(csfM0>0)
		// {
		//     outputSignal += csfSignal*csfM0;
		// }

		float t1Signal = 1/(wmM0/m_t1WM + gmM0/m_t1GM + csfM0/m_t1CSF);

		float t2StarSignal = 1/(wmM0/m_t2StarWM + gmM0/m_t2StarGM + csfM0/m_t2StarCSF);

		float protonDensitySignal = wmM0*m_protonDensityWM + gmM0*m_protonDensityGM + csfM0*m_protonDensityCSF;	    
	    
		outputSignal=this->signalSimulation(t1Signal,t2StarSignal,protonDensitySignal,1);	    

		outputIt.Set(outputSignal);
	    }
	}

	typedef typename itk::AdditiveGaussianNoiseImageFilter <TOutputImage,TOutputImage> AdditiveGaussianNoiseFilterType;
	typename AdditiveGaussianNoiseFilterType::Pointer gaussianNoiseFilter = AdditiveGaussianNoiseFilterType::New();

	if(m_NoiseLevel>0)
	{
	    gaussianNoiseFilter->SetInput(this->GetOutput());
	    gaussianNoiseFilter->SetMean(0);

	    float wmSignal = this->signalSimulation(m_t1WM,m_t2StarWM,m_protonDensityWM, 1);
	    
	    gaussianNoiseFilter->SetStandardDeviation(wmSignal*m_NoiseLevel/100);

	    try
	    {
		gaussianNoiseFilter->Update();
	    }
	    catch(itk::ExceptionObject &e)
	    {
		std::cerr << e << std::endl;
		exit(EXIT_FAILURE);
	    }
	    
	    itk::ImageRegionIterator<TOutputImage> gaussianNoiseIt(gaussianNoiseFilter->GetOutput(), gaussianNoiseFilter->GetOutput()->GetLargestPossibleRegion());

	    for(gaussianNoiseIt.GoToBegin(),outputIt.GoToBegin() ;
		!gaussianNoiseIt.IsAtEnd() ;
		++gaussianNoiseIt,++outputIt)
	    {
		float signal = fabs(gaussianNoiseIt.Get());

		outputIt.Set(signal);
	    }
	}
    }
    
    /* Run the simulation from T1 Map */
    template <class TInputImage, class TOutputImage>
    void
    MPRAGESimulationImageFilter<TInputImage, TOutputImage>
    ::SimulationFromT1Map()
    {
	itk::ImageRegionConstIterator<TInputImage> t1It(this->GetT1Map(), this->GetT1Map()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> outputIt(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());
	
	/* Simulate the effect of B1+ inhomogeneity on the signal */
	if(m_UseB1PlusMap)
	{
	    itk::ImageRegionConstIterator<TInputImage> b1PlusIt(this->GetB1PlusMap(), this->GetB1PlusMap()->GetLargestPossibleRegion());

	    if(m_UseB1MinusMap)
	    {
		itk::ImageRegionConstIterator<TInputImage> b1MinusIt(this->GetB1MinusMap(), this->GetB1MinusMap()->GetLargestPossibleRegion());
		
		/* Simulations */
		for(t1It.GoToBegin(),b1PlusIt.GoToBegin(),b1MinusIt.GoToBegin(),outputIt.GoToBegin() ;
		    !t1It.IsAtEnd() ;
		    ++t1It,++b1PlusIt,++b1MinusIt,++outputIt)
		{
		    float t1Value = ((float) t1It.Get())/1000;
		    float pdValue = 1;

		    if(m_Field==3) // Estimate the PD values as in Volz et al. "Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities"
		    {
			float k1 = 0.858;
			float k2 = 522;
			pdValue = t1Value/(k1*t1Value+k2);
		    }
		    
		    float b1PlusValue = ((float) b1PlusIt.Get())/1000;
		    float b1MinusValue = ((float) b1MinusIt.Get())/1000;
		    
		    float signal = this->signalSimulation(t1Value,100000,pdValue,b1PlusValue);
		    signal *= b1MinusValue;

		    outputIt.Set(signal);
		}
	    }
	    else
	    {
		/* Simulations */
		for(t1It.GoToBegin(),b1PlusIt.GoToBegin(),outputIt.GoToBegin() ;
		    !t1It.IsAtEnd() ;
		    ++t1It,++b1PlusIt,++outputIt)
		{
		    float t1Value = ((float) t1It.Get())/1000;
		    float pdValue = 1;

		    if(m_Field==3) // Estimate the PD values as in Volz et al. "Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities"
		    {
			float k1 = 0.858;
			float k2 = 522;
			pdValue = t1Value/(k1*t1Value+k2);
		    }		    

		    float b1PlusValue = ((float) b1PlusIt.Get())/1000;
		    
		    float signal = this->signalSimulation(t1Value,100000,pdValue,b1PlusValue);		   

		    outputIt.Set(signal);
		}
	    }
	}
	else if(m_UseB1MinusMap)
	{
	    itk::ImageRegionConstIterator<TInputImage> b1MinusIt(this->GetB1MinusMap(), this->GetB1MinusMap()->GetLargestPossibleRegion());
		
	    /* Simulations */
	    for(t1It.GoToBegin(),b1MinusIt.GoToBegin(),outputIt.GoToBegin() ;
		!t1It.IsAtEnd() ;
		++t1It,++b1MinusIt,++outputIt)
	    {
		float t1Value = ((float) t1It.Get())/1000;
		float pdValue = 1;

		if(m_Field==3) // Estimate the PD values as in Volz et al. "Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities"
		{
		    float k1 = 0.858;
		    float k2 = 522;
		    pdValue = t1Value/(k1*t1Value+k2);
		}
		
		float b1MinusValue = ((float) b1MinusIt.Get())/1000;
		    
		float signal = this->signalSimulation(t1Value,100000,pdValue,1);
		signal *= b1MinusValue;

		outputIt.Set(signal);
	    }
	}
	/* The effect of B1+ inhomogeneity is not taken into account to simulate the signal */
	else
	{
	    /* Simulations */
	    for(t1It.GoToBegin(),outputIt.GoToBegin() ;
		!t1It.IsAtEnd() ;
		++t1It,++outputIt)
	    {
		float t1Value = ((float) t1It.Get())/1000;
		float pdValue = 1;

		if(m_Field==3) // Estimate the PD values as in Volz et al. "Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities"
		{
		    float k1 = 0.858;
		    float k2 = 522;
		    pdValue = t1Value/(k1*t1Value+k2);
		}
		    
		float signal = this->signalSimulation(t1Value,100000,pdValue,1);
		outputIt.Set(signal);
	    }
	}
    }
    
    /* Run the filter */
    template <class TInputImage, class TOutputImage>
    void
    MPRAGESimulationImageFilter<TInputImage, TOutputImage>
    ::GenerateData()
    {
	/* Check the coherence in sequence parameters */
	this->CheckSequenceParameters();
	
	/* Memory allocation of the filters output and image iterator declaration */
	this->AllocateOutputs();

	if(this->GetT1Map()==nullptr)
	{
	    this->SimulationFromSegmentations();
	}
	else
	{
	    this->SimulationFromT1Map();
	}

    }
}  // end namespace

#endif
