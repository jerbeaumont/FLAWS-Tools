#ifndef __FLAWSProcessingImageFilter_hxx
#define __FLAWSProcessingImageFilter_hxx

#include <math.h>

#include "flawsProcessing.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <fstream>
#include <iostream>

namespace flaws
{
    /* Constructor */
    template <class TInputImage, class TOutputImage>
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::FLAWSProcessingImageFilter():
	m_DenoisingCoeff(0),
	m_Alpha1(0),
	m_Alpha2(0),	
	m_InversionTime1(0),
	m_InversionTime2(0),
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
	m_Alpha1Rad(0),
	m_Alpha2Rad(0),
	m_invEff(1),
	m_Alpha1b(0),
	m_Alpha2b(0),	
	m_InversionTime1b(0),
	m_InversionTime2b(0),
	m_EchoSpacingb(0),
	m_RepetitionTimeb(0),
	m_NumberOfSlicesb(0),
	m_PartialFourierb(0),
	m_timeGreb(0),
	m_timeGreBeforeb(0),
	m_timeGreAfterb(0),
	m_nbExb(0),
	m_nbExBeforeb(0),
	m_nbExAfterb(0),
	m_Alpha1Radb(0),
	m_Alpha2Radb(0)
    {	
	this->SetNumberOfRequiredInputs(3);
	this->SetNumberOfRequiredOutputs(13);

	this->SetNthOutput(1,this->MakeOutput(1));
	this->SetNthOutput(2,this->MakeOutput(2));
	this->SetNthOutput(3,this->MakeOutput(3));
	this->SetNthOutput(4,this->MakeOutput(4));
	this->SetNthOutput(5,this->MakeOutput(5));
	this->SetNthOutput(6,this->MakeOutput(6));
	this->SetNthOutput(7,this->MakeOutput(7));
	this->SetNthOutput(8,this->MakeOutput(8));
	this->SetNthOutput(9,this->MakeOutput(9));
	this->SetNthOutput(10,this->MakeOutput(10));
	this->SetNthOutput(11,this->MakeOutput(11));
	this->SetNthOutput(12,this->MakeOutput(12));

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
    FLAWSProcessingImageFilter<TInputImage,TOutputImage>
    ::CheckSequenceParameters()
    {	
	if(m_Alpha1==0 || m_Alpha2==0 || m_InversionTime1==0 || m_InversionTime2==0 || m_EchoTime==0 || m_EchoSpacing==0 || m_RepetitionTime==0 || m_NumberOfSlices==0 || m_PartialFourier==0)
	{
	    std::cerr << "Error: one of the following parameters has not been set: lambda, alpha1, alpha2, inversion time 1, inversion time 2, echo time, echo spacing, repetition time, number of slices, partial Fourier" << std::endl;
	    exit(EXIT_FAILURE);
	}
	
	if(m_Alpha1 > 15)
	{
	    std::cout << "Warning: uncommon flip angle value: alpha1 = " << m_Alpha1 << std::endl; 
	}

	if(m_Alpha2 > 15)
	{
	    std::cout << "Warning: uncommon flip angle value: alpha2 = " << m_Alpha2 << std::endl; 
	}

	if( !(m_PartialFourier==6 || m_PartialFourier==7 || m_PartialFourier==8) )
	{
	    std::cerr << "Error: partial Fourier must be set to {6,7,8}." << std::endl;
	    exit(EXIT_FAILURE);
	}

	if( ((double) m_NumberOfSlices)*m_PartialFourier/8 != std::floor(((double) m_NumberOfSlices)*m_PartialFourier/8) )
	{
	    std::cerr << "Error: number of slices * partial fourier must be an integer." << std::endl;
	    exit(EXIT_FAILURE);
	}

	m_nbEx = ((double) m_PartialFourier)/8 * m_NumberOfSlices;
	m_nbExBefore = ((double) m_PartialFourier-4)/8 * m_NumberOfSlices;
	m_nbExAfter = ((double) m_NumberOfSlices)/2;
	
	m_timeGre = m_EchoSpacing*m_nbEx;
	m_timeGreBefore = m_EchoSpacing*m_nbExBefore;
	m_timeGreAfter = m_EchoSpacing*m_nbExAfter;

	if(m_InversionTime1+m_timeGre >= m_InversionTime2)
	{
	    std::cerr << "Error: invalid parameter combination: inversion time 1 + GRE acquisition time >= inversion time 2." << std::endl;
	    exit(EXIT_FAILURE);
	}

	if(m_InversionTime2+m_timeGreAfter >= m_RepetitionTime)
	{
	    std::cerr << "Error: invalid parameter combination: inversion time 2 + GRE acquisition time/2 >= repetition time." << std::endl;
	    exit(EXIT_FAILURE);
	}

	m_Alpha1Rad = ((double) m_Alpha1) * G_PI/180;
	m_Alpha2Rad = ((double) m_Alpha2) * G_PI/180;

	if(this->GetB1PlusMap())
	{
	    if(m_Alpha1b==0 || m_Alpha2b==0 || m_InversionTime1b==0 || m_InversionTime2b==0 || m_EchoSpacingb==0 || m_RepetitionTimeb==0 || m_NumberOfSlicesb==0 || m_PartialFourierb==0)
	    {
		std::cerr << "Error: one of the following parameters has not been set: alpha1b, alpha2b, inversion time 1b, inversion time 2b, echo spacingb, repetition timeb, number of slicesb, partial Fourierb" << std::endl;
		exit(EXIT_FAILURE);
	    }

	    if(m_Alpha1b > 15)
	    {
		std::cout << "Warning: uncommon flip angle value: alpha1 b = " << m_Alpha1 << std::endl; 
	    }

	    if(m_Alpha2b > 15)
	    {
		std::cout << "Warning: uncommon flip angle value: alpha2 b = " << m_Alpha1 << std::endl; 
	    }

	    if( !(m_PartialFourierb==6 || m_PartialFourierb==7 || m_PartialFourierb==8) )
	    {
		std::cerr << "Error: partial Fourier b must be set to {6,7,8}." << std::endl;
		exit(EXIT_FAILURE);
	    }

	    if( ((double) m_NumberOfSlicesb)*m_PartialFourierb/8 != std::floor(((double) m_NumberOfSlicesb)*m_PartialFourierb/8) )
	    {
		std::cerr << "Error: number of slices b * partial fourier b bust be an integer." << std::endl;
		exit(EXIT_FAILURE);
	    }

	    m_nbExb = ((double) m_PartialFourierb)/8 * m_NumberOfSlicesb;
	    m_nbExBeforeb = ((double) m_PartialFourierb-4)/8 * m_NumberOfSlicesb;
	    m_nbExAfterb = ((double) m_NumberOfSlicesb)/2;
	
	    m_timeGreb = m_EchoSpacingb*m_nbExb;
	    m_timeGreBeforeb = m_EchoSpacingb*m_nbExBeforeb;
	    m_timeGreAfterb = m_EchoSpacingb*m_nbExAfterb;

	    if(m_InversionTime1b+m_timeGreb >= m_InversionTime2b)
	    {
		std::cerr << "Error: invalid parameter combination: inversion time 1 b + GRE acquisition time b >= inversion time 2 b." << std::endl;
		exit(EXIT_FAILURE);
	    }

	    if(m_InversionTime2b+m_timeGreAfterb >= m_RepetitionTimeb)
	    {
		std::cerr << "Error: invalid parameter combination: inversion time 2 b + GRE acquisition time b/2 >= repetition time." << std::endl;
		exit(EXIT_FAILURE);
	    }

	    m_Alpha1Radb = ((double) m_Alpha1b) * G_PI/180;
	    m_Alpha2Radb = ((double) m_Alpha2b) * G_PI/180;
	}
    }
    
    /** Simulation of the FLAWS sequence to generate the T1 lookup table */
    template <class TInputImage, class TOutputImage>
    void
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::T1Simulation(double _t1, double _b1, double* _signals)
    {
	double E1 = exp(-m_EchoSpacing/_t1);

	double TA = m_InversionTime1 - m_timeGreBefore;
	double TB = m_InversionTime2 - m_InversionTime1-m_timeGre;
	double TC = m_RepetitionTime - 2*m_timeGre - TA - TB;

	double EA = exp(-TA/_t1);
	double EB = exp(-TB/_t1);
	double EC = exp(-TC/_t1);
	
	double cosAlpha1E1 = cos(m_Alpha1Rad*_b1)*E1;
	double cosAlpha2E1 = cos(m_Alpha2Rad*_b1)*E1;
	
	double sinAlpha1 = sin(m_Alpha1Rad*_b1);
	double sinAlpha2 = sin(m_Alpha2Rad*_b1);
	
	double mzSteadyStateDen = 1+m_invEff*pow(cosAlpha1E1*cosAlpha2E1,m_nbEx)*EA*EB*EC;
	
	double mzSteadyStateNum = (1-EA);

	mzSteadyStateNum=mzSteadyStateNum*pow(cosAlpha1E1,m_nbEx) + (1-E1)*(1-pow(cosAlpha1E1,m_nbEx))/(1-cosAlpha1E1);
	mzSteadyStateNum=mzSteadyStateNum*EB+(1-EB);

	mzSteadyStateNum=mzSteadyStateNum*pow(cosAlpha2E1,m_nbEx) + (1-E1)*(1-pow(cosAlpha2E1,m_nbEx))/(1-cosAlpha2E1);
	mzSteadyStateNum=mzSteadyStateNum*EC+(1-EC);

	double mzSteadyState = mzSteadyStateNum/mzSteadyStateDen;
	
	double temp = (-m_invEff*mzSteadyState*EA+(1-EA))*pow(cosAlpha1E1,m_nbExBefore)+(1-E1)*(1-pow(cosAlpha1E1,m_nbExBefore))/(1-cosAlpha1E1);

	*(_signals) = temp*sinAlpha1;
	
	temp = temp*pow(cosAlpha1E1,m_nbExAfter)+(1-E1)*(1-pow(cosAlpha1E1,m_nbExAfter))/(1-cosAlpha1E1);
	temp = (temp*EB+(1-EB))*pow(cosAlpha2E1,m_nbExBefore)+(1-E1)*(1-pow(cosAlpha2E1,m_nbExBefore))/(1-cosAlpha2E1);

	*(_signals+1) = temp*sinAlpha2;	
    }

    /** Simulation of the SA2RAGE sequence to generate the B1 lookup table */
    template <class TInputImage, class TOutputImage>
    double
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::B1Simulation(double _t1, double _b1)
    {
	double invEff = -cos(_b1*G_PI/2);
	
	double E1 = exp(-m_EchoSpacingb/_t1);

	double TA = m_InversionTime1b - m_timeGreBeforeb;
	double TB = m_InversionTime2b - m_InversionTime1b-m_timeGreb;
	double TC = m_RepetitionTimeb - 2*m_timeGreb - TA - TB;

	double EA = exp(-TA/_t1);
	double EB = exp(-TB/_t1);
	double EC = exp(-TC/_t1);
	
	double cosAlpha1E1 = cos(m_Alpha1Radb*_b1)*E1;
	double cosAlpha2E1 = cos(m_Alpha2Radb*_b1)*E1;
	
	double sinAlpha1 = sin(m_Alpha1Radb*_b1);
	double sinAlpha2 = sin(m_Alpha2Radb*_b1);
	
	double mzSteadyStateDen = 1+invEff*pow(cosAlpha1E1*cosAlpha2E1,m_nbExb)*EA*EB*EC;
	
	double mzSteadyStateNum = (1-EA);

	mzSteadyStateNum=mzSteadyStateNum*pow(cosAlpha1E1,m_nbExb) + (1-E1)*(1-pow(cosAlpha1E1,m_nbExb))/(1-cosAlpha1E1);
	mzSteadyStateNum=mzSteadyStateNum*EB+(1-EB);

	mzSteadyStateNum=mzSteadyStateNum*pow(cosAlpha2E1,m_nbExb) + (1-E1)*(1-pow(cosAlpha2E1,m_nbExb))/(1-cosAlpha2E1);
	mzSteadyStateNum=mzSteadyStateNum*EC+(1-EC);

	double mzSteadyState = mzSteadyStateNum/mzSteadyStateDen;
	
	double temp = (-invEff*mzSteadyState*EA+(1-EA))*pow(cosAlpha1E1,m_nbExBeforeb)+(1-E1)*(1-pow(cosAlpha1E1,m_nbExBeforeb))/(1-cosAlpha1E1);

	double signalGre1 = temp*sinAlpha1;
	
	temp = temp*pow(cosAlpha1E1,m_nbExAfterb)+(1-E1)*(1-pow(cosAlpha1E1,m_nbExAfterb))/(1-cosAlpha1E1);
	temp = (temp*EB+(1-EB))*pow(cosAlpha2E1,m_nbExBeforeb)+(1-E1)*(1-pow(cosAlpha2E1,m_nbExBeforeb))/(1-cosAlpha2E1);

	double signalGre2 = temp*sinAlpha2;
	
	return signalGre1/signalGre2;
    }
    
    /** Estimation of the T1 from FLAWS-HC and FLAWS lookup table */
    template <class TInputImage, class TOutputImage>
    unsigned long
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::T1Estimation(unsigned long _minT1, unsigned long _t1Range, unsigned long _midT1Range, double* _t1LookupTable, double _combiValue)
    {
	unsigned long lowLookupIndex = _minT1;
	unsigned long highLookupIndex = _t1Range+_minT1-1;
	unsigned long lookupIndex = _midT1Range;

	unsigned long t1Estimate=0;

	if(_t1LookupTable[lowLookupIndex]<_t1LookupTable[highLookupIndex])
	{
	    while( (highLookupIndex-lowLookupIndex) > 1 )
	    {
		if(_combiValue <= _t1LookupTable[lowLookupIndex])
		{
		    return lowLookupIndex+1;
		}
		else if(_combiValue >= _t1LookupTable[highLookupIndex])
		{
		    return highLookupIndex+1;
		}	    
		else if(_combiValue < _t1LookupTable[lookupIndex])
		{
		    highLookupIndex=lookupIndex;
		}
		else
		{
		    lowLookupIndex=lookupIndex;
		}

		lookupIndex=lowLookupIndex+(highLookupIndex-lowLookupIndex)/2;
	    }
	}
	else
	{
	    while( (highLookupIndex-lowLookupIndex) > 1 )
	    {
		if(_combiValue >= _t1LookupTable[lowLookupIndex])
		{
		    return lowLookupIndex+1;
		}
		else if(_combiValue <= _t1LookupTable[highLookupIndex])
		{
		    return highLookupIndex+1;
		}	    
		else if(_combiValue < _t1LookupTable[lookupIndex])
		{
		    lowLookupIndex=lookupIndex;
		}
		else
		{
		    highLookupIndex=lookupIndex;
		}

		lookupIndex=lowLookupIndex+(highLookupIndex-lowLookupIndex)/2;
	    }
	}
	
	if(fabs(_combiValue-_t1LookupTable[highLookupIndex]) < fabs(_combiValue-_t1LookupTable[lookupIndex]))
	{
	    t1Estimate=highLookupIndex+1;
	}
	else if(fabs(_combiValue-_t1LookupTable[lowLookupIndex]) < fabs(_combiValue-_t1LookupTable[lookupIndex]))
	{
	    t1Estimate=lowLookupIndex+1;
	}
	else
	{
	    t1Estimate=lookupIndex+1;
	}

	return t1Estimate;
    }

    /** Estimation of the T1 from FLAWS-HC and FLAWS lookup table */
    template <class TInputImage, class TOutputImage>
    unsigned long
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::T1Estimation(unsigned long _minT1, unsigned long _t1Range, unsigned long _midT1Range, unsigned long _b1Range, double* _t1B1LookupTable, double _hcValue, unsigned long _b1Value)
    {
	unsigned long int lowT1 = _minT1;
	unsigned long int highT1 = _t1Range+_minT1-1;
	unsigned long int currentT1 = _midT1Range;

	unsigned long int lowLookupIndex = lowT1*_b1Range+_b1Value-1;
	unsigned long int highLookupIndex = highT1*_b1Range+_b1Value-1;
	unsigned long int lookupIndex = currentT1*_b1Range+_b1Value-1;
	
	unsigned long int t1Estimate=0;
	
	if(_t1B1LookupTable[lowLookupIndex]<_t1B1LookupTable[highLookupIndex])
	{
	    while( (highT1-lowT1) > 1 )
	    {
		if(_hcValue <= _t1B1LookupTable[lowLookupIndex])
		{		    
		    return lowT1+1;
		}
		else if(_hcValue >= _t1B1LookupTable[highLookupIndex])
		{
		    return highT1+1;
		}	    
		else if(_hcValue < _t1B1LookupTable[lookupIndex])
		{
		    highT1=currentT1;
		    highLookupIndex=lookupIndex;
		}
		else
		{
		    lowT1=currentT1;
		    lowLookupIndex=lookupIndex;
		}

		currentT1 = lowT1+(highT1-lowT1)/2;
		lookupIndex = currentT1*_b1Range+_b1Value-1;
	    }
	}
	else
	{
	    while( (highT1-lowT1) > 1 )
	    {
		if(_hcValue >= _t1B1LookupTable[lowLookupIndex])
		{
		    return lowT1+1;
		}
		else if(_hcValue <= _t1B1LookupTable[highLookupIndex])
		{
		    return highT1+1;
		}	    
		else if(_hcValue < _t1B1LookupTable[lookupIndex])
		{
		    lowT1=currentT1;
		    lowLookupIndex=lookupIndex;
		}
		else
		{
		    highT1=currentT1;
		    highLookupIndex=lookupIndex;
		}

		currentT1 = lowT1+(highT1-lowT1)/2;
		lookupIndex = currentT1*_b1Range+_b1Value-1;
	    }
	}

	if(fabs(_hcValue-_t1B1LookupTable[highLookupIndex]) < fabs(_hcValue-_t1B1LookupTable[lookupIndex]))
	{
	    t1Estimate=highT1+1;
	}
	else if(fabs(_hcValue-_t1B1LookupTable[lowLookupIndex]) < fabs(_hcValue-_t1B1LookupTable[lookupIndex]))
	{
	    t1Estimate=lowT1+1;
	}
	else
	{
	    t1Estimate=currentT1+1;
	}

	return t1Estimate;
    }    
    
    /** Estimation of the B1+ from the SA2RAGE division image and lookup tables */
    template <class TInputImage, class TOutputImage>
    unsigned long
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::B1Estimation(unsigned long _b1Range, unsigned long _midB1Range, double* _b1T1LookupTable, double _sa2rageValue, unsigned long _t1Value)
    {
	unsigned long int t1Index = _t1Value-1;
	
	unsigned long int lowB1 = 0;
	unsigned long int highB1 = _b1Range-1;
	unsigned long int currentB1 = _midB1Range;
	
	unsigned long int lowLookupIndex = t1Index*_b1Range+lowB1;
	unsigned long int highLookupIndex = t1Index*_b1Range+highB1;
	unsigned long int lookupIndex = t1Index*_b1Range+currentB1;

	unsigned long int b1Estimate=0;

	if(_b1T1LookupTable[lowLookupIndex]<_b1T1LookupTable[highLookupIndex])
	{
	    while( (highB1-lowB1) > 1 )
	    {
		if(_sa2rageValue >= _b1T1LookupTable[highLookupIndex]) 
		{
		    return highB1+1;
		}
		else if(_sa2rageValue <= _b1T1LookupTable[lowLookupIndex])
		{		    
		    return lowB1+1;
		}	    
		else if(_sa2rageValue > _b1T1LookupTable[lookupIndex])
		{
		    lowB1=currentB1;
		    lowLookupIndex=lookupIndex;
		}
		else
		{
		    highB1=currentB1;
		    highLookupIndex=lookupIndex;
		}

		currentB1 = lowB1+(highB1-lowB1)/2;
		lookupIndex = t1Index*_b1Range+currentB1;		
	    }
	}
	else
	{
	    while( (highB1-lowB1) > 1 )
	    {
		if(_sa2rageValue >= _b1T1LookupTable[lowLookupIndex]) 
		{
		    return lowB1+1;
		}
		else if(_sa2rageValue <= _b1T1LookupTable[highLookupIndex])
		{
		    return highB1+1;
		}	    
		else if(_sa2rageValue > _b1T1LookupTable[lookupIndex])
		{
		    highB1 = currentB1;
		    highLookupIndex = lookupIndex;
		}
		else
		{
		    lowB1 = currentB1;
		    lowLookupIndex=lookupIndex;
		}

		currentB1 = lowB1+(highB1-lowB1)/2;
		lookupIndex = t1Index*_b1Range+currentB1;		
	    }
	}
	
	if(fabs(_sa2rageValue-_b1T1LookupTable[highLookupIndex]) < fabs(_sa2rageValue-_b1T1LookupTable[lookupIndex]))
	{
	    b1Estimate=highB1+1;
	}
	else if(fabs(_sa2rageValue-_b1T1LookupTable[lowLookupIndex]) < fabs(_sa2rageValue-_b1T1LookupTable[lookupIndex]))
	{
	    b1Estimate=lowB1+1;
	}
	else
	{
	    b1Estimate=currentB1+1;
	}

	return b1Estimate;
    }    
    
    /* T1 Estimation */
    /* Get back the FLAWS1 and FLAWS2 signals signs */
    /* Generate the combined FLAWS images */
    template <class TInputImage, class TOutputImage>
    void
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::T1Mapping()
    {
	/* Iterator declaration */
	itk::ImageRegionConstIterator<TInputImage> inv1It(this->GetInv1(), this->GetInv1()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> inv2It(this->GetInv2(), this->GetInv2()->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<TInputImage> uniIt(this->GetUNI(), this->GetUNI()->GetLargestPossibleRegion());

	itk::ImageRegionIterator<TOutputImage> inv1SignedIt(this->GetInv1Signed(), this->GetInv1Signed()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> inv2SignedIt(this->GetInv2Signed(), this->GetInv2Signed()->GetLargestPossibleRegion());

	itk::ImageRegionIterator<TOutputImage> flawsHCIt(this->GetFLAWSHC(), this->GetFLAWSHC()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> flawsHCOIt(this->GetFLAWSHCO(), this->GetFLAWSHCO()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> flawsHCDenIt(this->GetFLAWSHCDen(), this->GetFLAWSHCDen()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> flawsHCODenIt(this->GetFLAWSHCODen(), this->GetFLAWSHCODen()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> minIt(this->GetMin(), this->GetMin()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> mipIt(this->GetMip(), this->GetMip()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> mipDenIt(this->GetMipDen(), this->GetMipDen()->GetLargestPossibleRegion());
	
	itk::ImageRegionIterator<TOutputImage> t1MapIt(this->GetT1Map(), this->GetT1Map()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> flawsHCB1PlusCorrIt(this->GetFLAWSHCB1PlusCorr(),this->GetFLAWSHCB1PlusCorr()->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TOutputImage> flawsHCOB1PlusCorrIt(this->GetFLAWSHCOB1PlusCorr(),this->GetFLAWSHCOB1PlusCorr()->GetLargestPossibleRegion());

        /* LookupTable generation and T1 estimation */
	unsigned long b1Range = 2000;
	unsigned long midB1Range = b1Range/2;
	unsigned long t1Range = 5000;
	unsigned long midT1Range = t1Range/2;

	double *t1LookupTable = new double[t1Range*2];
	double *t1HCLookupTable = new double[t1Range];
	double *t1UNILookupTable = new double[t1Range];
		
	double *t1B1LookupTable = new double[t1Range*2*b1Range];
	double *t1B1HCLookupTable = new double[t1Range*b1Range];
	double *t1B1UNILookupTable = new double[t1Range*b1Range];
	double *b1T1LookupTable = new double[b1Range*t1Range];	

	if(this->GetB1PlusMap())
	{
	    unsigned long *t1Zero = new unsigned long[b1Range];
	    unsigned long *t1Min = new unsigned long[b1Range];
	    unsigned long *t1Max = new unsigned long[b1Range];

	    unsigned long *t1RangeBefore = new unsigned long[b1Range];	       
	    unsigned long *t1RangeBetween = new unsigned long[b1Range];
	    unsigned long *t1RangeAfter = new unsigned long[b1Range];
		
	    unsigned long *midT1RangeBefore = new unsigned long[b1Range];
	    unsigned long *midT1RangeBetween = new unsigned long[b1Range];
	    unsigned long *midT1RangeAfter = new unsigned long[b1Range];		

	    double min=-1,max=1;
	    double distMin=100,distMax=100;
		
	    for(unsigned long b1=1 ; b1<=b1Range ; b1++)
	    {
		for(unsigned long t1=1 ; t1<=t1Range ; t1++)
		{
		    this->T1Simulation(((double) t1)/1000, ((double) b1)/1000, t1B1LookupTable+(((t1-1)*b1Range+b1-1)*2));
		    
		    b1T1LookupTable[(t1-1)*b1Range+b1-1]=this->B1Simulation(((double) t1)/1000, ((double) b1)/1000);

		    double flaws1 = t1B1LookupTable[((t1-1)*b1Range+b1-1)*2];
		    double flaws2 = t1B1LookupTable[((t1-1)*b1Range+b1-1)*2+1];

		    double hcValue = (fabs(flaws1)-fabs(flaws2))/(fabs(flaws1)+fabs(flaws2));		    
		    t1B1HCLookupTable[(t1-1)*b1Range+b1-1] = hcValue;
			
		    if(fabs(hcValue-min)<distMin)
		    {
			distMin = fabs(hcValue-min);
			t1Min[b1-1] = t1-1;
		    }
		    else if(fabs(hcValue-max)<distMax)
		    {
			distMax = fabs(hcValue-max);
			t1Max[b1-1] = t1-1;
		    }

		    double uniValue = (flaws1*flaws2)/(flaws1*flaws1+flaws2*flaws2);
		    t1B1UNILookupTable[(t1-1)*b1Range+b1-1] = uniValue;
		}		    
		
		distMin=100;
		distMax=100;

		t1RangeBefore[b1-1]=t1Min[b1-1];
		t1RangeBetween[b1-1]=t1Max[b1-1]-t1Min[b1-1];
		t1RangeAfter[b1-1]=t1Range-t1Max[b1-1];

		midT1RangeBefore[b1-1]=t1RangeBefore[b1-1]/2;
		midT1RangeBetween[b1-1]=t1RangeBetween[b1-1]/2+t1Min[b1-1];
		midT1RangeAfter[b1-1]=t1RangeAfter[b1-1]/2+t1Max[b1-1];
	    }
	    
	    /* Iterator declaration */
	    itk::ImageRegionConstIterator<TInputImage> b1MapIt(this->GetB1PlusMap(), this->GetB1PlusMap()->GetLargestPossibleRegion());
	    
	    /* t1 estimation and combination image correction using the B1+ map */
	    for(inv1It.GoToBegin(),inv2It.GoToBegin(),uniIt.GoToBegin(),inv1SignedIt.GoToBegin(),inv2SignedIt.GoToBegin(),flawsHCIt.GoToBegin(),flawsHCOIt.GoToBegin(),flawsHCDenIt.GoToBegin(),flawsHCODenIt.GoToBegin(),minIt.GoToBegin(),mipIt.GoToBegin(),mipDenIt.GoToBegin(),b1MapIt.GoToBegin(),t1MapIt.GoToBegin(),flawsHCB1PlusCorrIt.GoToBegin(),flawsHCOB1PlusCorrIt.GoToBegin() ;
		!flawsHCIt.IsAtEnd() ;
		++inv1It,++inv2It,++uniIt,++inv1SignedIt,++inv2SignedIt,++flawsHCIt,++flawsHCOIt,++flawsHCDenIt,++flawsHCODenIt,++minIt,++mipIt,++mipDenIt,++b1MapIt,++t1MapIt,++flawsHCB1PlusCorrIt,++flawsHCOB1PlusCorrIt)
	    {
		/* Get the image values and generate the combination images */	
		double inv1Value = inv1It.Get();
		double inv2Value = inv2It.Get();

		double sumValue = inv1Value+inv2Value;

		double inv1SignedValue = inv1Value;
		double inv2SignedValue = inv2Value;
		
		double uniValue = uniIt.Get();
		uniValue = uniValue/4095-0.5; // Get UNI back in its original value range

		unsigned long currentT1Value = 1500;
		unsigned long currentB1Value = 0;
		double sa2rageValue = b1T1LookupTable[(1500-1)*b1Range+((unsigned int) b1MapIt.Get())-1];

		double hcValue = (inv1Value-inv2Value)/sumValue;

		if(inv1Value==0 && inv2Value==0) // Avoid Nans
		{
		    hcValue=-1;
		}
		else if(hcValue > 1) //Constrains background noise values within the [-1;1] range
		{
		    hcValue=1;
		}
		else if(hcValue < -1)
		{
		    hcValue=-1;
		}

		double hcoValue = -1*hcValue;
		
		double minValue = inv1Value;

		if(inv2Value < minValue)
		{
		    minValue=inv2Value;
		}

		/* Compute the full range combination images according to the FLAWS1 and FLAWS2 signals signs */
		double mipValue,mipDenValue,hcFullRange,hcoFullRange,hcDenValue,hcoDenValue;
		
		if(uniValue<0)
		{
		    inv1SignedValue*=-1;

		    hcFullRange = hcValue;
		    hcoFullRange = hcoValue;

		    hcDenValue = (inv1Value-inv2Value-4*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);
		    hcoDenValue = (inv2Value-inv1Value-4*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);

		    mipValue = minValue/sumValue;
		    mipDenValue = (minValue-2*m_DenoisingCoeff)/(sumValue+2*m_DenoisingCoeff);
		}
		else if((uniValue >= 0) && (inv1Value > inv2Value))
		{
		    inv1SignedValue*=-1;
		    inv2SignedValue*=-1;

		    hcFullRange = 2 - hcValue;
		    hcoFullRange = -2 - hcoValue;
		    
		    hcDenValue = 2 - (inv1Value-inv2Value+4*m_DenoisingCoeff)/(sumValue+0.5*m_DenoisingCoeff);
		    hcoDenValue = -2 - (inv2Value-inv1Value+2*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);

		    mipValue = -minValue/sumValue;
		    mipDenValue = (-minValue-2*m_DenoisingCoeff)/(sumValue+2*m_DenoisingCoeff);
		}
		else // Case uniValue >=0 and inv1Value <= inv2Value
		{
		    hcFullRange = -2 - hcValue;
		    hcoFullRange = 2 - hcoValue;

		    hcDenValue = -2 - (inv1Value-inv2Value+2*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);		    
		    hcoDenValue = 2 - (inv2Value-inv1Value+4*m_DenoisingCoeff)/(sumValue+0.5*m_DenoisingCoeff);

		    mipValue = -minValue/sumValue;
		    mipDenValue = (-minValue-2*m_DenoisingCoeff)/(sumValue+2*m_DenoisingCoeff);
		}
		
		/* Perform the T1 mapping */
		if(inv1SignedValue>0 && inv2SignedValue >0) // T1 estimation before the -1 in FLAWS-HC
		{
		    for(size_t i=0 ; i<4 ; i++)
		    {
			currentB1Value = this->B1Estimation(b1Range,midB1Range,b1T1LookupTable,sa2rageValue,currentT1Value);			

			currentT1Value = this->T1Estimation(0,t1RangeBefore[currentB1Value-1],midT1RangeBefore[currentB1Value-1],b1Range,t1B1HCLookupTable,hcValue,currentB1Value);

			double uniError = fabs(uniValue - t1B1UNILookupTable[(currentT1Value-1)*b1Range+currentB1Value-1]);
			double hcError = fabs(hcValue - t1B1HCLookupTable[(currentT1Value-1)*b1Range+currentB1Value-1]);
		    }
		}
		else if(inv1SignedValue<=0 && inv2SignedValue >=0) // T1 estimation between the -1 and the 1 in FLAWS-HC
		{
		    for(size_t i=0 ; i<4 ; i++)
		    {
			currentB1Value = this->B1Estimation(b1Range,midB1Range,b1T1LookupTable,sa2rageValue,currentT1Value);			

			currentT1Value = this->T1Estimation(t1Min[currentB1Value-1],t1RangeBetween[currentB1Value-1],midT1RangeBetween[currentB1Value-1],b1Range,t1B1HCLookupTable,hcValue,currentB1Value);
		    }
		}
		else if(inv1SignedValue<0 && inv2SignedValue <0) // T1 estimation after the 1 in FLAWS-HC
		{
		    for(size_t i=0 ; i<4 ; i++)
		    {
			currentB1Value = this->B1Estimation(b1Range,midB1Range,b1T1LookupTable,sa2rageValue,currentT1Value);

			currentT1Value = this->T1Estimation(t1Max[currentB1Value-1],t1RangeAfter[currentB1Value-1],midT1RangeAfter[currentB1Value-1],b1Range,t1B1HCLookupTable,hcValue,currentB1Value);

			double uniError = fabs(uniValue - t1B1UNILookupTable[(currentT1Value-1)*b1Range+currentB1Value-1]);
			double hcError = fabs(hcValue - t1B1HCLookupTable[(currentT1Value-1)*b1Range+currentB1Value-1]);
		    }
		}
		else
		{
		    std::cerr << "Error: impossible signal signs combinations for FLAWS. Possible issues: the input images do not correspond to FLAWS images or the MP2RAGE UNI image was not computed online." << std::endl;
		    exit(EXIT_FAILURE);		    
		}

		/* Constrains background noise values within the intensity range of the combined images */
		if(hcDenValue > 1) //Constrains background noise values within the [-2;1] range
		{
		    hcDenValue=1;
		}
		else if(hcDenValue < -2)
		{
		    hcDenValue=-2;
		}

		if(hcoDenValue > 1) //Constrains background noise values within the [-2;1] range
		{
		    hcoDenValue=1;
		}
		else if(hcoDenValue < -2)
		{
		    hcoDenValue=-2;
		}
	    
		if(mipValue > 0.5) //Constrains background noise values within the [-0.5;0.5] range
		{
		    mipValue=0.5;
		}
		else if(mipValue < -0.5)
		{
		    mipValue=-0.5;
		}

		if(mipDenValue > 0.5) //Constrains background noise values within the [-0.5;0.5] range
		{
		    mipDenValue=0.5;
		}
		else if(mipDenValue < -0.5)
		{
		    mipDenValue=-0.5;
		}

		/* Set image values */
		inv1SignedIt.Set(inv1SignedValue);
		inv2SignedIt.Set(inv2SignedValue);	       

		flawsHCIt.Set(hcFullRange);
		flawsHCOIt.Set(hcoFullRange);

		flawsHCDenIt.Set(hcDenValue);
		flawsHCODenIt.Set(hcoDenValue);

		minIt.Set(minValue);
		mipIt.Set(mipValue);
		mipDenIt.Set(mipDenValue);
		
		t1MapIt.Set(currentT1Value);

		/* Generate the bias corrected hc and hco images */		
		double biasCorrectedHCValue = t1B1HCLookupTable[((currentT1Value-1)*b1Range)+999]; // Plus 999 to get a signal with a B1+ of 1.
		
		if(currentT1Value < t1Min[currentB1Value-1]) // Reconstruct the image in its full intensity range
		{
		    biasCorrectedHCValue=-2-biasCorrectedHCValue;
		}
		else if(currentT1Value > t1Max[currentB1Value-1])
		{
		    biasCorrectedHCValue=2-biasCorrectedHCValue;
		}
		
		flawsHCB1PlusCorrIt.Set(biasCorrectedHCValue);
		flawsHCOB1PlusCorrIt.Set(-1*biasCorrectedHCValue);
	    }
	}
	else
	{	   
	    /* lookupTable generation */
	    for(unsigned long t1=1 ; t1<=t1Range ; t1++)
	    {
		this->T1Simulation(((double) t1)/1000, 1, t1LookupTable+(t1-1)*2);
	    }
	    
	    double min=-1,max=1;
	    double distMin=100,distMax=100;
	    unsigned long t1Min,t1Max;

	    for(unsigned long t1=1 ; t1<=t1Range ; t1++)
	    {
		double flaws1 = t1LookupTable[(t1-1)*2];
		double flaws2 = t1LookupTable[(t1-1)*2+1];
		    
		double hcValue = (fabs(flaws1)-fabs(flaws2))/(fabs(flaws1)+fabs(flaws2));

		t1HCLookupTable[t1-1] = hcValue;

		if(fabs(hcValue-min)<distMin)
		{
		    distMin = fabs(hcValue-min);
		    t1Min = t1;
		}
		else if(fabs(hcValue-max)<distMax)
		{
		    distMax = fabs(hcValue-max);
		    t1Max = t1;
		}

		double uniValue = (flaws1*flaws2)/(flaws1*flaws1+flaws2*flaws2);
		t1UNILookupTable[t1-1] = uniValue;
	    }		

	    unsigned long t1RangeBefore = t1Min;
	    unsigned long t1RangeBetween = t1Max-t1Min;
	    unsigned long t1RangeAfter = t1Range-t1Max;
		
	    unsigned long midT1RangeBefore = t1RangeBefore/2;
	    unsigned long midT1RangeBetween = t1RangeBetween/2+t1Min;
	    unsigned long midT1RangeAfter = t1RangeAfter/2+t1Max;
		
	    /* t1 estimation */
	    for(inv1It.GoToBegin(),inv2It.GoToBegin(),flawsHCIt.GoToBegin(),flawsHCOIt.GoToBegin(),flawsHCDenIt.GoToBegin(),flawsHCODenIt.GoToBegin(),minIt.GoToBegin(),mipIt.GoToBegin(),mipDenIt.GoToBegin(),t1MapIt.GoToBegin(),inv1SignedIt.GoToBegin(),inv2SignedIt.GoToBegin(),uniIt.GoToBegin() ;
		!flawsHCIt.IsAtEnd() ;
		++inv1It,++inv2It,++flawsHCIt,++flawsHCOIt,++flawsHCDenIt,++flawsHCODenIt,++minIt,++mipIt,++mipDenIt,++t1MapIt,++inv1SignedIt,++inv2SignedIt,++uniIt)
	    {
		/* Get the image values and generate the combination images */	
		double inv1Value = inv1It.Get();
		double inv2Value = inv2It.Get();

		double sumValue = inv1Value+inv2Value;

		double inv1SignedValue = inv1Value;
		double inv2SignedValue = inv2Value;
		
		double uniValue = uniIt.Get();
		uniValue = uniValue/4095-0.5; // Get UNI back in its original value range

		double hcValue = (inv1Value-inv2Value)/sumValue;

		if(inv1Value==0 && inv2Value==0) // Avoid Nans
		{
		    hcValue=-1;
		}
		else if(hcValue > 1) //Constrains background noise values within the [-1;1] range
		{
		    hcValue=1;
		}
		else if(hcValue < -1)
		{
		    hcValue=-1;
		}

		double hcoValue = -1*hcValue;

		
		double minValue = inv1Value;

		if(inv2Value < minValue)
		{
		    minValue=inv2Value;
		}

		/* Compute the full range combination images according to the FLAWS1 and FLAWS2 signals signs */
		double mipValue,mipDenValue,hcFullRange,hcoFullRange,hcDenValue,hcoDenValue;
		
		if(uniValue<0)
		{
		    inv1SignedValue*=-1;

		    hcFullRange = hcValue;
		    hcoFullRange = hcoValue;

		    hcDenValue = (inv1Value-inv2Value-4*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);
		    hcoDenValue = (inv2Value-inv1Value-4*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);

		    mipValue = minValue/sumValue;
		    mipDenValue = (minValue-2*m_DenoisingCoeff)/(sumValue+2*m_DenoisingCoeff);	       
		}
		else if((uniValue >= 0) && (inv1Value > inv2Value))
		{
		    inv1SignedValue*=-1;
		    inv2SignedValue*=-1;

		    hcFullRange = 2 - hcValue;
		    hcoFullRange = -2 - hcoValue;

		    hcDenValue = 2 - (inv1Value-inv2Value+4*m_DenoisingCoeff)/(sumValue+0.5*m_DenoisingCoeff);
		    hcoDenValue = -2 - (inv2Value-inv1Value+2*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);		   

		    mipValue = -minValue/sumValue;
		    mipDenValue = (-minValue-2*m_DenoisingCoeff)/(sumValue+2*m_DenoisingCoeff);
		}
		else // Case uniValue >=0 and inv1Value <= inv2Value
		{
		    hcFullRange = -2 - hcValue;
		    hcoFullRange = 2 - hcoValue;

		    hcDenValue = -2 - (inv1Value-inv2Value+2*m_DenoisingCoeff)/(sumValue+m_DenoisingCoeff);		    
		    hcoDenValue = 2 - (inv2Value-inv1Value+4*m_DenoisingCoeff)/(sumValue+0.5*m_DenoisingCoeff);

		    mipValue = -minValue/sumValue;
		    mipDenValue = (-minValue-2*m_DenoisingCoeff)/(sumValue+2*m_DenoisingCoeff);
		}

		/* Perform T1 mapping */
		unsigned long currentT1Value=0;

		if(inv1SignedValue>=0 && inv2SignedValue >=0) // T1 estimation before the -1 in the combination image
		{

		    currentT1Value = this->T1Estimation(0,t1RangeBefore,midT1RangeBefore,t1HCLookupTable,hcValue);

		    double uniError = fabs(uniValue - t1UNILookupTable[currentT1Value-1]);
		    double hcError = fabs(hcValue - t1HCLookupTable[currentT1Value-1]);
		}
		else if(inv1SignedValue<0 && inv2SignedValue >=0)
		{
		    currentT1Value = this->T1Estimation(t1Min,t1RangeBetween,midT1RangeBetween,t1HCLookupTable,hcValue);
		}
		else if(inv1SignedValue<0 && inv2SignedValue <0)
		{
		    currentT1Value = this->T1Estimation(t1Max,t1RangeAfter,midT1RangeAfter,t1HCLookupTable,hcValue);

		    double uniError = fabs(uniValue - t1UNILookupTable[currentT1Value-1]);
		    double hcError = fabs(hcValue - t1HCLookupTable[currentT1Value-1]);
		}
		else
		{
		    std::cerr << "Error: impossible signal signs combinations for FLAWS. Possible issues: the input images do not correspond to FLAWS images or the phase unwrapping failed." << std::endl;
		    exit(EXIT_FAILURE);
		}
		
		/* Constrains background noise values within the intensity range of the combined images */
		if(hcDenValue > 1) //Constrains background noise values within the [-2;1] range
		{
		    hcDenValue=1;
		}
		else if(hcDenValue < -2)
		{
		    hcDenValue=-2;
		}

		if(hcoDenValue > 1) //Constrains background noise values within the [-2;1] range
		{
		    hcoDenValue=1;
		}
		else if(hcoDenValue < -2)
		{
		    hcoDenValue=-2;
		}
	    
		if(mipValue > 0.5) //Constrains background noise values within the [-0.5;0.5] range
		{
		    mipValue=0.5;
		}
		else if(mipValue < -0.5)
		{
		    mipValue=-0.5;
		}

		if(mipDenValue > 0.5) //Constrains background noise values within the [-0.5;0.5] range
		{
		    mipDenValue=0.5;
		}
		else if(mipDenValue < -0.5)
		{
		    mipDenValue=-0.5;
		}

		/* Set image values */
		inv1SignedIt.Set(inv1SignedValue);
		inv2SignedIt.Set(inv2SignedValue);	       

		flawsHCIt.Set(hcFullRange);
		flawsHCOIt.Set(hcoFullRange);

		flawsHCDenIt.Set(hcDenValue);
		flawsHCODenIt.Set(hcoDenValue);

		minIt.Set(minValue);
		mipIt.Set(mipValue);
		mipDenIt.Set(mipDenValue);
		
		t1MapIt.Set(currentT1Value);
	    }
	}
    }    
    
    /* Run the filter */
    template <class TInputImage, class TOutputImage>
    void
    FLAWSProcessingImageFilter<TInputImage, TOutputImage>
    ::GenerateData()
    {
	/* Check the coherence in sequence parameters */
	this->CheckSequenceParameters();

	/* Memory allocation of the filter outputs */
	this->AllocateOutputs();

	typename TOutputImage::Pointer t1Map = this->GetT1Map();
	t1Map->SetOrigin(this->GetInput()->GetOrigin());
	t1Map->SetDirection(this->GetInput()->GetDirection());
	t1Map->SetSpacing(this->GetInput()->GetSpacing());
	t1Map->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	t1Map->Allocate();
	
	typename TOutputImage::Pointer flawsHC = this->GetFLAWSHC();
	flawsHC->SetOrigin(this->GetInput()->GetOrigin());
	flawsHC->SetDirection(this->GetInput()->GetDirection());
	flawsHC->SetSpacing(this->GetInput()->GetSpacing());
	flawsHC->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	flawsHC->Allocate();

	typename TOutputImage::Pointer flawsHCO = this->GetFLAWSHCO();
	flawsHCO->SetOrigin(this->GetInput()->GetOrigin());
	flawsHCO->SetDirection(this->GetInput()->GetDirection());
	flawsHCO->SetSpacing(this->GetInput()->GetSpacing());
	flawsHCO->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	flawsHCO->Allocate();

	typename TOutputImage::Pointer flawsHCB1PlusCorr = this->GetFLAWSHCB1PlusCorr();
	flawsHCB1PlusCorr->SetOrigin(this->GetInput()->GetOrigin());
	flawsHCB1PlusCorr->SetDirection(this->GetInput()->GetDirection());
	flawsHCB1PlusCorr->SetSpacing(this->GetInput()->GetSpacing());
	flawsHCB1PlusCorr->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	flawsHCB1PlusCorr->Allocate();

	typename TOutputImage::Pointer flawsHCOB1PlusCorr = this->GetFLAWSHCOB1PlusCorr();
	flawsHCOB1PlusCorr->SetOrigin(this->GetInput()->GetOrigin());
	flawsHCOB1PlusCorr->SetDirection(this->GetInput()->GetDirection());
	flawsHCOB1PlusCorr->SetSpacing(this->GetInput()->GetSpacing());
	flawsHCOB1PlusCorr->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	flawsHCOB1PlusCorr->Allocate();

	typename TOutputImage::Pointer flawsHCDen = this->GetFLAWSHCDen();
	flawsHCDen->SetOrigin(this->GetInput()->GetOrigin());
	flawsHCDen->SetDirection(this->GetInput()->GetDirection());
	flawsHCDen->SetSpacing(this->GetInput()->GetSpacing());
	flawsHCDen->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	flawsHCDen->Allocate();

	typename TOutputImage::Pointer flawsHCODen = this->GetFLAWSHCODen();
	flawsHCODen->SetOrigin(this->GetInput()->GetOrigin());
	flawsHCODen->SetDirection(this->GetInput()->GetDirection());
	flawsHCODen->SetSpacing(this->GetInput()->GetSpacing());
	flawsHCODen->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	flawsHCODen->Allocate();

	typename TOutputImage::Pointer inv1Signed = this->GetInv1Signed();
	inv1Signed->SetOrigin(this->GetInput()->GetOrigin());
	inv1Signed->SetDirection(this->GetInput()->GetDirection());
	inv1Signed->SetSpacing(this->GetInput()->GetSpacing());
	inv1Signed->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	inv1Signed->Allocate();

	typename TOutputImage::Pointer inv2Signed = this->GetInv2Signed();
	inv2Signed->SetOrigin(this->GetInput()->GetOrigin());
	inv2Signed->SetDirection(this->GetInput()->GetDirection());
	inv2Signed->SetSpacing(this->GetInput()->GetSpacing());
	inv2Signed->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	inv2Signed->Allocate();

	typename TOutputImage::Pointer min = this->GetMin();
	min->SetOrigin(this->GetInput()->GetOrigin());
	min->SetDirection(this->GetInput()->GetDirection());
	min->SetSpacing(this->GetInput()->GetSpacing());
	min->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	min ->Allocate();

	typename TOutputImage::Pointer mip = this->GetMip();
	mip->SetOrigin(this->GetInput()->GetOrigin());
	mip->SetDirection(this->GetInput()->GetDirection());
	mip->SetSpacing(this->GetInput()->GetSpacing());
	mip->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	mip ->Allocate();

	typename TOutputImage::Pointer mipDen = this->GetMipDen();
	mipDen->SetOrigin(this->GetInput()->GetOrigin());
	mipDen->SetDirection(this->GetInput()->GetDirection());
	mipDen->SetSpacing(this->GetInput()->GetSpacing());
	mipDen->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	mipDen ->Allocate();

	/* Perform T1 mapping and image combination */
	this->T1Mapping();
    }
}  // end namespace

#endif
