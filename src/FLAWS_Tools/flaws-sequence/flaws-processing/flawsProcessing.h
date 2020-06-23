#ifndef __FLAWSProcessingImageFilter_h
#define __FLAWSProcessingImageFilter_h

#define G_PI 3.14159265358979323846 /* value of pi */

#include "itkImageToImageFilter.h"

#include "itkImage.h"

namespace flaws
{

/** \class FLAWSProcessingImageFilter
 *  \brief Generate the FLAWS combination images and T1 maps from FLAWS1, FLAWS2 and UNI
 *
 *  \author Jeremy Beaumont
 *  \version 0.2
 *  \date Wed 16 January 2019 12:32:00
 *  \ingroup flaws
 */
    template <class TInputImage, class TOutputImage>
	class FLAWSProcessingImageFilter:
	public itk::ImageToImageFilter<TInputImage, TOutputImage>
    {
    public:
	/** Standard class typedefs */
	typedef FLAWSProcessingImageFilter                     Self;
	typedef itk::ImageToImageFilter<TInputImage,TOutputImage> Superclass;
	typedef itk::SmartPointer<Self>                           Pointer;
	typedef itk::SmartPointer<const Self>                     ConstPointer;

	/** Method for creation through the object factory */
	itkNewMacro(Self);

	/** Run-time type information */
	itkTypeMacro(FLAWSProcessingImageFilter, itk::ImageToImageFilter);

	/** Inv1 */
	void SetInv1(const TInputImage *Inv1)
	{
	    this->SetInput(Inv1);
	}

	const TInputImage* GetInv1() const
	{
	    return this->GetInput();
	}
	
	/** Inv2 */
	void SetInv2(const TInputImage *Inv2)
	{
	    this->SetNthInput(1, const_cast<TInputImage *>(Inv2) );
	}

	const TInputImage* GetInv2() const
	{
	    return static_cast<const TInputImage*>(this->GetInput(1));
	}

	/** MP2RAGE UNI */
	void SetUNI(const TInputImage *UNI)
	{
	    this->SetNthInput(2, const_cast<TInputImage *>(UNI) );
	}

	const TInputImage* GetUNI() const
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
	
	/* T1Map */
	TOutputImage* GetT1Map()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(1) );
	}

	/* FLAWS HC */
	TOutputImage* GetFLAWSHC()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(2) );
	}

	/* FLAWS HCO */
	TOutputImage* GetFLAWSHCO()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(3) );
	}
	
	/* FLAWS HC B1 Plus corrected */
	TOutputImage* GetFLAWSHCB1PlusCorr()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(4) );
	}

	/* FLAWS HCO B1 Plus corrected */
	TOutputImage* GetFLAWSHCOB1PlusCorr()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(5) );
	}

	/* FLAWS HC Den*/
	TOutputImage* GetFLAWSHCDen()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(6) );
	}

	/* FLAWS HCO Den*/
	TOutputImage* GetFLAWSHCODen()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(7) );
	}

	/* Inversion 1 unwrapped phase image */
	TOutputImage* GetInv1Signed()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(8) );
	}

	/* Inversion 2 unwrapped phase image */
	TOutputImage* GetInv2Signed()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(9) );
	}

	/* FLAWS min image */
	TOutputImage* GetMin()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(10) );
	}

	/* FLAWS mip image */
	TOutputImage* GetMip()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(11) );
	}

	/* FLAWS mip-den image */
	TOutputImage* GetMipDen()
	{
	    return dynamic_cast<TOutputImage*>( this->itk::ProcessObject::GetOutput(12) );
	}

	
	/** Denoising coefficient */
	itkSetMacro(DenoisingCoeff, unsigned short);
	itkGetConstMacro(DenoisingCoeff, const unsigned short);

	/** FLAWS sequence information: Alpha1 (Degrees) */
	itkSetMacro(Alpha1, unsigned char);
	itkGetConstMacro(Alpha1, const unsigned char);

	/** FLAWS sequence information: Alpha2 (Degrees) */
	itkSetMacro(Alpha2, unsigned char);
	itkGetConstMacro(Alpha2, const unsigned char);

	/** FLAWS sequence information: Inversion time 1 (s) */
	itkSetMacro(InversionTime1, double);
	itkGetConstMacro(InversionTime1, const double);

	/** FLAWS sequence information: Inversion time 2 (s) */
	itkSetMacro(InversionTime2, double);
	itkGetConstMacro(InversionTime2, const double);
	
	/** FLAWS sequence information: Echo time (s) */
	itkSetMacro(EchoTime, double);
	itkGetConstMacro(EchoTime, const double);

	/** FLAWS sequence information: Echo spacing (s) */
	itkSetMacro(EchoSpacing, double);
	itkGetConstMacro(EchoSpacing, const double);
	
	/** FLAWS sequence information: Sequence repetition time (s) */
	itkSetMacro(RepetitionTime, double);
	itkGetConstMacro(RepetitionTime, const double);

	/** FLAWS sequence information: Number of slices */
	itkSetMacro(NumberOfSlices, unsigned short);
	itkGetConstMacro(NumberOfSlices, const unsigned short);
      
	/** FLAWS sequence information: Partial Fourier (6,7,8) */
	itkSetMacro(PartialFourier, unsigned char);
	itkGetConstMacro(PartialFourier, const unsigned char);

	/** SA2RAGE sequence information: Alpha1 (Degrees) */
	itkSetMacro(Alpha1b, unsigned char);
	itkGetConstMacro(Alpha1b, const unsigned char);

	/** SA2RAGE sequence information: Alpha2 (Degrees) */
	itkSetMacro(Alpha2b, unsigned char);
	itkGetConstMacro(Alpha2b, const unsigned char);

	/** SA2RAGE sequence information: Inversion time 1 (s) */
	itkSetMacro(InversionTime1b, double);
	itkGetConstMacro(InversionTime1b, const double);

	/** SA2RAGE sequence information: Inversion time 2 (s) */
	itkSetMacro(InversionTime2b, double);
	itkGetConstMacro(InversionTime2b, const double);
	
	/** SA2RAGE sequence information: Echo spacing (s) */
	itkSetMacro(EchoSpacingb, double);
	itkGetConstMacro(EchoSpacingb, const double);
	
	/** SA2RAGE sequence information: Sequence repetition time (s) */
	itkSetMacro(RepetitionTimeb, double);
	itkGetConstMacro(RepetitionTimeb, const double);

	/** SA2RAGE sequence information: Number of slices */
	itkSetMacro(NumberOfSlicesb, unsigned short);
	itkGetConstMacro(NumberOfSlicesb, const unsigned short);
      
	/** SA2RAGE sequence information: Partial Fourier (6,7,8) */
	itkSetMacro(PartialFourierb, unsigned char);
	itkGetConstMacro(PartialFourierb, const unsigned char);
	
    protected:
	/** Constructor */
	FLAWSProcessingImageFilter();
  
	/** Destructor */
	virtual ~FLAWSProcessingImageFilter() {};

	/** Check the coherence in sequence parameters */
	void CheckSequenceParameters();

	/** Simulation of the FLAWS sequence to generate the T1 lookup table */
	void T1Simulation(double,double,double*);

	/** Simulation of the SA2RAGE sequence to generate the B1 lookup table */
	double B1Simulation(double,double);

	/** Estimation of the T1 from FLAWS-HC and FLAWS lookup table */
	unsigned long T1Estimation(unsigned long,unsigned long,unsigned long,double*,double);
	    
	/** Estimation of the T1 from FLAWS-HC and FLAWS lookup table */
	unsigned long T1Estimation(unsigned long,unsigned long,unsigned long,unsigned long,double*,double,unsigned long);

	/** Estimation of the B1+ from the SA2RAGE division image and lookup tables */
	unsigned long B1Estimation(unsigned long,unsigned long,double*,double,unsigned long);

	/* T1 Estimation */
	void T1Mapping();
		  
	/** Run the filter */
	virtual void GenerateData();

    private:
	FLAWSProcessingImageFilter(const Self&);  // purposely not implemented
	void operator=(const Self&);  // purposely not implemented

	/* Sequence parameters  */
	unsigned char m_PartialFourier,m_Alpha1,m_Alpha2,m_PartialFourierb,m_Alpha1b,m_Alpha2b;
	unsigned short m_NumberOfSlices,m_nbEx,m_nbExBefore,m_nbExAfter,m_NumberOfSlicesb,m_nbExb,m_nbExBeforeb,m_nbExAfterb,m_DenoisingCoeff;
	double m_InversionTime1,m_InversionTime2,m_RepetitionTime,m_EchoTime,m_EchoSpacing,m_Alpha1Rad,m_Alpha2Rad,m_timeGre,m_timeGreBefore,m_timeGreAfter,m_invEff,m_InversionTime1b,m_InversionTime2b,m_RepetitionTimeb,m_EchoSpacingb,m_Alpha1Radb,m_Alpha2Radb,m_timeGreb,m_timeGreBeforeb,m_timeGreAfterb;
    };  // end class

}  // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "flawsProcessing.hxx"
#endif

#endif
