#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <flawsProcessing.h>

#include <tclap/CmdLine.h>

int main(int argc, const char** argv)
{
    // Parsing arguments
    TCLAP::CmdLine cmd("FLAWS Tools - FLAWS Processing", ' ',FLAWS_TOOLS_VERSION);

    TCLAP::ValueArg<std::string> inv1Arg("","inv1","Input INV1 image",true,"","INV1",cmd);
    TCLAP::ValueArg<std::string> inv2Arg("","inv2","Input INV2 image",true,"","INV2",cmd);

    TCLAP::ValueArg<std::string> uniArg("","uni","Input UNI image",true,"","UNI",cmd);

    TCLAP::ValueArg<std::string> b1PlusArg("","b1Plus","B1 Plus Map",false,"","B1PLUS",cmd);

    TCLAP::ValueArg<unsigned short> denoisingCoeffArg("","denoising","Denoising coefficient (Default: 20)",false,20,"DENOISINGCOEFF",cmd);
    
    TCLAP::ValueArg<std::string> outputArg("o","output","Output folder.",true,"","OUTPUT",cmd);

    // FLAWS parameters
    TCLAP::ValueArg<unsigned short> alpha1Arg("","alpha1","First flip angle (degrees).",true,0,"ALPHA1",cmd);
    TCLAP::ValueArg<unsigned short> alpha2Arg("","alpha2","Second flip angle (degrees).",true,0,"ALPHA2",cmd);

    TCLAP::ValueArg<double> ti1Arg("","ti1","First inversion time (sec).",true,0,"TI1",cmd);
    TCLAP::ValueArg<double> ti2Arg("","ti2","Second inversion time (sec).",true,0,"TI2",cmd);

    TCLAP::ValueArg<double> teArg("","te","Echo time (sec).",true,0,"TE",cmd);
    TCLAP::ValueArg<double> esArg("","es","Echo spacing (sec).",true,0,"ES",cmd);
    TCLAP::ValueArg<double> trArg("","tr","Repetition time (sec).",true,0,"TR",cmd);

    TCLAP::ValueArg<unsigned short> nbSlicesArg("","nbSlices","Number of slices.",true,0,"NBSLICES",cmd);

    std::vector<unsigned int> allowedPF;
    allowedPF.push_back(6);
    allowedPF.push_back(7);
    allowedPF.push_back(8);
    
    TCLAP::ValuesConstraint<unsigned int> allowedPartialFourier(allowedPF);

    TCLAP::ValueArg<unsigned int> partialFourierArg("","partialFourier","Set to 6 if 6/8 partial fourier is used (same idea for 7/8, Default: 8/8).",true,8,&allowedPartialFourier);
    cmd.add(partialFourierArg);

    // SA2RAGE Parameters
    TCLAP::ValueArg<unsigned short> alpha1bArg("","alpha1b","SA2RAGE: First flip angle (degrees).",false,0,"ALPHA1B",cmd);
    TCLAP::ValueArg<unsigned short> alpha2bArg("","alpha2b","SA2RAGE: Second flip angle (degrees).",false,0,"ALPHA2B",cmd);

    TCLAP::ValueArg<double> ti1bArg("","ti1b","SA2RAGE: First inversion time (sec).",false,0,"TI1B",cmd);
    TCLAP::ValueArg<double> ti2bArg("","ti2b","SA2RAGE: Second inversion time (sec).",false,0,"TI2B",cmd);

    TCLAP::ValueArg<double> esbArg("","esb","SA2RAGE: Echo spacing (sec).",false,0,"ESB",cmd);
    TCLAP::ValueArg<double> trbArg("","trb","SA2RAGE: Repetition time (sec).",false,0,"TRB",cmd);
    
    TCLAP::ValueArg<unsigned short> nbSlicesbArg("","nbSlicesb","SA2RAGE: Number of slices.",false,0,"NBSLICESB",cmd);

    std::vector<unsigned int> allowedPFb;
    allowedPFb.push_back(6);
    allowedPFb.push_back(7);
    allowedPFb.push_back(8);
    
    TCLAP::ValuesConstraint<unsigned int> allowedPartialFourierb(allowedPFb);

    TCLAP::ValueArg<unsigned int> partialFourierbArg("","partialFourierb","SA2RAGE: Set to 6 if 6/8 partial fourier is used (same idea for 7/8, Default: 8/8).",false,8,&allowedPartialFourierb);
    cmd.add(partialFourierbArg);

    try
    {
	cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
	std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
	return EXIT_FAILURE;
    }

    if(b1PlusArg.isSet() && !(alpha1bArg.isSet() && alpha2bArg.isSet() && ti1bArg.isSet() && ti2bArg.isSet() && esbArg.isSet() && trbArg.isSet() && nbSlicesbArg.isSet() && partialFourierbArg.isSet()))
    {
	std::cerr << "IOError: The SA2RAGE arguments must be provided when a B1 Plus map is provided." << std::endl;
	return EXIT_FAILURE;
    }
    
    // Constants
    const unsigned int Dimension = 3;

    // Typedefs
    typedef itk::Image<double,Dimension> DoubleImageType;
    typedef itk::ImageFileReader <DoubleImageType> DoubleReaderType;
    typedef itk::ImageFileWriter <DoubleImageType> DoubleWriterType;

    typedef flaws::FLAWSProcessingImageFilter <DoubleImageType,DoubleImageType> FLAWSProcessingFilterType;
    
    //Read images
    DoubleReaderType::Pointer inv1Reader = DoubleReaderType::New();
    inv1Reader->SetFileName(inv1Arg.getValue());

    DoubleReaderType::Pointer inv2Reader = DoubleReaderType::New();
    inv2Reader->SetFileName(inv2Arg.getValue());

    DoubleReaderType::Pointer uniReader = DoubleReaderType::New();
    uniReader->SetFileName(uniArg.getValue());
    
    DoubleReaderType::Pointer b1PlusReader = DoubleReaderType::New();
    if(b1PlusArg.isSet())
    {
	b1PlusReader->SetFileName(b1PlusArg.getValue());
    }
    
    try
    {
	inv1Reader->Update();
	inv2Reader->Update();
	uniReader->Update();

	if(b1PlusArg.isSet())
	{
	    b1PlusReader->Update();
	}	
    }
    catch(itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return EXIT_FAILURE;
    }

    FLAWSProcessingFilterType::Pointer flawsProcessor = FLAWSProcessingFilterType::New();

    flawsProcessor->SetInv1(inv1Reader->GetOutput());
    flawsProcessor->SetInv2(inv2Reader->GetOutput());
    flawsProcessor->SetUNI(uniReader->GetOutput());

    flawsProcessor->SetDenoisingCoeff(denoisingCoeffArg.getValue());
    
    flawsProcessor->SetAlpha1(alpha1Arg.getValue());
    flawsProcessor->SetAlpha2(alpha2Arg.getValue());

    flawsProcessor->SetInversionTime1(ti1Arg.getValue());
    flawsProcessor->SetInversionTime2(ti2Arg.getValue());

    flawsProcessor->SetEchoTime(teArg.getValue());
    flawsProcessor->SetEchoSpacing(esArg.getValue());
    flawsProcessor->SetRepetitionTime(trArg.getValue());

    flawsProcessor->SetNumberOfSlices(nbSlicesArg.getValue());
    flawsProcessor->SetPartialFourier(partialFourierArg.getValue());

    if(b1PlusArg.isSet())
    {
	flawsProcessor->SetB1PlusMap(b1PlusReader->GetOutput());

	flawsProcessor->SetAlpha1b(alpha1bArg.getValue());
	flawsProcessor->SetAlpha2b(alpha2bArg.getValue());

	flawsProcessor->SetInversionTime1b(ti1bArg.getValue());
	flawsProcessor->SetInversionTime2b(ti2bArg.getValue());

	flawsProcessor->SetEchoSpacingb(esbArg.getValue());
	flawsProcessor->SetRepetitionTimeb(trbArg.getValue());

	flawsProcessor->SetNumberOfSlicesb(nbSlicesbArg.getValue());
	flawsProcessor->SetPartialFourierb(partialFourierbArg.getValue());
    }
	
    try
    {
    	flawsProcessor->Update();
    }
    catch (itk::ExceptionObject &e)
    {
    	std::cerr << e << std::endl;
    	return EXIT_FAILURE;
    }

    std::string t1MapPath=outputArg.getValue()+"/FLAWS_T1MAP.nii.gz";
    std::string flawsHCPath=outputArg.getValue()+"/FLAWS_HC.nii.gz";
    std::string flawsHCOPath=outputArg.getValue()+"/FLAWS_HCO.nii.gz";
    std::string flawsHCB1PlusCorrPath=outputArg.getValue()+"/FLAWS_HC_B1PlusCorr.nii.gz";
    std::string flawsHCOB1PlusCorrPath=outputArg.getValue()+"/FLAWS_HCO_B1PlusCorr.nii.gz";
    std::string flawsHCDenPath=outputArg.getValue()+"/FLAWS_HC-DEN.nii.gz";
    std::string flawsHCODenPath=outputArg.getValue()+"/FLAWS_HCO-DEN.nii.gz";
    
    std::string inv1SignedPath=outputArg.getValue()+"/FLAWS_INV1_Signed.nii.gz";
    std::string inv2SignedPath=outputArg.getValue()+"/FLAWS_INV2_Signed.nii.gz";

    std::string minPath=outputArg.getValue()+"/FLAWS_MIN.nii.gz";
    std::string mipPath=outputArg.getValue()+"/FLAWS_MIP.nii.gz";
    std::string mipDenPath=outputArg.getValue()+"/FLAWS_MIP-DEN.nii.gz";
    
    DoubleWriterType::Pointer t1MapWriter = DoubleWriterType::New();	
    t1MapWriter->SetFileName(t1MapPath);		
    t1MapWriter->SetUseCompression(true);
    t1MapWriter->SetInput(flawsProcessor->GetT1Map());

    DoubleWriterType::Pointer flawsHCWriter = DoubleWriterType::New();	
    flawsHCWriter->SetFileName(flawsHCPath);		
    flawsHCWriter->SetUseCompression(true);
    flawsHCWriter->SetInput(flawsProcessor->GetFLAWSHC());

    DoubleWriterType::Pointer flawsHCOWriter = DoubleWriterType::New();	
    flawsHCOWriter->SetFileName(flawsHCOPath);		
    flawsHCOWriter->SetUseCompression(true);
    flawsHCOWriter->SetInput(flawsProcessor->GetFLAWSHCO());

    DoubleWriterType::Pointer flawsHCB1PlusCorrWriter = DoubleWriterType::New();	
    flawsHCB1PlusCorrWriter->SetFileName(flawsHCB1PlusCorrPath);		
    flawsHCB1PlusCorrWriter->SetUseCompression(true);
    flawsHCB1PlusCorrWriter->SetInput(flawsProcessor->GetFLAWSHCB1PlusCorr());

    DoubleWriterType::Pointer flawsHCOB1PlusCorrWriter = DoubleWriterType::New();	
    flawsHCOB1PlusCorrWriter->SetFileName(flawsHCOB1PlusCorrPath);		
    flawsHCOB1PlusCorrWriter->SetUseCompression(true);
    flawsHCOB1PlusCorrWriter->SetInput(flawsProcessor->GetFLAWSHCOB1PlusCorr());

    DoubleWriterType::Pointer flawsHCDenWriter = DoubleWriterType::New();	
    flawsHCDenWriter->SetFileName(flawsHCDenPath);		
    flawsHCDenWriter->SetUseCompression(true);
    flawsHCDenWriter->SetInput(flawsProcessor->GetFLAWSHCDen());

    DoubleWriterType::Pointer flawsHCODenWriter = DoubleWriterType::New();	
    flawsHCODenWriter->SetFileName(flawsHCODenPath);		
    flawsHCODenWriter->SetUseCompression(true);
    flawsHCODenWriter->SetInput(flawsProcessor->GetFLAWSHCODen());
    
    DoubleWriterType::Pointer inv1SignedWriter = DoubleWriterType::New();	
    inv1SignedWriter->SetFileName(inv1SignedPath);		
    inv1SignedWriter->SetUseCompression(true);
    inv1SignedWriter->SetInput(flawsProcessor->GetInv1Signed());

    DoubleWriterType::Pointer inv2SignedWriter = DoubleWriterType::New();	
    inv2SignedWriter->SetFileName(inv2SignedPath);		
    inv2SignedWriter->SetUseCompression(true);
    inv2SignedWriter->SetInput(flawsProcessor->GetInv2Signed());

    DoubleWriterType::Pointer minWriter = DoubleWriterType::New();	
    minWriter->SetFileName(minPath);		
    minWriter->SetUseCompression(true);
    minWriter->SetInput(flawsProcessor->GetMin());

    DoubleWriterType::Pointer mipWriter = DoubleWriterType::New();	
    mipWriter->SetFileName(mipPath);		
    mipWriter->SetUseCompression(true);
    mipWriter->SetInput(flawsProcessor->GetMip());

    DoubleWriterType::Pointer mipDenWriter = DoubleWriterType::New();	
    mipDenWriter->SetFileName(mipDenPath);		
    mipDenWriter->SetUseCompression(true);
    mipDenWriter->SetInput(flawsProcessor->GetMipDen());

    try
    {
	t1MapWriter->Update();
	flawsHCWriter->Update();
	flawsHCOWriter->Update();

	if(b1PlusArg.isSet())
	{
	    flawsHCB1PlusCorrWriter->Update();
	    flawsHCOB1PlusCorrWriter->Update();
	}
	
	flawsHCDenWriter->Update();
	flawsHCODenWriter->Update();

	inv1SignedWriter->Update();
	inv2SignedWriter->Update();

	minWriter->Update();
	mipWriter->Update();
	mipDenWriter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
    	std::cerr << e << std::endl;
    	return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
