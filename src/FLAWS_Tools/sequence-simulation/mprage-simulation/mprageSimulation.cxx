#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <mprageSimulation.h>

#include <tclap/CmdLine.h>

int main(int argc, const char** argv)
{
    // Parsing arguments
    TCLAP::CmdLine cmd("FLAWS Tools - MPRAGE Signal simulation", ' ',FLAWS_TOOLS_VERSION);

    TCLAP::ValueArg<std::string> wmArg("","wmMap","Input WM map",false,"","WMMAP",cmd);
    TCLAP::ValueArg<std::string> gmArg("","gmMap","Input GM map",false,"","GMMAP",cmd);
    TCLAP::ValueArg<std::string> csfArg("","csfMap","Input CSF map",false,"","CSFMAP",cmd);

    TCLAP::ValueArg<std::string> t1Arg("","t1Map","Input T1 map",false,"","T1MAP",cmd);

    TCLAP::ValueArg<std::string> b1PlusArg("","b1PlusMap","B1+ Map",false,"","B1PLUSMAP",cmd);
    TCLAP::ValueArg<std::string> b1MinusArg("","b1MinusMap","B1- Map",false,"","B1MINUSMAP",cmd);
    
    TCLAP::ValueArg<std::string> outputArg("o","output","Output image.",true,"","OUTPUT",cmd);

    std::vector<float> allowedField;
    allowedField.push_back(1.5);
    allowedField.push_back(3);
    allowedField.push_back(7);

    TCLAP::ValuesConstraint<float> allowedFieldValues(allowedField);

    TCLAP::ValueArg<float> fieldArg("","field","Field used to simulate the images (enabled: 1.5, 3 and 7T).",true,0,&allowedFieldValues);
    cmd.add(fieldArg);

    TCLAP::ValueArg<unsigned int> alphaArg("","alpha","Flip angle (degrees).",true,0,"ALPHA",cmd);

    TCLAP::ValueArg<float> tiArg("","ti","Inversion time (sec).",true,0,"TI",cmd);

    TCLAP::ValueArg<float> teArg("","te","Echo time (sec).",true,0,"TE",cmd);
    TCLAP::ValueArg<float> esArg("","es","Echo spacing (sec).",true,0,"ES",cmd);
    TCLAP::ValueArg<float> trArg("","tr","Repetition time (sec).",true,0,"TR",cmd);

    TCLAP::ValueArg<unsigned int> nbSlicesArg("","nbSlices","Number of slices.",true,0,"NBSLICES",cmd);

    std::vector<unsigned int> allowedPF;
    allowedPF.push_back(6);
    allowedPF.push_back(7);
    allowedPF.push_back(8);
    
    TCLAP::ValuesConstraint<unsigned int> allowedPartialFourier(allowedPF);

    TCLAP::ValueArg<unsigned int> partialFourierArg("","partialFourier","Set to 6 if 6/8 slice partial fourier is used (same idea for 7/8, Default: 8/8).",true,8,&allowedPartialFourier);
    cmd.add(partialFourierArg);

    TCLAP::ValueArg<unsigned int> noiseLevelArg("","noiseLevel","Percentage of noise level (default: 0%). This is only available for simulations from pv maps",false,0,"NOISELEVEL",cmd);    
	
    try
    {
	cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
	std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
	return EXIT_FAILURE;
    }

    // Check if the simulation can be run.
    if((!wmArg.isSet() || !gmArg.isSet() || !csfArg.isSet()) && !t1Arg.isSet())
    {
	std::cerr << "Error: If no T1 Map provided, wm, gm and csf segmentations needed." << std::endl;
	exit(EXIT_FAILURE);
    }
    
    // Constants
    const unsigned int Dimension = 3;

    // Typedefs
    typedef itk::Image<float,Dimension> FloatImageType;
    
    typedef itk::ImageFileReader <FloatImageType> FloatReaderType;
    typedef itk::ImageFileWriter <FloatImageType> FloatWriterType;

    typedef flaws::MPRAGESimulationImageFilter <FloatImageType,FloatImageType> MPRAGESimulationFilterType;
    
    //Read images
    FloatReaderType::Pointer wmReader = FloatReaderType::New();
    wmReader->SetFileName(wmArg.getValue());

    FloatReaderType::Pointer gmReader = FloatReaderType::New();
    gmReader->SetFileName(gmArg.getValue());

    FloatReaderType::Pointer csfReader = FloatReaderType::New();
    csfReader->SetFileName(csfArg.getValue());

    FloatReaderType::Pointer t1Reader = FloatReaderType::New();
    t1Reader->SetFileName(t1Arg.getValue());

    FloatReaderType::Pointer b1PlusReader = FloatReaderType::New();
    b1PlusReader->SetFileName(b1PlusArg.getValue());

    FloatReaderType::Pointer b1MinusReader = FloatReaderType::New();
    b1MinusReader->SetFileName(b1MinusArg.getValue());
    
    try
    {
	if(!t1Arg.isSet())
	{
	    wmReader->Update();
	    gmReader->Update();
	    csfReader->Update();
	}
	else
	{
	    t1Reader->Update();
	}
	
	if(b1PlusArg.isSet())
	{
	    b1PlusReader->Update();
	}

	if(b1MinusArg.isSet())
	{
	    b1MinusReader->Update();
	}
    }
    catch(itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return EXIT_FAILURE;
    }
    
    MPRAGESimulationFilterType::Pointer mprageSimulator = MPRAGESimulationFilterType::New();

    if(!t1Arg.isSet())
    {
	mprageSimulator->SetWMMap(wmReader->GetOutput());
	mprageSimulator->SetGMMap(gmReader->GetOutput());
	mprageSimulator->SetCSFMap(csfReader->GetOutput());
    }
    else
    {
	mprageSimulator->SetWMMap(t1Reader->GetOutput()); // Avoid error due to the absence of input image
	mprageSimulator->SetT1Map(t1Reader->GetOutput());
    }
	
    if(b1PlusArg.isSet())
    {
	mprageSimulator->SetUseB1PlusMap(true);
	mprageSimulator->SetB1PlusMap(b1PlusReader->GetOutput());
    }

    if(b1MinusArg.isSet())
    {
	mprageSimulator->SetUseB1MinusMap(true);
	mprageSimulator->SetB1MinusMap(b1MinusReader->GetOutput());
    }
    
    mprageSimulator->SetField(fieldArg.getValue());
    
    mprageSimulator->SetAlpha(alphaArg.getValue());

    mprageSimulator->SetInversionTime(tiArg.getValue());

    mprageSimulator->SetEchoTime(teArg.getValue());
    mprageSimulator->SetEchoSpacing(esArg.getValue());
    mprageSimulator->SetRepetitionTime(trArg.getValue());

    mprageSimulator->SetNumberOfSlices(nbSlicesArg.getValue());
    mprageSimulator->SetPartialFourier(partialFourierArg.getValue());

    mprageSimulator->SetNoiseLevel(noiseLevelArg.getValue());

    try
    {
    	mprageSimulator->Update();
    }
    catch (itk::ExceptionObject &e)
    {
    	std::cerr << e << std::endl;
    	return EXIT_FAILURE;
    }
    
    FloatWriterType::Pointer outputWriter = FloatWriterType::New();	
    outputWriter->SetFileName(outputArg.getValue());		
    outputWriter->SetUseCompression(true);
    outputWriter->SetInput(mprageSimulator->GetOutput());
  
    try
    {
	outputWriter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
    	std::cerr << e << std::endl;
    	return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
