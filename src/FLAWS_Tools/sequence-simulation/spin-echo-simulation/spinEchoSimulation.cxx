#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <spinEchoSimulation.h>

#include <tclap/CmdLine.h>

int main(int argc, const char** argv)
{
    // Parsing arguments
    TCLAP::CmdLine cmd("FLAWS Tools - Spin echo simulation", ' ',FLAWS_TOOLS_VERSION);

    TCLAP::ValueArg<std::string> wmArg("","wmMap","Input WM map",true,"","WMMAP",cmd);
    TCLAP::ValueArg<std::string> gmArg("","gmMap","Input GM map",true,"","GMMAP",cmd);
    TCLAP::ValueArg<std::string> csfArg("","csfMap","Input CSF map",true,"","CSFMAP",cmd);
    TCLAP::ValueArg<std::string> fatArg("","fatMap","Input fat map",false,"","FATMAP",cmd);
    
    TCLAP::ValueArg<std::string> outputArg("o","output","Output image.",true,"","OUTPUT",cmd);

    std::vector<float> allowedField;
    allowedField.push_back(1.5);
    allowedField.push_back(3);
    allowedField.push_back(7);

    TCLAP::ValuesConstraint<float> allowedFieldValues(allowedField);

    TCLAP::ValueArg<float> fieldArg("","field","Field used to simulate the images (enabled: 1.5, 3 and 7T).",true,0,&allowedFieldValues);
    cmd.add(fieldArg);

    TCLAP::ValueArg<float> teArg("","te","Echo time (sec).",true,0,"TE",cmd);
    TCLAP::ValueArg<float> trArg("","tr","Repetition time (sec).",true,0,"TR",cmd);   
	
    try
    {
	cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
	std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
	return EXIT_FAILURE;
    }
    
    // Constants
    const unsigned int Dimension = 3;

    // Typedefs
    typedef itk::Image<float,Dimension> FloatImageType;
    
    typedef itk::ImageFileReader <FloatImageType> FloatReaderType;
    typedef itk::ImageFileWriter <FloatImageType> FloatWriterType;

    typedef flaws::SpinEchoSimulationImageFilter <FloatImageType,FloatImageType> SpinEchoSimulationFilterType;
    
    //Read images
    FloatReaderType::Pointer wmReader = FloatReaderType::New();
    wmReader->SetFileName(wmArg.getValue());

    FloatReaderType::Pointer gmReader = FloatReaderType::New();
    gmReader->SetFileName(gmArg.getValue());

    FloatReaderType::Pointer csfReader = FloatReaderType::New();
    csfReader->SetFileName(csfArg.getValue());

    FloatReaderType::Pointer fatReader = FloatReaderType::New();
    fatReader->SetFileName(fatArg.getValue());
    
    try
    {
	wmReader->Update();
	gmReader->Update();
	csfReader->Update();

	if(fatArg.isSet())
	{
	  fatReader->Update();
	} 
    }
    catch(itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return EXIT_FAILURE;
    }

    SpinEchoSimulationFilterType::Pointer spinEchoSimulator = SpinEchoSimulationFilterType::New();

    spinEchoSimulator->SetWMMap(wmReader->GetOutput());
    spinEchoSimulator->SetGMMap(gmReader->GetOutput());
    spinEchoSimulator->SetCSFMap(csfReader->GetOutput());

    if(fatArg.isSet())
    {
      spinEchoSimulator->SetSimuFat(true);
      spinEchoSimulator->SetFatMap(fatReader->GetOutput());
    } 
    
    spinEchoSimulator->SetField(fieldArg.getValue());

    spinEchoSimulator->SetEchoTime(teArg.getValue());
    spinEchoSimulator->SetRepetitionTime(trArg.getValue());

    try
    {
    	spinEchoSimulator->Update();
    }
    catch (itk::ExceptionObject &e)
    {
    	std::cerr << e << std::endl;
    	return EXIT_FAILURE;
    }
        
    FloatWriterType::Pointer spinEchoWriter = FloatWriterType::New();	
    spinEchoWriter->SetFileName(outputArg.getValue());		
    spinEchoWriter->SetUseCompression(true);
    spinEchoWriter->SetInput(spinEchoSimulator->GetOutput());
  
    try
    {
    	spinEchoWriter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
    	std::cerr << e << std::endl;
    	return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
