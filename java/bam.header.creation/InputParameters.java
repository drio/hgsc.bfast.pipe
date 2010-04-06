import java.util.ArrayList;

/**
 * @author niravs
 * Class that encapsulates parameters passed at the command line
 */
public class InputParameters
{
  String refType = "";            // Reference Type
  String inputFile = null;        // BAM / SAM file to read
  String outputFile = null;       // BAM / SAM file to generate
  Reference reference = null;     // Reference instance
  private ReferenceReader refList = null;
  
  /**
   * Class constructor - parse and validate command line arguments
   * @param cmdLine
   * @throws Exception
   */
  public InputParameters(String[] cmdLine) throws Exception
  {
    String [] paramValue;
 
    for(int i = 0; i < cmdLine.length; i++)
    {
      paramValue = cmdLine[i].split("=", 2);

      if(paramValue[0].compareToIgnoreCase("type") == 0)
      {
        refType = paramValue[1];
      }
      else
      if(paramValue[0].compareToIgnoreCase("I") == 0)
      {
        inputFile = paramValue[1];
      }
      else
      if(paramValue[0].compareToIgnoreCase("O") == 0)
      {
        outputFile = paramValue[1];
      }
    }
    validateParams();
  }
  
  /**
   * Private helper method to validate the parameters.
   */
  private void validateParams() throws Exception
  {
    boolean errorFound = false;
    refList = new ReferenceReader();
    
    if(inputFile == null)
    {
      System.err.println("Error : Parameter \"I\" (InputFile) MUST be specified");
      errorFound = true;
    }
    if(outputFile == null)
    {
      System.err.println("Error : Parameter \"O\" (OutputFile) MUST be specified");
      errorFound = true;
    }
    if(refType == null || refType.isEmpty())
    {
      System.err.println("Error : Parameter \"Type\" MUST be specified");
      errorFound = true;
    }
    
    reference = refList.getReference(refType);
    if(reference == null)
    {
      System.err.println("Did not find a reference for : " + refType);
      errorFound = true;
    }
    if(errorFound == true)
    {
      printUsage();
      System.exit(-2);
    }
    else
    {
      System.out.println("    Input File  : " + inputFile);
      System.out.println("    Output File : " + outputFile);
    }
  }
  
  /**
   * Helper method to display usage information
   */
  private void printUsage()
  {
    System.err.println();
    System.err.println("Usage:");
    System.err.println("RegenerateBAM.jar I=InputFile O=OutputFile Type=SampleType");
    System.err.println("    I     : Input BAM/SAM file to read");
    System.err.println("    O     : Ouptut BAM/SAM file to generate");
    System.err.println("    Type  : Type of Sample");
    
    ArrayList<String> refNames = refList.getAllReferenceNames();
    
    if(refNames.size() > 0)
    {
      System.out.println("            Allowed Values are : ");
      
      for(int i = 0; i < refNames.size(); i++)
      {
        System.out.print(refNames.get(i) + " ");
        if(i > 0 && i % 3 == 0)
        {
          System.out.println("");
        }
      }
    }
  }
  
  /**
   * For debugging purpose
   */
  void showParameters()
  {
    System.out.println("Parameters   : ");
    System.out.println("Input File   : " + inputFile);
    System.out.println("Output File  : " + outputFile);
  }
}
