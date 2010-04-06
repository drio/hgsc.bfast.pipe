/**
 * Class that encapsulates the "main" method
 */
import java.io.*;
import java.util.*;

/**
 * @author niravs
 *
 */
public class Driver
{

 /**
  * @param args
  */
  public static void main(String[] args) 
  {
    if(args.length != 3)
    {
      printUsage();
      System.exit(-1);
    }
    try
    {
      InputFiles in = new InputFiles(args[1]);
      ArrayList<String>csFastaFiles = in.getCSFastaFiles();
    
      if(csFastaFiles.size() < 1)
      {
        throw new Exception("Could not find any CSFasta file within " + args[1]);
      }
      ValidateBAM val = null;
      
      if(csFastaFiles.size() > 1)
      {
        val = new ValidateBAM(args[0], csFastaFiles.get(0), csFastaFiles.get(1),
                              args[2]);
      }
      else
      {
        val = new ValidateBAM(args[0], csFastaFiles.get(0), "", args[2]);
      }
      boolean isValid = val.vaildate();
    
      if(isValid == true)
      {
        System.exit(0);
      }
      else
      {
        System.exit(1);
      }
    }
    catch(Exception e)
    {
      System.err.println(e.getMessage());
      e.printStackTrace();
      System.exit(-2);
    }
  }
 
  /**
   * Private helper method to show usage information
   */
  private static void printUsage()
  {
    System.err.println("Usage:");
    System.err.println("java -jar jarFileName BAMFile InputDir LogFile");
    System.err.println("  BAMFile  - BAM/SAM file generated");
    System.err.println("  InputDir - Parent directory of csfasta files");
    System.err.println("             (Files used to generate BAMFile)");
    System.err.println("  LogFile  - File to write output to");
  }
}
