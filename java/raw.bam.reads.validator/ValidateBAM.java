/**
 * This class validates that the BAM file has the same number of reads
 * as the input file (namely CSFASTA).
 */

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import java.io.*;

/**
 * @author niravs
 *
 */
public class ValidateBAM
{
  private File bamFile = null;        // BAM/SAM file handle
  private File inputFileRead1 = null; // Input file (Read 1/ fragment) handle
  private File inputFileRead2 = null; // Input file (Read 2) handle
  private String bamName;             // Names of corresponding files
  private String inputRead1Name;
  private String inputRead2Name;
  private String outputFileName;
  
  /**
   * Class Constructor
   * @param bamFile
   * @param inputRead1
   * @param inputRead2
   * @throws Exception
   */
  public ValidateBAM(String bamFile, String inputRead1, String inputRead2,
                     String outputFile) throws Exception
  {
    this.bamName = bamFile;
    this.inputRead1Name = inputRead1;
    this.inputRead2Name = inputRead2;
    this.outputFileName = outputFile;
  
    inputFileRead1 = validateFile(inputRead1Name);
    this.bamFile = validateFile(bamName);
        
    if(inputRead2 != null && !inputRead2.isEmpty())
    {
      inputFileRead2 = validateFile(inputRead2Name);
    }
  }
  
  /**
   * Validate the BAM File with Input CSFasta file and show output to display
   * and a log file
   * @throws Exception
   */
  public boolean vaildate() throws Exception
  {
    BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
    boolean isValid = true;
 
    long numReadsFile1 = 0;
    long numReadsFile2 = 0;
    long numReadsBAM = 0;
    long startTime = 0;
    long stopTime = 0;
    
    System.out.println("Calculating Reads in " + inputRead1Name);
    
    startTime = System.currentTimeMillis();
    numReadsFile1 = countReadsCsFasta(inputFileRead1);
    stopTime = System.currentTimeMillis();
    System.out.println("Time taken : " + (stopTime - startTime)/1000.0 + "sec");
    
    if(inputRead2Name != null && !inputRead2Name.isEmpty())
    {
      System.out.println("Calculating Reads in " + inputRead2Name);
      startTime = System.currentTimeMillis();
      numReadsFile2 = countReadsCsFasta(inputFileRead2);
      stopTime = System.currentTimeMillis();
      System.out.println("Time taken : " + (stopTime - startTime)/1000.0 + "sec");
    }
    
    System.out.println("Calculating Reads in " + bamName);
    
    startTime = System.currentTimeMillis();
    numReadsBAM = countReadsInBAM();
    stopTime = System.currentTimeMillis();
    System.out.println("Time taken : " + (stopTime - startTime)/1000.0 + "sec");
    
    if(numReadsBAM != (numReadsFile1 + numReadsFile2))
    {
      isValid = false;
    }
    if(isValid == true)
    {
      writer.write("OK");
      System.out.println("");
      System.out.println("");
      System.out.println("OK");
    }
    else
    {
      writer.write("ERROR");
      System.out.println("ERROR");
    }
    writer.newLine();
    
    writer.write("CsFasta File : " + inputRead1Name);
    writer.newLine();
    System.out.println("CsFasta File : " + inputRead1Name);
    writer.write("Num Reads in CsFasta : " + numReadsFile1);
    writer.newLine();
    System.out.println("Num Reads in CsFasta : " + numReadsFile1);
    
    if(inputRead2Name != null && !inputRead2Name.isEmpty())
    {
      writer.write("CsFasta File : " + inputRead2Name);
      writer.newLine();
      writer.write("Num Reads in CsFasta : " + numReadsFile2);
      writer.newLine();
      System.out.println("CsFasta File : " + inputRead2Name);
      System.out.println("Num Reads in CsFasta : " + numReadsFile2);
    }
    writer.write("BAM File : " + bamName);
    System.out.println("BAM File : " + bamName);
    writer.newLine();
    writer.write("Num Reads in BAM File : " + numReadsBAM);
    writer.newLine();
    System.out.println("Num Reads in BAM File : " + numReadsBAM);
    
    writer.close();

    return isValid;
  }
  
  /**
   * Verify that the specified file is existing and can be read
   * @param fileName - File to check for read permission
   * @throws - Exception if file does not exist or cannot be read
   */
  private File validateFile(String fileName) throws Exception
  {
    File file = new File(fileName);
    
    if(!file.exists() || !file.canRead())
      throw new Exception(fileName + " cannot be opened for reading");
    return file;
  }
  
  /**
   * Method to calculate the number of reads in a BAM/SAM file
   * @return - Number of reads in file
   */
  private long countReadsInBAM()
  {
    long numReads = 0;
    
    SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
    SAMFileReader reader = new SAMFileReader(bamFile);
    
    for (final SAMRecord samRecord : reader)
    {
      numReads++;
      if(numReads % 1000000 == 0)
      {
        System.err.print("\r" + numReads);
      }
    }
    reader.close();
    return numReads;
  }
  
  /**
   * Method to calculate the number of reads in input CSFASTA file
   * @return - Number of reads in file
   */
  private long countReadsCsFasta(File inputFile) throws Exception
  {
    long numReads = 0;
    String line;
    
    BufferedReader reader = new BufferedReader(new FileReader(inputFile));
    
    while((line = reader.readLine()) != null)
    {
      if(line.startsWith(">"))
      {
        numReads++;
        if(numReads % 1000000 == 0)
        {
          System.err.print("\r" + numReads);
        }
      }
    }
    reader.close();
    return numReads;
  }
}
