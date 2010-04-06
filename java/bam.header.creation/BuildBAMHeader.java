//package RegenerateBAMHeader;
import java.io.*;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import java.util.List;

/**
 * Class to re-generate BAM header from specified command line parameters
 * @author niravs
 *
 */
public class BuildBAMHeader
{
  private InputParameters inputParam; // InputParameters reference
  private SAMFileHeader header;       // Header of SAM/BAM file
  private SAMFileReader reader;       // To read file whose header needs to be fixed
    
  /**
   * Class constructor
   * @param inputParam - Instance of InputParameters class
   */
  public BuildBAMHeader(InputParameters inputParam)
  {
    this.inputParam = inputParam;
    SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
    reader = new SAMFileReader(new File(inputParam.inputFile));
  }
  
  /**
   * Helper method to regenerate header by adding sequence dictionary
   * and possibly add RG tag
   */
  public void regenerateHeader()
  {
    try
    {
      System.out.println("Building header");
      BuildSeqDict seqDict = new BuildSeqDict(inputParam);
      header = reader.getFileHeader();
      
      List<SAMReadGroupRecord> listRG = header.getReadGroups();
      
      // Current behavior is to abort if read group is not present in header.
      if(listRG == null || listRG.size() <= 0)
      {
        throw new Exception("Specified BAM/SAM file : " + inputParam.inputFile + " does not have RG records");
      }
      header.setSequenceDictionary(seqDict.makeSequenceDictionary());
      System.out.println("Header built");
      buildOutputFile();
    }
    catch(Exception e)
    {
      reader.close();
      System.err.println(e.getMessage());
      e.printStackTrace();
    }
  }
  
  /**
   * Private helper method to generate output file
   */
  private void buildOutputFile()
  {
    File f = new File(inputParam.outputFile);
    final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header,true, f);
    
    System.out.println("Copying reads");
    long idx = 0;
    for (final SAMRecord samRecord : reader)
    {
      writer.addAlignment(samRecord);
      idx++;
      
      if(idx % 1000000 == 0)
      {
        System.err.print("\r" + idx);
      }
    }
    writer.close();
    reader.close();
  }
}
