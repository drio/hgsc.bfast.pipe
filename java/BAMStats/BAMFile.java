/**
 * Package encapsulating the project that reads a BAM file and generates 
 * statistics on reads of the BAM file
 */
package BAMStats;

import java.io.*;
import java.io.FileWriter;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;


/**
 * The class that encapsulates the functions on a BAM file
 *
 */
public class BAMFile {
  private String inputFileName;	    // Input Filename
  private ReadType rType;			// Read-type
  private Statistics fileStats;     // To generate stats on the file
  
  /**
   * Class constructor that gets the name for the BAM file
   * @param inputFileName - String representing full path of BAM file 
   * @param rType - Type of read (i.e. read 1 or 2 in paired read)
   */
  public BAMFile(String inputFileName, ReadType rType)
  {
    this.inputFileName = inputFileName;
    this.rType = rType;
  }
  
  /**
   * Public helper method to compute statistics on the BAM file
   * @param sequencerType - Type of sequencer, "solid" or empty string
   */
  public void computeStatistics(String sequencerType)
  {
    SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
	
    final SAMFileReader fileReader = new SAMFileReader(new File(inputFileName));
    CloseableIterator<SAMRecord> recordList = fileReader.iterator();

    fileStats = new Statistics(rType);
    
    SAMRecord rec = null;
    
    long startTime = System.currentTimeMillis();
    while(recordList.hasNext())
    {
    	rec = null;
    	rec = recordList.next();
    	
    	if(rec != null)
    	{
    	  fileStats.processRead(rec);
    	  rec = null;
    	}
    }
    long stopTime = System.currentTimeMillis();
    
    recordList.close();
    fileReader.close();
    
    System.out.println("");
    System.out.println("BAM/SAM File : " + inputFileName);
	
    if(rType == ReadType.Read1)
    {
      if(sequencerType.isEmpty())
      {
        System.out.println("Read Type    : Read 1");
      }
      else
      {
        System.out.println("Read Type    : F3");
      }
    }
    else
    if(rType == ReadType.Read2)
    {
      if(sequencerType.isEmpty())
      {
        System.out.println("Read Type    : Read 2");
      }
      else
      {
    	System.out.println("Read Type    : R3");  
      }
    }
    else
    {
      System.out.println("Read Type    : Fragment");
    }
    reportStats();
    System.out.format("%nRunning Time      : %.3f sec%n%n", (stopTime - startTime)/1000.0);
  }
  
  /**
   * Log map quality distribution
   */
  private void reportStats()
  {
    fileStats.showStats();
    long[] mapQualDist = fileStats.getMapQualityDistribution();
    long[] nonUniqQualDist = fileStats.getNonUniqueMapQualityDistribution();
    
    String mapQualLog = inputFileName.substring(0, inputFileName.length() - 4) + "_log.txt";
    String mapNonUniqQualLog = inputFileName.substring(0, inputFileName.length() - 4) + "_nonunique_log.txt";
    
    System.out.println("");
    System.out.println("Logging Quality Distrubution of Mapped Reads to : " + mapQualLog);
    logMapQualDistData(mapQualLog, mapQualDist);
    
    System.out.println("Logging Quality Distrubution of Mapped Reads With Multiple Alignments  to : " + mapNonUniqQualLog);
    logMapQualDistData(mapNonUniqQualLog, nonUniqQualDist);
  }
  
  /**
   * Helper method to log distribution data to file
   * @param fileName - Log file name
   * @param data     - Map quality distribution data 
   */
  private void logMapQualDistData(String fileName, long[] data)
  {
    try
    {
      FileWriter fw = new FileWriter(fileName);
      
      for(int i = 0; i < data.length; i++)
      {
        fw.write(i + " " + data[i] + "\r\n");
      }
      fw.close();
    }
    catch(Exception e)
    {
      System.err.println(e.getMessage());
      e.printStackTrace();
    }
  }
}
