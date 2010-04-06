/**
 * Class that encapsulates the main function
 */
package BAMStats;

/**
 * @author niravs
 * Class to drive this application.
 */
  public class Driver
  {
	/**
	 * Main method - entry point of the application
	 * @param args - Command line parameter strings
	 */
    public static void main(String[] args)
    {
      if(args.length < 2 || args.length > 3)
      {
        printUsage();
        System.exit(-1);
      }
      try
      {
    	if(!args[1].equals("1") && !args[1].equals("2") && !args[1].equals("3") ||
    	   ((args.length == 3) && (args[2].compareToIgnoreCase("solid") != 0)))
        {
    	  printUsage();
    	  System.exit(-1);
        }
    	
    	ReadType rType = ReadType.Read1;
    	if(args[1].equals("2"))
    	{
          rType = ReadType.Read2;
    	}
    	else
    	if(args[1].equals("3"))
    	{
          rType = ReadType.Fragment;
    	}
    	
        BAMFile bf = new BAMFile(args[0], rType);
        if(args.length < 3)
        {
          bf.computeStatistics("");
        }
        else
        {
          bf.computeStatistics("solid");
        }
      }
      catch(Exception e)
      {
        System.out.println("Exception encountered : " + e.getMessage());
        e.printStackTrace();
        System.exit(-2);
      }
    }
	
    /**
     * Helper method to explain command line parameters
     */
    private static void printUsage()
    {
      System.out.println("Usage:");
	  System.out.println("java -jar BAMStats.jar FileName ReadType Sequencer");
	  System.out.println("\tFileName             - BAM/SAM file to be analyzed");
	  System.out.println("\tReadType             - 1:Read1, 2:Read2, 3:Fragment");
      System.out.println("\tSequencer {Optional} - Please specify \"solid\" for");
      System.out.println("\t                     - solid reads, omit for other types");
    }
}
