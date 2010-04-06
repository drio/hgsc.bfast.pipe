//package RegenerateBAMHeader;
/**
 * Class that acts as entry point of the application
 * @author niravs
 */

public class Driver
{
  public static void main(String args[])
  {
    try
    {
      InputParameters ip = new InputParameters(args);
      BuildBAMHeader bh = new BuildBAMHeader(ip);
      bh.regenerateHeader();
    }
    catch(Exception e)
    {
      System.err.println(e.getMessage());
      e.printStackTrace();
    }
  }
}
