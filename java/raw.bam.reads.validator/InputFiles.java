/**
 * Class to search for input files
 * It uses breadth-first search to find .csfasta files within the specified
 * directory
 */
import java.io.File;
import java.util.*;
/**
 * @author niravs
 * Class to search and return csfasta files
 */
public class InputFiles
{
  private String dirPath = null;                // Directory path
  private ArrayList<String> targetFiles = null; // List of csfasta files
  
  /**
   * Class constructor
   * @param path - Parent directory to search for csfasta files
   */
  public InputFiles(String path)
  {
    this.dirPath = path;
    targetFiles = new ArrayList<String>();
  }

  /**
   * Method to return all csfasta files using BFS search
   * This method uses breadth-first search to find csfasta files within
   * the specified parent directory
   * @return - ArrayList<String> of csfasta files 
   */
  public ArrayList<String> getCSFastaFiles()
  {
    ArrayList<String> childQueue = new ArrayList<String>();
    String children[] = new File(dirPath).list();
    
    for(int i = 0; i < children.length; i++)
    {
      File temp = new File(dirPath + File.separator + children[i]);
      if(temp.isDirectory() && temp.canRead())
      {
        childQueue.add(temp.getAbsolutePath());
      }
      else
      if(temp.exists() && temp.isFile() &&
         children[i].toLowerCase().endsWith(".csfasta"))
      {
        targetFiles.add(temp.getAbsolutePath());
      }
    }
    
    while(!childQueue.isEmpty())
    {
      File temp = new File(childQueue.remove(0));
      if(temp.isDirectory() && temp.canRead())
      {
        children = temp.list();
        for(int i = 0; i < children.length; i++)
        {
          childQueue.add(temp.getAbsolutePath() + File.separator + children[i]);
        }
        children = null;
      }
      else
      if(temp.exists() && temp.isFile() &&
         temp.getAbsolutePath().toLowerCase().endsWith(".csfasta"))
      {
        targetFiles.add(temp.getAbsolutePath());
      }
    }
    return targetFiles;
  }
}

