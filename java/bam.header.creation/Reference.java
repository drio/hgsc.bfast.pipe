/**
 * Class to encapsulate a reference
 */

/**
 * @author niravs
 *
 */
public class Reference
{
  private String name;           // Name of the reference
  private String path;           // Path to fasta file
  private String genomeAssembly; // Value for GA tag
  private String species;        // Value for SP tag
  
  /**
   * Class Constructor
   * @param name - Reference name
   * @param path - Path to reference fasta file
   * @param ga   - Genome Assembly, value for GA tag
   * @param sp   - Species, value for SP tag 
   */
  public Reference(String name, String path, String ga, String sp)
  {
    this.name = name;
    this.path = path;
    this.genomeAssembly = ga;
    this.species = sp;
  }
  
  public String getName()
  {
    return name;
  }
  
  public String getPath()
  {
    return path;
  }
  
  public String getGenomeAssembly()
  {
    return genomeAssembly;
  }
  
  public String getSpecies()
  {
    return species;
  }
}
