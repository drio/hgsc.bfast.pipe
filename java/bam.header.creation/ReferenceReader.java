/**
 * Class to encapsulate reading reference XML and generate an
 * instance of reference class
 */

import java.io.*;
import java.util.*;
import org.w3c.dom.Document;
import org.w3c.dom.*;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;

/**
 * @author niravs
 * Class that reads ReferenceList.xml and build "Reference" class instance
 */
public class ReferenceReader
{
  // List of reference objects built from XML file
  private ArrayList<Reference> referenceList = null;
  
  /**
   * Class Constructor
   * Reads the ReferenceList.xml packed as a resource with the JAR file.
   * @throws Exception
   */
  public ReferenceReader() throws Exception
  {
    //File f = new File(referenceFile);
    InputStream inStream = ReferenceReader.class.getResourceAsStream("/ReferenceList.xml");
    if(inStream == null)
    {
      throw new Exception("Error : Could not load ReferenceList.xml");
    }
    referenceList = new ArrayList<Reference>();
    parseXMLFile(inStream);
  }
  
  /**
   * Private Helper method to read XML file and build List of Reference
   * objects
   * @param xmlFile - Instance of File
   * @throws Exception
   */
  private void parseXMLFile(InputStream xmlFile) throws Exception
  {
	  System.out.println("Reading Reference File");
    DocumentBuilderFactory docBF = DocumentBuilderFactory.newInstance();
    DocumentBuilder docBuilder = docBF.newDocumentBuilder();
    Document doc = docBuilder.parse(xmlFile);
    
    doc.getDocumentElement().normalize();
    NodeList refList = doc.getElementsByTagName("Reference");
    
    for(int i = 0; i < refList.getLength(); i++)
    {
      Node n = refList.item(i);
      
      if(n.getNodeType() == Node.ELEMENT_NODE)
      {
        Element ref = (Element)n;
        String refName = ref.getAttribute("Name");
        String path = ref.getAttribute("Path");
        String ga = ref.getAttribute("GenomeAssembly");
        String species = ref.getAttribute("Species");
        System.out.println("Found Reference : " + refName);
        Reference r = new Reference(refName, path, ga, species);
        referenceList.add(r);
        r = null;
      }
    }
  }
  
  /**
   * Return the Reference object corresponding to the name
   * If the specified name does not match any in the reference list,
   * return null.
   * @param refName - "Name" of the reference
   * @return - Reference object or null
   */
  public Reference getReference(String refName)
  {
    Reference temp = null;
    
    for(int i = 0; i < referenceList.size(); i++)
    {
      temp = referenceList.get(i);
      
      if(temp.getName().equalsIgnoreCase(refName))
      {
        return temp;
      }
    }
    return null;
  }
  
  /**
   * Return the list of all reference names that could be used
   * as values for cmd line parameter "Type"
   * @return - String of reference names
   */
  public ArrayList<String>getAllReferenceNames()
  {
    ArrayList<String> result = new ArrayList<String>();
    
    for(int i = 0; i < referenceList.size(); i++)
    {
      result.add(referenceList.get(i).getName());
    }
    return result;
  }
}
