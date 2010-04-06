import java.io.File;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.StringUtil;

/**
 * @author niravs
 * Class to build sequence dictionary
 */
public class BuildSeqDict
{
  private String referenceName;          // reference path
  private MessageDigest md5 = null;      // MD5
  private InputParameters inputParam;    // InputParameters instance
  private String URI;                    // URI of reference
  private String genomeAssembly;         // The AS tag
  private String species;                // The species (SP) tag
  
  /**
   * Class constructor - Generate sequence dictionary based on specified
   * "type" input parameter
   * @param inputParam - Instance of InputParameters class
   * @throws Exception
   */
  public BuildSeqDict(InputParameters inputParam) throws Exception
  {
    try 
    {
      md5 = MessageDigest.getInstance("MD5");
      URI = null;
      genomeAssembly = null;
      species = null;
    }
    catch(NoSuchAlgorithmException e)
    {
      throw new Exception("Cannot create MessageDigest MD5 instance"); 
    }
    this.inputParam = inputParam;
    setReference();
  }
  
  /**
   * Private helper method to select correct reference path
   */
  private void setReference() throws Exception
  {
    if(inputParam.reference != null)
    {
      referenceName = inputParam.reference.getPath();
      species = inputParam.reference.getSpecies();
      genomeAssembly = inputParam.reference.getGenomeAssembly();
    }
    else
    {
      throw new Exception("Unknown reference specified");
    }
  }
  
  /**
   * Read all the sequences from the given reference file, and convert into SAMSequenceRecords
   * @param referenceFile fasta or fasta.gz
   * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
   * @throws Exception 
   */
  public SAMSequenceDictionary makeSequenceDictionary() throws Exception
  {
    final File referenceFile = new File(referenceName);
    final ReferenceSequenceFile refSeqFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, true);

    ReferenceSequence refSeq;
      
    if(URI == null)
    {
      URI = "file:" + referenceFile.getAbsolutePath();
    }
      
    final List<SAMSequenceRecord> ret = new ArrayList<SAMSequenceRecord>();
    final Set<String> sequenceNames = new HashSet<String>();
      
    for (int numSequences = 0; (refSeq = refSeqFile.nextSequence()) != null; ++numSequences)
    {
      if (sequenceNames.contains(refSeq.getName()))
      {
        throw new Exception("Sequence name appears more than once in reference: " + refSeq.getName());
      }
      sequenceNames.add(refSeq.getName());
      ret.add(makeSequenceRecord(refSeq));
    }
    return new SAMSequenceDictionary(ret);
  }
  
  /**
   * Private helper method to generate sequence dictionary
   * @param refSeq
   * @return
   */
  private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) 
  {
    final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName(), refSeq.length());

    // Compute MD5 of upcased bases
    final byte[] bases = refSeq.getBases();
    for (int i = 0; i < bases.length; ++i) 
    {
      bases[i] = StringUtil.toUpperCase(bases[i]);
    }

    ret.setAttribute(SAMSequenceRecord.MD5_TAG, md5Hash(bases));
    if (genomeAssembly != null)
    {
      ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, genomeAssembly);
    }
      
    ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
    if (species != null)
    {
      ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, species);
    }
    return ret;
  }
  
  /**
   * Generate MD5 hash
   * @param bytes
   * @return
   */
  private String md5Hash(final byte[] bytes) 
  {
    md5.reset();
    md5.update(bytes);
    String s = new BigInteger(1, md5.digest()).toString(16);
    if (s.length() != 32)
    {
      final String zeros = "00000000000000000000000000000000";
      s = zeros.substring(0, 32 - s.length()) + s;
    }
    return s;
  }
}
