# RNA-Seq Scripts

 ** filterBwaMappings.sh **
 
 Input: <SAM FILE>
      SAM file with BWA alignments sorted by read name.
 Output: <SAM FILE>
      One alignment per Read. For multi-mappings, select the best alignment score (field AS).
 
 
 ** filterGTF.sh **
 
 Input: <GTF File> <function>  
        <function> : median; max; min;
 Output: Three column file
      <gene id> <transcript length> <transcript_id>
      Obtain the minimum, median or maximum transcript length.

** selectGtf.sh **
Input: <GTF File>
Output: Three column file
      <gene id> <transcript length> <transcript_id>

** statistics_bwa_single_end.sh **
Input: <SAM FILE>
      SAM file with BWA alignments.
Output: Mapping Statistics
