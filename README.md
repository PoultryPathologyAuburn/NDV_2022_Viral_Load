# NDV_2022_Viral_Load

### Workflow
script: <a href="NDV_viralload.sh">NDV_viralload.sh</a>
<ol>
  <li>Indexing reference genomes: </li> <ul>
    <li>bwa index</li> 
    <li>samtools faidx</li>
    <li>Picard CreateSequenceDictionary</li>
    <li>NCBI reference genomes: AF077761.1.fna (LaSota 1999), ASM283408.1.fna (Chinese ducks 2018), ASM478661.1.fna (Ireland chickens 2023)</li></ul>
  
  <li>Trimming with Trimmomatic</li>
  
  <li>Alignment against reference: </li><ul>
    <li>bwa mem</li>
    <li>samtools sort</li></ul>
    
  <li>Statistics: </li><ul>
    <li>samtools flagstat</li>
    <li>samtools depth</li>
    <li>HTSeq (viral load)</li>
    <li>viral load/total reads</li></ul>
</ol>
