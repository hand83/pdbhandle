# pdbhandle
Content: lib_pdbhandle.py

Library to handle pdb protein structure files in python.
Uses Modeller 9.19 (https://salilab.org/modeller) and TM-align (2016/05/21, https://zhanglab.ccmb.med.umich.edu/TM-align/) programs.

Classes:
  pdb_data
    Extract sequence information from a pdb file.
    Properties:
      PDBID - protein structure identifier
      CHAIN - protein chain identifier
      PDBPATH - pdb file path
      ALLAMINO - aminoacid information with spatial and sequence positions
      FULLSEQ - complete protein sequence
      PDBSEQ - sequence with visible aminoacids only, "-" dashes indicate missing amioacids
      PDBSEQ_HTM - sequence with visible aminoacids only, "." dots indicate heteroatoms
     Functions:
      genModelSeq(ul) - fills up missing loops in the structure up to ul loop size limit (amioacids length)

Functions:
  GenModel(pdb_data, [destination], [ul], [pir])
    Creates model using a pdb_data object as template.
    The function uses Modeller program.
    Writes model to file if destination is provided. Otherwise returns the content of the model in pdb format.
    Fills and models up missig loops in the template up to ul size limit (aminoacids length).
    Custom alignment file can be provided in pir format. [pir] stands for the path of the file.
    Custom alignment file should contain the pdbid, chain and sequence of the template pdb in pir format.
    Any kind of target sequence can be specified in the pir file to be modelled on the template structure.
    
  RotateStruct(pdb_data, mtx, [allchain], [destination])
    Rotates the pdb coordinates of the given pdb_data objects according to the mtx rotation matrix.
    The matrix is a 4x3 matrix as [translation, x, y, z][x, y z] values.
    By default, only the selected chain will be rotated (specified by the CHAIN property). Set allchain = True to rotate all chains.
    By default, the function returns the structure content in pdb format. Specifying the destination lets the fuction to write it into a file.
    
  alignStruct(pdb_data, pdb_data_ref, [fasta], [destination])
    Aligns a pdb_data object to an other pdb_data reference object.
    This is a rotation function combined with a spatial alignment.
    The function uses TM-Align program.
    Sequence-wise alignment can be specified in fasta format.
    By default, the function returns the structure content in pdb format. Specifying the destination lets the fuction to write it into a file.
    
  GetRMSD(pdb_data, pdb_data, [fasta])
    Aligns a pdb_data objecto to an other pdb_data object and returns the Root Mean Squared Deviation of alpha-carbon atoms and the TM-Align score.
    The function uses TM-Align program.
    Sequence-wise alignment can be specified in fasta format.
    
  get_uniprot(pdbid_chain)
    A function that maps a protein chain to its Uniprot identifier.
    The server used by the function is based on the PDB/UniProt Mapping Query the database of Dr. Andrew C.R. Martin, University College London: http://www.bioinf.org.uk/pdbsws/index.html
    Reference: Martin, Andrew C. R. (2005) Mapping PDB chains to UniProtKB entries, Bioinformatics, 21:4297-4301
    
