#! /usr/bin/python

import os, gzip, shutil, urllib2
from modeller import *
from modeller.automodel import *
from subprocess import Popen, PIPE



AA2A = {
    "ALA" : "A",
    "CYS" : "C",
    "ASP" : "D",
    "GLU" : "E",
    "PHE" : "F",
    "GLY" : "G",
    "HIS" : "H",
    "ILE" : "I",
    "LYS" : "K",
    "LEU" : "L",
    "MET" : "M",
    "ASN" : "N",
    "PRO" : "P",
    "GLN" : "Q",
    "ARG" : "R",
    "SER" : "S",
    "THR" : "T",
    "VAL" : "V",
    "TRP" : "W",
    "TYR" : "Y",
    # alternative aminoacids:
    "MSE" : "M",
    "MEX" : "C",
    "ABU" : "C",
    "CSO" : "C",
    "ALS" : "A",
    "FME" : "M",
    "MLY" : "K"
}



def isBreak(A, B):
# segCondition
# condition for segmenting sequence of positions
# segmented if A and B numbers are not consecutive
    return A["pos"] + 1 != B["pos"]



def isNewSegment(A, B):
# segCondition
# condition for segmenting aminoacid sequences
# segmented at "X-" or "-X" boundaries where X is an aminoacid
    boundary = A + B
    return len(boundary.strip("-")) == 1



def segmentArray(A, segCondition):
# SPLIT ARRAY AT THE SPECIFIED BOUNDARY CONDITION
# segCondition is a function to evaluate boundary
    bp = []
    en = 0
    for i in range(len(A) - 1):
        if segCondition(A[i], A[i + 1]):
            st = en
            en = i + 1
            bp.append( {"b" : st, "e" : en} )
    st = en
    bp.append( {"b" : st, "e" : len(A)} )
    return bp



def wrapSeq(sequence, W = 60):
# wrap aminoacid sequences to W line width
    L = len(sequence) / W
    wseq = ""
    for i in range(L):
        wseq += sequence[i * W : (i + 1) * W] + "\n"
    if len(sequence) % W > 0:
        wseq += sequence[L * W :]
    return wseq



class pdb_data:
# CLASS FOR PDB DATA TO ACCESS SEQUENCES, MISSING RESIDUES AND C-ALPHA POSITIONS
    # properties:
    # PDBPATH: location of the pdb file
    # PDBID: PDB identifier
    # CHAIN: 
    # 

    def __init__(self, pdbid = "unknown", chain = " ", pdbpath = ""):
        # pdbid: PDB identifier (4 alphanumeric characters)
        # chain: protein chain (1 letter or space if not given)
        # pdbpath: custom pdb file to read from

        if os.path.exists(pdbpath):
            # read from custom pdb file if exists
            self.PDBPATH = pdbpath
            with open(self.PDBPATH, "r") as f:
                c = f.readlines()
        elif "PDBPATH" in os.environ:
            # access pdb file from pdb database if accessible
            # in this case pdb database path should be accessible from environment as PDBPATH
            self.PDBPATH = os.environ["PDBPATH"] + "/" + pdbid[1:3] + "/pdb" + pdbid + ".ent.gz"
            if not os.path.exists(self.PDBPATH):
                raise Exception("MissingPDBFile", pdbid, chain, self.PDBPATH)
            f = gzip.open(self.PDBPATH, "r")
            c = f.readlines()
            f.close()
        else:
            raise Exception("MissingPDBFile", pdbid, chain, pdbpath)

        # MISSING RESIDUES
        # listed within REMARKS section if there are
        MissingAA = []
        beg = -1
        end = -1
        i = 0
        while beg == -1 and i < len(c):
            line = c[i].split()
            if len(line) == 4 and line[0] == "REMARK" and line[2] == "MISSING" and line[3] == "RESIDUES":
                beg = i + 6
            i += 1
        if beg != -1:
            i = beg
            while end == -1 and i < len(c):
                line = c[i].split()
                if len(line) == 2 and line[0] == "REMARK":
                    end = i
                i += 1
        for i in range(beg, end):
            # only selected chains [A-Za-z]
            if c[i][19] == chain:
                elem = {
                    "pos" : int(c[i][21:26]), # aminoacid position
                    "AA" : AA2A[c[i][15:18]] # aminoacid 3 letter code
                }
                if elem not in MissingAA:
                    MissingAA.append(elem)

        # STRUCTURE RESIDUES WITH C-ALPHA POSITIONS
        # filtering for chain atom, heteroatom and terminator records
        StructAA = []
        block = filter(lambda l: l[:6].strip() in ["ATOM", "HETATM", "TER"] and l[21] == chain and l[17:20] in AA2A, c)
        aminoacids = []
        positions = []
        poscont = []
        ishetatm = []
        i = 0
        while i < len(block) and block[i][:3] != "TER":
            pos = int(block[i][22:26])
            if pos not in positions:
                positions.append(pos)
                aminoacids.append(AA2A[block[i][17:20]])
                poscont.append([block[i]])
                if block[i][:6] == "HETATM":
                    ishetatm.append(True)
                else:
                    ishetatm.append(False)
            else:
                poscont[-1].append(block[i])
            i += 1
        for i in range(len(positions)):
            ca = filter(lambda l: l[13:16].strip() == "CA", poscont[i])
            if len(ca) == 0:
                cax, cay, caz = "", "", ""
            else:
                cax, cay, caz = float(ca[0][30:38]), float(ca[0][38:46]), float(ca[0][46:54])
            StructAA.append({
                "pos" : positions[i], # aminoacid position
                "AA" : aminoacids[i], # aminoacid 3 letter code
                "CAX" : cax, # c-alpha atom x-coordinate
                "CAY" : cay, # c-alpha atom y-coordinate
                "CAZ" : caz, # c-alpha atom z-coordinate
                "hetatm" : ishetatm[i] # True if the atom is heteroatom, usually refers to an alternative aminoacid
            })
            
        # COMBING MISSING AND STRUCTURE PARTS
        # mathch structure aminoacid numbering with missing aminoacid indices
        # structure aminoacid numbering might be not consecutive or successive

        # segment structure and missing amioacid lists where they are not consecutive
        S = segmentArray(StructAA, isBreak)
        M = []
        if len(MissingAA) > 0:
            M = segmentArray(MissingAA, isBreak)
        # order of segments with type label, fill first with structure segments
        segOrd = map(lambda x: {"seg" : x, "type" : "S"}, S)
        # insert mismatch segments into segOrd at certain positions
        for Mseg in M:
            # insert if the a structure segment is consecutive to a mismatch segment according to position numbers
            i = 0
            nextMatch = False
            while i < len(segOrd) and not nextMatch:
                if segOrd[i]["type"] == "S" and MissingAA[ Mseg["e"] - 1 ]["pos"] + 1 == StructAA[ segOrd[i]["seg"]["b"] ]["pos"]:
                    nextMatch = True
                else:
                    i += 1
            if nextMatch:
                segOrd = segOrd[:i] + [{"seg" : Mseg, "type" : "M"}] + segOrd[i:]
            else:
            # insert if the mismatch segment is consecutive to a structure segment according to position numbers
                i = 0
                prevMatch = False
                while i < len(segOrd) and not prevMatch:
                    if segOrd[i]["type"] == "S" and MissingAA[ Mseg["b"] ]["pos"] - 1 == StructAA[ segOrd[i]["seg"]["e"] - 1 ]["pos"]:
                        prevMatch = True
                    else:
                        i += 1
                if prevMatch:
                    if i < len(segOrd) - 1:
                        segOrd = segOrd[:i + 1] + [{"seg" : Mseg, "type": "M"}] + segOrd[i + 1:]
                    else:
                    # append to the end if this is the last structure segment
                        segOrd += [{"seg" : Mseg, "type" : "M"}]
                else:
                # append to the end if there is no preceeding or consecuting structure segment
                    segOrd += [{"seg" : Mseg, "type" : "M"}]

        # create a merged list of residues according to the segment order
        Merged = []
        for seg in segOrd:
            if seg["type"] == "M":
                Merged += MissingAA[ seg["seg"]["b"] : seg["seg"]["e"] ]
            else:
                Merged += StructAA[ seg["seg"]["b"] : seg["seg"]["e"] ]
        self.ALLAMINO = Merged

        # GETTING SEQUENCES
        self.FULLSEQ = "".join(map(lambda x: x["AA"], Merged)) # ALL RESIDUES
        self.PDBSEQ = "".join(map(lambda x: x["AA"] if len(x) == 6 else "-", Merged)) # ONLY VISIBLE (non-missing) RESIDUES
        self.PDBSEQ_HTM = map(lambda x: x["AA"] if len(x) == 6 else "-", Merged)
        for i in range(len(self.PDBSEQ_HTM)):
            if len(Merged[i]) > 3 and Merged[i]["hetatm"]:
                self.PDBSEQ_HTM[i] = "."
        self.PDBSEQ_HTM = "".join(self.PDBSEQ_HTM) # visible residues, alternative residues replaced with a "."


    # pdb_data class function
    def genModelSeq(self, loop_threshold):
        # MODEL SEQUENCE GENERATION
        # segments of missing residues will be filled if they are shorter than the loop_threshold
        # chain break is inserted if the missing segment is longer than the threshold
        MS = segmentArray(self.PDBSEQ, isNewSegment)
        ModelSeq = ""
        ModelSeq += self.PDBSEQ[ MS[0]["b"] : MS[0]["e"] ]
        if len(MS) > 2:
            for i in range(1, len(MS) - 1):
                if self.PDBSEQ[ MS[i]["b"] ] != "-":
                    ModelSeq += self.PDBSEQ[ MS[i]["b"] : MS[i]["e"] ]
                elif MS[i]["e"] - MS[i]["b"] < loop_threshold:
                # fill with residues below threshold
                    ModelSeq += self.FULLSEQ[ MS[i]["b"] : MS[i]["e"] ]
                else:
                # insert break above threshold
                    ModelSeq += "/" + self.PDBSEQ[ MS[i]["b"] + 1 : MS[i]["e"] ]
        if len(MS) > 1:
            ModelSeq += self.PDBSEQ[ MS[-1]["b"] : MS[-1]["e"] ]
        return ModelSeq


# generate model
def GenModel(pdb, destination = "protein.pdb", ulen = 12, CustomPir = ""):
# Generates model from pdb file
# Alternative aminoacids with heteroatoms are replaced with standard aminoacids
# Missing segments will be modelled up to the length threshold, ulen
# Modeling loops over 11-12 aminoacids is not reliable
# CustomPir: pir format alignment file can contain any target sequence aligned to template pdb
# REQUIRES MODELLER PROGRAM: https://salilab.org/modeller/
# compatible release: Modeller 9.19, released Jul. 25th, 2017, Linux 64-bit
# modeller path should be accessible from environment as MODELLER_PATH
    modelSeq = pdb.genModelSeq(ulen)
    if pdb.CHAIN == " ":
        chain = "A"
    tmpName = pdb.PDBID + "_" + chain
    # model name from custom alignment
    if os.path.exists(CustomPir):
        with open(CustomPir, "r") as f:
            PIRAL = f.readlines()
        modelName = filter(lambda x: x[:4] == ">P1;" and x[4:-1] != tmpName)[0][4:-1]
    # model name if alignment is not available
    else:
        modelName = os.path.basename(destination)
        if modelName[-4:] == ".pdb":
            modelName = modelName[:-4]
    # first and last position of structure aminoacids is required
    pdbpos = [ x["pos"] for x in pdb.ALLAMINO if len(x) == 6 ]
    # Set environment variables for Modeller
    if "MODELLER_PATH" not in os.environ:
        raise Exception("MissingEnvVariable", "MODELLER_PATH")
    modeller_intel = os.environ["MODELLER_PATH"] + "/lib/x86_64-intel8"
    modeller_modlib = os.environ["MODELLER_PATH"] + "/modlib"
    CENV = os.environ.copy()
    if "LD_LIBRARY_PATH" in CENV:
        CENV["LD_LIBRARY_PATH"] = CENV["LD_LIBRARY_PATH"] + ":" + modeller_intel
    else:
        CENV["LD_LIBRARY_PATH"] = modeller_intel
    if "PYTHONPATH" in CENV:
        CENV["PYTHONPATH"] = CENV["PYTHONPATH"] + ":" + modeller_intel + ":" + modeller_modlib
    else:
        CENV["PYTHONPATH"] = modeller_intel + ":" + modeller_modlib
    # Modeller works in current directory, therefore we have to navigate to destination directory to run Modeller
    CWD = os.getcwd()
    os.chdir(os.path.dirname(destination))
    # model-template alignment in pir format
    pir_path = modelName + ".pir"
    # if specified alignment file is provided
    if os.path.exists(CustomPir):
        if os.path.abspath(CustomPir) != os.path.abspath(pir_path):
            shutil.copy(CustomPir, pir_path)
    # if alignment file is not provided
    else:
        PIRAL = (
        ">P1;" + modelName + "\n" +
        "sequence:::::::::\n" +
        wrapSeq(modelSeq) + "*\n" +
        "\n" +
        ">P1;" + tmpName + "\n" +
        "structureX:{0}:{1:d}:{2}:{3:d}:{2}::::\n".format(pdb.PDBID, pdbpos[0], chain, pdbpos[-1]) +
        wrapSeq(pdb.PDBSEQ_HTM) + "*\n"
        )
    with open(pir_path, "w") as f:
        f.write(PIRAL)
    # accessible template pdb file
    template_path = tmpName + ".pdb"
    if pdb.PDBPATH[-3:] == ".gz":
        with open(template_path, "w") as fx:
            fz = gzip.open(pdb.PDBPATH, "r")
            shutil.copyfileobj(fz, fx)
            fz.close()
    else:
        if os.path.abspath(pdb.PDBPATH) != os.path.abspath(template_path):
            shutil.copy(pdb.PDBPATH, template_path)
    # run Modeller
    log.verbose()
    model_env = environ()
    model_env.io.atom_files_directory = ["./"]
    model_env.io.hetatm = True
    amod = automodel(model_env, alnfile = pir_path, knowns = tmpName, sequence = modelName)
    amod.starting_model = 1
    amod.ending_model = 1
    amod.make()
    model = filter(lambda x: x["failure"] is None, amod.outputs)
    if len(model) > 0:
        model = model[0]
    else:
        raise Exception("NoModelOutput", pir_path)
    # Keep resulted pdb file only
    model_pdb = modelName + ".pdb"
    dircontent = os.listdir("./")
    for f in dircontent:
        if f[len(target_name):] in [".rsr", ".sch", ".ini", ".V99990001", ".D00000001"]:
            os.remove(f)
    os.rename(model["name"], model_pdb)
    if os.path.abspath(pdb.PDBPATH) != os.path.abspath(template_path):
        os.remove(template_path)
    if os.path.exists(CustomPir) and os.path.abspath(CustomPir) != os.path.(pir_path):
        os.remove(pir_path)
    # go back to current directory
    os.chdir(CWD)
    os.rename(os.path.dirname(destination) + "/" + model_pdb, destination)
    if not os.path.exists(destination):
        raise Exception("ModelNotGenerated", pdb.PDBID, pdb.CHAIN)



def SelectChain(pdbcont, chain):
    # filters out the selected chain in the pdb file content
    header = filter(lambda x: x[:6].strip() in ["REMARK", "EXPDTA"], c)
    prot = filter(lambda x: x[:6].strip() in ["ATOM", "HETATM", "TER"] and x[21] == chain, c)
    footer = filter(lambda x: x[:6].strip() == "END", c)
    return header + prot + footer



def RotateStruct(pdb, rotmat, allchain = False, destination = ""):
# Rotate pdb file with a given rotation matrix
# rotation matrix is a 3x4 list with rotation coordinates [X, Y, Z][translation, x, y, z]
# Set allchain True to rotate all chains in the pdb file
# by default, pdb.CHAIN will be extracted and rotated only
    if pdb.PDBPATH[-3:] == ".gz":
        zf = gzip.open(pdb.PDBPATH, "r")
        c = f.readlines()
        zf.close()
    else:    
        with open(pdb.PDBPATH) as f:
            c = f.readlines()
    header = filter(lambda x: x[:6].strip() in ["REMARK", "EXPDTA"], c)
    if allchain:
        prot = filter(lambda x: x[:6].strip() in ["ATOM", "HETATM", "TER"], c)
    else:
        prot = filter(lambda x: x[:6].strip() in ["ATOM", "HETATM", "TER"] and x[21] == pdb.CHAIN, c)
    footer = filter(lambda x: x[:6].strip() == "END", c)
    for i in range(len(prot)):
        if prot[i][:6] in ["ATOM", "HETATM"]:
            x, y, z = float(prot[i][30:38]), float(prot[i][38:46]), float(prot[i][46:54])
            tx = rotmat[0][0] + rotmat[0][1] * x + rotmat[0][2] * y + rotmat[0][3] * z
            ty = rotmat[1][0] + rotmat[1][1] * x + rotmat[1][2] * y + rotmat[1][3] * z
            tz = rotmat[2][0] + rotmat[2][1] * x + rotmat[3][2] * y + rotmat[2][3] * z
            prot[i] = prot[i][:30] + "{0: 8.3f}{1: 8.3f}{2: 8.3f}".format(tx, ty, tz) + prot[i][54:]
    # write file if destination is given
    if len(destination) > 0:
        with open(destination, "w") as f:
            for l in header:
                f.write(l)
            for l in prot:
                f.write(prot)
            for l in footer:
                f.write(footer)
        return destination
    # return rotated file content if destination is not given
    else:
        return header + prot + footer
    


def alignStruct(pdb, pdbref, fasta_path = "", destination = "")
# align a structure to a given reference pdb structure
# uses TM-align fortran version: https://zhanglab.ccmb.med.umich.edu/TM-align/
# latest TM-align version implemented for: 2016/05/21
# TM-align executable should be accessible from the environment as TMALIGN_PATH
# fasta sequence alignment file can be provided as alignment restraint
    # accessible target protein
    if pdb.CHAIN == " ":
        trg_chain = "A"
    trg = pdb.PDBID + "_" + trg_chain
    trg_path = trg + "_trg.pdb"
    if pdb.PDBPATH[-3:] == ".gz":
        fz = gzip.open(pdb.PDBPATH, "r")
        c = fz.readlines()
        fz.close()
    else:
        with open(pdb.PDBPATH, "r") as f:
            c = f.readlines()
    c = SelectChain(c, pdb.CHAIN)
    with open(trg_path, "w") as f:
        for l in c:
            write(l)
    # accessible template protein
    if pdbref.CHAIN = " ":
        tmp_chain = "A" 
    tmp = pdbref.PDBID + "_" + tmp_chain
    tmp_path = tmp + "_tmp.pdb"
    if pdbref.PDBPATH[-3:] == ".gz":
        tmp_path = tmp + ".pdb"
        with open(tmp_path, "w") as fx:
            fz = gzip.open(pdbref.PDBPATH, "r")
            c = fz.readlines()
            fz.close()
    else:
        with open(pdbref.PDBPATH) as f:
            c = f.readlines()
    c = SelectChain(c, pdbref.CHAIN)
    with open(tmp_path, "w") as f:
        for l in c:
            write(l)
    # rotation matrix file
    rotmat_path = trg + "__" + tmp + ".mat"
    # run TM-align
    if "TMALIGN_PATH" not in os.environ:
        raise Exception("MissingEnvVariable", "TMALIGN_PATH")
    if len(fasta_path) > 0:
        call(os.environ["TMALIGN_PATH"] + " " + trg_path + " " + tmp_path + " -I " + fasta_path + " -m " + rotmat_path, shell = True)
    else:
        call(os.environ["TMALIGN_PATH"] + " " + trg_path + " " + tmp_path + " -m " + rotmat_path, shell = True)
    # read rotation matrix
    with open(rotmat_path, "r") as f:
        c = f.readlines()
    m = [[], [], []]
    m[0] = [float(c[2][5:20]), float(c[2][20:35]), float(c[2][35:50]), float(c[2][50:])]
    m[1] = [float(c[3][5:20]), float(c[3][20:35]), float(c[3][35:50]), float(c[3][50:])]
    m[2] = [float(c[4][5:20]), float(c[4][20:35]), float(c[4][35:50]), float(c[4][50:])]
    # rotate target pdb:
    if len(destination) > 0:
        result = RotateStruct(pdb, m, destination = destination)
    else:
        result = RotateStruct(pdb, m)
    if os.path.abspath(pdb.PDBPATH) != trg_path:
        os.remove(trg_path)
    if os.path.abspath(pdbref.PDBPATH) != tmp_path:
        os.remove(tmp_path)
    os.remove(rotmat_path)
    return result



# get rmsd (pdbid, chain or pdbpath)
def GetRMSD(pdb, pdbref, fasta_path = "")
# determine RMSD and TM-score between two proteins
# uses TM-align fortran version: https://zhanglab.ccmb.med.umich.edu/TM-align/
# latest TM-align version implemented for: 2016/05/21
# TM-align executable should be accessible from the environment as TMALIGN_PATH
# fasta sequence alignment file can be provided as alignment restraint
    # accessible target protein
    if pdb.CHAIN == " ":
        trg_chain = "A"
    trg = pdb.PDBID + "_" + trg_chain
    trg_path = trg + "_trg.pdb"
    if pdb.PDBPATH[-3:] == ".gz":
        fz = gzip.open(pdb.PDBPATH, "r")
        c = fz.readlines()
        fz.close()
    else:
        with open(pdb.PDBPATH, "r") as f:
            c = f.readlines()
    c = SelectChain(c, pdb.CHAIN)
    with open(trg_path, "w") as f:
        for l in c:
            write(l)
    # accessible template protein
    if pdbref.CHAIN = " ":
        tmp_chain = "A" 
    tmp = pdbref.PDBID + "_" + tmp_chain
    tmp_path = tmp + "_tmp.pdb"
    if pdbref.PDBPATH[-3:] == ".gz":
        tmp_path = tmp + ".pdb"
        with open(tmp_path, "w") as fx:
            fz = gzip.open(pdbref.PDBPATH, "r")
            c = fz.readlines()
            fz.close()
    else:
        with open(pdbref.PDBPATH) as f:
            c = f.readlines()
    c = SelectChain(c, pdbref.CHAIN)
    with open(tmp_path, "w") as f:
        for l in c:
            write(l)
    # run TM-align
    if "TMALIGN_PATH" not in os.environ:
        raise Exception("MissingEnvVariable", "TMALIGN_PATH")
    if len(fasta_path) > 0:
        tmout = Popen(os.environ["TMALIGN_PATH"] + " " + tmp_path + " " + trg_path + " -I " + fasta_path, shell = True, bufsize = 1, stdout = PIPE).stdout.readlines()
    else:
        tmout = Popen(os.environ["TMALIGN_PATH"] + " " + tmp_path + " " + trg_path, shell = True, bufsize = 1, stdout = PIPE).stdout.readlines()
    if len(tmout) < 19 or len(tmout[17].split()) < 5 or len(tmout[19].split()) < 2:
        for l in tmout:
            print l[:-1]
        raise Exception("WrongTMAlignOutput", tmp_path, trg_path)
    return {
        "RMSD" : float(tmout[17].split()[4][:-1]),
        "TMscore" : float(tmout[19].split()[1])
    }



def get_uniprot(pdbid, chain):
    # Fetch UniProt accession numbers from PDB/UniProt Mapping webserver
    # Query url:
    url = "http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&id=" + pdbid + "&chain=" + chain
    out = urllib2.urlopen(url)
    res = out.read()
    res = res.split("\n")[:-1]

    # outfile = pdbid + "_" + chain + ".up"
    # call("wget " + "\"" + url + "\" -q -O " + outfile, shell = True)
    # if not os.path.exists(outfile):
    #     raise Exception("OutPutNotGenerated", outfile)
    # with open(outfile, "r") as f:
    #     res = f.readlines()
    # os.remove(outfile)

    if len(res) > 0 and len(res) % 7 > 0:
        raise Exception("WrongOutput", url)
    AC = []
    for i in range(len(res) / 7):
        r = res[7 * i + 2].split()
        if r[0] != "AC:" or len(r) != 2:
            print r
            raise Exception("MissingUniprotResult")
        else:
            AC.append(r[1])

    return AC

