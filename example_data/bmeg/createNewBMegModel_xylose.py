#!/usr/bin/env python
## 
## @file    createExampleSBML.py
## @brief   Creates example SBML models presented in the SBML specification.
## @author  Akiya Jouraku
## @author  Michael Hucka
## @author  Sarah Keating
## 
## <!--------------------------------------------------------------------------
## This sample program is distributed under a different license than the rest
## of libSBML.  This program uses the open-source MIT license, as follows:
##
## Copyright (c) 2013-2014 by the California Institute of Technology
## (California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
## and the University of Heidelberg (Germany), with support from the National
## Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Neither the name of the California Institute of Technology (Caltech), nor
## of the European Bioinformatics Institute (EMBL-EBI), nor of the University
## of Heidelberg, nor the names of any contributors, may be used to endorse
## or promote products derived from this software without specific prior
## written permission.
## ------------------------------------------------------------------------ -->
## 


import sys
import os.path
from libsbml import *
import getopt

#
# These variables are used in writeExampleSBML when writing an SBML
# document.  They are handed to libSBML functions in order to include
# the program information into comments within the SBML file.
#
ProgramName = "createExampleModels";
ProgramVersion = "1.0.0";

#
# The SBML Level and Version of the example SBML models.
#
Level = 2;
Version = 4;




#===============================================================================
#
#
# Functions for creating the Example SBML documents.
#
#
#===============================================================================



# 
# 
# Creates an SBML model represented in "7.2 Example involving units"
# in the SBML Level 2 Version 4 Specification.
# 
# 


#===============================================================================
#
#
# Helper functions for writing/validating the given SBML documents.
# 
#
#===============================================================================

#
# 
# Validates the given SBMLDocument.
#
#  This function is based on validateSBML.py implemented by
#  Sarah Keating, Ben Bornstein, and Michael Hucka.
#
#

def validateExampleSBML(sbmlDoc):
    if (sbmlDoc == None):
        print("validateExampleSBML: given a None SBML Document");
        return False;

    consistencyMessages = "";
    validationMessages = "";
    noProblems = True;
    numCheckFailures = 0;
    numConsistencyErrors = 0;
    numConsistencyWarnings = 0;
    numValidationErrors = 0;
    numValidationWarnings = 0;

    # LibSBML 3.3 is lenient when generating models from scratch using the
    # API for creating objects.  Once the whole model is done and before it
    # gets written out, it's important to check that the whole model is in
    # fact complete, consistent and valid.

    numCheckFailures = sbmlDoc.checkInternalConsistency();
    if (numCheckFailures > 0):
        noProblems = False;
        for i in range(0, numCheckFailures):
            sbmlErr = sbmlDoc.getError(i);
            if (sbmlErr.isFatal() or sbmlErr.isError()):
                numConsistencyErrors = 1 + numConsistencyErrors;
            else:
                numConsistencyWarnings = 1+numConsistencyWarnings;
        sbmlDoc.printErrors();

    # If the internal checks fail, it makes little sense to attempt
    # further validation, because the model may be too compromised to
    # be properly interpreted.

    if (numConsistencyErrors > 0):
        consistencyMessages = consistencyMessages  + "Further validation aborted.";
    else:
        numCheckFailures = sbmlDoc.checkConsistency();
        if (numCheckFailures > 0):
            noProblems = False;
            for i in range(0, numCheckFailures):
                sbmlErr = sbmlDoc.getError(i);
                if (sbmlErr.isFatal() or sbmlErr.isError()):
                    numValidationErrors = 1+numValidationErrors;
                else:
                    numValidationWarnings = 1+numValidationWarnings;
            sbmlDoc.printErrors();

    if (noProblems):
        return True;
    else:
        if (numConsistencyErrors > 0):
			tmp = ""
			if (numConsistencyErrors > 1):
				tmp = "s";
			print("ERROR: encountered " + str(numConsistencyErrors) + " consistency error" + tmp + " in model '" + sbmlDoc.getModel().getId() + "'.");
		
        if (numConsistencyWarnings > 0):
			tmp = ""
			if (numConsistencyWarnings > 1):
				tmp = "s"
			print("Notice: encountered " + str(numConsistencyWarnings)
            + " consistency warning" + tmp
            + " in model '" + sbmlDoc.getModel().getId() + "'.");
        print(consistencyMessages);

        if (numValidationErrors > 0):
			tmp = ""
			if (numValidationErrors > 1):
				tmp = "s"
			print("ERROR: encountered " + str(numValidationErrors) + " validation error" + (tmp)
            + " in model '" + sbmlDoc.getModel().getId() + "'.");
        if (numValidationWarnings > 0):
			tmp = ""
			if (numValidationWarnings > 1):
				tmp = "s"
			print("Notice: encountered " + str(numValidationWarnings)
            + " validation warning" + (tmp)
            + " in model '" + sbmlDoc.getModel().getId() + "'.");
        print(validationMessages);

        return (numConsistencyErrors == 0 and numValidationErrors == 0);


# 
# 
# Writes the given SBMLDocument to the given file.
# 
# 
def writeExampleSBML(sbmlDoc, filename):
    result = writeSBML(sbmlDoc, filename);
    if (result == 1):
        print("Wrote file '" + filename + "'");
        return True;
    else:
        print("Failed to write '" + filename + "'");
        return False;


#===============================================================================
#
# Main routine
#
#  Creates SBML models represented in "Example models expressed in XML using
#  SBML" in Section 7 of the SBML Level 2 Version 4 specification(*). 
#
#   (*) The specification document is available at the following URL:
#       http://sbml.org/Documents/Specifications
#
#===============================================================================
#


def createFunctions(model):

    fdef = model.createFunctionDefinition();
    fdef.setId("f");
    mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
      <lambda>
        <bvar><ci> x </ci></bvar>
        <apply>
          <times/>
          <ci> x </ci>
          <cn> 2 </cn>
        </apply>
      </lambda>
    </math>
	"""
    astMath = readMathMLFromString(mathXMLString);
    fdef.setMath(astMath);

    fdef = model.createFunctionDefinition();
    fdef.setId("hill_func");
    mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
      <lambda>
        <bvar><ci> x </ci></bvar>
        <bvar><ci> K </ci></bvar>
        <bvar><ci> n </ci></bvar>
        <apply>
            <divide/>
            <cn> 1 </cn>
            <apply>
                <plus/>
                <cn>1</cn>
                <apply>
                    <power/>
                    <apply>
                        <divide/>
                        <ci> x </ci>
                        <ci> K </ci>
                    </apply>
                    <ci>n</ci>
                </apply>
            </apply>
        </apply>
      </lambda>
    </math>
	"""
    astMath = readMathMLFromString(mathXMLString);
    fdef.setMath(astMath);

    fdef = model.createFunctionDefinition();
    fdef.setId("a_prime");
    mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
      <lambda>
        <bvar><ci> a </ci></bvar>
        <bvar><ci> Km </ci></bvar>
        <apply>
            <divide/>
            <ci> a  </ci>
            <ci> Km </ci>
        </apply>
      </lambda>
    </math>
	"""
    astMath = readMathMLFromString(mathXMLString);
    fdef.setMath(astMath);

    fdef = model.createFunctionDefinition();
    fdef.setId("one_plus_a_prime");
    mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
      <lambda>
        <bvar><ci> a </ci></bvar>
        <bvar><ci> Km </ci></bvar>
        <apply>
            <plus/>
            <cn> 1 </cn>
            <apply>
                <ci> a_prime </ci>
                <ci> a  </ci>
                <ci> Km </ci>
            </apply>
        </apply>
      </lambda>
    </math>
	"""
    astMath = readMathMLFromString(mathXMLString);
    fdef.setMath(astMath);

    return


def createSpecies(model, compName):
    #---------------------------------------------------------------------------
    #
    # Creates Species objects inside the Model object.
    #
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # (Species1) Creates a Species object ("S1")
    #---------------------------------------------------------------------------



    mat_gfp_name = "Mat_GFP"
    print "creating species: " + mat_gfp_name
    mat_gfp = model.createSpecies();
    mat_gfp.setId(mat_gfp_name);
    # Sets the "compartment" attribute of the Species object to identify the
    # compartnet in which the Species object located.
    mat_gfp.setCompartment(compName);
    # Sets the "initialConcentration" attribute of the Species object.
    #
    #  The units of this Species object is determined by two attributes of this
    #  Species object ("substanceUnits" and "hasOnlySubstanceUnits") and the
    #  "spatialDimension" attribute of the Compartment object ("cytosol") in which
    #  this species object located.
    #  Since the default values are used for "substanceUnits" (substance (mole))
    #  and "hasOnlySubstanceUnits" (False) and the value of "spatialDimension" (3)
    #  is greater than 0, the units of this Species object is  mole/litre .
    #
    mat_gfp.setInitialConcentration(0);


    #mat_rfp_name = "Mat_RFP"
    #print "creating species: " + mat_rfp_name
    #mat_rfp = model.createSpecies();
    #mat_rfp.setId(mat_rfp_name);
    #mat_rfp.setCompartment(compName);
    #mat_rfp.setInitialConcentration(0);


    gfp_name = "GFP"
    print "creating species: " + gfp_name
    gfp = model.createSpecies();
    gfp.setId(gfp_name);
    gfp.setCompartment(compName);
    gfp.setInitialConcentration(0);

    #rfp_name = "RFP"
    #print "creating species: " + rfp_name
    #rfp = model.createSpecies();
    #rfp.setId(rfp_name);
    #rfp.setCompartment(compName);
    #rfp.setInitialConcentration(0);

    #rfp_name = "RFP_i1"
    #print "creating species: " + rfp_name
    #rfp = model.createSpecies();
    #rfp.setId(rfp_name);
    #rfp.setCompartment(compName);
    #rfp.setInitialConcentration(0);

    #rfp_name = "RFP_i2"
    #print "creating species: " + rfp_name
    #rfp = model.createSpecies();
    #rfp.setId(rfp_name);
    #rfp.setCompartment(compName);
    #rfp.setInitialConcentration(0);

    # MRNA_GFP
    mrna_gfp_name = "mRNA_GFP"
    print "creating species: " + mrna_gfp_name
    mrna_gfp = model.createSpecies();
    mrna_gfp.setId(mrna_gfp_name);
    mrna_gfp.setCompartment(compName);
    mrna_gfp.setInitialConcentration(0);

    # MRNA_RFP
    #mrna_rfp_name = "mRNA_RFP"
    #print "creating species: " + mrna_rfp_name
    #mrna_rfp = model.createSpecies();
    #mrna_rfp.setId(mrna_rfp_name);
    #mrna_rfp.setCompartment(compName);
    #mrna_rfp.setInitialConcentration(0);

    # prom_GFP (i.e. DNA conc)
    prom_gfp_name = "prom_GFP"
    print "creating species: " + prom_gfp_name
    prom_gfp = model.createSpecies();
    prom_gfp.setId(prom_gfp_name);
    prom_gfp.setCompartment(compName);
    prom_gfp.setInitialConcentration(1);

    # prom_RFP (i.e. DNA conc)
    #prom_rfp_name = "prom_RFP"
    #print "creating species: " + prom_rfp_name
    #prom_rfp = model.createSpecies();
    #prom_rfp.setId(prom_rfp_name);
    #prom_rfp.setCompartment(compName);
    #prom_rfp.setInitialConcentration(1);

    # NTP
    ntp_name = "NTP"
    print "creating species: " + ntp_name
    ntp = model.createSpecies();
    ntp.setId(ntp_name);
    ntp.setCompartment(compName);
    ntp.setInitialConcentration(5000);

    # NMP
    nmp_name = "NMP"
    print "creating species: " + nmp_name
    nmp = model.createSpecies();
    nmp.setId(nmp_name);
    nmp.setCompartment(compName);
    nmp.setInitialConcentration(0);

    # AA
    aa_name = "AA"
    print "creating species: " + aa_name
    aa = model.createSpecies();
    aa.setId(aa_name);
    aa.setCompartment(compName);
    aa.setInitialConcentration(30000);

    # Ribo
    ribo_name = "Ribo"
    print "creating species: " + ribo_name
    ribo = model.createSpecies();
    ribo.setId(ribo_name);
    ribo.setCompartment(compName);
    ribo.setInitialConcentration(10);

    # RNAP
    rnap_name = "RNAP"
    print "creating species: " + rnap_name
    rnap = model.createSpecies();
    rnap.setId(rnap_name);
    rnap.setCompartment(compName);
    rnap.setInitialConcentration(1);

    #---------------------------------------------------------------------------
    #
    # Hybrid species
    #
    #---------------------------------------------------------------------------

    # RNAP-prom_GFP
    rnap_prom_gfp_name = "RNAP_prom_GFP"
    print "creating species: " + rnap_prom_gfp_name
    rnap_prom_gfp = model.createSpecies();
    rnap_prom_gfp.setId(rnap_prom_gfp_name);
    rnap_prom_gfp.setCompartment(compName);
    rnap_prom_gfp.setInitialConcentration(0);

    #
    #rnap_prom_rfp_name = "RNAP_prom_RFP"
    #print "creating species: " + rnap_prom_rfp_name
    #rnap_prom_rfp = model.createSpecies();
    #rnap_prom_rfp.setId(rnap_prom_rfp_name);
    #rnap_prom_rfp.setCompartment(compName);
    #rnap_prom_rfp.setInitialConcentration(0);


    # RNAP-prom_GFP
    rnap_prom_gfp_name = "RNAP_prom_GFP_elong"
    print "creating species: " + rnap_prom_gfp_name
    rnap_prom_gfp = model.createSpecies();
    rnap_prom_gfp.setId(rnap_prom_gfp_name);
    rnap_prom_gfp.setCompartment(compName);
    rnap_prom_gfp.setInitialConcentration(0);

    #
    #rnap_prom_rfp_name = "RNAP_prom_RFP_elong"
    #print "creating species: " + rnap_prom_rfp_name
    #rnap_prom_rfp = model.createSpecies();
    #rnap_prom_rfp.setId(rnap_prom_rfp_name);
    #rnap_prom_rfp.setCompartment(compName);
    #rnap_prom_rfp.setInitialConcentration(0);

    #
    ribo_mrna_gfp_name = "Ribo_mRNA_GFP"
    print "creating species: " + ribo_mrna_gfp_name
    ribo_mrna_gfp = model.createSpecies();
    ribo_mrna_gfp.setId(ribo_mrna_gfp_name);
    ribo_mrna_gfp.setCompartment(compName);
    ribo_mrna_gfp.setInitialConcentration(0);

    #
    #ribo_mrna_rfp_name = "Ribo_mRNA_RFP"
    #print "creating species: " + ribo_mrna_rfp_name
    #ribo_mrna_rfp= model.createSpecies();
    #ribo_mrna_rfp.setId(ribo_mrna_rfp_name);
    #ribo_mrna_rfp.setCompartment(compName);
    #ribo_mrna_rfp.setInitialConcentration(0);

    #
    ribo_mrna_gfp_name = "Ribo_mRNA_GFP_elong"
    print "creating species: " + ribo_mrna_gfp_name
    ribo_mrna_gfp = model.createSpecies();
    ribo_mrna_gfp.setId(ribo_mrna_gfp_name);
    ribo_mrna_gfp.setCompartment(compName);
    ribo_mrna_gfp.setInitialConcentration(0);

    #
    #ribo_mrna_rfp_name = "Ribo_mRNA_RFP_elong"
    #print "creating species: " + ribo_mrna_rfp_name
    #ribo_mrna_rfp= model.createSpecies();
    #ribo_mrna_rfp.setId(ribo_mrna_rfp_name);
    #ribo_mrna_rfp.setCompartment(compName);
    #ribo_mrna_rfp.setInitialConcentration(0);




    # secondary energy source
    second_nrg_name = "IInd_nrg_src"
    print "creating species: " + second_nrg_name
    second_nrg= model.createSpecies();
    second_nrg.setId(second_nrg_name);
    second_nrg.setCompartment(compName);
    second_nrg.setInitialConcentration(5000);

    return


def createParams(model):
    #---------------------------------------------------------------------------
    #
    # Creates a global Parameter object inside the Model object.
    #
    #---------------------------------------------------------------------------

    # assume different translation rates
    para = model.createParameter();
    para.setId("tl_rate_gfp");
    para.setValue(10000);

    # assume different translation rates
    #para = model.createParameter();
    #para.setId("tl_rate_rfp");
    #para.setValue(10000);

    # assume different translation rates
    para = model.createParameter();
    para.setId("tl_escape_rate_gfp");
    para.setValue(100);

    # assume different translation rates
    #para = model.createParameter();
    #para.setId("tl_escape_rate_rfp");
    #para.setValue(100);

    # assume same transcription rate for both promoters
    para = model.createParameter();
    para.setId("tx_rate");
    para.setValue(1000);

    # assume same transcription escape rate for both promoters
    para = model.createParameter();
    para.setId("tx_escape_rate_gfp");
    para.setValue(100);

    # assume same transcription escape rate for both promoters
    #para = model.createParameter();
    #para.setId("tx_escape_rate_rfp");
    #para.setValue(100);

    para = model.createParameter();
    para.setId("rna_deg_rate_gfp");
    para.setValue(0.2);

    #para = model.createParameter();
    #para.setId("rna_deg_rate_rfp");
    #para.setValue(0.2);


    # GFP maturation rate
    para = model.createParameter();
    para.setId("gfp_mat_rate");
    para.setValue(0.06); # min-1  = 1E-3 s-1

    # RFP maturation rate
    #para = model.createParameter();
    #para.setId("rfp_mat_rate1");
    #para.setValue(0.01); # min-1  = 1E-3 s-1

    # RFP maturation rate
    #para = model.createParameter();
    #para.setId("rfp_mat_rate2");
    #para.setValue(0.01); # min-1  = 1E-3 s-1

    # RFP maturation rate
    #para = model.createParameter();
    #para.setId("rfp_mat_rate3");
    #para.setValue(0.01); # min-1  = 1E-3 s-1

    # NTP deg rate
    para = model.createParameter();
    para.setId("ntp_deg_rate");
    para.setValue(0.01);

    # NMP deg rate
    para = model.createParameter();
    para.setId("nmp_deg_rate");
    para.setValue(0.01);

    #  Km_NTP_RNAP
    para = model.createParameter();
    para.setId("Km_ntp_rnap");
    para.setValue(0.1);

    #  Km_NTP_Ribo
    para = model.createParameter();
    para.setId("Km_ntp_ribo");
    para.setValue(10);

    #  Km_AA_Ribo
    para = model.createParameter();
    para.setId("Km_aa_ribo");
    para.setValue(1000);

    #  kon_prom RNAP - assume the same for both promoters
    para = model.createParameter();
    para.setId("kon_prom");
    para.setValue(1);

    #  koff_prom RNAP - assume the same for both promoters
    para = model.createParameter();
    para.setId("koff_prom_gfp");
    para.setValue(4);


    #  koff_prom RNAP - assume the same for both promoters
    #para = model.createParameter();
    #para.setId("koff_prom_rfp");
    #para.setValue(4);


    #  kon_mRNA_GFP Ribo
    para = model.createParameter();
    para.setId("kon_mrna_gfp");
    para.setValue(400.0);

    #  kon_mRNA_RFP Ribo
    #para = model.createParameter();
    #para.setId("kon_mrna_rfp");
    #para.setValue(400.0);

    #  koff_mRNA_GFP Ribo
    para = model.createParameter();
    para.setId("koff_mrna_gfp");
    para.setValue(100.0);

    #  koff_mRNA_RFP Ribo
    #para = model.createParameter();
    #para.setId("koff_mrna_rfp");
    #para.setValue(100.0);

    # length of GFP transcript
    para = model.createParameter();
    para.setId("mrna_gfp_len");
    para.setValue(1150);

    # length of RFP transcript
    #para = model.createParameter();
    #para.setId("mrna_rfp_len");
    #para.setValue(1102);

    # length of GFP protein
    para = model.createParameter();
    para.setId("prot_gfp_len");
    para.setValue(245);

    # length of RFP protein
    #para = model.createParameter();
    #para.setId("prot_rfp_len");
    #para.setValue(236);

    # secondary energy regen rate
    para = model.createParameter();
    para.setId("regen_rate");
    para.setValue(60);

    # secondary energy regen rate
    para = model.createParameter();
    para.setId("Km_IInd_nrg_src");
    para.setValue(800);

    #Km_nmp_regen
    para = model.createParameter();
    para.setId("Km_nmp_regen");
    para.setValue(15);


    #
    para = model.createParameter();
    para.setId("xylr_conc");
    para.setValue(0);

    #
    para = model.createParameter();
    para.setId("xylose_conc");
    para.setValue(0);

    # repressor binding to DNA
    para = model.createParameter();
    para.setId("Kd");
    para.setValue(1);

    # xylose binding to xylr
    para = model.createParameter();
    para.setId("Kx");
    para.setValue(1);

    # stoichiometry of XylR binding
    para = model.createParameter();
    para.setId("n");
    para.setValue(1);

    # stoichiometry of xylose binding
    para = model.createParameter();
    para.setId("m");
    para.setValue(1);

    # ribosome degradation rate
    para = model.createParameter();
    para.setId("ribo_deg_rate");
    para.setValue(0.01);

    return


def createGFPMaturationRxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("gfp_mat_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("GFP");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("Mat_GFP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("gfp_mat_rate");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);

    kl.setMath(astTimes1)

    return


# OK
def createRFPMaturationRxn1(model, compName):

    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("rfp_mat_rxn1")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("RFP");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("RFP_i1");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("rfp_mat_rate1");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);

    kl.setMath(astTimes1)

    return

# OK
def createRFPMaturationRxn2(model, compName):

    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("rfp_mat_rxn2")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("RFP_i1");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("RFP_i2");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("rfp_mat_rate2");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("RFP_i1");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);

    kl.setMath(astTimes1)

    return


# OK
def createRFPMaturationRxn3(model, compName):

    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("rfp_mat_rxn3")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("RFP_i2");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("Mat_RFP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("rfp_mat_rate3");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("RFP_i2");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);

    kl.setMath(astTimes1)

    return


def createNTPdegRxn(model, compName):

    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("ntp_deg_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("NTP");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("NMP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("ntp_deg_rate");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("NTP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);

    kl.setMath(astTimes1)

    return


def createNMPdegRxn(model, compName):

    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("nmp_deg_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("NMP");


    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("nmp_deg_rate");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("NMP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);

    kl.setMath(astTimes1)

    return


def createRNAGFPDegRxn(model, compName):


    # create mRNA degradation
    rna_deg_rxn = model.createReaction();
    rna_deg_rxn.setId("mrna_gfp_deg_rxn")
    rna_deg_rxn.setReversible(False);

    react = rna_deg_rxn.createReactant();
    react.setSpecies("mRNA_GFP");
    react.setStoichiometry(1)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("NMP");
    stoichmath = ASTNode(AST_NAME);
    stoichmath.setName("mrna_gfp_len");
    sm = spr.createStoichiometryMath()
    sm.setMath(stoichmath)

    kl = rna_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("rna_deg_rate_gfp");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("mRNA_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createRNARFPDegRxn(model, compName):


    # create mRNA degradation
    rna_deg_rxn = model.createReaction();
    rna_deg_rxn.setId("mrna_rfp_deg_rxn")
    rna_deg_rxn.setReversible(False);

    react = rna_deg_rxn.createReactant();
    react.setSpecies("mRNA_RFP");
    react.setStoichiometry(1)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("NMP");
    stoichmath = ASTNode(AST_NAME);
    stoichmath.setName("mrna_rfp_len");
    sm = spr.createStoichiometryMath()
    sm.setMath(stoichmath)

    kl = rna_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("rna_deg_rate_rfp");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("mRNA_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createRibo_mRNA_GFPDegRxn(model, compName):

    # create mRNA degradation
    rna_deg_rxn = model.createReaction();
    rna_deg_rxn.setId("ribo_mrna_gfp_deg_rxn")
    rna_deg_rxn.setReversible(False);

    react = rna_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_GFP");
    react.setStoichiometry(1)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("NMP");
    stoichmath = ASTNode(AST_NAME);
    stoichmath.setName("mrna_gfp_len");
    sm = spr.createStoichiometryMath()
    sm.setMath(stoichmath)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("Ribo");

    kl = rna_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("rna_deg_rate_gfp");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)

def createRibo_mRNA_RFPDegRxn(model, compName):

    # create mRNA degradation
    rna_deg_rxn = model.createReaction();
    rna_deg_rxn.setId("ribo_mrna_rfp_deg_rxn")
    rna_deg_rxn.setReversible(False);

    react = rna_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_RFP");
    react.setStoichiometry(1)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("NMP");
    stoichmath = ASTNode(AST_NAME);
    stoichmath.setName("mrna_rfp_len");
    sm = spr.createStoichiometryMath()
    sm.setMath(stoichmath)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("Ribo");

    kl = rna_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("rna_deg_rate_rfp");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createRibo_mRNA_GFP_elongDegRxn(model, compName):

    # create mRNA degradation
    rna_deg_rxn = model.createReaction();
    rna_deg_rxn.setId("ribo_mrna_gfp_elong_deg_rxn")
    rna_deg_rxn.setReversible(False);

    react = rna_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_GFP_elong");
    react.setStoichiometry(1)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("Ribo");

    kl = rna_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("rna_deg_rate_gfp");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_GFP_elong");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)

def createRibo_mRNA_RFP_elongDegRxn(model, compName):

    # create mRNA degradation
    rna_deg_rxn = model.createReaction();
    rna_deg_rxn.setId("ribo_mrna_rfp_elong_deg_rxn")
    rna_deg_rxn.setReversible(False);

    react = rna_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_RFP_elong");
    react.setStoichiometry(1)

    spr = rna_deg_rxn.createProduct();
    spr.setSpecies("Ribo");

    kl = rna_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("rna_deg_rate_rfp");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_RFP_elong");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createProm_GFP_RNAP_on_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("prom_gfp_rnap_on_rxn")
    gfp_mat_rxn.setReversible(False);

    react1 = gfp_mat_rxn.createReactant();
    react1.setSpecies("prom_GFP");

    react2 = gfp_mat_rxn.createReactant();
    react2.setSpecies("RNAP");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("RNAP_prom_GFP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("kon_prom");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("prom_GFP");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("RNAP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return





def createProm_RFP_RNAP_on_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("prom_rfp_rnap_on_rxn")
    gfp_mat_rxn.setReversible(False);

    react1 = gfp_mat_rxn.createReactant();
    react1.setSpecies("prom_RFP");

    react2 = gfp_mat_rxn.createReactant();
    react2.setSpecies("RNAP");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("RNAP_prom_RFP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("kon_prom");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("prom_RFP");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("RNAP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return


def createProm_GFP_RNAP_off_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("prom_gfp_rnap_off_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("RNAP_prom_GFP");

    prod1 = gfp_mat_rxn.createProduct();
    prod1.setSpecies("prom_GFP");

    prod2 = gfp_mat_rxn.createProduct();
    prod2.setSpecies("RNAP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("koff_prom_gfp");


    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("RNAP_prom_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return

def createProm_RFP_RNAP_off_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("prom_rfp_rnap_off_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("RNAP_prom_RFP");

    prod1 = gfp_mat_rxn.createProduct();
    prod1.setSpecies("prom_RFP");

    prod2 = gfp_mat_rxn.createProduct();
    prod2.setSpecies("RNAP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("koff_prom_rfp");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("RNAP_prom_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return



def createmRNA_GFP_Ribo_on_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("mrna_gfp_ribo_on_rxn")
    gfp_mat_rxn.setReversible(False);

    react1 = gfp_mat_rxn.createReactant();
    react1.setSpecies("mRNA_GFP");

    react2 = gfp_mat_rxn.createReactant();
    react2.setSpecies("Ribo");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("Ribo_mRNA_GFP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("kon_mrna_gfp");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("Ribo");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("mRNA_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return

def createmRNA_RFP_Ribo_on_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("mrna_rfp_ribo_on_rxn")
    gfp_mat_rxn.setReversible(False);

    react1 = gfp_mat_rxn.createReactant();
    react1.setSpecies("mRNA_RFP");

    react2 = gfp_mat_rxn.createReactant();
    react2.setSpecies("Ribo");

    prod = gfp_mat_rxn.createProduct();
    prod.setSpecies("Ribo_mRNA_RFP");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("kon_mrna_rfp");

    astGFP = ASTNode(AST_NAME);
    astGFP.setName("Ribo");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("mRNA_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astGFP);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return


def createmRNA_GFP_Ribo_off_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("mrna_gfp_ribo_off_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_GFP");

    prod1 = gfp_mat_rxn.createProduct();
    prod1.setSpecies("mRNA_GFP");

    prod2 = gfp_mat_rxn.createProduct();
    prod2.setSpecies("Ribo");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("koff_mrna_gfp");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("Ribo_mRNA_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return

def createmRNA_RFP_Ribo_off_Rxn(model, compName):


    gfp_mat_rxn = model.createReaction();
    gfp_mat_rxn.setId("mrna_rfp_ribo_off_rxn")
    gfp_mat_rxn.setReversible(False);

    react = gfp_mat_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_RFP");

    prod1 = gfp_mat_rxn.createProduct();
    prod1.setSpecies("mRNA_RFP");

    prod2 = gfp_mat_rxn.createProduct();
    prod2.setSpecies("Ribo");

    kl = gfp_mat_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astmat_rate = ASTNode(AST_NAME);
    astmat_rate.setName("koff_mrna_rfp");

    astRNAP = ASTNode(AST_NAME);
    astRNAP.setName("Ribo_mRNA_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astmat_rate);
    astTimes1.addChild(astRNAP);

    kl.setMath(astTimes1)

    return

def createTX_GFP_Rxn(model, compName):


    tx_rxn = model.createReaction();
    tx_rxn.setId("tx_gfp_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("RNAP");

    pr2 = tx_rxn.createProduct();
    pr2.setSpecies("mRNA_GFP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("NTP");
    stoichmath = ASTNode(AST_NAME);
    stoichmath.setName("mrna_gfp_len");
    sm = react1.createStoichiometryMath()
    sm.setMath(stoichmath)

    react2 = tx_rxn.createReactant()
    react2.setSpecies("RNAP_prom_GFP_elong")

    kl = tx_rxn.createKineticLaw();

    # Need to redeclare as otherwise segmentation faults - seems nodes can not be reused
    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("tx_rate");

    astKm_mrna_cap = ASTNode(AST_NAME);
    astKm_mrna_cap.setName("Km_ntp_rnap")

    astMrna_cap = ASTNode(AST_NAME)
    astMrna_cap.setName("NTP")

    astDNA = ASTNode(AST_NAME);
    astDNA.setName("RNAP_prom_GFP_elong");



    astMRna_cap_prime =  ASTNode(AST_FUNCTION);
    astMRna_cap_prime.setName("a_prime")
    astMRna_cap_prime.addChild(astMrna_cap)
    astMRna_cap_prime.addChild(astKm_mrna_cap)


    # create again as can not reuse ast objects
    astKm_mrna_cap = ASTNode(AST_NAME);
    astKm_mrna_cap.setName("Km_ntp_rnap")
    astMrna_cap = ASTNode(AST_NAME)
    astMrna_cap.setName("NTP")

    astMRna_cap_op_prime =  ASTNode(AST_FUNCTION);
    astMRna_cap_op_prime.setName("one_plus_a_prime")
    astMRna_cap_op_prime.addChild(astMrna_cap)
    astMRna_cap_op_prime.addChild(astKm_mrna_cap)

    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("mrna_gfp_len");

    astTimes_denom = ASTNode(AST_TIMES);
    astTimes_denom.addChild(astGFPlen)
    astTimes_denom.addChild(astMRna_cap_op_prime)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astMRna_cap_prime)
    astDivide.addChild(astTimes_denom)

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDNA)
    astTimes1.addChild(astDivide);

    kl.setMath(astTimes1)

    return


def createTX_RFP_Rxn(model, compName):


    tx_rxn = model.createReaction();
    tx_rxn.setId("tx_rfp_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("RNAP");

    pr2 = tx_rxn.createProduct();
    pr2.setSpecies("mRNA_RFP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("NTP");
    stoichmath = ASTNode(AST_NAME);
    stoichmath.setName("mrna_rfp_len");
    sm = react1.createStoichiometryMath()
    sm.setMath(stoichmath)

    react2 = tx_rxn.createReactant()
    react2.setSpecies("RNAP_prom_RFP_elong")

    kl = tx_rxn.createKineticLaw();

    # Need to redeclare as otherwise segmentation faults - seems nodes can not be reused
    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("tx_rate");

    astKm_mrna_cap = ASTNode(AST_NAME);
    astKm_mrna_cap.setName("Km_ntp_rnap")

    astMrna_cap = ASTNode(AST_NAME)
    astMrna_cap.setName("NTP")

    astDNA = ASTNode(AST_NAME);
    astDNA.setName("RNAP_prom_RFP_elong");



    astMRna_cap_prime =  ASTNode(AST_FUNCTION);
    astMRna_cap_prime.setName("a_prime")
    astMRna_cap_prime.addChild(astMrna_cap)
    astMRna_cap_prime.addChild(astKm_mrna_cap)


    # create again as can not reuse ast objects
    astKm_mrna_cap = ASTNode(AST_NAME);
    astKm_mrna_cap.setName("Km_ntp_rnap")
    astMrna_cap = ASTNode(AST_NAME)
    astMrna_cap.setName("NTP")

    astMRna_cap_op_prime =  ASTNode(AST_FUNCTION);
    astMRna_cap_op_prime.setName("one_plus_a_prime")
    astMRna_cap_op_prime.addChild(astMrna_cap)
    astMRna_cap_op_prime.addChild(astKm_mrna_cap)


    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("mrna_rfp_len");

    astTimes_denom = ASTNode(AST_TIMES);
    astTimes_denom.addChild(astGFPlen)
    astTimes_denom.addChild(astMRna_cap_op_prime)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astMRna_cap_prime)
    astDivide.addChild(astTimes_denom)

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDNA)
    astTimes1.addChild(astDivide);

    kl.setMath(astTimes1)

    return


def createTL_GFPRxn(model, compName):
    # Create translation reaction

    tl_rxn = model.createReaction();
    tl_rxn.setId("tl_gfp_rxn");
    tl_rxn.setReversible(False);


    react1 = tl_rxn.createReactant();
    react1.setSpecies("NTP");
    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_gfp_len");
    astReal = ASTNode(AST_REAL);
    astReal.setName("const_mul");
    astReal.setValue(2)
    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astReal);
    astTimes1.addChild(astGFPlen);
    sm = react1.createStoichiometryMath()
    sm.setMath(astTimes1)

    react2 = tl_rxn.createReactant();
    react2.setSpecies("AA");
    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_gfp_len");
    sm = react2.createStoichiometryMath()
    sm.setMath(astGFPlen)

    react3 = tl_rxn.createReactant()
    react3.setSpecies("Ribo_mRNA_GFP_elong")

    pr1 = tl_rxn.createProduct();
    pr1.setSpecies("GFP");

    pr2 = tl_rxn.createProduct();
    pr2.setSpecies("NMP");
    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_gfp_len");
    astReal = ASTNode(AST_REAL);
    astReal.setName("const_mul");
    astReal.setValue(2)
    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astReal);
    astTimes1.addChild(astGFPlen);
    sm = pr2.createStoichiometryMath()
    sm.setMath(astTimes1)

    pr3 = tl_rxn.createProduct();
    pr3.setSpecies("Ribo");


    kl = tl_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTL_rate = ASTNode(AST_NAME);
    astTL_rate.setName("tl_rate_gfp");

    astRNA = ASTNode(AST_NAME);
    astRNA.setName("Ribo_mRNA_GFP_elong");


    astAA = ASTNode(AST_NAME);
    astAA.setName("AA");
    astKmAA = ASTNode(AST_NAME);
    astKmAA.setName("Km_aa_ribo");
    astAA_prime =  ASTNode(AST_FUNCTION);
    astAA_prime.setName("a_prime")
    astAA_prime.addChild(astAA)
    astAA_prime.addChild(astKmAA)

    astAA = ASTNode(AST_NAME);
    astAA.setName("AA");
    astKmAA = ASTNode(AST_NAME);
    astKmAA.setName("Km_aa_ribo");
    astAA_op_prime =  ASTNode(AST_FUNCTION);
    astAA_op_prime.setName("one_plus_a_prime")
    astAA_op_prime.addChild(astAA)
    astAA_op_prime.addChild(astKmAA)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_prime =  ASTNode(AST_FUNCTION);
    astNTP_prime.setName("a_prime")
    astNTP_prime.addChild(astNTP)
    astNTP_prime.addChild(astKmNTP)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_op_prime =  ASTNode(AST_FUNCTION);
    astNTP_op_prime.setName("one_plus_a_prime")
    astNTP_op_prime.addChild(astNTP)
    astNTP_op_prime.addChild(astKmNTP)


    astTimes_num = ASTNode(AST_TIMES);
    astTimes_num.addChild(astAA_prime)
    astTimes_num.addChild(astNTP_prime)

    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_gfp_len");

    astTimes_denom = ASTNode(AST_TIMES);
    astTimes_denom.addChild(astGFPlen)
    astTimes_denom.addChild(astAA_op_prime)
    astTimes_denom.addChild(astNTP_op_prime)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astTimes_num)
    astDivide.addChild(astTimes_denom)

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTL_rate);
    astTimes1.addChild(astDivide)
    astTimes1.addChild(astRNA);
    #astTimes1.addChild(astProtCap)

    kl.setMath(astTimes1)
    return


def createTL_RFPRxn(model, compName):
    # Create translation reaction

    tl_rxn = model.createReaction();
    tl_rxn.setId("tl_rfp_rxn");
    tl_rxn.setReversible(False);


    react1 = tl_rxn.createReactant();
    react1.setSpecies("NTP");
    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_rfp_len");
    astReal = ASTNode(AST_REAL);
    astReal.setName("const_mul");
    astReal.setValue(2)
    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astReal);
    astTimes1.addChild(astGFPlen);
    sm = react1.createStoichiometryMath()
    sm.setMath(astTimes1)

    react2 = tl_rxn.createReactant();
    react2.setSpecies("AA");
    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_rfp_len");
    sm = react2.createStoichiometryMath()
    sm.setMath(astGFPlen)

    react3 = tl_rxn.createReactant()
    react3.setSpecies("Ribo_mRNA_RFP_elong")

    pr1 = tl_rxn.createProduct();
    pr1.setSpecies("RFP");

    pr2 = tl_rxn.createProduct();
    pr2.setSpecies("NMP");
    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_rfp_len");
    astReal = ASTNode(AST_REAL);
    astReal.setName("const_mul");
    astReal.setValue(2)
    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astReal);
    astTimes1.addChild(astGFPlen);
    sm = pr2.createStoichiometryMath()
    sm.setMath(astTimes1)

    pr3 = tl_rxn.createProduct();
    pr3.setSpecies("Ribo");

    kl = tl_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTL_rate = ASTNode(AST_NAME);
    astTL_rate.setName("tl_rate_rfp");

    astRNA = ASTNode(AST_NAME);
    astRNA.setName("Ribo_mRNA_RFP_elong");


    astAA = ASTNode(AST_NAME);
    astAA.setName("AA");
    astKmAA = ASTNode(AST_NAME);
    astKmAA.setName("Km_aa_ribo");
    astAA_prime =  ASTNode(AST_FUNCTION);
    astAA_prime.setName("a_prime")
    astAA_prime.addChild(astAA)
    astAA_prime.addChild(astKmAA)

    astAA = ASTNode(AST_NAME);
    astAA.setName("AA");
    astKmAA = ASTNode(AST_NAME);
    astKmAA.setName("Km_aa_ribo");
    astAA_op_prime =  ASTNode(AST_FUNCTION);
    astAA_op_prime.setName("one_plus_a_prime")
    astAA_op_prime.addChild(astAA)
    astAA_op_prime.addChild(astKmAA)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_prime =  ASTNode(AST_FUNCTION);
    astNTP_prime.setName("a_prime")
    astNTP_prime.addChild(astNTP)
    astNTP_prime.addChild(astKmNTP)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_op_prime =  ASTNode(AST_FUNCTION);
    astNTP_op_prime.setName("one_plus_a_prime")
    astNTP_op_prime.addChild(astNTP)
    astNTP_op_prime.addChild(astKmNTP)

    astTimes_num = ASTNode(AST_TIMES);
    astTimes_num.addChild(astAA_prime)
    astTimes_num.addChild(astNTP_prime)

    astGFPlen = ASTNode(AST_NAME);
    astGFPlen.setName("prot_rfp_len");

    astTimes_denom = ASTNode(AST_TIMES);
    astTimes_denom.addChild(astGFPlen)
    astTimes_denom.addChild(astAA_op_prime)
    astTimes_denom.addChild(astNTP_op_prime)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astTimes_num)
    astDivide.addChild(astTimes_denom)

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTL_rate);
    astTimes1.addChild(astDivide)
    astTimes1.addChild(astRNA);
    #astTimes1.addChild(astProtCap)

    kl.setMath(astTimes1)
    return


def createTX_GFP_Escape_Rxn(model, compName):

    tx_rxn = model.createReaction();
    tx_rxn.setId("tx_gfp_escape_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("RNAP_prom_GFP_elong");

    pr2 = tx_rxn.createProduct();
    pr2.setSpecies("prom_GFP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("RNAP_prom_GFP");

    mod = tx_rxn.createModifier()
    mod.setSpecies("NTP")

    kl = tx_rxn.createKineticLaw();

    #

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_rnap");
    astNTP_prime = ASTNode(AST_FUNCTION);
    astNTP_prime.setName("a_prime")
    astNTP_prime.addChild(astNTP)
    astNTP_prime.addChild(astKmNTP)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_rnap");
    astNTP_op_prime = ASTNode(AST_FUNCTION);
    astNTP_op_prime.setName("one_plus_a_prime")
    astNTP_op_prime.addChild(astNTP)
    astNTP_op_prime.addChild(astKmNTP)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astNTP_prime)
    astDivide.addChild(astNTP_op_prime)

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("tx_escape_rate_gfp");


    astDNA = ASTNode(AST_NAME);
    astDNA.setName("RNAP_prom_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDNA)
    astTimes1.addChild(astDivide)

    kl.setMath(astTimes1)

    return

def createTX_RFP_Escape_Rxn(model, compName):

    tx_rxn = model.createReaction();
    tx_rxn.setId("tx_rfp_escape_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("RNAP_prom_RFP_elong");

    pr2 = tx_rxn.createProduct();
    pr2.setSpecies("prom_RFP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("RNAP_prom_RFP");

    mod = tx_rxn.createModifier()
    mod.setSpecies("NTP")

    kl = tx_rxn.createKineticLaw();

    #

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_rnap");
    astNTP_prime = ASTNode(AST_FUNCTION);
    astNTP_prime.setName("a_prime")
    astNTP_prime.addChild(astNTP)
    astNTP_prime.addChild(astKmNTP)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_rnap");
    astNTP_op_prime = ASTNode(AST_FUNCTION);
    astNTP_op_prime.setName("one_plus_a_prime")
    astNTP_op_prime.addChild(astNTP)
    astNTP_op_prime.addChild(astKmNTP)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astNTP_prime)
    astDivide.addChild(astNTP_op_prime)

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("tx_escape_rate_rfp");


    astDNA = ASTNode(AST_NAME);
    astDNA.setName("RNAP_prom_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDNA)
    astTimes1.addChild(astDivide)

    kl.setMath(astTimes1)

    return


def createTL_GFP_Escape_Rxn(model, compName):

    tx_rxn = model.createReaction();
    tx_rxn.setId("tl_gfp_escape_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("Ribo_mRNA_GFP_elong");

    pr2 = tx_rxn.createProduct();
    pr2.setSpecies("mRNA_GFP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("Ribo_mRNA_GFP");

    mod = tx_rxn.createModifier()
    mod.setSpecies("NTP")

    kl = tx_rxn.createKineticLaw();

    #

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_prime = ASTNode(AST_FUNCTION);
    astNTP_prime.setName("a_prime")
    astNTP_prime.addChild(astNTP)
    astNTP_prime.addChild(astKmNTP)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_op_prime = ASTNode(AST_FUNCTION);
    astNTP_op_prime.setName("one_plus_a_prime")
    astNTP_op_prime.addChild(astNTP)
    astNTP_op_prime.addChild(astKmNTP)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astNTP_prime)
    astDivide.addChild(astNTP_op_prime)

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("tl_escape_rate_gfp");


    astDNA = ASTNode(AST_NAME);
    astDNA.setName("Ribo_mRNA_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDNA)
    astTimes1.addChild(astDivide)

    kl.setMath(astTimes1)

    return

def createTL_RFP_Escape_Rxn(model, compName):

    tx_rxn = model.createReaction();
    tx_rxn.setId("tl_rfp_escape_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("Ribo_mRNA_RFP_elong");

    pr2 = tx_rxn.createProduct();
    pr2.setSpecies("mRNA_RFP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("Ribo_mRNA_RFP");

    mod = tx_rxn.createModifier()
    mod.setSpecies("NTP")

    kl = tx_rxn.createKineticLaw();

    #
    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_prime = ASTNode(AST_FUNCTION);
    astNTP_prime.setName("a_prime")
    astNTP_prime.addChild(astNTP)
    astNTP_prime.addChild(astKmNTP)

    astNTP = ASTNode(AST_NAME);
    astNTP.setName("NTP");
    astKmNTP = ASTNode(AST_NAME);
    astKmNTP.setName("Km_ntp_ribo");
    astNTP_op_prime = ASTNode(AST_FUNCTION);
    astNTP_op_prime.setName("one_plus_a_prime")
    astNTP_op_prime.addChild(astNTP)
    astNTP_op_prime.addChild(astKmNTP)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astNTP_prime)
    astDivide.addChild(astNTP_op_prime)

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("tl_escape_rate_rfp");


    astDNA = ASTNode(AST_NAME);
    astDNA.setName("Ribo_mRNA_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDNA)
    astTimes1.addChild(astDivide)

    kl.setMath(astTimes1)

    return


def createRegen_Rxn(model, compName):


    tx_rxn = model.createReaction();
    tx_rxn.setId("regen_rxn");
    tx_rxn.setReversible(False);

    pr1 = tx_rxn.createProduct();
    pr1.setSpecies("NTP");

    react1 = tx_rxn.createReactant();
    react1.setSpecies("NMP");

    react2 = tx_rxn.createReactant();
    react2.setSpecies("IInd_nrg_src");


    kl = tx_rxn.createKineticLaw();

    # Need to redeclare as otherwise segmentation faults - seems nodes can not be reused
    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astTX_rate = ASTNode(AST_NAME);
    astTX_rate.setName("regen_rate");



    astKm_nmp = ASTNode(AST_NAME);
    astKm_nmp.setName("Km_nmp_regen")
    astNMP = ASTNode(AST_NAME)
    astNMP.setName("NMP")
    astNMP_prime =  ASTNode(AST_FUNCTION);
    astNMP_prime.setName("a_prime")
    astNMP_prime.addChild(astNMP)
    astNMP_prime.addChild(astKm_nmp)


    astKm_nmp = ASTNode(AST_NAME);
    astKm_nmp.setName("Km_nmp_regen")
    astNMP = ASTNode(AST_NAME)
    astNMP.setName("NMP")
    astNMP_op_prime =  ASTNode(AST_FUNCTION);
    astNMP_op_prime.setName("one_plus_a_prime")
    astNMP_op_prime.addChild(astNMP)
    astNMP_op_prime.addChild(astKm_nmp)

    astKm_sec = ASTNode(AST_NAME);
    astKm_sec.setName("Km_IInd_nrg_src")
    astSec = ASTNode(AST_NAME)
    astSec.setName("IInd_nrg_src")
    astSec_prime =  ASTNode(AST_FUNCTION);
    astSec_prime.setName("a_prime")
    astSec_prime.addChild(astSec)
    astSec_prime.addChild(astKm_sec)

    astKm_sec = ASTNode(AST_NAME);
    astKm_sec.setName("Km_IInd_nrg_src")
    astSec = ASTNode(AST_NAME)
    astSec.setName("IInd_nrg_src")
    astSec_op_prime =  ASTNode(AST_FUNCTION);
    astSec_op_prime.setName("one_plus_a_prime")
    astSec_op_prime.addChild(astSec)
    astSec_op_prime.addChild(astKm_sec)


    astTimes_num = ASTNode(AST_TIMES);
    astTimes_num.addChild(astNMP_prime)
    astTimes_num.addChild(astSec_prime)

    astTimes_denom = ASTNode(AST_TIMES);
    astTimes_denom.addChild(astNMP_op_prime)
    astTimes_denom.addChild(astSec_op_prime)

    astDivide = ASTNode(AST_DIVIDE);
    astDivide.addChild(astTimes_num)
    astDivide.addChild(astTimes_denom)

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astTX_rate);
    astTimes1.addChild(astDivide);

    kl.setMath(astTimes1)

    return


def createRibo_DegRxn(model, compName):

    # create mRNA degradation
    ribo_deg_rxn = model.createReaction();
    ribo_deg_rxn.setId("ribo_deg_rxn")
    ribo_deg_rxn.setReversible(False);

    react = ribo_deg_rxn.createReactant();
    react.setSpecies("Ribo");
    react.setStoichiometry(1)


    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("NMP");
    #stoichmath = ASTNode(AST_NAME);
    #stoichmath.setName("mrna_rfp_len");
    #sm = spr.createStoichiometryMath()
    #sm.setMath(stoichmath)

    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("Ribo");

    kl = ribo_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("ribo_deg_rate");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createRibo_elong_Deg_GFP_Rxn(model, compName):

    # create mRNA degradation
    ribo_deg_rxn = model.createReaction();
    ribo_deg_rxn.setId("ribo_elong_deg_gfp_rxn")
    ribo_deg_rxn.setReversible(False);

    react = ribo_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_GFP_elong");
    react.setStoichiometry(1)


    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("NMP");
    #stoichmath = ASTNode(AST_NAME);
    #stoichmath.setName("mrna_rfp_len");
    #sm = spr.createStoichiometryMath()
    #sm.setMath(stoichmath)

    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("Ribo");

    kl = ribo_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("ribo_deg_rate");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_GFP_elong");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)

def createRibo_elong_Deg_RFP_Rxn(model, compName):

    # create mRNA degradation
    ribo_deg_rxn = model.createReaction();
    ribo_deg_rxn.setId("ribo_elong_deg_rfp_rxn")
    ribo_deg_rxn.setReversible(False);

    react = ribo_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_RFP_elong");
    react.setStoichiometry(1)


    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("NMP");
    #stoichmath = ASTNode(AST_NAME);
    #stoichmath.setName("mrna_rfp_len");
    #sm = spr.createStoichiometryMath()
    #sm.setMath(stoichmath)

    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("Ribo");

    kl = ribo_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("ribo_deg_rate");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_RFP_elong");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createRibo_Deg_GFP_Rxn(model, compName):

    # create mRNA degradation
    ribo_deg_rxn = model.createReaction();
    ribo_deg_rxn.setId("ribo_deg_gfp_rxn")
    ribo_deg_rxn.setReversible(False);

    react = ribo_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_GFP");
    react.setStoichiometry(1)


    spr = ribo_deg_rxn.createProduct();
    spr.setSpecies("mRNA_GFP");
    #stoichmath = ASTNode(AST_NAME);
    #stoichmath.setName("mrna_rfp_len");
    #sm = spr.createStoichiometryMath()
    #sm.setMath(stoichmath)

    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("Ribo");

    kl = ribo_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("ribo_deg_rate");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_GFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)


def createRibo_Deg_RFP_Rxn(model, compName):

    # create mRNA degradation
    ribo_deg_rxn = model.createReaction();
    ribo_deg_rxn.setId("ribo_deg_rfp_rxn")
    ribo_deg_rxn.setReversible(False);

    react = ribo_deg_rxn.createReactant();
    react.setSpecies("Ribo_mRNA_RFP");
    react.setStoichiometry(1)


    spr = ribo_deg_rxn.createProduct();
    spr.setSpecies("mRNA_RFP");
    #stoichmath = ASTNode(AST_NAME);
    #stoichmath.setName("mrna_rfp_len");
    #sm = spr.createStoichiometryMath()
    #sm.setMath(stoichmath)

    #spr = ribo_deg_rxn.createProduct();
    #spr.setSpecies("Ribo");

    kl = ribo_deg_rxn.createKineticLaw();

    astTube = ASTNode(AST_NAME);
    astTube.setName(compName);

    astdeg_rate = ASTNode(AST_NAME);
    astdeg_rate.setName("ribo_deg_rate");

    astMRNA = ASTNode(AST_NAME);
    astMRNA.setName("Ribo_mRNA_RFP");

    astTimes1 = ASTNode(AST_TIMES);
    astTimes1.addChild(astTube);
    astTimes1.addChild(astdeg_rate);
    astTimes1.addChild(astMRNA);

    kl.setMath(astTimes1)

def createModel():
    level = Level;
    version = Version;

    #---------------------------------------------------------------------------
    #
    # Creates an SBMLDocument object
    #
    #---------------------------------------------------------------------------

    sbmlDoc = SBMLDocument(level, version);

    #---------------------------------------------------------------------------
    #
    # Creates a Model object inside the SBMLDocument object.
    #
    #---------------------------------------------------------------------------

    model = sbmlDoc.createModel();
    model.setId("BMegModel");

    #old_createFunctions(model)


    #---------------------------------------------------------------------------
    #
    # Creates a Compartment object inside the Model object.
    #
    #---------------------------------------------------------------------------
    compName = "TestTube";
    # Creates a Compartment object ("compartmentOne")
    comp = model.createCompartment();
    comp.setId(compName);
    # Sets the "size" attribute of the Compartment object.
    #
    #   The units of this Compartment object is the default SBML
    #   units of volume (litre), and thus we don't have to explicitly invoke
    #   setUnits("litre") function to set the default units.
    #
    comp.setSize(1);

    createFunctions(model);
    createSpecies(model, compName);
    createParams(model);

    createGFPMaturationRxn(model, compName);
    #createRFPMaturationRxn1(model, compName);
    #createRFPMaturationRxn2(model, compName);
    #createRFPMaturationRxn3(model, compName);
    createRNAGFPDegRxn(model, compName);
    #createRNARFPDegRxn(model, compName);
    createNTPdegRxn(model, compName)
    createNMPdegRxn(model, compName)
    createProm_GFP_RNAP_on_Rxn(model, compName)
    #createProm_RFP_RNAP_on_Rxn(model, compName)
    createProm_GFP_RNAP_off_Rxn(model, compName)
    #createProm_RFP_RNAP_off_Rxn(model, compName)
    createmRNA_GFP_Ribo_on_Rxn(model, compName)
    #createmRNA_RFP_Ribo_on_Rxn(model, compName)
    createmRNA_GFP_Ribo_off_Rxn(model, compName)
    #createmRNA_RFP_Ribo_off_Rxn(model, compName)
    createTX_GFP_Rxn(model, compName)
    #createTX_RFP_Rxn(model, compName)
    createTL_GFPRxn(model, compName)
    #createTL_RFPRxn(model, compName)
    createRegen_Rxn(model, compName)
    createTX_GFP_Escape_Rxn(model, compName)
    #createTX_RFP_Escape_Rxn(model, compName)
    createTL_GFP_Escape_Rxn(model, compName)
    #createTL_RFP_Escape_Rxn(model, compName)
    createRibo_mRNA_GFPDegRxn(model, compName)
    #createRibo_mRNA_RFPDegRxn(model, compName)
    createRibo_mRNA_GFP_elongDegRxn(model, compName)
    #createRibo_mRNA_RFP_elongDegRxn(model, compName)

    createRibo_DegRxn(model, compName)
    createRibo_elong_Deg_GFP_Rxn(model, compName)
    #createRibo_elong_Deg_RFP_Rxn(model, compName)
    createRibo_Deg_GFP_Rxn(model, compName)
    #createRibo_Deg_RFP_Rxn(model, compName)

    return sbmlDoc

def main (args):
    sbmlDoc = None;
    SBMLok = False;
    default_root = "bmegModel_xylose"


    letters = 'r:'
    opts, params = getopt.getopt(sys.argv[1:], letters)
    for o, p in opts:
        if o == "-r":
            default_root = p




    sbmlDoc = createModel();
    SBMLok = validateExampleSBML(sbmlDoc);
    #writeExampleSBML(sbmlDoc, "bmegModel.xml");
    if (SBMLok):
        print "writing"
        writeExampleSBML(sbmlDoc, default_root+".xml");
    if (not SBMLok):
        print "ERROR"
        return 1;

    # get rid of functions
    config = ConversionProperties()
    if config != None:
        config.addOption('expandFunctionDefinitions')

    status = sbmlDoc.convert(config)
    if status != LIBSBML_OPERATION_SUCCESS:
        # Handle error somehow.
        print('Error: conversion failed due to the following:')
        sbmlDoc.printErrors()
        return 1

    SBMLok = validateExampleSBML(sbmlDoc);
    if (SBMLok):
        print "writing"
        writeExampleSBML(sbmlDoc, default_root+"_nofuncs.xml");
    if (not SBMLok):
        print "ERROR"
        return 1;

    return 0

  
if __name__ == '__main__':
  main(sys.argv)  
