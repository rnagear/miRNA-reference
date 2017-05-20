## Experimentally-validated microRNA:target gene data ##

These are flat text, tab-delimited files containing experimentally-validated microRNA:target gene interactions.
Included files:

  * `complete_set_validated_interactions.txt`
    - all experimentally-validated interaction data collected for this study
    - columns: miRID (microRNA ID), GeneName, Type ([pos]itive or [neg]ative interaction 
      between the miRNA and the gene), and Source (the type od experimental validation)
        * VALID: publication databases (miRTarBase and TarBase)
        * PROEO: validated from high-throughput proteomics data of protein changes on miRNA overexpression
        * EXPRS: validated from gene-expression changes in array data on miRNA overexpression
        * OEKD: gene expression changes measured by array from overexpression and knockdown of miRNA in singe set of experiments
        * CLIP: cross-linking and immunoprecipitation data, where the identity of the interacting miRNA is independently known
        
  * `train_data.txt`
    - data used to train the machine learning models to classify miRNA:target interactions
    - columns: Type (as above), the remaining columns (diana, mirmap, mirtarget, mirza_U, paccmit, targetscan_C, targetscan_NC)
      correspond to the prediction scores of the input algorithms (see associated publication for further details). The \_U suffix
      indicates the union set of all MIRZA predictions, tne \_C and \_NC indicate the TargetScan conserved, and non-conserved
      interactions respectively.
      
  * `test_data.txt`
    - data unseen during the training of the models used to evaluate performance of the resulting integrated dataset, and for
      the other input algorithms (Diana microT-CDS, miRMap, miRTarget, MIRZA, PACCMIT-CDS, and TargetScan)
    - columns: GeneName, miRID (as above), Type (as above), Source (as above)
