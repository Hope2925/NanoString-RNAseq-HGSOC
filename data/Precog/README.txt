* Bonome_GSE26712: 
	* 11,979 genes & 185 samples
	* Meta: DSS_Status and DSS_Time
		* Split between 90 optimal and 95 suboptimal
		* If from GEO: tissue (Normal or HGSOC), surgery outcome (Optimal or Suboptimal), status (AWD,NED,DOD), survival years
	* Advanced stage, high-grade papillary serous ovarian cancer
* Crijnns_GSE13876:
	* Obtained at primary surgery prior to chemotherapy
	* 37,632 genes & 415 samples
	* Meta: OS_Status and OS_Time
		* If from GEO: patient id, status of 0/1, fumnd, age, sample nr, 
* Mok_GSE18521:
	* 17,788 genes & 76 samples
	* Meta: OS_Status and OS_Time
		* If from GEO: Sample Title (some tumor, normal, Cell line), tissue, treatment, tumor stage, tumor grade, surv data)
* Yoshishara_GSE17260:
	* 41,000 genes & 110 samples
	* Meta: OS_Status and OS_Time
		* If from GEO: tumor grade, Stage, cytoreductive surgery (optimal or not optimal), progression-free survival (m), recurrence (1), overall survival (m), death (1), disease state
	* Also downloaded GSE17260_clinical_information.txt: sample ID, stage, cytoreduce survery, Progression-free survival (M), recurrence (1), Overall survival (M), death (1)
* Yoshishara_GSE32062:
	* 41,093 genes & 260 samples
	* Meta: OS_Status and OS_Time
		* If from GEO: tissue, grading, Stage, surgery status (optimal or suboptimal), taxane, platinum, pfs (m), rec (1), os (m), death (1)
	* Also downloaded GSE32062_clinical_information.txt: sample name, Disease Histology, Grading, stage, Surgery status, Taxane, Platinum, PFS (M), Rec(1), OS (M), Death (1)
		* Note that some do in fact overlap with GSE17260 but not all
