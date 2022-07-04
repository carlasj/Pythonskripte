import os

#os.system fuehrt einen Befehl direkt aus, damit kannst du mehrere Befehle hintereinanderschreiben. Das ist mega praktisch wenn du mehrere Proben hast

#Sam sortieren und nach bam konvertieren und die bam Datei in den Ordner bam_file_folder unter dem Namen sample_name speichern



def Gliom_Prozess (fastq_datei, sample_name):
	#os.system("cat " + fastq_datei +"*.fastq.gz > /vol/mf-tumor/Gliom/" + sample_name +".fastq.gz")
	
	#os.system("/vol/mf-tumor/tools/Porechop/./porechop-runner.py -i /vol/mf-tumor/Gliom/" + sample_name +".fastq.gz -o /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz")
	#os.system("minimap2 -a -x map-ont -t 14 /vol/mf-tumor/reference/HG38.canon.fasta /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz | samtools view -b | samtools sort  > /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +".bam")
	#os.system("samtools sort -@ 16 /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-canon.sam -o /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-canon-sort.sam")
	#os.system("samtools view -bh -@ 16 /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-canon-sort.sam -o /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + ".bam")
	#os.system("samtools index /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + ".bam")
	#os.system("samtools bedcov -Q 1 /vol/mf-tumor/intervals/Gliom-Mutationen.bed /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +".bam > /vol/mf-tumor/Gliom/Gliom-Amplikons/" + sample_name +"-Amplikons.bed")
	os.system("/homes/carlasj/.conda/envs/condaverzeichnis/bin/longshot --max_cov 100000 --bam /vol/mf-tumor/2022_02_WGS_LSK/" + sample_name + ".bam --ref /vol/mf-tumor/reference/Homo_sapiens_assembly38.fasta --out /vol/mf-tumor/Gliom/" + sample_name + "-2-longshot2.vcf")
	#os.system("/vol/mf-tumor/tools/Sniffles-master/bin/sniffles-core-1.0.12 sniffles -s 1 -m /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +".bam -v /vol/mf-tumor/Gliom/" + sample_name + "-sniffles.vcf")
	os.system("java -Xmx8g -jar /homes/carlasj/snpEff/snpEff.jar GRCh38.99 /vol/mf-tumor/Gliom/" + sample_name + "-2-longshot2.vcf > /vol/mf-tumor/Gliom/" + sample_name + "-2-snpEff.vcf")
	os.system("cat /vol/mf-tumor/Gliom/" + sample_name + '-2-snpEff.vcf | java -jar /homes/carlasj/snpEff/SnpSift.jar filter  "( QUAL >= 50 )" > /vol/mf-tumor/Gliom/' +sample_name +"-2-snpEff-filtered.vcf")
	os.system("java -jar /homes/carlasj/snpEff/SnpSift.jar filter \"(ANN[*].EFFECT has 'missense_variant') \" /vol/mf-tumor/Gliom/" + sample_name + "-2-snpEff-filtered.vcf > /vol/mf-tumor/Gliom/Gliom-Longshot-snpEff/" + sample_name + "-2-variants.vcf")
	os.system("java -jar /homes/carlasj/snpEff/SnpSift.jar extractFields /vol/mf-tumor/Gliom/Gliom-Longshot-snpEff/" + sample_name + "-2-variants.vcf ANN[0].GENE ANN[0].EFFECT ANN[0].HGVS_P ANN[0].HGVS_C CHROM POS REF ALT QUAL DP AC > /vol/mf-tumor/Gliom/Gliom-Longshot-snpEff/" + sample_name + "-2-snpEff-alle_Mutationen.vcf")


	#os.system("samtools bedcov -Q 1 /vol/mf-tumor/intervals/hg38-halbe-Chromosomen.bed /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +".bam > /vol/mf-tumor/Gliom/Gliom-Amplikons/" + sample_name +"-hg38-halbe-Chromosomen.bed")
#def Samtools_ngmlr_bedcov(sample_name):
	#os.system("/vol/mf-tumor/tools/ngmlr-0.2.7/ ngmlr -t 4 -r /vol/mf-tumor/reference/HG38.canon.fasta -q /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz -o /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +"-ngmlr.sam -x ont")
	#os.system("samtools sort -@16 /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-ngmlr.sam -o /vol/mf-tumor/Gliom/Gliom-Mappings/" +sample_name + "-ngmlr-sort.sam")
	#os.system("samtools view -bh -@ 16 /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-ngmlr-sort.sam -o /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-ngmlr.bam")
	#os.system("samtools index /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + "-ngmlr.bam")
	#os.system("samtools bedcov -Q -1 /vol/mf-tumor/intervals/hg38-halbe-Chromosomen.bed /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +"-ngmlr.bam > /vol/mf-tumor/Gliom/Gliom-Amplikons/" + sample_name +"-ngmlr-halbe-Chromosomen.bed")
	#os.system("bedtools intersect -c -a /vol/mf-tumor/intervals/hg38-halbe-Chromosomen.bed -b /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +"-ngmlr.bam > /vol/mf-tumor/Gliom/Gliom-Amplikons/" + sample_name +"-ngmlr-halbe-Chromosomen-intersect.bed")

def Gliom_Prozess_MGMT (fastq_datei, sample_name):
	os.system("cat " + fastq_datei +"*.fastq.gz > /vol/mf-tumor/Gliom/" + sample_name +".fastq.gz")
	
	#os.system("/vol/mf-tumor/tools/Porechop/./porechop-runner.py -i /vol/mf-tumor/Gliom/" + sample_name +".fastq.gz -o /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz")
	#os.system("minimap2 -a -x map-ont -t 16 /vol/mf-tumor/Gliom/MGMT_Promotor/Reference.fasta /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz | samtools view -b | samtools sort   > /vol/mf-tumor/Gliom/MGMT_Promotor/" + sample_name +".bam")
	#os.system("samtools sort -@ 16 /vol/mf-tumor/Gliom/MGMT_Promotor/" + sample_name + "-canon.sam -o /vol/mf-tumor/Gliom/MGMT_Promotor/" + sample_name + "-canon-sort.sam")
	#os.system("samtools view -bh -@ 16 /vol/mf-tumor/Gliom/MGMT_Promotor/" + sample_name + "-canon-sort.sam -o /vol/mf-tumor/Gliom/MGMT_Promotor/" + sample_name + ".bam")
	#os.system("samtools index /vol/mf-tumor/Gliom/MGMT_Promotor/" + sample_name + ".bam")
	
	#os.system("samtools mpileup /vol/mf-tumor/Gliom/MGMT_Promoter/" + sample_name + ".bam  -f /vol/mf-tumor/Gliom/MGMT_Promoter/Reference.fasta -l /vol/mf-tumor/Gliom/MGMT_Promoter/MGMT_CpG.bed > /vol/mf-tumor/Gliom/MGMT_Promoter/" +sample_name + "-Pileup.vcf")

def Gliom_Prozess_wgs (fastq_datei, sample_name):
	#os.system("cat " + fastq_datei +"*.fastq.gz > /vol/mf-tumor/Gliom/" + sample_name +".fastq.gz")
	
	#os.system("/vol/mf-tumor/tools/Porechop/./porechop-runner.py -i /vol/mf-tumor/Gliom/" + sample_name +".fastq.gz -o /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz")
	#os.system("minimap2 -a -x map-ont -t 14 /vol/mf-tumor/reference/HG38.canon.fasta /vol/mf-tumor/Gliom/" + sample_name +"-porechop.fastq.gz | samtools view -b | samtools sort  > /vol/mf-tumor/Gliom/wgs/" + sample_name +".bam")

	#os.system("samtools index /vol/mf-tumor/Gliom/wgs/" + sample_name + ".bam")
	#os.system("rm /vol/mf-tumor/Gliom/" + sample_name + "-porechop.fastq.gz")
	#os.system("rm /vol/mf-tumor/Gliom/" + sample_name + ".fastq.gz")
	os.system("samtools bedcov -Q 1 /vol/mf-tumor/intervals/hg38-chromosomband2.bed /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name +".bam > /vol/mf-tumor/Gliom/Gliom-Amplikons/" + sample_name +"-hg38-chromosomband.bed")
	#os.system("samtools bedcov -Q 1 /vol/mf-tumor/intervals/hg38-halbe-Chromosomen.bed /vol/mf-tumor/Gliom/wgs/" + sample_name +".bam > /vol/mf-tumor/Gliom/wgs/" + sample_name +"-hg38-halbe-Chromosomen.bed")
	#os.system("samtools bedcov -Q 1 /vol/mf-tumor/intervals/wgs_calling_small_regions.hg38.fragment.10kB.bed /vol/mf-tumor/Gliom/wgs/" + sample_name +".bam > /vol/mf-tumor/Gliom/wgs/" + sample_name + "-10kB.bed")

def Zelllinien (fastq_datei, sample_name):
	#os.system("samtools view -bh " + sample_name + ".bam 
	os.system("bcftools mpileup -d 1000000 -f /vol/mf-tumor/reference/Homo_sapiens_assembly38.fasta /vol/mf-tumor/Gliom/Gliom-Mappings/" + sample_name + ".bam | bcftools call -mv -Ob -o /vol/mf-tumor/Gliom/" + sample_name + "-calls.bcf")
	os.system("bcftools convert -Ov /vol/mf-tumor/Gliom/" + sample_name + "-calls.bcf -o /vol/mf-tumor/Gliom/" + sample_name + "-calls.vcf")
	
	os.system("java -Xmx8g -jar /homes/carlasj/snpEff/snpEff.jar GRCh38.99 /vol/mf-tumor/Gliom/" + sample_name + "-calls.vcf > /vol/mf-tumor/Gliom/" + sample_name + "-calls-snpEff.vcf")
	os.system("cat /vol/mf-tumor/Gliom/" + sample_name + '-calls-snpEff.vcf | java -jar /homes/carlasj/snpEff/SnpSift.jar filter  "( QUAL >= 50 )" > /vol/mf-tumor/Gliom/' +sample_name +"-calls-snpEff-filtered.vcf")
	os.system("java -jar /homes/carlasj/snpEff/SnpSift.jar filter \"(ANN[*].EFFECT has 'missense_variant') \" /vol/mf-tumor/Gliom/" + sample_name + "-calls-snpEff-filtered.vcf > /vol/mf-tumor/Gliom/Gliom-Longshot-snpEff/" + sample_name + "-calls-variants.vcf")
	os.system("java -jar /homes/carlasj/snpEff/SnpSift.jar extractFields /vol/mf-tumor/Gliom/" + sample_name + "-calls-snpEff-filtered.vcf ANN[0].GENE ANN[0].EFFECT ANN[0].HGVS_P ANN[0].HGVS_C CHROM POS REF ALT QUAL DP AC > /vol/mf-tumor/Gliom/Gliom-Longshot-snpEff/" + sample_name + "-calls-Mutationen.vcf")



#def longshot ("bam_file", "sample_name", "bam_file_folder"):
	#os.system("cd /homes/carlasj/.conda/envs/condaverzeichnis/bin")
	#os.system("longshot --bam " + bam_file + "--ref /vol/mf-tumor/reference/Homo_sapiens_assembly38.fasta --out /vol/mf-tumor/Gliom" +sample_name+ "-longshot.vcf")
	
#def snpEff ():
	#os.system("/homes/carlasj/snpEff java -Xmx8g -jar snpEff.jar GRCh38.99 /vol/mf-tumor/Gliom/test-GliomRun-longshot.vcf > /vol/mf-tumor/Gliom/Test-Gliom-snpEff.vcf")



#Gliom_Prozess_MGMT("/vol/mf-tumor/FASTQ/Run00987_MIN106_RBK004/fastq_pass/barcode01/", "Run987-1")
#Gliom_Prozess_MGMT("/vol/mf-tumor/FASTQ/Run00987_MIN106_RBK004/fastq_pass/barcode02/", "Run987-2")
#Gliom_Prozess_MGMT("/vol/mf-tumor/FASTQ/Run00987_MIN106_RBK004/fastq_pass/barcode03/", "Run987-3")



Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-1")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-2")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-3")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-4")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-5")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-6")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-7")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-8")
Gliom_Prozess_wgs("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode04/", "Run822-9")
#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode05/", "Gliom_12")
#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run01026_MIN106_RBK004/fastq_pass/barcode06/", "Gliom_16")
#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run00982_MIN106_RBK004/fastq_pass/barcode02/", "Run982-02")




#Gliom_Prozess("/vol/mf-tumor/fastq_pass/barcode03/", "Gliom_15")


#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run00982_MIN106_RBK004/fastq_pass/barcode02/", "Run982-2")

#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run00968_MIN106_RBK004/no_sample/20220330_1306_X1_FAS48140_e3227cce/fastq_pass/barcode02/", "Run968-2")


#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run00937_MIN106_RBK004/fastq_pass/barcode04/", "Run937-4")
#Gliom_Prozess("/vol/mf-tumor/FASTQ/Run00937_MIN106_RBK004/fastq_pass/barcode05/", "Run937-5")
#Gliom_Prozess_MGMT("/vol/mf-tumor/Gliom/fastq_pass/834-barcode04/", "Run834-4")

#Samtools_ngmlr_bedcov("Run822-1")
#Samtools_ngmlr_bedcov("Run822-2")
#Samtools_ngmlr_bedcov("Run822-3")
#Samtools_ngmlr_bedcov("Run822-4")
#Samtools_ngmlr_bedcov("Run822-5")
#Samtools_ngmlr_bedcov("Run822-6")
#Samtools_ngmlr_bedcov("Run822-7")
#Samtools_ngmlr_bedcov("Run822-9")

