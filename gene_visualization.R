library(cBioPortalData) #cbioportal api ����� ���� ���̺귯��
library(AnVIL) #bioconductor ����� ���� ���̺귯��
library(g3viz) #gene mutation data�� �ð�ȭ�ϱ����� ���̺귯��
library(maftools) #maf������ �ٷ�� ���� ���̺귯�� 
library(TCGAbiolinks) #GDC �����Ϳ� ���� �м��� ���� bioconductor���̺귯��
library(dplyr) #������ó���� ���� ���̺귯��

maf_mutect2 <- GDCquery_Maf("ACC", pipelines = "mutect2") %>% subset(Sequencer == "Illumina HiSeq 2000") %>% read.maf 
#TCGA �� ������ ���̽����� ���ϴ� �Ͽ� ���� ���� ������
getGeneSummary(maf_mutect2) #�� ������ ��� ���

plotmafSummary(maf = maf_mutect2 , rmOutlier = TRUE, addStat = 'median', dashboard = TRUE) #�ϵ����� �ð�ȭ

oncoplot(maf = maf_mutect2 , top = 10, removeNonMutated = TRUE) #���� 10�� mutation gene�� ���� oncoplot

titv = titv(maf = maf_mutect2 , plot = FALSE, useSyn = TRUE) #mutation data�� SNP data  
plotTiTv(res = titv) #SNP data �ð�ȭ

OncogenicPathways(maf = maf_mutect2 ) #data�� Oncogenic Signaling Pathways �ð�ȭ
PlotOncogenicPathways(maf = maf_mutect2 , pathways = "RTK-RAS") #������� �����ڴ� ���������� ���� �����ڴ� �Ķ������� �ð�ȭ

cbio <- cBioPortal() 
cbio

cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)

getStudies(api, buildReport = FALSE)


studies <- getStudies(cbio)
studies #cbioportal ���� �� ������ �ҷ�����

cat("cancer types:",length(unique(studies$cancerTypeId))) #cancer type�� ���� ���� ���

cancerTypeCounts <- 
  studies %>%
  group_by(cancerTypeId)%>%
  summarise(totalSamples=sum(allSampleCount)) %>%
  arrange(desc(totalSamples)) %>%
  top_n(20)

cancerTypeCounts #cancer type ���� sample ���� ���

cancerTypeCounts <- cancerTypeCounts %>% arrange(cancerTypeCounts$totalSamples)
par(mar=c(4,6,1,1))
barplot(cancerTypeCounts$totalSamples,names=cancerTypeCounts$cancerTypeId,
        main = "Number of samples by Top 20 Primary Sites",
        horiz = TRUE,
        las =1) #cancer type�� sample ���� �ð�ȭ

mutation.dat <- g3viz::getMutationsFromCbioportal("acc_tcga", "TP53") #���ϴ� �� ������ ���õ� gene mutation data �ҷ�����
g3Lollipop(mutation.dat, gene.symbol = "TP53") #gene muation data �ð�ȭ