library(cBioPortalData) #cbioportal api 사용을 위한 라이브러리
library(AnVIL) #bioconductor 사용을 위한 라이브러리
library(g3viz) #gene mutation data를 시각화하기위한 라이브러리
library(maftools) #maf파일을 다루기 위한 라이브러리 
library(TCGAbiolinks) #GDC 데이터와 통합 분석을 위한 bioconductor라이브러리
library(dplyr) #데이터처리를 위한 라이브러리

maf_mutect2 <- GDCquery_Maf("ACC", pipelines = "mutect2") %>% subset(Sequencer == "Illumina HiSeq 2000") %>% read.maf 
#TCGA 암 데이터 베이스에서 원하는 암에 대한 파일 가져옴
getGeneSummary(maf_mutect2) #암 데이터 요약 출력

plotmafSummary(maf = maf_mutect2 , rmOutlier = TRUE, addStat = 'median', dashboard = TRUE) #암데이터 시각화

oncoplot(maf = maf_mutect2 , top = 10, removeNonMutated = TRUE) #상위 10개 mutation gene에 대한 oncoplot

titv = titv(maf = maf_mutect2 , plot = FALSE, useSyn = TRUE) #mutation data의 SNP data  
plotTiTv(res = titv) #SNP data 시각화

OncogenicPathways(maf = maf_mutect2 ) #data의 Oncogenic Signaling Pathways 시각화
PlotOncogenicPathways(maf = maf_mutect2 , pathways = "RTK-RAS") #종양억제 유전자는 빨간색으로 종양 유전자는 파란색으로 시각화

cbio <- cBioPortal() 
cbio

cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)

getStudies(api, buildReport = FALSE)


studies <- getStudies(cbio)
studies #cbioportal 연결 및 데이터 불러오기

cat("cancer types:",length(unique(studies$cancerTypeId))) #cancer type의 종류 갯수 출력

cancerTypeCounts <- 
  studies %>%
  group_by(cancerTypeId)%>%
  summarise(totalSamples=sum(allSampleCount)) %>%
  arrange(desc(totalSamples)) %>%
  top_n(20)

cancerTypeCounts #cancer type 별로 sample 갯수 출력

cancerTypeCounts <- cancerTypeCounts %>% arrange(cancerTypeCounts$totalSamples)
par(mar=c(4,6,1,1))
barplot(cancerTypeCounts$totalSamples,names=cancerTypeCounts$cancerTypeId,
        main = "Number of samples by Top 20 Primary Sites",
        horiz = TRUE,
        las =1) #cancer type의 sample 갯수 시각화

mutation.dat <- g3viz::getMutationsFromCbioportal("acc_tcga", "TP53") #원하는 암 종류와 관련된 gene mutation data 불러오기
g3Lollipop(mutation.dat, gene.symbol = "TP53") #gene muation data 시각화
