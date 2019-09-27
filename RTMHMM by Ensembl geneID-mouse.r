#GeneID <- Ensembl gene IDs.csv, should be a csv file address str


RTMHMM <- function(GeneID,b="~/RTMHMM/AA Seq.fa",c="~/RTMHMM/Results Of TMHMM.txt",d="~/RTMHMM/TMHMM ProteinID.csv",e="~/RTMHMM/TMHMM GeneID.csv"){

  library(biomaRt)
  require(reshape2)
  library(RSelenium)


  dir.create("~/RTMHMM")

  # input1
  EnsemblGeneID <- read.csv(file=GeneID, header=F, stringsAsFactors = FALSE)
  # define biomart object
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  # query biomart
  results <- apply(EnsemblGeneID, 2, function(x){
    getBM(attributes = "ensembl_peptide_id",
          filters = "ensembl_gene_id", values = x,
          mart = mart)})
  results <- melt(results)
  results <- results[,1]
  protein <- getSequence(id= results,
                         type="ensembl_peptide_id",
                         seqType="peptide",
                         mart= mart)
  #output1
  exportFASTA(protein, file = "~/RTMHMM/AA Seq.fa")

  remDr <- remoteDriver(
    remoteServerAddr = "localhost",
    port = 4444,
    browserName = "chrome")
  remDr$open()
  remDr$navigate("http://www.cbs.dtu.dk/services/TMHMM/")
  webElem <- remDr$findElement(using ='xpath','/html/body/table/tbody/tr/td[2]/table[2]/tbody/tr/td/form/p[3]/input[3]')
  webElem$clickElement()

  webElem <- remDr$findElement(using = 'xpath','/html/body/table/tbody/tr/td[2]/table[2]/tbody/tr/td/form/p[1]/input')
  # input2
  path_upload <- normalizePath("~/RTMHMM/AA Seq.fa")
  webElem$sendKeysToElement(list(path_upload))
  webElem <- remDr$findElement(using = 'xpath','/html/body/table/tbody/tr/td[2]/table[2]/tbody/tr/td/form/p[5]/input[1]')
  webElem$clickElement()

  webElem <-NULL
  while(is.null(webElem)){
    webElem <- tryCatch({remDr$findElement(using = 'xpath', value = '/html/body/a')},
                        error = function(e){NULL})
    Sys.sleep(10)
  }

  webElem <- remDr$findElement(using = 'xpath','/html/body/pre')
  Results_TMHMM <- webElem$getElementText()[[1]]
  Results_TMHMM <- melt(Results_TMHMM)

  #output2
  write.table(Results_TMHMM, file=c, row.names = F, col.names = F, quote=F)

  # input1
  my_data <- read.table(c, header = F)
  colname_sp <- c("ensembl_peptide_id", "Len","ExpAA",	"First60",	"PredHel", "Topology")

  names(my_data) <- colname_sp
  results <- grep("PredHel=0", my_data$PredHel)
  data <- my_data[-results,]

  geneID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_peptide_id","description"),
                  filters = "ensembl_peptide_id", values = data$ensembl_peptide_id,
                  mart = mart)

  MPs <- merge(data, geneID, by = "ensembl_peptide_id")
  #output3
  write.csv(MPs, file=d,quote=F)

  geneID1 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                   filters = "ensembl_peptide_id", values = data$ensembl_peptide_id,
                   mart = mart)
  #output4
  write.csv(geneID1, file=e,quote=F)

  print("Done, output 4 files in ~/RTMHMM")

}
