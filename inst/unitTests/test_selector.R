test_selector <- function(){
    library(oligo)
    library(oligoData)
    library(RUnit)
    pkgs <- c('pd.huex.1.0.st.v2', 'pd.hugene.1.0.st.v1')
    getRandom <- function(dt){
        set.seed(1)
        i <- sort(sample(nrow(dt), 100))
        d0 <- dt[i,]
        rownames(d0) <- NULL
        d0
    }
    fld <- c('fid', 'fsetid', 'type')
    out <- c('fid', 'man_fsetid')
    data(affyExonFS)
    coreExon <- getProbeInfo(affyExonFS, field=fld, target='core')[out]
    coreExon <- getRandom(coreExon)
    psetExon <- getProbeInfo(affyExonFS, field=fld, target='probeset')[out]
    psetExon <- getRandom(psetExon)
    agenExon <- getProbeInfo(affyExonFS, field=fld, target='probeset', subset= type == 'control->bgp->antigenomic')[out]
    agenExon <- getRandom(agenExon)

    data(affyGeneFS)
    coreGene <- getProbeInfo(affyGeneFS, field=fld, target='core')[out]
    coreGene <- getRandom(coreGene)
    psetGene <- getProbeInfo(affyGeneFS, field=fld, target='probeset')[out]
    psetGene <- getRandom(psetGene)
    agenGene <- getProbeInfo(affyGeneFS, field=fld, target='probeset', subset= type == 'control->bgp->antigenomic')[out]
    agenGene <- getRandom(agenGene)
    core <- list(coreExon, coreGene)
    pset <- list(psetExon, psetGene)
    agen <- list(agenExon, agenGene)
    names(core) <- names(pset) <- names(agen) <- pkgs
    load('fids_ref0.rda')
    checkEquals(core, core0) & checkEquals(pset, pset0) & checkEquals(agen, agen0)
}
