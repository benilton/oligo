test_selector <- function(){
    library(oligo)
    library(oligoData)
    library(RUnit)
    pkgs <- c('pd.huex.1.0.st.v2', 'pd.hugene.1.0.st.v1')
    message('Loading reference data')
    load(system.file('unitTests', 'fids_ref0.rda', package='oligo'))
    fld <- c('fid', 'fsetid', 'type')
    out <- c('fid', 'man_fsetid')
    message('Loading sample dataset: Exon')
    data(affyExonFS)
    message('Getting probe info: core')
    coreExon <- getProbeInfo(affyExonFS, field=fld, target='core')[out][icore0[[pkgs[1]]],]
    message('Getting probe info: probeset')
    psetExon <- getProbeInfo(affyExonFS, field=fld, target='probeset')[out][ipset0[[pkgs[1]]],]
    message('Getting probe info: antigenomic backgroung probes')
    agenExon <- getProbeInfo(affyExonFS, field=fld, target='probeset', subset= type == 'control->bgp->antigenomic')[out][iagen0[[pkgs[1]]],]
    rownames(coreExon) <- rownames(psetExon) <- rownames(agenExon) <- NULL
    message('Loading sample dataset: Gene')
    data(affyGeneFS)
    message('Getting probe info: core')
    coreGene <- getProbeInfo(affyGeneFS, field=fld, target='core')[out][icore0[[pkgs[2]]],]
    message('Getting probe info: probeset')
    psetGene <- getProbeInfo(affyGeneFS, field=fld, target='probeset')[out][ipset0[[pkgs[2]]],]
    message('Getting probe info: antigenomic backgroung probes')
    agenGene <- getProbeInfo(affyGeneFS, field=fld, target='probeset', subset= type == 'control->bgp->antigenomic')[out][iagen0[[pkgs[2]]],]
    rownames(coreGene) <- rownames(psetGene) <- rownames(agenGene) <- NULL
    core <- list(coreExon, coreGene)
    pset <- list(psetExon, psetGene)
    agen <- list(agenExon, agenGene)
    names(core) <- names(pset) <- names(agen) <- pkgs
    checkEquals(core, core0) & checkEquals(pset, pset0) & checkEquals(agen, agen0)
}
