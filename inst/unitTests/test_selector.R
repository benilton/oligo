test_selector <- function(){
    library(oligo)
    library(oligoData)
    library(RUnit)
    pkgs <- c('pd.huex.1.0.st.v2', 'pd.hugene.1.0.st.v1')
    load(system.file('unitTests', 'fids_ref0.rda', package='oligo'))
    fld <- c('fid', 'fsetid', 'type')
    out <- c('fid', 'man_fsetid')
    data(affyExonFS)
    coreExon <- getProbeInfo(affyExonFS, field=fld, target='core')[out][icore0[[pkgs[1]]],]
    psetExon <- getProbeInfo(affyExonFS, field=fld, target='probeset')[out][ipset0[[pkgs[1]]],]
    agenExon <- getProbeInfo(affyExonFS, field=fld, target='probeset', subset= type == 'control->bgp->antigenomic')[out][iagen0[[pkgs[1]]],]

    data(affyGeneFS)
    coreGene <- getProbeInfo(affyGeneFS, field=fld, target='core')[out][icore0[[pkgs[2]]],]
    psetGene <- getProbeInfo(affyGeneFS, field=fld, target='probeset')[out][ipset0[[pkgs[2]]],]
    agenGene <- getProbeInfo(affyGeneFS, field=fld, target='probeset', subset= type == 'control->bgp->antigenomic')[out][iagen0[[pkgs[2]]],]
    core <- list(coreExon, coreGene)
    pset <- list(psetExon, psetGene)
    agen <- list(agenExon, agenGene)
    names(core) <- names(pset) <- names(agen) <- pkgs
    checkEquals(core, core0) & checkEquals(pset, pset0) & checkEquals(agen, agen0)
}
