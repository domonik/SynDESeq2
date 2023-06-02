remotes::install_github("YuLab-SMU/clusterProfiler")
remotes::install_github("YuLab-SMU/createKEGGdb")
library(createKEGGdb)
species <- c(snakemake@config["keggOrgID"])
if (require("KEGG.db")) {
    remove.packages("KEGG.db")
}
create_kegg_db2 <- function(species, path) {
    packagedir <- tempfile() # tempdir() maynot empty

    ## skeleton
    prepare_pkg_skeleton(packagedir)

    ## sqlite
    sqlite_path <- paste(packagedir, "inst", "extdata", sep=.Platform$file.sep)
    prepare_kegg_db(species, sqlite_path)

    ## build pkg
    pkgbuild::build(packagedir, dest_path = path)
}

createKEGGdb::create_kegg_db2(species, snakemake@params[["localKEGGdb"]])
install.packages(file.path(snakemake@params[["localKEGGdb"]] , "/KEGG.db_1.0.tar.gz"), repos=NULL, type="source")