setMethod("show", "PDInfo", function(object) {
    cat("Class........:", class(object), "\n")
    cat("Manufacturer.:", manufacturer(object), "\n")
    cat("Genome Build.:", genomeBuild(object), "\n")
    cat("Chip Geometry:", geometry(object)[1], "rows x ", geometry(object)[2], "columns\n")
})

setMethod("show", "platformDesign", function(object) {
    callNextMethod()
    cat("Array type:", object@type, "\n")
    cat("Number of rows:", nrow(object), "\n")
    cat("Number of columns:", ncol(object), "\n")
    cat("Number of features:", nProbes(object), "\n")
})
