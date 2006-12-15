setMethod("show", "PDInfo", function(object) {
    cat("Class:", class(object), "\n")
    cat("Manufacturer:", object@manufacturer, "\n")
    cat("Genome Build:", object@genomebuild, "\n")
})

setMethod("show", "platformDesign", function(object) {
    callNextMethod()
    cat("Array type:", object@type, "\n")
    cat("Number of columns:", object@ncol, "\n")
    cat("Number of rows:", object@nrow, "\n")
    cat("Number of features:", nProbes(object), "\n")
})
