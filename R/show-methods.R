setMethod("show", "PDInfo", function(object) {
    cat("Class:", class(object), "\n")
    cat("Manufacturer:", manufacturer(object), "\n")
    cat("Genome Build:", genomeBuild(object), "\n")
})

setMethod("show", "platformDesign", function(object) {
    callNextMethod()
    cat("Array type:", type(object), "\n")
    cat("Number of rows:", nrow(object), "\n")
    cat("Number of columns:", ncol(object), "\n")
    cat("Number of features:", nProbes(object), "\n")
})
