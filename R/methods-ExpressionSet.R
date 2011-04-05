## Methods for ExpressionSet's

setMethod("boxplot", signature(x="ExpressionSet"),
          function(x, which, transfo=identity, nsample=10000, ...){
            stopifnot(is.function(transfo))
            if(!missing(which))
                warning("Argument 'which' ignored (not meaningful for ExpressionSet)")
            if (nrow(x) > nsample){
              idx <- sort(sample(nrow(x), nsample))
            }else{
              idx <- 1:nrow(x)
            }
            toPlot <- transfo(exprs(x[idx,]))
            toPlot <- as.data.frame(toPlot)
            dots <- list(...)
            dots[["x"]] <- toPlot
            rm(toPlot)
            if (is.null(dots[["col"]]))
              dots[["col"]] <- darkColors(ncol(x))
            if (is.null(dots[["range"]]))
              dots[["range"]] <- 0
            if (is.null(dots[["main"]]))
                dots[["main"]] <- "exprs"
            do.call("boxplot", dots)
          })

setMethod("hist", "ExpressionSet",
          function(x, transfo=identity, nsample=10000, ...){
            stopifnot(is.function(transfo))
            if (nrow(x) > nsample){
              idx <- sort(sample(nrow(x), nsample))
            }else{
              idx <- 1:nrow(x)
            }
            tmp <- transfo(exprs(x[idx,]))
            res <- matDensity(tmp)

            dots <- list(...)
            if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "density"
            if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "log-intensity"
            if (is.null(dots[["xlim"]])) dots[["xlim"]] <- range(res[["x"]])
            if (is.null(dots[["ylim"]])) dots[["ylim"]] <- range(res[["y"]])
            if (is.null(dots[["col"]]))
              dots[["col"]] <- darkColors(ncol(x))
            if (is.null(dots[["type"]])) dots[["type"]] <- "l"

            dots[["x"]] <- res[["x"]]
            dots[["y"]] <- res[["y"]]
            do.call("matplot", dots)

            invisible(res)
          })


setMethod("MAplot", "ExpressionSet",
          function(object, what=exprs, transfo=identity, groups, refSamples, which,
                   pch=".", summaryFun=rowMedians, plotFun=smoothScatter,
                   main="vs pseudo-median reference chip", pairs=FALSE, ...){
              stopifnot(is.function(what))
              maplot(x=what(object), transfo=transfo, groups=groups,
                     refSamples=refSamples, which=which, pch=pch,
                     summaryFun=summaryFun, main=main, pairs=pairs, ...)
          })
