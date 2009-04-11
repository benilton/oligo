setMethod("getX", "DBPDInfo",
          function(object, type){
            stopifnot(!missing(type))
            sql <- "SELECT fid, x FROM"
            if (type == "pm"){
              sql <- paste(sql, "pmfeature")
            }else if(type == "bg"){
              sql <- paste(sql, "bgfeature")
            }else{
              stop("Method not implemented for type ", type)
            }
            res <- dbGetQuery(db(object), sql)
            res[order(res[["fid"]]), "x"]
          })

setMethod("getY", "DBPDInfo",
          function(object, type){
            stopifnot(!missing(type))
            sql <- "SELECT fid, y FROM"
            if (type == "pm"){
              sql <- paste(sql, "pmfeature")
            }else if(type == "bg"){
              sql <- paste(sql, "bgfeature")
            }else{
              stop("Method not implemented for type ", type)
            }
            res <- dbGetQuery(db(object), sql)
            res[order(res[["fid"]]), "y"]
          })
