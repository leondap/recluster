biodecrypt<-function (mat, id, alpha = NULL, ratio = 2.5, buffer = 90000, polygon = NULL, 
    checkdist = T, minimum = 7, plot = T, map = NULL, xlim = NULL, 
    ylim = NULL, main = NULL) 
{
    res <- NULL
    res$type <- "sep"
    borders <- NULL
    taxa <- length(which(unique(id)>0))
    colnames(mat) <- c("Long", "Lat")
	usenondupl<-which(!(duplicated(cbind(mat,id))))
	mat<-mat[usenondupl,]
	id<-id[usenondupl]
	matpunti<-st_as_sf(as.data.frame(mat),coords = c("Long","Lat"), crs = 4326)
    distances <- matrix(0, nrow(mat), taxa)
    distances2 <- matrix(0, nrow(mat), taxa)
    if (is.null(alpha)) {
        alpha = rep(10, taxa)
    }
	if (length(alpha)==1) {
	alph<-alpha
        alpha = rep(alph, taxa)
    }

    if (is.null(xlim)) {
        xlim <- c(min(mat[, 1]), max(mat[, 1]))
    }
    if (is.null(ylim)) {
        ylim <- c(min(mat[, 2]), max(mat[, 2]))
    }
    if (plot) {
        plot(cbind(xlim, ylim), type = "n", main = main)
        if (!is.null(map)) {
            plot(map, add = T)
        }
    }
    vectab <- matrix(0, nrow(mat), taxa + 1)
    hulls <- list()
    hullpl <- NULL
    areas <- NULL
    oldw <- getOption("warn")
    options(warn = -1)
    for (spec in 1:taxa) {
#spec<-1
        taxsp <- which(id == spec)
        hullas <- matpunti[taxsp, ]
        if (nrow(hullas) >= minimum) {
		hullasp<-mat[taxsp, ]
            hull <- ahull(hullasp, alpha = alpha[spec])
            hull2 <- ah2sf(hull)
            hullspat<- st_as_sf(hull2)
            if (!(is.null(polygon))) {
                hullspat <- st_intersection(hullspat, polygon)
            }
        }
        if (nrow(hullas) < minimum) {
            hull<-st_convex_hull(st_make_valid(st_union(hullas)))
		hullspat<-hull
		if (!(is.null(polygon))) {
                hullspat <- st_intersection(hull, polygon)
            }
        }
        hulls[[spec]] <- hullspat
        hullpl[[spec]] <- hull
        if (plot) {
            plot(hull, add = T)
		points(hullas)

        }
        areas[spec] <- sum(st_area(hulls[[spec]]))
	pointin <- unlist(st_contains(hulls[[spec]], matpunti))
        vectab[pointin, spec] <- 1
	dist_matrix <- st_distance(matpunti , hullspat)
		distances[,spec]  <- apply(dist_matrix, 1, min)
        
    }

    vectab[, ncol(vectab)] <- rowSums(vectab[, 1:taxa])
    id2 <- id
    uncertain1 <- which(vectab[, ncol(vectab)] > 1 & id == 0)
	#uncertain1 are those within two hulls
    uncertain2 <- which(vectab[, ncol(vectab)] == 0 & id == 0)
	#uncertain2 are those outside any hull
    inside <- which(vectab[, ncol(vectab)] == 1 & id == 0)
	#inside are those within a single hull
    if (length(uncertain2) > 1) {
        distancesunc <- distances[uncertain2, ]
        order <- matrix(NA, length(uncertain2), taxa)
        for (unc2 in 1:length(uncertain2)) {
            wh <- uncertain2[unc2]
            order[unc2, 1:taxa] <- c(1:taxa)[order(distancesunc[unc2, 
                1:taxa])]
        }
        attribution <- order[, 1]
        for (k in 1:length(uncertain2)) {
            dist <- distancesunc[k, ]
            ordereddist <- dist[order(dist)]
            if (ordereddist[2] > buffer & (ordereddist[2]/ordereddist[1]) > 
                ratio) {
                id2[uncertain2[k]] <- attribution[k]
            }
        }
    }
    if (length(uncertain2) == 1) {
        distancesunc <- distances[uncertain2, ]
        order <- matrix(NA, 1, taxa)
        order[1, 1:taxa] <- c(1:taxa)[order(distancesunc[1:taxa])]
        attribution <- order[, 1]
        ordereddist <- distancesunc[order(distancesunc)]
        if (ordereddist[2] > buffer & (ordereddist[2]/ordereddist[1]) > 
            ratio) {
            id2[uncertain2] <- attribution
        }
    }
    if (length(inside) > 1) {
        distancesunc <- distances[inside, ]
        for (k in 1:length(inside)) {
            attr <- which(vectab[inside[k], 1:taxa] == 1)
            diste <- distancesunc[k, ]
            diste <- diste[-attr]
            if (min(diste) > buffer) {
                id2[inside[k]] <- attr
            }
        }
    }
    if (length(inside) == 1) {
        distancesunc <- distances[inside, ]
        attr <- which(vectab[inside, 1:taxa] == 1)
        diste <- distancesunc[-attr]
        if (min(diste) > buffer) {
            id2[inside] <- attr
        }
    }
    if (checkdist) {
        check <- which(id2 > 0 & id == 0)
        if (length(check > 0)) {
            for (ch in 1:length(check)) {
			#ch<-1
                attrp <- id2[check[ch]]
			use<-id[which(id > 0)]
                dist1 <- st_distance(matpunti[check[ch],], matpunti[which(id > 0),])
			mini<-which.min(dist1)
                 if (attrp != use[mini]) {
                  id2[check[ch]] <- 0
                }
            }
        }
    }
    intersect <- matrix(NA, taxa, taxa)
    sympatry <- intersect
    for (k in 1:(taxa - 1)) {
        for (c in (k+1):taxa) {
#k<-1
#c<-2
            tryintersect <- try(st_intersection(hulls[[k]], hulls[[c]]), silent = TRUE)
            if (inherits(tryintersect, "try-error")) {
                intersect[k, c] <- 0
                sympatry[k, c] <- 0
                intersect[c, k] <- intersect[k, c]
                sympatry[c, k] <- sympatry[k, c]
            }
            else {
                if (!is.null(tryintersect)) {
		intersect[k, c] <- sum(st_area(st_make_valid(tryintersect)))
                  sympatry[k, c] <- intersect[k, c]/(areas[c] + 
                    areas[k] - intersect[k, c])
                  intersect[c, k] <- intersect[k, c]
                  sympatry[c, k] <- sympatry[k, c]
                }
                if (is.null(inter)) {
                  intersect[k, c] <- 0
                  sympatry[k, c] <- 0
                  intersect[c, k] <- intersect[k, c]
                  sympatry[c, k] <- sympatry[k, c]
                }
            }
        }
    }
    options(warn = oldw)
    res$areas <- areas
    res$intersections <- intersect
    res$sympatry <- sympatry
    res$NUR <- (length(which(id2 == 0))/length(which(id == 0))) * 
        100
    res$table <- cbind(mat, id2, id)
    res$hulls <- hulls
    res$hullpl <- hullpl
    return(res)
    if (plot) {
        points(matpunti, col = id2, cex = 0.5)
    }
}
