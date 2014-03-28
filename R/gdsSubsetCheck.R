gdsSubsetCheck <- function(parent.gds,
                           sub.gds,
                           sample.include=NULL,
                           snp.include=NULL,
                           sub.storage=NULL,
                           verbose=TRUE) {
  
  # this assumes that sample.id is the only 1D sample variable in the GDS
  
  gds <- openfn.gds(parent.gds)
  gds.sub <- openfn.gds(sub.gds)
  
  # check sampleID
  sampID.parent <- read.gdsn(index.gdsn(gds, "sample.id"))
  sampID.sub <- read.gdsn(index.gdsn(gds.sub, "sample.id"))
  chk <- all(is.element(sampID.sub, sampID.parent))
  if (!chk) {
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("sample.id in sub GDS is not a subset of parent GDS")
  }
  if (is.null(sample.include)){
    sample.include <- sampID.parent
  }
  chk <- setequal(sample.include, sampID.sub)
  if (!chk) {
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("samples in sub GDS are not the same as sample.include")
  }
  
  # check snp variables
  snpID.parent <- read.gdsn(index.gdsn(gds, "snp.id"))
  snpID.sub <- read.gdsn(index.gdsn(gds.sub, "snp.id"))
  chk <- all(is.element(snpID.sub, snpID.parent))
  if (!chk){
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("snp.id in sub GDS is not a subset of parent GDS")
  }
  if (is.null(snp.include)){
    snp.include <- snpID.parent
  }
  chk <- setequal(snp.include, snpID.sub)
  if (!chk) {
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("snps in sub GDS are not the same as snp.include")
  }
  snpsel <- snpID.parent %in% snp.include
  
  variables <- c("snp.id", "snp.position", "snp.chromosome")
  chk <- all(variables %in% ls.gdsn(gds))
  if (!chk) {
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("parent GDS does not have snp.id, snp.position, and snp.chromosome")
  }
  chk <- all(variables %in% ls.gdsn(gds.sub))
  if (!chk) {
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("sub GDS does not have snp.id, snp.position, and snp.chromosome")
  }
  
  for (variable in variables){
    node.parent <- index.gdsn(gds, variable)
    node.sub <- index.gdsn(gds.sub, variable)
    
    vals.parent <- read.gdsn(node.parent)
    vals.sub <- read.gdsn(node.sub)
    
    chk <- allequal(vals.parent[snpsel], vals.sub)
    if (!chk) {
      closefn.gds(gds)
      closefn.gds(gds.sub)
      stop("snps in sub GDS are not the same as snp.include")
    }
    
    # TO DO: check attributes here.
    if (variable == "snp.chromosome"){
      attributes.parent <- get.attr.gdsn(node.parent)
      attributes.sub <- get.attr.gdsn(node.sub)
      chk <- setequal(names(attributes.parent), names(attributes.sub))
      if (!chk) {
        closefn.gds(gds)
        closefn.gds(gds.sub)
        stop(paste("sub GDS has different attributes than parent GDS for", variable))
      }
      
      for (attribute in names(attributes.parent)){
        chk <- allequal(attributes.parent[[attribute]], attributes.sub[[attribute]])
        if (!chk) {
          closefn.gds(gds)
          closefn.gds(gds.sub)
          stop(paste("sub GDS has different attribute values than parent GDS for", variable))
        }
      }
    }
    
  }
  
  # check other variables
  node.names <- ls.gdsn(gds)
  node.names <- node.names[!(node.names %in% c("sample.id", "snp.id", "snp.position", "snp.chromosome", "snp.allele", "snp.rs.id"))]
  # check that they exist in the sub file
  chk <- all(variables %in% ls.gdsn(gds.sub))
  if (!chk) {
    closefn.gds(gds)
    closefn.gds(gds.sub)
    stop("sub GDS does not include the same variables as parent GDS")
  }
  
  # compare the variables node by node
  for (node.name in node.names){
    
    # node in parent file
    node.parent <- index.gdsn(gds, node.name)
    
    # check dimensions
    desc <- objdesp.gdsn(node.parent)
    
    # only check on 2-d variables
    if (length(desc$dim) != 2) next
    if (verbose) message(paste("working on", node.name))
    
    # get the subset node
    node.sub <- index.gdsn(gds.sub, node.name)
    
    # check attributes
    attributes.parent <- get.attr.gdsn(node.parent)
    attributes.sub <- get.attr.gdsn(node.sub)
    
    chk <- setequal(names(attributes.parent), names(attributes.sub))
    if (!chk) {
      closefn.gds(gds)
      closefn.gds(gds.sub)
      stop(paste("sub GDS has different attributes than parent GDS for", variable))
    }
    # check attribute values
    for (attribute in names(attributes.parent)[]) {
      # these don't have values
      if (attribute %in% c("snp.order", "sample.order", "missing.value")) next
      if (attribute == "storage" & !is.null(sub.storage)) {
        chk(attributes.sub[[attribute]] == sub.storage)
        if (!chk) {
          closefn.gds(gds)
          closefn.gds(gds.sub)
          stop(paste("sub GDS has incorrect attribute for", attribute))
        }
      } else{
        chk <- allequal(attributes.parent[[attribute]], attributes.sub[[attribute]])
        if (!chk) {
          closefn.gds(gds)
          closefn.gds(gds.sub)
          stop(paste("sub GDS has different attribute values than parent GDS for", variable))
        }
      }
    }
    
    if ("missing.value" %in% names(attributes.parent)) {
      parent.miss <- attributes.parent[["missing.value"]]
    } else {
      parent.miss <- NA
    }
    
    if ("missing.value" %in% names(attributes.sub)) {
      sub.miss <- attributes.sub[["missing.value"]]
    } else {
      sub.miss <- NA
    }
    
    # check data sample by sample
    for (i in 1:length(sampID.sub)){
      
      sample <- sampID.sub[i]
      i.sub <- which(sampID.sub %in% sample)
      i.parent <- which(sampID.parent %in% sample)
      
      if ("sample.order" %in% names(attributes.sub)) dimType <- "scan,snp" else dimType <- "snp,scan"
      
      if (dimType == "scan,snp"){
        start.parent <- c(i.parent, 1)
        start.sub <- c(i.sub, 1)
        count <- c(1, -1)
      } else {
        start.parent <- c(1, i.parent)
        start.sub <- c(1, i.sub)
        count <- c(-1, 1)
      }
      
      vals.parent <- read.gdsn(node.parent, start=start.parent, count=count)[snpsel]
      vals.sub <- read.gdsn(node.sub, start=start.sub, count=count)
      
      # set missing values
      vals.parent[which(vals.parent == parent.miss)] <- NA
      vals.sub[which(vals.sub == sub.miss)] <- NA
      
      if (length(vals.parent) != length(vals.sub)) {
        closefn.gds(gds)
        closefn.gds(gds.sub)
        stop(paste("lengths of variable", node.name, "are not the same"))
      }
      
      chk <- allequal(vals.parent, vals.sub)
      if (!chk) {
        closefn.gds(gds)
        closefn.gds(gds.sub)
        stop(paste("values of variable", node.name, "are not the same."))
      }
    }
  }
  
  closefn.gds(gds)
  closefn.gds(gds.sub)
  message("All variables match.")
}