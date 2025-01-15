# Functions for calculating and plotting R data


calculateVals <- function(vct.a, vct.b, gene_len) {
	# get sample_a solo values
	out.a <- c(vct.a, get_stall_ratio(vct.a), get_pause_ratio(vct.a,gene_len), get_unloading_ratio(vct.a,gene_len))
	
	# get sample_b solo values
	out.b <- c(vct.b, get_stall_ratio(vct.b), get_pause_ratio(vct.b,gene_len), get_unloading_ratio(vct.b,gene_len))
	
	# get comparison l2fc values
	j <- 1
	comp.vector <- c()
	while(j <= length(out.a)){
		cur.l2fc <- log2(out.b[k] + 0.1) - log2(out.a[j] + 0.1)
		comp.vector <- append(comp.vector, cur.l2fc)
		j <- j + 1
	}
	return(c(out.a, out.b, comp.vector))
}

get_stall_ratio <- function(x) {
  ret_val = 0
  prom_len <- .GlobalEnv$tssr_start - .GlobalEnv$pro_start + 1
  tssr_len <- .GlobalEnv$gene_start - .GlobalEnv$tssr_start + 1
  if (x[2] != 0) {
    ret_val <- ((x[1]+0.1)/prom_len) / ((x[2]+0.1)/tssr_len)
  }
  return(ret_val)
}

get_pause_ratio <- function(x, glen) {
  ret_val = 0
  tssr_len <- .GlobalEnv$gene_start - .GlobalEnv$tssr_start + 1
  if (x[3] != 0) {
    ret_val <- ((x[2]+0.1)/tssr_len) / ((x[3]+0.1)/(glen+1))
  }
  return(ret_val)
}

get_unloading_ratio <- function(x, glen) {
  ret_val = 0
  tesr_len <- .GlobalEnv$tesr_end
  if (x[4] != 0) {
    ret_val <- (x[3]/(glen+ 1)) / (x[4]/tesr_len)
  }
  return(ret_val)
}