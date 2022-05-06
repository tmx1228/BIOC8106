
# This function counts Kmers in a given string
#
# @param str A character vector
# @param k The k, or length of the kmer, you wish to count
# @returns This function should return a table of kmer counts,
# 		   indexed by kmer name


countKmers = function(str, k) {
  kmers=c()
  for(i in 1:(nchar(str)-k+1)){
    iter <- substring(str,i,i+k-1)
    if(iter %in% names(kmers)){
      kmers[iter] = kmers[iter] + 1
    }else{
      kmers[iter] = 1
    }
  }
  return(kmers)
}

