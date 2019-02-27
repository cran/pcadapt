/******************************************************************************/

#include <pcadapt/bed-acc.h>

/******************************************************************************/

inline size_t ceil4(size_t n) {
  return (n + 3) / 4;
}

/******************************************************************************/

bed::bed(std::string path, int n, int p) : 
  n(n), p(p), n_byte(ceil4(n)) {
  
  std::error_code error;
  mio::ummap_source magic_number;
  magic_number.map(path, 0, 3, error);
  if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());
  
  if (!(magic_number[0] == '\x6C' && magic_number[1] == '\x1B')) 
    Rcpp::stop("File is not a binary PED file.");
  
  /* Check mode: 00000001 indicates the default variant-major mode (i.e.
   list all samples for first variant, all samples for second variant,
   etc), 00000000 indicates the unsupported sample-major mode (i.e. list
   all variants for the first sample, list all variants for the second
   sample, etc */
  if (magic_number[2] != '\x01') Rcpp::stop("Sample-major mode is not supported.");
  
  // Map after this magic number
  this->ro_ummap.map(path, 3, mio::map_entire_file, error);
  if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());
  
  // Check if given dimensions match the file
  if ((this->n_byte * this->p) != this->ro_ummap.size())
    Rcpp::stop("n or p does not match the dimensions of the file.");
}

/******************************************************************************/

// [[Rcpp::export]]
SEXP bedXPtr(std::string path, int n, int p) {
  
  // http://gallery.rcpp.org/articles/intro-to-exceptions/
  try {
    /* Create a pointer to a bed object and wrap it as an external pointer
     http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf */
    XPtr<bed> ptr(new bed(path, n, p), true);
    // Return the external pointer to the R side
    return ptr;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
}

/******************************************************************************/
