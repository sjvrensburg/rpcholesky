# Using the rpcholesky C++ API

To use the rpcholesky C++ API in your package:

1. Add this to your DESCRIPTION file:
```
LinkingTo: Rcpp, RcppArmadillo, rpcholesky
```

2. Include the headers in your C++ code:
```cpp
// For using all exported functions
#include <rpcholesky_all.h>
```

3. Access the functions from the rpcholesky namespace:
```cpp
// For inline functions
arma::uvec samples = rpcholesky::sample_with_replacement_optimized(n, b, weights);

// For exported functions
Rcpp::List result = rpcholesky::accelerated_rpcholesky(A, k, b_in, stoptol, verbose);
```

See the examples directory for complete usage examples.
