# This script finalizes the package after compilation
# Run this after successful compilation and installation

# Create directory if it doesn't exist
if (!dir.exists("tools/templates")) {
  dir.create("tools/templates", recursive = TRUE)
}

# Write a README with updated information
cat('# Using the rpcholesky C++ API

To use the rpcholesky C++ API in your package:

1. Add this to your DESCRIPTION file:
```
LinkingTo: Rcpp, RcppArmadillo, rpcholesky
```

2. Include the headers in your C++ code:
```cpp
// For using all exported functions
#include <rpcholesky.h>
```

3. Access the functions from the rpcholesky namespace:
```cpp
// For inline functions
arma::uvec samples = rpcholesky::sample_with_replacement_optimized(n, b, weights);

// For exported functions
Rcpp::List result = rpcholesky::accelerated_rpcholesky(A, k, b_in, stoptol, verbose);
```

See the examples directory for complete usage examples.
', file = "tools/templates/README_API.md")

# After a successful installation, copy this to inst/include/README.md
message("Package compilation and finalization complete.")
message("After successful installation, you may want to run:")
message('file.copy("tools/templates/README_API.md", "inst/include/README.md", overwrite = TRUE)')
