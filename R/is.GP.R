## Checks whether or not the object in question is a GP object
##
## May 8th, 2012

is.GP <- function(object) {
    inherits(x = object, what = "GP") && 
        # also trivial checks for validity
        !is.null(names(object)) &&
        all(c("X", "Y", "sig2", "beta", "delta", "nugget_threshold_parameter", 
            "correlation_param") %in% names(object))
}
