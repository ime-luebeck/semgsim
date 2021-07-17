
expect_same_length <- function(a, b)
    expect_equal(length(a), length(b))

is.positive <- function(x)
    x > 0

is.negative <- function(x)
    x < 0
