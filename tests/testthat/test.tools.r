test_that("rotation around the x axis works", {

    vec <- c(1, 3, 0)

    angles <- c(pi/2, 0, 0)

    res <- vec %>% rotate3D_lh(angles)
    res %>% expect_equal(c(1, 0, -3))

    res2 <- res %>% rotate3D_lh(angles)
    res2 %>% expect_equal(c(1, -3, 0))

    res3 <- res2 %>% rotate3D_lh(angles)
    res3 %>% expect_equal(c(1, 0, 3))

    res4 <- res3 %>% rotate3D_lh(angles)
    res4 %>% expect_equal(vec)

    res5 <- vec %>% rotate3D_lh(2*angles)
    res5 %>% expect_equal(res2)
})

test_that("rotation around y axis works", {

    vec <- c(1, 2, 3)
    angles <- c(0, pi/2, 0)

    res <- vec %>% rotate3D_lh(angles)
    res %>% expect_equal(c(-3, 2, 1))

    res2 <- res %>% rotate3D_lh(angles)
    res2 %>% expect_equal(c(-1, 2, -3))

    res3 <- res2 %>% rotate3D_lh(2*angles)
    res3 %>% expect_equal(vec)
})

test_that("rotation around z axis works", {

    vec <- c(1, 2, 3)
    angles <- c(0, 0, pi/2)

    res <- vec %>% rotate3D_lh(angles)
    res %>% expect_equal(c(2, -1, 3))

    res2 <- res %>% rotate3D_lh(angles)
    res2 %>% expect_equal(c(-1, -2, 3))

    res3 <- res2 %>% rotate3D_lh(2*angles)
    res3 %>% expect_equal(vec)
})

test_that("order of application is x-y-z", {

    set.seed(2)

    vec <- rnorm(3)

    angles <- runif(3, 0, 2*pi)

    expect_equal(vec %>% rotate3D_lh(angles),
                 vec %>%
                   rotate3D_lh(c(angles[1], 0, 0)) %>%
                   rotate3D_lh(c(0, angles[2], 0)) %>%
                   rotate3D_lh(c(0, 0, angles[3])))
})

test_that("negative angles work aswell", {

    set.seed(1)

    vec <- rnorm(3)
    angles <- runif(3, 0, 2*pi)

    expect_equal(vec,
                 vec %>%
                   rotate3D_lh(angles) %>%
                   rotate3D_lh(c(0, 0, -angles[3])) %>%
                   rotate3D_lh(c(0, -angles[2], 0)) %>%
                   rotate3D_lh(c(-angles[1], 0, 0)))
})

test_that("working with multiple vectors to rotate3D_lh at once works", {

    set.seed(1)
    
    vecs <- matrix(rnorm(90), nrow = 3)

    angles <- runif(3, 0, 2*pi)
    
    vecs.rotated <- vecs %>% rotate3D_lh(angles)
    vecs.rotated.one.by.one <- vecs %>% split(col(vecs)) %>%
        lapply(function(vec) vec %>% rotate3D_lh(angles))

    expect_equal(vecs.rotated %>% as.vector,
                 vecs.rotated.one.by.one %>%
                 do.call(cbind, .) %>%
                 as.vector)

    vecs.arr <- vecs %>% array(dim = c(3, 30))

    vecs.arr.rotated <- vecs.arr %>% rotate3D_lh(angles)
    vecs.arr.rotated.one.by.one <- vecs.arr %>% split(col(vecs.arr)) %>%
        lapply(function(vec) vec %>% rotate3D_lh(angles))
    
    expect_equal(vecs.arr.rotated %>% as.vector,
                 vecs.arr.rotated.one.by.one %>% do.call(cbind, .) %>%
                 as.vector)
})

test_that("Rotation does not change vector length", {

    set.seed(1)
    
    vecs <- matrix(10 * rnorm(90), nrow = 3)

    angles <- runif(3, 0, 2*pi)

    vecs.rotated <- vecs %>% rotate3D_lh(angles)
    expect_equal(vecs^2 %>% colSums %>% sqrt,
                 vecs.rotated^2 %>% colSums %>% sqrt)
})
