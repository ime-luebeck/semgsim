context("Integration test")

## Actual test
test_that("Integration works", {
  if (!skip_integration) {
    timing_sequential <- system.time({
    
        res_sequential <-
            simulate_muscle("../../config/integration.r",
                          deterministic = TRUE,
                          num_cores = 1)
    })
    
    
    timing_parallel <- system.time({
        res_parallel <- simulate_muscle("../../config/integration.r",
                                      deterministic = TRUE,
                                      num_cores = 4)
    })
    
    
    potential_seq <-
        res_sequential$surface_potentials %>% transform(potential_seq = potential)
    potential_par <-
        res_parallel$surface_potentials %>% transform(potential_par = potential)
    
    results_combined <-
        potential_seq[, c("electrode", "time", "potential_seq")] %>%
              merge(potential_par[, c("electrode", "time", "potential_par")]) %>%
                  dplyr::mutate(potentials_are_equal = (potential_seq == potential_par))
    
    expect_true(all(results_combined$potentials_are_equal))
    
    
    ## Teardown
    
    if (!file.exists("integration_test_timings.txt"))
        file.create("integration_test_timings.txt")
    
    file_connection <- file("integration_test_timings.txt")
    writeLines(c("Sequential timings",
                 paste(timing_sequential),
                 "\n",
                 "Parallel timings",
                 paste(timing_parallel)),
               file_connection)
    close(file_connection)
  } else
    skip("Skipping integration test as requested by skip_integration=TRUE.")
})