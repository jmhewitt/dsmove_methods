context('GpsLik.cpp C++ methods')

test_that('Validating basic ll inputs and outputs', {
  
  set.seed(2021)
  
  n = 10
  
  obs_lons = runif(n = n, min = -180, max = 180)
  obs_lats = runif(n = n, min = -90, max = 90)
  
  semi_majors = rep(250, n)
  semi_minors = rep(250, n)
  
  orientations = runif(n = n, min = 0, max = 180)
  
  na_ind = 2
  obs_lons[na_ind] = NA
  obs_lats[na_ind] = NA
  semi_majors[na_ind] = NA
  semi_minors[na_ind] = NA
  orientations[na_ind] = NA
  
  thresh = .05
  
  # missing observations return flat likelihood
  expect_equal(
    GpsLikEval(obs_lons = obs_lons, obs_lats = obs_lats, 
               semi_majors = semi_majors, semi_minors = semi_minors, 
               orientations = orientations, alpha = thresh, 
               test_lon = runif(n = 1, min = -180, max = 180), 
               test_lat = runif(n = 1, min = -90, max = 90), 
               ind = na_ind - 1),
    0
  )
  
  # likelihood is infinite for distant observations
  expect_equal(
    GpsLikEval(obs_lons = obs_lons, obs_lats = obs_lats, 
               semi_majors = semi_majors, semi_minors = semi_minors, 
               orientations = orientations, alpha = thresh, 
               test_lon = obs_lons[1] + 1, 
               test_lat = obs_lats[1],
               ind = 0),
    -Inf
  )
  
  # likelihood is the integration constant for exact observations
  expect_equal(
    GpsLikEval(obs_lons = obs_lons, obs_lats = obs_lats, 
               semi_majors = semi_majors, semi_minors = semi_minors, 
               orientations = orientations, alpha = thresh, 
               test_lon = obs_lons[1], 
               test_lat = obs_lats[1],
               ind = 0),
    -12.1363584271863
  )
  
  # likelihood is slightly off from the integration constant for close points
  expect_equal(
    GpsLikEval(obs_lons = obs_lons, obs_lats = obs_lats, 
               semi_majors = semi_majors, semi_minors = semi_minors, 
               orientations = orientations, alpha = thresh, 
               test_lon = obs_lons[1] + .01, 
               test_lat = obs_lats[1],
               ind = 0),
    -12.2815062581422
  )
})