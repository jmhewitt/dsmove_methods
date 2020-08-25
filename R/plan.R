lapply(list.files("./R/subplans", full.names = TRUE, recursive = TRUE), source)

the_plan = bind_plans(
  simulation_plan
)