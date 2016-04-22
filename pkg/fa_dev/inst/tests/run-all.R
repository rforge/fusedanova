library(testthat)
library(fusedanova)

test_package("fusedanova")

# list of tests :
# - comparaison with cluster path OK
# - test predict with lambda in fa and after -> OK
# - test the cross validation ... -> to do
# - comp with flsa OK for no split + other when maxflow is coded -> to do (cannot compare unless force 0 and 1 weights in W)
# - comp with genlasso ? --> to do
