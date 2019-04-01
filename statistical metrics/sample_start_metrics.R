weights_vector <- c(scan(file.choose()))
coordinates_vector <- matrix(scan(file.choose()), ncol=2, byrow=TRUE)
print('weighted mean')
w.mean(coordinates_vector, weights_vector)
w.var(coordinates_vector, weights_vector)
w.ad(coordinates_vector, weights_vector)


