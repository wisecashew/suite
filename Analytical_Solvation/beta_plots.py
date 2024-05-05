from scipy.stats import beta
import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
	a, b, loc, scale = 8.382758059994785, 1435239.3443798767, -2.5967567335988404, 5720842.862537522
	x = np.linspace(0, 200, 1000)
	rv = beta(a, b, loc, scale)
	pdf = rv.pdf(x)
	for i in range(len(x)):
		print(f"x: {x[i]}, pd: {pdf[i]}")
	