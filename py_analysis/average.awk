#!/bin/awk -f

BEGIN {
    sum = 0
    count = 0
    rcount = 100
}

NR>=rcount {
		x = $2 + 0
		sum += x
		sumsq += x*x
		count++
}

END {
	if (count > 0) {
		average  = (sum / count)
		variance = (sumsq / count) - (mean*mean)
		print "Average  is " average
		print "Variance is " variance
		} 
		else {
		print "No entries to process from the "rcount"th entry onwards."
	}
}

