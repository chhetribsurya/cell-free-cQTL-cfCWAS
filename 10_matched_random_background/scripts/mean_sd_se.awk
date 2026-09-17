BEGIN {
	n = 0
	s = 0
	ss = 0
}
{
	n++
	s += $1
	ss += $1 * $1
}
END {
	if (n > 0) {
		mean = s / n
		if (n > 1) {
			var = (ss - s * s / n) / (n - 1)
			if (var < 0) var = 0
			sd = sqrt(var)
			se = sd / sqrt(n)
		} else {
			sd = 0
			se = 0
		}
		printf "%s %s %s %s\n", mean, sd, se, n
	} else {
		print "NA NA NA 0"
	}
}
