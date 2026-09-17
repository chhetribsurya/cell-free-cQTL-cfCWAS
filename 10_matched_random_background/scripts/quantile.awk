# Empirical quantile for a pre-sorted numeric stream.
# Usage: sort -n file | awk -v p=0.025 -f scripts/quantile.awk
# Uses R type-7 interpolation: h = 1 + (n-1)*p
{
	a[NR] = $1
}
END {
	if (NR == 0) {
		print "NA"
		exit
	}
	if (p < 0) p = 0
	if (p > 1) p = 1
	h = 1 + (NR - 1) * p
	i = int(h)
	f = h - i
	if (i < 1) {
		print a[1]
		exit
	}
	if (i >= NR) {
		print a[NR]
		exit
	}
	print a[i] + f * (a[i + 1] - a[i])
}
