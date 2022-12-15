def TrimmedMean(UpCoverage):
	nonzero_count = 0
	for i in range(1,50):
		if (UpCoverage[-i]>0):
			nonzero_count += 1

	total = 0
	count=0
	for cov in UpCoverage:
		if(cov>0):
			total += cov
			count += 1
	trimMean = 0
	if nonzero_count > 0 and count >20:
	#if count >20:
		trimMean = total/count;
	return trimMean
