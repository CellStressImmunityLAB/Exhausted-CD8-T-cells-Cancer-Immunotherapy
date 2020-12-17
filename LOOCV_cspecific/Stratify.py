import numpy as np
import itertools

# ---------------------------------------------------------------------------------------------------
# -------------------------------------- Function Definitions ---------------------------------------
# ---------------------------------------------------------------------------------------------------


def stratify(n_splits, response, nframe):
	# As input, takes a matrix with observations in rows, columns as features
	# One feature must be listed as the response variable (e.g. LOG_IC50)
	mygenerator = []
	n_samples = len(nframe)
	mysorted = nframe.sort_values(response, axis=0, ascending=True, inplace=False, kind='mergesort')

	# Get fold lengths
	n_splits = n_splits
	fold_sizes = (n_samples // n_splits) * np.ones(n_splits, dtype=np.int)  # Equal partitioning
	fold_sizes[:n_samples % n_splits] += 1  # increase by one until all remaining data samples are filled in

	lister = []
	for fold_size in fold_sizes:  # Create the folds
		#  with zeroes
		lister.append([0] * fold_size)

	for arun in range(0, min(fold_sizes)):
		bucket = np.random.choice(range(0, n_splits), n_splits, replace=False)

		for element in range(len(bucket)):
			lister[element][arun] = bucket[element] + arun * n_splits

	# Always add to the end, starting from the first
	# Is always the last item anyway or min fold size would just go higher
	mydifference = len(mysorted) - n_splits * min(fold_sizes)

	if mydifference > 0:
		bucket = np.random.choice(range(0, mydifference), mydifference, replace=False)
		for something in range(len(bucket)):
			lister[something][-1] = bucket[something] + min(fold_sizes) * n_splits
	else:
		pass

	listerset = set(tuple(sublist) for sublist in lister)

	for thingy in itertools.combinations(lister, n_splits - 1):  # Itertools will shuffle! (first no longer longest)
		setthingy = set(
			tuple(subsection) for subsection in thingy)  # We need a set for the most efficient difference call
		thingy = sorted([item for subl in thingy for item in subl])  # Flatten & sort training set
		# Set difference (test) + flatten
		diff = sorted([myitem for subset in listerset.difference(setthingy) for myitem in subset])
		mygenerator.append((thingy, diff))  # Create generator object

	return mygenerator
