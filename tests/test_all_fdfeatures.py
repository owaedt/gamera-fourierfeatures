
from gamera.core import *
from gamera.toolkits.fourierfeatures import single
from gamera.toolkits.fourierfeatures import broken
init_gamera()
import math

import sys
def euclid(a,b):
	sum = 0.0
	for i in range(len(a)):
		sum += (a[i] - b[i]) * (a[i] - b[i])
	return math.sqrt(sum)


def calculate_fdfeatures(imgfile, prefix="fd"):
	img = load_image(imgfile)
	img = img .to_onebit()

#	gtresults = eval(open("gtresults.py", "r").read())
	results = {}
	for f in img.get_feature_functions()[0]:
		#print f
		if f[0].startswith(prefix):
			ff = f[1]
			results[f[0]] = list(ff()(img))

#			print gtresults[f[0]] == results[f[0]]

	return results


def writeToFile(features, filename):
	f = open(filename, "w")
	f.write(str(features))
	f.close()

def readFromFile(filename):
	f = open(filename, "r")
	features = eval(f.read())
	f.close()
	return features


def compare_features(gt, test, testtogtmap={}):

	correct = 0
	for f, v in test.items():
		gtf = f
		if f in testtogtmap:
			gtf = testtogtmap[f]

		if gtf is None:
			print("skipping")
			correct += 1
		elif gtf not in gt:
			print("GT-Data for %s missing" % f)
		else:
			ec = euclid(gt[gtf], test[f])
			if abs(ec) < 1e-10:
				correct += 1

				print(f, gtf, ec, "passed")
			else:
				print(gt[gtf])
				print(test[f])
				print(f, gtf, ec, "======= FAILEd =======")

	print("checked %d items and got %d correct." % (len(test), correct))
	assert(correct == len(test))
	return correct == len(test)

testtogtmap = {"fdsingle_complex_position": "fdsingle_abs_maxnorm",
		"fdsingle_complex_position_r1": "fdsingle_abs",
		"fdsingle_complex_position_phase_r1": "fdsingle_dimov",
		"fdsingle_complex_position_phase": "fdsingle_dimov_maxnorm",
		"fdsingle_real_position": "fdsingle_shridhar_maxnorm",
		"fdsingle_centroid_distance": "fdsingle_centroid_distance_abs",
		"fdsingle_centroid_distance_phase": "fdsingle_centroid_distance_maxnorm",
		"fdbroken_a": "fdbroken_a_abs_maxnorm",
		"fdbroken_a_phase": "fdbroken_a_maxnorm",
		"fdbroken_a_phase_s1": "fdbroken_a",
		"fdbroken_c_phase": "fdbroken_c_maxnorm",
		"fdbroken_c_phase_s1": "fdbroken_c",
		"fdbroken_c": "fdbroken_c_abs_maxnorm",
		"fourierDescriptorBrokenShape": "fdbroken_a_abs_maxnorm",
			"fdsingle_dimov_original": None}


gtfeatures = readFromFile("./testdata-z-08.dat")

def test_fdsingle():
	testfeatures = calculate_fdfeatures("./z-08.png", "fdsingle")

	res = compare_features(gtfeatures, testfeatures, testtogtmap)
	assert(res == True)

def test_fdbroken():
	testfeatures = calculate_fdfeatures("./z-08.png", "fdbroken")

	res = compare_features(gtfeatures, testfeatures, testtogtmap)
	assert(res == True)


