import numpy as np
from sklearn import svm
from sklearn import cross_validation
from sklearn import datasets

# Load data
concentrations = np.loadtxt(open("concentration_healthy_asthma.txt","rb"), delimiter=",", skiprows=0)

# Split data
X_train, X_test, y_train, y_test = cross_validation.train_test_split(concentrations[:,0:28],concentrations[:,29], test_size=0.4, random_state=0)

# Train a fancy SVM
clf = svm.SVC().fit(X_train, y_train)
score=clf.score(X_test, y_test)


loo = cross_validation.LeaveOneOjg,ut(140)
len(loo)