# Learning-To-Rank-

Implementation of Linear Regression Model Using Batch and Stochastic Models in Matlab
=====================================================================================
For real world dataset, I have used the Microsoft LETOR 4.0 Dataset. LETOR is a package of
benchmark data sets for research on Learning To Rank released by Microsoft Research Asia.

The latest version, 4.0, can be found at

http://research.microsoft.com/en-us/um/beijing/projects/letor/letor4dataset.aspx

It contains 8 datasets for four ranking settings derived from the two query sets and the Gov2 web
page collection. For this project, download MQ2007. There are three versions for each dataset:
\NULL", \MIN" and \QueryLevelNorm". In this project, only the \QueryLevelNorm" version
\Querylevelnorm.txt" will be used. The entire dataset consists of 69623 query-document
pairs(rows), each having 46 features. Here are two sample rows from the MQ2008 dataset.

The meaning of each column are as follows.

1. The first column is the relevance label of the for the corresponding feature values. It takes value 0, 1 or 2. This is the
objective output y we expect our linear regression to give.

2. The second column qid is the query id. which is not used in the model.

3. The following 46 columns are the features. They are the 46-dimensional input vector x
for our linear regression model. All the features are normalized to fall in the interval of
[0; 1].

OUTPUT:

1. The linear regression model gives the root mean square values trainPer1,trainPer2 and validPer1,validPer2 of the errors from the training set and the validation set.

2. The initial weights w01,w02[batch size X 1] for the stochastic regression model and the the weight differences dw1,dw2 [batch size X training set size] between each weight update.

3. The learning rate [1 X training set size] which stores updated values for each learning rate hyper parameter value in each iteration.
