# CurveSketching Package

This is a simple and buggy SageMath package for Curve Sketching done as a project of the course [MAT007 Introduction to computer algebra and mathematical programming] lectured by Johannes Schmitt at the Institute of Mathematics, University of Zurich.
It contains four methods to get the roots of a function (using built-in Sage 9.2 solve(), roots(), find_root() and a sloppy self implemented Newton's method). Use the roots() method over the solve() method for big polynomials, find_root is fast as well but can't handle many functions.