function [B] = unique_no_sort_rows(A)

%UNIQUE_NO_SORT_ROWS Set unique rows unsorted.
%   B = UNIQUE_NO_SORT_ROWS(A) for the array A returns a vector of the unique 
%   rows of A in the order that they appear in A, i.e. B is unsorted.  
%
%   Michael Thomas Petralia
%   Harvard: August 31, 2009


% indexing the first occurence of each unique row in A
[unique_sorted,i,j] = unique(A,'rows','first');
% adding the index to the matrix of unique sorted rows
indexed = [i unique_sorted];
% sorting the rows based on the index (i.e. unsorting the rows)
unsorted = sortrows(indexed,1);
% removing the index column
B = unsorted(:,2:end);


