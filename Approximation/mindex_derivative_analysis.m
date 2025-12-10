function[a_structure, dimension_jumps] = mindex_derivative_analysis(a)
% [a_structure, dimension_jumps] = mindex_derivative_analysis(a)
%
% Analyzes entries in a and provides indexing structure to facilitate fast
% polynomial differentiation. The output matrix a_structure has arows =
% max(a(:))+1 rows, and each column is a set of row indices in a corresponding
% to a single dimension-deflated multi-index.
%
% a_structure(:,i) contains arows entries, such that a(a_structure(:,i), :)
% contains all multi-indices with (d-1) coordinates the same, and the d'th
% coordinate ordered from 0,...,arows-1. 
%
% When a does not contain enough high-degree terms, the entries of
% a_structure(:,i) have entry N+1, where N = size(a,1). Thus, for this i the
% indexing a(a_structure(:,i),:) will fail, and one must pad a with an (N+1)st
% row.
%
% The columns of a_structure corresponding to a dimension-q deflation are given
% by dimension_jumps(q):(dimension_jumps(q+1)-1).
%
% This whole procedure essentially rearranges a into a format that is fast for
% matlab analysis.  In a compiled language, something simpler can be used.

[N,d] = size(a);
Kmax = max(a(:));

a_structure = zeros([Kmax+1 0]);
dimension_jumps = [1];

% for each dimension, we analyze remaining multi-indices
for q = 1:d

  deflated_a = a(:,setdiff(1:d,q));

  [sort_a, sortinds] = sortrows(deflated_a);

  jumps = [1; find(any(abs(diff(sort_a, 1, 1)) > 0, 2))+1; N+1];

  % sort_a(jumps(i):(jump(i+1)-1)) are identical rows

  M = length(jumps) - 1;
  dimension_jumps = [dimension_jumps; M+dimension_jumps(end)];
  a_size = size(a_structure, 2);
  a_structure = [a_structure zeros([Kmax+1 M])];

  for qq = 1:M

    % We don't check to make sure the degrees are distinct, but they should be
    qq_inds = sortinds(jumps(qq):(jumps(qq+1)-1));
    padding = Kmax+1 - numel(qq_inds);
    degrees = a(qq_inds,q);

    [degsort,deginds] = sort(degrees);
    if any(diff(degsort) > 1)
      warning('The input index set is not a lower set. Differentiation will not be accurate.');
    end
    a_structure(:,a_size+qq) = [sortinds(jumps(qq):(jumps(qq+1)-1)); (N+1)*ones([padding 1])];

  end

end
