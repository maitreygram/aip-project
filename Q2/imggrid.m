function [ m n ] = imggrid( dm, dn, sz )
% imggrid -- Returns rectilinear coordinate vectors
%
% SYNOPSIS
%   [ m n ] = imggrid( dm, dn, sz )
%
% DESCRIPTION
%   Given the sample spacings and the image size, this
%   function returns the row and column coordinate vectors
%   for the image. Both vectors are centered about zero.
%
% INPUT
%   dm     Spacing between rows.
%   dn     Spacing between columns.
%   sz     2x1 vector of the image size: [ Nrows Ncols ].
%
% OUTPUT
%   m      sz(1) x 1 row coordinate vector.
%   n      1 x sz(2) column coordinate vector.

M = sz(1);
N = sz(2);

if( mod( M, 2 ) == 0 )
  m = dm * ( ceil( -M/2 ) : floor( M/2 ) - 1 )';
else
  m = dm * ( ceil( -M/2 ) : floor( M/2 ) )';
end

if( mod( N, 2 ) == 0 )
  n = dn * ( ceil( -N/2 ) : floor( N/2 ) - 1 );
else
  n = dn * ( ceil( -N/2 ) : floor( N/2 ) );
end