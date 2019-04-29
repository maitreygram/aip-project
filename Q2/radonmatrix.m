function [ R rho theta ] = radonmatrix( drho, dtheta, M, N )
% radonmatrix - Discrete Radon Trasnform matrix
%
% SYNOPSIS
%   [ R rho theta ] = radonmatrix( drho, dtheta, M, N )
%
% DESCRIPTION
%   Returns a matrix representation of a Discrete Radon
%   Transform (DRT).
%
% INPUT
%   drho     Radial spacing the the DRT.
%   dtheta   Angular spacing of the DRT (rad).
%   M        Number of rows in the image.
%   N        Number of columns in the image.
%
% OUTPUT
%   R        LP x MN DRT matrix. The values of the L and
%            P will depend on the radial and angular spacings.
%   rho      Vector of radial sample locations.
%   theta    Vector of angular sample locations (rad).
%

% For each angle, we define a set of rays parameterized
% by rho. We then find the pixels on the MxN grid that
% are closest to each line. The elements in R corresponding
% to those pixels are given the value of 1.

% The maximum extent of the region of support. It's for
% rho = 0 and theta = pi/4, the line that runs caddy-corner.
W = sqrt( M^2 + N^2 );

rho = -W/2 : drho : W/2;
theta = 0 : dtheta : 180 - dtheta;

L = length( rho );
P = length( theta );

R = false( L*P, M*N );

% Define a meshgrid w/ (0,0) in the middle that
% we can use a standard coordinate system.
[ mimg nimg ] = imggrid( 1, 1, [ M N ] );

% We loop over each angle and define all of the lines.
% We then just figure out which indices each line goes
% through and put a 1 there.
for ii = 1 : P

  phi = theta(ii) * pi/180;

  % The equaiton is rho = m * sin(phi) + n * cos(phi).
  % We either define a vector for m and solve for n
  % or vice versa. We chose which one based on angle
  % so that we never g4et close to dividing by zero.
  if( phi >= pi/4 && phi <= 3*pi/4 )

    t =  -W : min( 1/sqrt(2), 1/abs(cot(phi)) ) : +W;
    T = length( t );
    

    rhom = repmat( rho(:), 1, T );
    tn = repmat( t(:)', L, 1 );
    mline = ( rhom - tn * cos(phi) ) ./ sin(phi);

    for jj = 1 : L
      p = round( tn(jj,:) - min( nimg ) ) + 1;
      q = round( mline(jj,:) - min( mimg ) ) + 1;  
      inds = p >= 1 & p <= N & q >= 1 & q <= M;
      R( (ii-1)*L + jj, unique( sub2ind( [ M N ], q(inds), p(inds) ) ) ) = 1;
    end

  else

    t =  -W : min( 1/sqrt(2), 1/abs(tan(phi)) ) : +W;
    T = length( t );

    rhon = repmat( rho(:)', T, 1 );    
    tm = repmat( t(:), 1, L );
    nline = ( rhon - tm * sin(phi) ) ./ cos(phi);

    for jj = 1 : L
      p = round( nline(:,jj) - min( nimg ) ) + 1;
      q = round( tm(:,jj) - min( mimg ) ) + 1;  
      inds = p >= 1 & p <= N & q >= 1 & q <= M;
      R( (ii-1)*L + jj, unique( sub2ind( [ M N ], q(inds), p(inds) ) ) ) = 1;
    end

  end

end

R = double( sparse( R ) );

return;