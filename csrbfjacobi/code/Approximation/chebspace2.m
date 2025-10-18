function [ x ] = chebspace2(a,b,n )
%
  x = zeros ( n, 1 );

  if ( n == 1 )

    x(1) = ( a + b ) / 2.0;

  else

    for i = 1 : n

      theta = ( n - i ) * pi / ( n - 1 );

      c = cos ( theta );

      if ( mod ( n, 2 ) == 1 )
        if ( 2 * i - 1 == n )
          c = 0.0;
        end
      end

      x(i) = ( ( 1.0 - c ) * a  ...
             + ( 1.0 + c ) * b ) ...
             /   2.0;

    end

  end

  return
end
