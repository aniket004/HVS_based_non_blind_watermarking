function [Y,Cu,Cv] = rgb2ycucv(R,G,B)

      [ nr nc ] = size( R );
      
      rgb(:,:,1) = double(R);
      rgb(:,:,2) = double(G);
      rgb(:,:,3) = double(B);
      
    T = [1/4    1/2     1/4;
          0     -1       1 ;
          1     -1       0];

    ycucv = zeros(nr,nc,3);
    
    for i = 1:nr
        for j = 1:nc
            x = zeros(3,1);
            for k = 1:3
                x(k,1) = rgb(i,j,k);
            end
            y = ( T * x );
            %y = round( 1000 * T * x );        
            for k = 1:3
                ycucv(i,j,k) = double( y(k,1) );
            end
        end
    end

    Y = ycucv(:,:,1);
    Cu = ycucv(:,:,2);
    Cv = ycucv(:,:,3);
    

end
