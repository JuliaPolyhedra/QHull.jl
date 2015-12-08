using Base.Test
import CHull

pts = [-1.0 0;
       1 1;
       3 0;
       3 3;
       2 2;
       0 3;
       -1 2]

hull = CHull.chull(pts)
@assert hull.vertices == [1, 3, 4, 6, 7]
@assert size(hull.points) == size(pts)
                        
