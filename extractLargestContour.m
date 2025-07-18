function [max_contour_x, max_contour_y] = extractLargestContour(C)
    index = 1;
    max_points = 0;
    max_contour_x = [];
    max_contour_y = [];
    
    while index < size(C, 2)
        nPoints = C(2, index);
        
        if nPoints > max_points
            max_points = nPoints;
            max_contour_x = C(1, index+1:index+nPoints);
            max_contour_y = C(2, index+1:index+nPoints);
        end
        
        index = index + nPoints + 1;
    end
end