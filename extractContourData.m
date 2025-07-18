function contourData = extractContourData(C)
    contourData = struct('level', {}, 'x', {}, 'y', {});
    index = 1;
    while index < size(C, 2)
        level = C(1, index);
        nPoints = C(2, index);
        x = C(1, index+1:index+nPoints);
        y = C(2, index+1:index+nPoints);
        
        % Check if the contour meets the criteria
        duration = max(y) - min(y);
        x_range = [min(x), max(x)];
        
        if (duration >= 1/360*30 && duration <= 1/360*90) && ...
           (x_range(1) < 180 && x_range(2) > 130)
            contourData(end+1).level = level;
            contourData(end).x = x;
            contourData(end).y = y;
        end
        
        index = index + nPoints + 1;
    end
end