function CannyEdgeDetector(file)

    close all;
    sigma = 2; % Gaussian filter sigma 1.4
    high_threshold_ratio = 0.2; % High threshold ratio
    low_threshold_ratio = 0.15; % Low threshold ratio
    
    slice = imread(file);
    [h,w] = size(slice);
    
    % Gaussian smoothing
    fslice = imgaussfilt(im2double(slice), sigma);
    % figure; imshow(fslice, [min(fslice(:)),max(fslice(:))]); title('Gaussian smoothing');
    
    % % Roberts mask
%     kx = [-1, 1];
%     ky = [-1; 1];
    
    % Prewitt mask/cross
    kx = [-1 0 +1; -1 0 +1; -1 0 +1];
    ky = [-1 -1 -1; 0 0 0; +1 +1 +1];

    % Horizontal gradient
    gx = conv2(fslice, kx, 'same');
    % figure; imshow(gx,[]); title('Horizontal gradient');
    
    % Vertical gradient
    gy = conv2(fslice, ky, 'same');
    % figure; imshow(gy,[]); title('Vertical gradient');
    
    % Magnitude
    mag = sqrt(gx.^2 + gy.^2);
    % figure; imshow(mag, [min(mag(:)),max(mag(:))]); title('Magnitude');
    
    % Angle
    angle = atan2d(gy, gx);
    % figure; imshow(angle, [min(angle(:)),max(angle(:))]); title('Angle');
    
    % N1: Thin the ridges using non-maximum suppression (using magnitude and angle)
    non_maximum_supression = zeros(h,w);
    % row
    for i = 2 : h - 1
        % column
        for j = 2 : w - 1
            % horizontal edge
            if (-22.5 <= angle(i,j) && angle(i,j) <= 22.5) || (157.5 <= angle(i,j) || angle(i,j) <= -157.5)
                if mag(i,j) < mag(i,j+1) || mag(i,j) < mag(i,j-1)
                    non_maximum_supression(i,j)= 0;
                else
                    non_maximum_supression(i,j)= mag(i,j);
                end
            % -45 edge
            elseif (22.5 <= angle(i,j) && angle(i,j) <= 67.5) || (-157.5 <= angle(i,j) && angle(i,j) <= -112.5)
                if mag(i,j) < mag(i+1,j+1) || mag(i,j) < mag(i-1,j-1)
                    non_maximum_supression(i,j)= 0;
                else
                    non_maximum_supression(i,j)= mag(i,j);
                end
            % vertical edge
            elseif (67.5 <= angle(i,j) && angle(i,j) <= 112.5) || (-112.5 <= angle(i,j) && angle(i,j) <= -67.5)
                if mag(i,j) < mag(i+1,j) || mag(i,j) < mag(i-1,j)
                    non_maximum_supression(i,j)= 0;
                else
                    non_maximum_supression(i,j)= mag(i,j);
                end
            % 45 edge
            elseif (112.5 <= angle(i,j) && angle(i,j) <= 157.5) || (-67.5 <= angle(i,j) && angle(i,j) <= -22.5)
                if mag(i,j) < mag(i+1,j-1) || mag(i,j) < mag(i-1,j+1)
                    non_maximum_supression(i,j)= 0;
                else
                    non_maximum_supression(i,j)= mag(i,j);
                end
            end
        end
    end
    
    % figure; imshow(non_maximum_supression, [min(non_maximum_supression(:)),max(non_maximum_supression(:))]); title('Non Maximum Suppression');
    
    % N2: Hysteresis thresholding (separation into strong and weak edge pixels)

    high_threshold = max(max(non_maximum_supression)) * high_threshold_ratio;
    low_threshold = high_threshold * low_threshold_ratio;

    strong_edges = zeros(2, h);
    weak_edges = zeros(2, w);

    strong_index = 1;
    weak_index = 1;

    % row
    for i = 2 : h - 1
        % col
        for j = 2 : w - 1
            % Strong edge
            if non_maximum_supression(i, j) >= high_threshold
                non_maximum_supression(i, j) = 1;
                strong_edges(strong_index, 1) = i;
                strong_edges(strong_index, 2) = j;
                strong_index = strong_index + 1;
            % Weak edge
            elseif non_maximum_supression(i, j) >= low_threshold
                weak_edges(weak_index, 1) = i;
                weak_edges(weak_index, 2) = j;
                weak_index = weak_index + 1;
            % No edge
            else
                non_maximum_supression(i, j) = 0;
            end
        end
    end

    figure; imshow(non_maximum_supression, [min(non_maximum_supression(:)),max(non_maximum_supression(:))]); title('Hysteresis thresholding');
    
    % N3: Form longer edges (edge-linking w/ 8-connectivity of weak pixels to strong pixels)
end
