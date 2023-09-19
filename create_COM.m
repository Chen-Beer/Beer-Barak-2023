function com = create_COM(event,electrodes,electrode_names_ordered,electrodes_names)
    E = length(electrodes);
    T = size(event,2);
    a = 10; 
    b = 12;

    labeled_im = logical(true([a b]));

    com = zeros(2,T);
    event_2D = zeros(a,b,T);
    for e = 1:E
        ei = electrodes(e);
        [r,c] = ind2sub([a b],find(contains(electrode_names_ordered,electrodes_names(ei))));
        event_2D(r,c,:) = squeeze(event(e,:));
    end
    for t = 1:T
        tmp = squeeze(event_2D(:,:,t));
        props = regionprops(labeled_im, tmp, 'Centroid', 'WeightedCentroid');
        com(:,t) = props.WeightedCentroid;
    end
    
end
