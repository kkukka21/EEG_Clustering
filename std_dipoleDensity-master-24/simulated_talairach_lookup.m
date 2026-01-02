function labels = simulated_talairach_lookup(talCoords, radius)
    % SIMULATED_TALAIRACH_LOOKUP Simulates a Talairach atlas lookup for probabilistic labels
    % Inputs:
    %   talCoords - Talairach coordinates [x, y, z] (e.g., [-57.8 -66.2 8.6])
    %   radius - Search radius in mm (e.g., confusion sphere diameter / 2)
    % Output:
    %   labels - Struct with anatomical labels and probabilities

    labels = struct();
    x = talCoords(1);
    y = talCoords(2);
    z = talCoords(3);

    % Adjust probabilities based on radius (larger radius increases uncertainty)
    % For simplicity, we'll scale probabilities inversely with radius
    probScale = max(0.5, 1 - radius / 20); % Radius > 20mm reduces probabilities

    % Determine hemisphere
    if x < 0 % Left hemisphere
        if y < -50 % Posterior regions
            if z < -20 % Inferior (e.g., cerebellum)
                labels.Cerebellum = 0.500 * probScale;
                labels.Cerebellar_Tonsil = 0.300 * probScale;
                labels.Declive = 0.150 * probScale;
            elseif z < 20 % Middle (e.g., occipital)
                labels.Occipital_Lobe = 0.450 * probScale;
                labels.Precuneus = 0.300 * probScale;
                labels.Cuneus = 0.150 * probScale;
                labels.Brodmann_area_19 = 0.400 * probScale;
                labels.Brodmann_area_18 = 0.200 * probScale;
            else % Superior (e.g., parietal)
                labels.Parietal_Lobe = 0.400 * probScale;
                labels.Precuneus = 0.300 * probScale;
                labels.Superior_Parietal_Lobule = 0.200 * probScale;
                labels.Brodmann_area_7 = 0.350 * probScale;
                labels.Brodmann_area_5 = 0.150 * probScale;
            end
        elseif y < 0 % Mid-posterior (e.g., temporal, parietal)
            if z < -10 % Inferior (e.g., temporal)
                labels.Temporal_Lobe = 0.450 * probScale;
                labels.Inferior_Temporal_Gyrus = 0.300 * probScale;
                labels.Fusiform_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_20 = 0.350 * probScale;
                labels.Brodmann_area_37 = 0.150 * probScale;
            elseif z < 30 % Middle (e.g., temporal-parietal junction)
                labels.Temporal_Lobe = 0.400 * probScale;
                labels.Superior_Temporal_Gyrus = 0.300 * probScale;
                labels.Parietal_Lobe = 0.200 * probScale;
                labels.Brodmann_area_22 = 0.350 * probScale;
                labels.Brodmann_area_39 = 0.150 * probScale;
            else % Superior (e.g., parietal)
                labels.Parietal_Lobe = 0.450 * probScale;
                labels.Inferior_Parietal_Lobule = 0.300 * probScale;
                labels.Supramarginal_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_40 = 0.350 * probScale;
                labels.Brodmann_area_7 = 0.150 * probScale;
            end
        else % Anterior (e.g., frontal)
            if z < -10 % Inferior (e.g., orbitofrontal)
                labels.Frontal_Lobe = 0.500 * probScale;
                labels.Inferior_Frontal_Gyrus = 0.300 * probScale;
                labels.Orbital_Gyrus = 0.150 * probScale;
                labels.Brodmann_area_47 = 0.350 * probScale;
                labels.Brodmann_area_11 = 0.200 * probScale;
            elseif z < 30 % Middle (e.g., lateral frontal)
                labels.Frontal_Lobe = 0.450 * probScale;
                labels.Inferior_Frontal_Gyrus = 0.300 * probScale;
                labels.Middle_Frontal_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_45 = 0.350 * probScale;
                labels.Brodmann_area_46 = 0.150 * probScale;
            else % Superior (e.g., superior frontal)
                labels.Frontal_Lobe = 0.400 * probScale;
                labels.Superior_Frontal_Gyrus = 0.300 * probScale;
                labels.Precentral_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_6 = 0.350 * probScale;
                labels.Brodmann_area_8 = 0.150 * probScale;
            end
        end
    else % Right hemisphere
        if y < -50 % Posterior regions
            if z < -20 % Inferior (e.g., cerebellum)
                labels.Cerebellum = 0.500 * probScale;
                labels.Cerebellar_Tonsil = 0.300 * probScale;
                labels.Declive = 0.150 * probScale;
            elseif z < 20 % Middle (e.g., occipital)
                labels.Occipital_Lobe = 0.400 * probScale;
                labels.Inferior_Occipital_Gyrus = 0.300 * probScale;
                labels.Lingual_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_19 = 0.350 * probScale;
                labels.Brodmann_area_18 = 0.150 * probScale;
            else % Superior (e.g., parietal)
                labels.Parietal_Lobe = 0.400 * probScale;
                labels.Precuneus = 0.300 * probScale;
                labels.Superior_Parietal_Lobule = 0.200 * probScale;
                labels.Brodmann_area_7 = 0.350 * probScale;
                labels.Brodmann_area_5 = 0.150 * probScale;
            end
        elseif y < 0 % Mid-posterior (e.g., temporal, parietal)
            if z < -10 % Inferior (e.g., temporal)
                labels.Temporal_Lobe = 0.450 * probScale;
                labels.Inferior_Temporal_Gyrus = 0.300 * probScale;
                labels.Fusiform_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_20 = 0.350 * probScale;
                labels.Brodmann_area_37 = 0.150 * probScale;
            elseif z < 30 % Middle (e.g., temporal-parietal junction)
                labels.Temporal_Lobe = 0.400 * probScale;
                labels.Superior_Temporal_Gyrus = 0.300 * probScale;
                labels.Parietal_Lobe = 0.200 * probScale;
                labels.Brodmann_area_22 = 0.350 * probScale;
                labels.Brodmann_area_39 = 0.150 * probScale;
            else % Superior (e.g., parietal)
                labels.Parietal_Lobe = 0.450 * probScale;
                labels.Inferior_Parietal_Lobule = 0.300 * probScale;
                labels.Supramarginal_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_40 = 0.350 * probScale;
                labels.Brodmann_area_7 = 0.150 * probScale;
            end
        else % Anterior (e.g., frontal)
            if z < -10 % Inferior (e.g., orbitofrontal)
                labels.Frontal_Lobe = 0.500 * probScale;
                labels.Inferior_Frontal_Gyrus = 0.300 * probScale;
                labels.Orbital_Gyrus = 0.150 * probScale;
                labels.Brodmann_area_47 = 0.350 * probScale;
                labels.Brodmann_area_11 = 0.200 * probScale;
            elseif z < 30 % Middle (e.g., lateral frontal)
                labels.Frontal_Lobe = 0.450 * probScale;
                labels.Inferior_Frontal_Gyrus = 0.300 * probScale;
                labels.Middle_Frontal_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_45 = 0.350 * probScale;
                labels.Brodmann_area_46 = 0.150 * probScale;
            else % Superior (e.g., superior frontal)
                labels.Frontal_Lobe = 0.400 * probScale;
                labels.Superior_Frontal_Gyrus = 0.300 * probScale;
                labels.Precentral_Gyrus = 0.200 * probScale;
                labels.Brodmann_area_6 = 0.350 * probScale;
                labels.Brodmann_area_8 = 0.150 * probScale;
            end
        end
    end

    % If no labels were assigned, return Unknown
    if isempty(fieldnames(labels))
        labels.Unknown = 1.0;
    end
end